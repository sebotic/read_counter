use rust_htslib::bam::{self, record::Aux, Read};
use rust_htslib::errors::Error as HtslibError; // This import is crucial
//use rayon::prelude::*;
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::process;

use ahash::AHashMap;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        print_usage(&args[0]);
        process::exit(1);
    }
    
    // --- Argument Parsing ---
    let mut input_path_str: Option<String> = None;
    let mut ref_fasta_path_str: Option<String> = None;
    let mut max_records: Option<usize> = None;

    let mut arg_iter = args.iter().skip(1);
    while let Some(arg) = arg_iter.next() {
        match arg.as_str() {
            "-n" | "--limit" => {
                if let Some(val_str) = arg_iter.next() {
                    match val_str.parse::<usize>() {
                        Ok(n) => max_records = Some(n),
                        Err(_) => {
                            eprintln!("Error: --limit value '{}' is not a valid positive integer.", val_str);
                            process::exit(1);
                        }
                    }
                } else {
                    eprintln!("Error: --limit flag requires a number.");
                    process::exit(1);
                }
            },
            _ if arg.starts_with('-') => {
                eprintln!("Error: Unknown flag '{}'", arg);
                print_usage(&args[0]);
                process::exit(1);
            }
            _ => { // Positional arguments
                if input_path_str.is_none() {
                    input_path_str = Some(arg.clone());
                } else if ref_fasta_path_str.is_none() {
                    ref_fasta_path_str = Some(arg.clone());
                } else {
                    eprintln!("Error: Too many positional arguments provided.");
                    print_usage(&args[0]);
                    process::exit(1);
                }
            }
        }
    }

    let input_path_str = input_path_str.ok_or_else(|| {
        eprintln!("Error: Missing required input BAM/CRAM file.");
        print_usage(&args[0]);
        "Missing input file".to_string()
    })?;
    
    // --- BAM/CRAM Reader Setup ---
    let input_path = Path::new(&input_path_str);
    let mut bam_reader = bam::Reader::from_path(input_path)
        .map_err(|e| format!("Error opening BAM/CRAM file '{}': {}", input_path.display(), e))?;

    let file_is_cram = input_path_str.ends_with(".cram") || input_path_str.ends_with(".crai");

    if file_is_cram {
        if let Some(ref_path_str) = ref_fasta_path_str {
            let ref_fasta_path = Path::new(&ref_path_str);
            if let Err(e) = bam_reader.set_reference(ref_fasta_path) {
                return Err(format!(
                    "Error setting reference FASTA '{}' for CRAM file '{}': {}. Ensure FASTA is valid and indexed.",
                    ref_fasta_path.display(),
                    input_path.display(),
                    e
                )
                .into());
            }
        } else {
            println!(
                "Info: No explicit reference FASTA provided for CRAM file '{}'. HTSlib will attempt automatic reference discovery.",
                input_path.display()
            );
        }
    } else if ref_fasta_path_str.is_some() {
        eprintln!(
            "Warning: Reference FASTA provided, but input file '{}' does not appear to be CRAM. The reference will be ignored.",
            input_path.display()
        );
    }
    
    if let Some(limit) = max_records {
        println!("Processing up to {} records from '{}'...", limit, input_path.display());
    } else {
        println!("Processing all records from '{}'...", input_path.display());
    }

    // // --- Core Processing Logic ---
    // let records_iterator: Box<dyn Iterator<Item = Result<bam::Record, HtslibError>> + Send> = 
    //     if let Some(limit) = max_records {
    //         Box::new(bam_reader.records().take(limit))
    //     } else {
    //         Box::new(bam_reader.records())
    //     };

    // --- Combined Phase: Read records and count barcodes directly ---
    println!("Reading records and counting barcodes...");
    let mut barcode_counts: AHashMap<String, usize> = AHashMap::new();
    
    let records_iterator = bam_reader.records();

    // Conditionally apply the limit
    let limited_iterator: Box<dyn Iterator<Item = Result<bam::Record, HtslibError>>> = 
        if let Some(limit) = max_records {
            Box::new(records_iterator.take(limit))
        } else {
            Box::new(records_iterator)
        };

    for record_result in limited_iterator {
        match record_result {
            Ok(record) => match record.aux(b"CB") {
                Ok(Aux::String(bc_str)) => {
                    *barcode_counts.entry(bc_str.to_string()).or_insert(0) += 1;
                },
                Err(HtslibError::BamAuxTagNotFound { .. }) => (), // Tag not found, do nothing
                _ => (), // Other tag types or errors, do nothing
            },
            Err(e) => eprintln!("Error reading BAM/CRAM record: {}. Skipping.", e),
        }
    }

    // --- Output Results (unchanged) ---
    let mut sorted_barcodes: Vec<(String, usize)> = barcode_counts.into_iter().collect();
    sorted_barcodes.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    
    let output_file = File::create("reads_per_barcode")?;
    let mut writer = BufWriter::new(output_file);
    let mut total_barcoded_reads = 0;
    for (barcode, count) in &sorted_barcodes {
        writeln!(writer, "{:>7} {}", count, barcode)?;
        total_barcoded_reads += count;
    }
    writer.flush()?;

    println!(
        "Finished processing. Found {} unique barcodes from a total of {} barcoded reads.",
        sorted_barcodes.len(),
        total_barcoded_reads
    );
    if let Some(limit) = max_records {
        println!("(Scanned a maximum of {} records).", limit);
    }
    println!("Results written to 'reads_per_barcode'");

    Ok(())
}

fn print_usage(program_name: &str) {
    eprintln!("A parallel BAM/CRAM barcode counter.");
    eprintln!("\nUsage:");
    eprintln!("  {} <input.bam_or_cram> [reference.fasta_if_cram] [--limit N | -n N]", program_name);
    eprintln!("\nArguments:");
    eprintln!("  <input.bam_or_cram>    Path to the input file.");
    eprintln!("  [reference.fasta_if_cram]  Optional path to the reference FASTA (required for CRAM).");
    eprintln!("\nOptions:");
    eprintln!("  -n, --limit <N>        Process only the first N records from the file.");
}
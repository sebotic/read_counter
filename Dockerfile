FROM rust:1.91

RUN apt-get update && apt-get install -y clang libclang-dev && rm -rf /var/lib/apt/lists/*

WORKDIR /usr/src/read_counter

RUN git --version

COPY . .

RUN cargo install --path .

CMD ["read_counter"]


# Use a newer Rust official image as builder
FROM rust:1.82-slim-bookworm AS builder

# Install required dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    g++ \
    make \
    cmake \
    pkg-config \
    libssl-dev \
    libclang-dev \
    clang \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /build

# Clone the repository and checkout the specific release
RUN git clone https://github.com/bguo068/ishare.git . && \
    git checkout v0.1.11

# Build the release binaries
RUN cargo build --release --bin gtencode && \
    cargo build --release --bin ibdutils && \
    cargo build --release --bin asibd

# Create final minimal image
FROM debian:bookworm-slim

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    libssl3 \
    zlib1g \
    libbz2-1.0 \
    liblzma5 \
    libcurl4 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Copy binaries from builder
COPY --from=builder /build/target/release/gtencode /usr/local/bin/
COPY --from=builder /build/target/release/ibdutils /usr/local/bin/
COPY --from=builder /build/target/release/asibd /usr/local/bin/

# Set working directory
WORKDIR /data

# Default command shows help for gtencode
# CMD ["gtencode", "--help"]

# ==============================================================================
# STAGE 1: High-Performance Build Environment
# ==============================================================================
FROM ubuntu:24.04 AS builder

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install essential compilation tools and OpenMP runtime dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

# Copy source code into the build container
COPY main.* .
COPY Makefile .

# Compile the headless engine with aggressive compiler optimizations
# -O3: Enables loop vectorization, unrolling, and aggressive optimizations
# -DHEADLESS: Decouples the engine from X11/OpenGL graphics libraries
# -fopenmp: Enables shared-memory multi-threading support
# -march=native / -mtune=native: Maximizes instruction set usage for target microarchitecture
RUN make core-release

# ==============================================================================
# STAGE 2: Minimalist Enterprise Runtime Environment
# ==============================================================================
FROM ubuntu:24.04 AS runtime

# Install only the light shared libraries required to execute OpenMP binaries
RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy the optimized binary from the builder stage (keeps final layer slim)
COPY --from=builder /build/pde_engine_headless .

# Create a non-privileged system user for secure container execution in AWS Batch
RUN useradd -m hpc_user && chown -R hpc_user:hpc_user /app
USER hpc_user

# Set strict OpenMP environment variables to control hardware thread affinity.
# These maximize L1/L2 cache locality and prevent core-hopping inside EC2.
ENV OMP_PLACES=cores
ENV OMP_PROC_BIND=close

# Default execution entrypoint
ENTRYPOINT ["./pde_engine_headless"]

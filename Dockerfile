# Stage 1: Build the Conda environment
FROM continuumio/miniconda3 AS builder

# Copy the environment YAML file into the container
COPY squidpipe_env.yaml /tmp/environment.yml

# Create the Conda environment and remove unnecessary files to save space
RUN conda env create -f /tmp/environment.yml && \
    conda clean --all --yes && \
    rm -rf /opt/conda/pkgs/* /root/.cache/pip

# Stage 2: Create the final lightweight image
FROM debian:bullseye-slim

# Copy only the Conda environment from the builder stage
COPY --from=builder /opt/conda/envs/squidpipe_env /opt/conda/envs/squidpipe_env

# Install minimal dependencies required to run Conda
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        bzip2 \
        ca-certificates \
        wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set the PATH to include the Conda environment
ENV PATH="/opt/conda/envs/squidpipe_env/bin:$PATH"

# Set the default shell to bash
SHELL ["/bin/bash", "-c"]

# Directly activate the environment in the CMD step
CMD ["conda activate squidpipe_env"]
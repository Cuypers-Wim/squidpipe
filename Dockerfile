# Stage 1: Build the Conda environment
FROM continuumio/miniconda3 AS builder

# Copy the environment YAML file into the container
COPY squidpipe_env.yaml /tmp/environment.yml

# Create the Conda environment and remove cache to save space
RUN conda env create -f /tmp/environment.yml && \
    conda clean --all --yes && \
    rm -rf /opt/conda/pkgs/* /root/.cache/pip

# Stage 2: Create the final lightweight image
FROM continuumio/miniconda3

# Copy the prepared Conda environment from the builder stage
COPY --from=builder /opt/conda /opt/conda

# Set the PATH to include the Conda environment
ENV PATH="/opt/conda/bin:$PATH"

# Activate the environment by default (replace 'your_env_name' with actual name)
SHELL ["conda", "run", "-n", "squidpipe_env", "/bin/bash", "-c"]

# Set a default command (optional, replace as needed)
CMD ["conda", "env list"]

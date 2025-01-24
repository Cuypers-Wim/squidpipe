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

# Initialize Conda and add environment activation to bashrc
RUN conda init && echo "conda activate viral_metagenomics_env" >> ~/.bashrc

# Use default shell (no need for exec bash in SHELL directive)
SHELL ["bash", "-c"]

# Start an interactive shell when container runs
CMD ["bash"]

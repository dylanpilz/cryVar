FROM mambaorg/micromamba:latest

# Set working directory
WORKDIR /workspace

# Copy environment file
COPY environment.yml /workspace/environment.yml

# Install dependencies using micromamba into base environment
RUN micromamba install -y -n base -f /workspace/environment.yml && \
    micromamba clean -afy && \
    rm -rf /opt/conda/pkgs/cache/* /tmp/* /var/tmp/*

# Set default command
CMD ["/bin/bash"]


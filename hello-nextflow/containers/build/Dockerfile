# Start from the micromamba base docker image
FROM mambaorg/micromamba:1.5.10-noble
# Copy the conda.yml file into the container
COPY --chown=$MAMBA_USER:$MAMBA_USER conda.yml /tmp/conda.yml
# Install various utilities for Nextflow to use and the packages in the conda.yml file
RUN micromamba install -y -n base -f /tmp/conda.yml \
    && micromamba install -y -n base conda-forge::procps-ng \
    && micromamba env export --name base --explicit > environment.lock \
    && echo ">> CONDA_LOCK_START" \
    && cat environment.lock \
    && echo "<< CONDA_LOCK_END" \
    && micromamba clean -a -y
# Run the container as the root user
USER root
# Set the PATH environment variable to include the micromamba installation directory
ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
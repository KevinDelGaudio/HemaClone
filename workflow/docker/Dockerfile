# CHIP workflow image: pinned tools plus the workflow code and default config.
# Build from the repository root (the build context is the whole repository,
# filtered by .dockerignore):
#   bash workflow/chip.sh build
# which runs: docker build -f workflow/docker/Dockerfile -t chip-workflow:<VERSION> .
# Behind an HTTP proxy, add:
#   bash workflow/chip.sh build --build-arg HTTPS_PROXY=$HTTPS_PROXY --build-arg HTTP_PROXY=$HTTP_PROXY

FROM mambaorg/micromamba:2.0.5

# 1. Tools. This layer is cached: changing only the code does not reinstall it.
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/docker/environment.yaml /tmp/environment.yaml
RUN micromamba install -y -n base -f /tmp/environment.yaml \
    && micromamba clean --all --yes

# 2. Workflow code and default config
USER root
COPY VERSION /opt/chip-workflow/VERSION
COPY workflow /opt/chip-workflow/workflow
COPY config/config.yaml config/multiqc_config.yaml /opt/chip-workflow/config/
COPY config/profile /opt/chip-workflow/config/profile

# The launcher copied out of the image defaults to this image's tag
RUN version="$(tr -d '[:space:]' < /opt/chip-workflow/VERSION)" \
    && sed -i "s|^DEFAULT_IMAGE=.*|DEFAULT_IMAGE=\"chip-workflow:${version}\"|" /opt/chip-workflow/workflow/chip.sh \
    && chmod -R a+rX /opt/chip-workflow
USER $MAMBA_USER

# Tools on PATH for any user ID (the container runs as the host user)
ENV PATH=/opt/conda/bin:$PATH \
    HOME=/tmp

WORKDIR /tmp

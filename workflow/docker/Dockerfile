# HemaClone image: pinned tools plus the workflow code and default config.
# Build from the repository root (the build context is the whole repository,
# filtered by .dockerignore):
#   bash workflow/hemaclone.sh build
# which runs: docker build -f workflow/docker/Dockerfile -t hemaclone:<VERSION> .
# Behind an HTTP proxy, add:
#   bash workflow/hemaclone.sh build --build-arg HTTPS_PROXY=$HTTPS_PROXY --build-arg HTTP_PROXY=$HTTP_PROXY

FROM mambaorg/micromamba:2.0.5

# 1. Tools. This layer is cached: changing only the code does not reinstall it.
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/docker/environment.yaml /tmp/environment.yaml
RUN micromamba install -y -n base -f /tmp/environment.yaml \
    && micromamba clean --all --yes

# 1b. somalier in its own environment (its htslib pin differs from the tools above)
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/docker/somalier.yaml /tmp/somalier.yaml
RUN micromamba create -y -f /tmp/somalier.yaml \
    && micromamba clean --all --yes

# 2. Workflow code and default config
USER root
COPY VERSION /opt/hemaclone/VERSION
COPY workflow /opt/hemaclone/workflow
COPY config/config.yaml config/multiqc_config.yaml /opt/hemaclone/config/
COPY config/profile /opt/hemaclone/config/profile

# The launcher copied out of the image defaults to this image's tag
RUN version="$(tr -d '[:space:]' < /opt/hemaclone/VERSION)" \
    && sed -i "s|^DEFAULT_IMAGE=.*|DEFAULT_IMAGE=\"hemaclone:${version}\"|" /opt/hemaclone/workflow/hemaclone.sh \
    && chmod -R a+rX /opt/hemaclone
USER $MAMBA_USER

# Tools on PATH for any user ID (the container runs as the host user)
ENV PATH=/opt/conda/bin:$PATH \
    HOME=/tmp

WORKDIR /tmp

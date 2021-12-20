FROM rocker/r-ver:3.5.2

LABEL maintainer "Diego Garrido-Martín <diego.garrido@crg.eu>" \
      version "1.1.0" \
      description "Docker image for sqtlseeker2-nf"

# Install needed tools
RUN apt-get update --fix-missing -qq && apt-get install -y -q \
       libcurl4-openssl-dev \
       libssl-dev \
       libssh2-1-dev \
       libxml2-dev \
       zlib1g-dev \
       procps \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

## Install R packages
RUN R -e 'install.packages(c("devtools", "optparse", "data.table"), repos="http://cloud.r-project.org/")'

## Install sQTLseekeR2
RUN R -e 'devtools::install_github("guigolab/sQTLseekeR2")'

## Change ENTRYPOINT
ENTRYPOINT /bin/bash

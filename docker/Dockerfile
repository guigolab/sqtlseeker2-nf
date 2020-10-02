FROM debian:stable

LABEL maintainer "Diego Garrido-Martín <diego.garrido@crg.eu>" \
      version "1.0.0" \
      description "Docker image for sqtlseeker2-nf"

# Install needed tools
RUN apt-get update --fix-missing -qq && apt-get install -y -q \
    libssl-dev \
    libcurl4-openssl-dev \
    libnlopt-dev \
    r-base \
    procps \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN R -e 'install.packages(c("devtools", "optparse"), repos="http://cloud.r-project.org/"); \
          source("https://bioconductor.org/biocLite.R"); \
          devtools::install_github("guigolab/sQTLseekeR2")'


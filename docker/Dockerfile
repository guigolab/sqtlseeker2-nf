FROM rocker/r-ver:3.5.2

RUN apt-get update --fix-missing -qq && apt-get install -y -q \
       libcurl4-openssl-dev \
       libssl-dev \
       libssh2-1-dev \ 
       libxml2-dev \
       zlib1g-dev \
       procps \
       python3 \
       python3-pip \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

## Install R packages
RUN R -e 'install.packages(c("devtools", "optparse", "data.table"), repos="http://cloud.r-project.org/")'

## Install python packages
RUN pip3 install numpy pandas scipy scikit-learn

## Install sQTLseekeR2.int
RUN R -e 'devtools::install_github("guigolab/sQTLseekeR2@interaction")'

## Change ENTRYPOINT
ENTRYPOINT /bin/bash


# Docker image for FAVITES-ABM-HIV-SD-Lite
FROM ubuntu:20.04

# Set up environment and install dependencies
RUN apt-get update && apt-get -y upgrade && \
    DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get install -y cmake libcurl4-openssl-dev libssl-dev libxml2-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev git parallel python3 python3-pip r-base-core wget && \

    # Install required Python packages
    pip3 install scipy && \
    pip3 install treesap && \
    pip3 install openpyxl && \
    pip3 install networkx && \
    pip3 install treeswift && \
    pip3 install niemads && \

    # Install required R packages
    R -e "install.packages('gtools')" && \
    R -e "install.packages('ensurer')" && \
    R -e "install.packages('truncnorm')" && \
    R -e "install.packages('assertthat')" && \
    R -e "install.packages('rlang')" && \
    R -e "install.packages('readxl')" && \
    R -e "install.packages('tidyverse')" && \
    R -e "install.packages('tictoc')" && \
    R -e "install.packages('fastRG')" && \
    R -e "install.packages('mice')" && \
    R -e "update.packages()" && \

    # Install abm_hiv-HRSA_SD
    wget -q https://github.com/mathematica-pub/abm_hiv/archive/refs/heads/HRSA_SD.zip && \
    unzip HRSA_SD.zip && rm HRSA_SD.zip && \
    mv abm_hiv-HRSA_SD /usr/local/bin/abm_hiv-HRSA_SD && \

    # Install CoaTran
    wget -qO- "https://github.com/niemasd/CoaTran/archive/refs/tags/0.0.3.tar.gz" | tar -zx && \
    cd CoaTran-* && \
    make && \
    mv coatran_* /usr/local/bin/ && \
    cd .. && \
    rm -rf CoaTran-* && \

    # Install FAVITES-ABM-HIV-SD-Lite in /home directory
    cd /home && \
    git clone https://github.com/niemasd/FAVITES-ABM-HIV-SD-Lite.git && \

    # Clean up
    rm -rf /root/.cache /tmp/*

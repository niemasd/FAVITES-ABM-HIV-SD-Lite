# Docker image for FAVITES-ABM-HIV-SD-Lite
FROM ubuntu:20.04

# Set up environment and install dependencies
RUN apt-get update && apt-get -y upgrade && \
    DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt install -y --no-install-recommends software-properties-common dirmngr wget && \
    wget -qO- "https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc" >> "/etc/apt/trusted.gpg.d/cran_ubuntu_key.asc" && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" && \
    DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get install -y --no-install-recommends cmake g++ gfortran git libcurl4-openssl-dev liblapack-dev libopenblas-dev libssl-dev libxml2-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev make parallel python3 python3-pip r-base-core && \

    # Install required Python packages
    pip3 install scipy && \
    pip3 install treesap && \
    pip3 install openpyxl && \
    pip3 install networkx && \
    pip3 install treeswift && \
    pip3 install niemads && \

    # Install required R packages
    R -e "withCallingHandlers(install.packages('devtools'), warning = function(w) stop(w))" && \
    R -e "withCallingHandlers(install.packages('gtools'), warning = function(w) stop(w))" && \
    R -e "withCallingHandlers(install.packages('truncnorm'), warning = function(w) stop(w))" && \
    R -e "withCallingHandlers(install.packages('assertthat'), warning = function(w) stop(w))" && \
    R -e "withCallingHandlers(install.packages('rlang'), warning = function(w) stop(w))" && \
    R -e "withCallingHandlers(install.packages('readxl'), warning = function(w) stop(w))" && \
    R -e "withCallingHandlers(install.packages('tidyverse'), warning = function(w) stop(w))" && \
    R -e "withCallingHandlers(install.packages('tictoc'), warning = function(w) stop(w))" && \
    R -e "withCallingHandlers(install.packages('fastRG'), warning = function(w) stop(w))" && \
    R -e "withCallingHandlers(install.packages('mice'), warning = function(w) stop(w))" && \
    R -e "update.packages()" && \
    R -e "require(devtools); withCallingHandlers(install_version('ensurer', version = '1.1', repos = 'http://cran.us.r-project.org'), warning = function(w) stop(w))" && \

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

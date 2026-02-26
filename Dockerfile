FROM bioconductor/bioconductor_docker:RELEASE_3_22

MAINTAINER Vinh Tran

ENV DEBIAN_FRONTEND noninteractive

# Install Ubuntu packages
RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
    pandoc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
	  libxml2-dev \
    curl \
    xtail \
    wget

# Download and install ShinyServer (latest version)
RUN wget --no-verbose https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb

# Install R dependencies (that cannot be install automatically)
RUN R -e "install.packages(c('bsplus','xfun'), dependencies=TRUE, repos='http://cran.rstudio.com/')"

# Install PhyloProfile package
RUN R -e 'BiocManager::install(c("ExperimentHub", "PhyloProfile", "PhyloProfileData"))'

# Copy configuration files into the Docker image
COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf
COPY /app /srv/shiny-server/

# Generate local taxonomy DB
RUN if [ ! -f /srv/shiny-server/data/preProcessedTaxonomy.txt ]; then \
        Rscript /srv/shiny-server/generate_taxonomy.R; \
    else \
        echo "preProcessedTaxonomy.txt already exists, skipping R script"; \
    fi
RUN chmod 777 /srv/shiny-server/data/*.txt

# Replace line in server.R during build
#RUN sed -i "s|homePath = c(wd='~/')|homePath <- c(taxDB='/srv/shiny-server/data')|" /srv/shiny-server/server.R

# Checking
RUN echo "=== Listing /srv/shiny-server ===" && \
    ls -la /srv/shiny-server && \
    echo "=== Listing data folder ===" && \
    ls -la /srv/shiny-server/data && \
    echo "===== Confirm homePath line in server.R =====" && \
    grep "homePath" /srv/shiny-server/server.R

# Make the ShinyApp available at port 80
EXPOSE 80

# Copy further configuration files into the Docker image
COPY shiny-server.sh /usr/bin/shiny-server.sh
RUN chmod +x /usr/bin/shiny-server.sh
ENTRYPOINT ["/usr/bin/shiny-server.sh"]

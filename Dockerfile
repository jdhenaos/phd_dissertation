# Use the official R image with version >4.4 as the base
FROM rocker/r-ver:4.2.3

# Install necessary dependencies for RStudio Server
RUN apt-get update && apt-get install -y \
    # Install basic utilities like bash, wget, curl, etc.
    bash \
    wget \
    curl \
    git \
    sudo \
    vim \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libpq-dev \
    gdebi-core \
    # Install RStudio Server dependencies
    libfontconfig1 \
    libx11-6 \
    libxtst6 \
    libxrender1 \
    && wget http://archive.ubuntu.com/ubuntu/pool/main/o/openssl/libssl1.1_1.1.1f-1ubuntu2_amd64.deb \
    && gdebi -n libssl1.1_1.1.1f-1ubuntu2_amd64.deb \
    && rm libssl1.1_1.1.1f-1ubuntu2_amd64.deb \
    # Install RStudio Server
    && wget https://download2.rstudio.org/server/focal/amd64/rstudio-server-2024.12.1-563-amd64.deb  \
    && gdebi -n rstudio-server-2024.12.1-563-amd64.deb \
    && rm rstudio-server-2024.12.1-563-amd64.deb \
    # Clean up to reduce image size
    && apt-get clean
 
RUN apt-get update && apt-get install -y \
  libcurl4-openssl-dev \
  libxml2-dev \
  libssl-dev \
  libz-dev libbz2-dev \
  liblzma-dev libpcre2-dev \
  build-essential \
  r-base-dev

RUN apt-get update && apt-get install -y libglpk40

# Install R packages commonly used in data science
RUN R -e "install.packages(c('devtools', 'ggplot2==3.5.1'), repos='https://cloud.r-project.org')"

# Set up RStudio server user
RUN useradd -m rstudio && echo "rstudio:rstudio" | chpasswd && adduser rstudio sudo

# Expose RStudio Server port (8787)
EXPOSE 8787

# Set working directory
WORKDIR /home/rstudio

# Start RStudio Server
CMD ["/usr/lib/rstudio-server/bin/rserver", "--server-daemonize", "0"]

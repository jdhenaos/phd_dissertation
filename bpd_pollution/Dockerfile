# Use the official R image with version >4.4 as the base
FROM rocker/r-ver:4.4.0

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

# Install Quarto CLI
RUN wget https://quarto.org/download/latest/quarto-linux-amd64.deb && \
    apt-get update && \
    apt-get install -y ./quarto-linux-amd64.deb && \
    rm quarto-linux-amd64.deb

# Install R packages commonly used in data science
RUN R -e "install.packages(c('languageserver','httpgd','devtools', 'car','rmarkdown','ggplot2==3.5.1'), repos='https://cloud.r-project.org')"

# Install a system library often needed for 'httpgd' (the graphic device we'll use)
RUN apt-get update && apt-get install -y \
    libfontconfig1-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /home/project


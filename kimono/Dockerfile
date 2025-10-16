# Use a stable R base
FROM rocker/r-ver:4.2.3

# -------------------------
# System dependencies
# -------------------------
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libgit2-dev \
    libssh2-1-dev \
    build-essential \
    git \
    wget \
    pandoc \
    pandoc-citeproc \
    && rm -rf /var/lib/apt/lists/*

# -------------------------
# Install R packages
# -------------------------
RUN R -e "install.packages(c('remotes', 'renv'))"
RUN R -e "remotes::install_github('benmarwick/rrtools')"

# -------------------------
# Set working directory
# -------------------------
WORKDIR /project

# -------------------------
# Copy renv lockfiles (optional, if present)
# -------------------------
#COPY renv.lock renv/activate.R /project/
#RUN if [ -f renv.lock ]; then R -e "renv::restore()"; fi

# -------------------------
# Install Quarto CLI (optional, if using Quarto)
# -------------------------
RUN wget https://github.com/quarto-dev/quarto-cli/releases/download/v1.5.57/quarto-1.5.57-linux-amd64.deb \
    && dpkg -i quarto-1.5.57-linux-amd64.deb \
    && rm quarto-1.5.57-linux-amd64.deb

# Default command
CMD ["R"]


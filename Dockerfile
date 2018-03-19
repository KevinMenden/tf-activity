FROM continuumio/anaconda3

# Install container-wide requrements gcc, pip, zlib, libssl, make, libncurses, fortran77, g++, R
RUN apt-get update && \
    apt-get install -y \
    git \
    gcc \
    g++ \
    make \
    perl \
    wget \
    build-essential \
    libbz2-dev \
    libcurl4-openssl-dev \
    #libgsl-dev \
    #libgsl2 \
    liblzma-dev \
    libncurses5-dev \
    libpcre3-dev \
    libreadline-dev \
    libssl-dev \
    curl \
    zip \
    unzip \
    zlib1g-dev \
    python-dev

# Install HOMER
RUN mkdir /opt/homer
RUN curl -fsSL http://homer.ucsd.edu/homer/configureHomer.pl -o /opt/homer/configureHomer.pl
RUN perl /opt/homer/configureHomer.pl -install

# Install MEME
RUN curl -fsSL meme-suite.org/meme-software/4.12.0/meme_4.12.0.tar.gz -o /opt/meme_4.12.0.tar.gz
RUN cd /opt/ && tar -xzf /opt/meme_4.12.0.tar.gz
RUN cd /opt/meme_4.12.0; ./configure --prefix=/opt/meme --enable-build-libxml2 --enable-build-libxslt; \
    make; make install

# Update conda
RUN conda update -n base conda

# Install biopython
RUN conda install -c bioconda biopython

# Install R
#RUN conda install -c r r

RUN apt-get update && apt-get install -y\
    gfortran

# Install bed tools
RUN conda install -c bioconda bedtools

RUN apt-get update && apt-get install -y \
     libcurl4-gnutls-dev

# Install R
ENV R_VERSION="R-3.4.3"
RUN curl -fsSL https://cran.r-project.org/src/base/R-3/${R_VERSION}.tar.gz -o /opt/${R_VERSION}.tar.gz && \
    tar xvzf /opt/${R_VERSION}.tar.gz -C /opt/ && \
    cd /opt/${R_VERSION};./configure --with-x=no;make;make install && \
rm /opt/${R_VERSION}.tar.gz

# TODO decide whether to keep this
# Install biomaRt - only necessary if this is to be included in pipeline ... (probably not though)
#RUN conda install -c bioconda bioconductor-biomart

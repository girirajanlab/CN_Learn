From rocker/r-ver:3.4.4

# Install basic LINUX tools and Java8
RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y \
    autoconf \
    gcc \
    git \
    make \
    ssh \
    wget \
    vim \
    build-essential \
    software-properties-common \
    ca-certificates \
    libssl-dev \
    libcurl4-openssl-dev \
    libatlas-base-dev \
    libmariadbclient-dev \
    libffi-dev \
    libxml2-dev \
    libncurses5-dev \
    libsm6 \
    libxrender1 \
    libfontconfig1 \
    libxt6 \
    libtcl8.6 \
    libtk8.6 \
    python3-pip \
    && add-apt-repository ppa:webupd8team/java -y \
    && apt-get update \
    && echo oracle-java8-installer shared/accepted-oracle-license-v1-1 select true | /usr/bin/debconf-set-selections \
    && apt-get install -y --allow-unauthenticated oracle-java8-installer \
    oracle-java8-set-default \
    && apt-get autoremove -y \
    && apt-get clean -y

# Install all the required R packages
RUN R -e "install.packages(c('Rcpp', 'fs', 'usethis'), repos='http://cran.rstudio.com/')" && \
    R -e "install.packages(c('pkgload', 'xml2', 'htmltools'), repos='http://cran.rstudio.com/')" && \
    R -e "install.packages('devtools', repos='http://cran.rstudio.com/')" && \
    R -e "install.packages(c('nnls', 'Hmisc', 'mgcv', 'plyr'), repos='http://cran.rstudio.com/')" && \
    R -e "install.packages(c('sqldf', 'matrixStats'), repos='http://cran.rstudio.com/')" && \
    R -e "install.packages('caret', repos='http://cran.rstudio.com/')" && \
    R -e "source('https://bioconductor.org/biocLite.R'); biocLite('Biostrings'); biocLite('rtracklayer'); \
          biocLite('GenomeInfoDb'); biocLite('IRanges'); biocLite('BSgenome'); biocLite('GenomicAlignments'); \
          biocLite('BiocParallel')" && \
    R -e "library(devtools); source('https://bioconductor.org/biocLite.R'); install_github('yuchaojiang/CODEX/package')"

# Download and install CN_Learn from Github
WORKDIR /opt/tools
RUN git clone --recursive https://github.com/girirajanlab/CN_Learn.git

# Install python
WORKDIR /opt/tools/CN_Learn/software
RUN wget https://www.python.org/ftp/python/3.7.3/Python-3.7.3.tgz && \
    tar xzf Python-3.7.3.tgz && \
    cd Python-3.7.3 && \
    ./configure && make && make install

# Install all the required python packages using PIP3
RUN pip3 install -U 'numpy==1.16.1' && \
    pip3 install -U 'Cython==0.27.3' && \
    pip3 install -U 'pandas==0.24.2' && \
    pip3 install -U 'scipy==1.2.1' && \
    pip3 install -U 'scikit-learn==0.20.3' && \
    pip3 install -U 'pydot==1.4.1'

# Install the tools required to run individual CNV callers
WORKDIR /opt/tools/CN_Learn/software
RUN tar -zxvf gatk-3.5.tar.gz && \
    tar -zxvf xhmm.tar.gz && \
    tar -zxvf clamms.tar.gz && \
    wget http://psychgen.u.hpc.mssm.edu/plinkseq_downloads/plinkseq-x86_64-latest.zip && \
    unzip plinkseq-x86_64-latest.zip && \
    cd plinkseq-0.10 && \
    wget http://psychgen.u.hpc.mssm.edu/plinkseq_resources/hg19/seqdb.hg19.gz && \
    gunzip seqdb.hg19.gz

# Install htslib
WORKDIR /opt/tools/CN_Learn/software
RUN wget -c https://github.com/samtools/htslib/archive/1.3.2.tar.gz && \
    tar -zxvf 1.3.2.tar.gz && \
    mv htslib-1.3.2 htslib && \
    cd htslib && \
    autoreconf && \
    ./configure && make && make install

# Install samtools
WORKDIR /opt/tools/CN_Learn/software
RUN wget -c https://github.com/samtools/samtools/archive/1.3.1.tar.gz && \
    tar -zxvf 1.3.1.tar.gz && \
    cd samtools-1.3.1 && \
    make && make install

# Install bedtools
WORKDIR /opt/tools/CN_Learn/software
RUN apt-get install -y python-pip bedtools && \
    wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz && \
    tar -zxvf bedtools-2.27.1.tar.gz && \
    cd bedtools2 && \
    make

# Install the Linux library used by KentUtils and cleanup the tar files
WORKDIR /opt/tools/CN_Learn/software
RUN apt install ./libpng12-0_1.2.54_amd64.deb && \
    rm gatk-3.5.tar.gz  && rm xhmm.tar.gz  && rm clamms.tar.gz && \
    rm Python-3.7.3.tgz && rm 1.3.2.tar.gz && rm 1.3.1.tar.gz  && \
    rm bedtools-2.27.1.tar.gz && rm plinkseq-x86_64-latest.zip


ENV CLAMMS_DIR=/opt/tools/CN_Learn/software/clamms/

WORKDIR /opt/tools

CMD ["/bin/bash"]

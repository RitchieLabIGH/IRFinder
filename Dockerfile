FROM ubuntu:bionic

LABEL version=v2.0.0


ENV LD_LIBRARY_PATH="/usr/local/lib/:$LD_LIBRARY_PATH"
ENV PYTHONNOUSERSITE="true"
ENV PATH="/Utils/bin/:${PATH}"

ARG DEBIAN_FRONTEND=noninteractive

### All the dependencies
RUN apt-get update && apt-get install -y gnupg2 && apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys \
        E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    && apt-get update \
    && apt-get install -y --no-install-recommends software-properties-common \
    && add-apt-repository -y 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' \
    && apt-get purge -y --autoremove software-properties-common && \
    apt-get update && \
    apt-get -y upgrade && \
    export DEBIAN_FRONTEND=noninteractive && \ 
    apt-get install -qy make build-essential libxml2-dev libcurl4-openssl-dev gcc bedtools samtools git gzip \
		zlib1g gawk libz-dev wget libboost-iostreams-dev python3.6 apt-transport-https software-properties-common \
		 python3-pip r-base && \
    apt-get clean && apt-get purge && \
	rm -rf /var/lib/apt/lists/* 
	
RUN	pip3 install -U --no-cache-dir numpy pandas \
    	scikit-learn scipy \
    	statsmodels

RUN mkdir -p /Utils/bin/ && \
    cd /Utils/ && \
    git clone https://github.com/alexdobin/STAR.git && \
    cd STAR/source && \
    make STAR && \
    ln -s /Utils/STAR/source/STAR /Utils/bin/STAR && \
	cd /Utils && \
	git clone https://github.com/comprna/SUPPA.git && \
	echo "#!/usr/bin/env python3" > /Utils/SUPPA/suppa.py.tmp && \
	cat /Utils/SUPPA/suppa.py >> /Utils/SUPPA/suppa.py.tmp && \
	mv /Utils/SUPPA/suppa.py.tmp /Utils/SUPPA/suppa.py && \
	chmod +x /Utils/SUPPA/suppa.py && \
	ln -s /Utils/SUPPA/suppa.py /Utils/bin/suppa.py

RUN cd /Utils/ && git clone https://github.com/lh3/minimap2 && \
	cd minimap2 && make && \
	ln -s /Utils/minimap2/minimap2 /Utils/bin/minimap2	

	 
RUN	Rscript -e 'options(warn=2); if (!requireNamespace("BiocManager", quietly = TRUE)) { install.packages("BiocManager",lib="/usr/lib/R/library/", force=TRUE) } ; BiocManager::install(c("BiocManager", "tximport", "readr", "RCurl", "DESeq2"), force=TRUE,ask=F,lib="/usr/lib/R/library/", quiet=F)' && \
	Rscript -e 'installed.packages()' | awk  'BEGIN {v=0} $1=="Version" {v=1; } v==1 && $1 == "DESeq2" { gsub("\"", ""); print $2;v=0 } '

ADD "https://www.random.org/cgi-bin/randbyte?nbytes=10&format=h" skipcache	

COPY ./bin /IRFinder/bin
COPY ./REF /IRFinder/REF
COPY ./src /IRFinder/src
COPY ./install.sh /IRFinder/
RUN    cd /IRFinder/ && \
	./install.sh

	 

    
ENTRYPOINT ["IRFinder"]


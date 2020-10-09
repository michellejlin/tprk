FROM ubuntu:bionic

## Install Julia
RUN set -eux; \
	apt-get update; \
	apt-get install -y --no-install-recommends \
		ca-certificates \
# ERROR: no download agent available; install curl, wget, or fetch
		curl \
	; \
	rm -rf /var/lib/apt/lists/*
ENV JULIA_PATH /usr/local/julia
ENV PATH $JULIA_PATH/bin:$PATH

# https://julialang.org/juliareleases.asc
# Julia (Binary signing key) <buildbot@julialang.org>
ENV JULIA_GPG 3673DF529D9049477F76B37566E3C7DC03D6E495

# https://julialang.org/downloads/
ENV JULIA_VERSION 0.6.4

RUN set -eux; \
	\
	savedAptMark="$(apt-mark showmanual)"; \
	if ! command -v gpg > /dev/null; then \
		apt-get update; \
		apt-get install -y --no-install-recommends \
			gnupg \
			dirmngr \
		; \
		rm -rf /var/lib/apt/lists/*; \
	fi; \
	\
# https://julialang.org/downloads/#julia-command-line-version
# https://julialang-s3.julialang.org/bin/checksums/julia-1.5.0-rc1.sha256
# this "case" statement is generated via "update.sh"
	dpkgArch="$(dpkg --print-architecture)"; \
	case "${dpkgArch##*-}" in \
#  amd64
# 		amd64) tarArch='x86_64'; dirArch='x64'; sha256='a4ea36aa86269116992393067e5afc182707cb4f26eac9fddda08e04a9c7b94d' ;; \
#  arm64v8
# 		arm64) tarArch='aarch64'; dirArch='aarch64'; sha256='7e9f3fac46264a2c861c542adce9d6b47276976dab40cbe19ca7ed2c97a82b66' ;; \
#  i386
# 		i386) tarArch='i686'; dirArch='x86'; sha256='ebc76bc879f722375e658a2c3cd43304ee8b05035fa46b8ad7c5c8eef1091a42' ;; \
#		*) echo >&2 "error: current architecture ($dpkgArch) does not have a corresponding Julia binary release"; exit 1 ;; \
	esac; \
	\
	folder="$(echo "$JULIA_VERSION" | cut -d. -f1-2)"; \
	#curl -fL -o julia.tar.gz.asc "https://julialang-s3.julialang.org/bin/linux/${dirArch}/${folder}/julia-${JULIA_VERSION}-linux-${tarArch}.tar.gz.asc"; \
	#curl -fL -o julia.tar.gz     "https://julialang-s3.julialang.org/bin/linux/x86/${folder}/julia-${JULIA_VERSION}-linux-${tarArch}.tar.gz"; \
	curl -fL -o julia.tar.gz     "https://julialang-s3.julialang.org/bin/linux/x64/0.6/julia-0.6.4-linux-x86_64.tar.gz"; \
	\
	#echo "${sha256} *julia.tar.gz" | sha256sum -c -; \
	\
	#export GNUPGHOME="$(mktemp -d)"; \
	#gpg --batch --keyserver ha.pool.sks-keyservers.net --recv-keys "$JULIA_GPG"; \
	#gpg --batch --verify julia.tar.gz.asc julia.tar.gz; \
	#command -v gpgconf > /dev/null && gpgconf --kill all; \
	#rm -rf "$GNUPGHOME" julia.tar.gz.asc; \
	#\
	mkdir "$JULIA_PATH"; \
	tar -xzf julia.tar.gz -C "$JULIA_PATH" --strip-components 1; \
	rm julia.tar.gz; \
	\
	apt-mark auto '.*' > /dev/null; \
	[ -z "$savedAptMark" ] || apt-mark manual $savedAptMark; \
	apt-get purge -y --auto-remove -o APT::AutoRemove::RecommendsImportant=false; \
	\
# smoke test
	julia --version

CMD ["julia"]


LABEL org.label-schema.license="GPL-2.0" \
      org.label-schema.vcs-url="https://github.com/rocker-org/r-apt" \
      org.label-schema.vendor="Rocker Project" \
      maintainer="Dirk Eddelbuettel <edd@debian.org>"

## Set a default user. Available via runtime flag `--user docker` 
## Add user to 'staff' group, granting them write privileges to /usr/local/lib/R/site.library
## User should also have & own a home directory (for rstudio or linked volumes to work properly). 
RUN useradd docker \
	&& mkdir /home/docker \
	&& chown docker:docker /home/docker \
	&& addgroup docker staff

RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		software-properties-common \
                ed \
		less \
		locales \
		vim-tiny \
		wget \
		ca-certificates \
        && add-apt-repository --enable-source --yes "ppa:marutter/rrutter3.5" \
	&& add-apt-repository --enable-source --yes "ppa:marutter/c2d4u3.5" 

## Configure default locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

## This was not needed before but we need it now
ENV DEBIAN_FRONTEND noninteractive

# Now install R and littler, and create a link for littler in /usr/local/bin
# Default CRAN repo is now set by R itself, and littler knows about it too
# r-cran-docopt is not currently in c2d4u so we install from source
RUN apt-get update \
        && apt-get install -y --no-install-recommends \
                 littler \
 		 r-base \
 		 r-base-dev \
 		 r-recommended \
  	&& ln -s /usr/lib/R/site-library/littler/examples/install.r /usr/local/bin/install.r \
 	&& ln -s /usr/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r \
 	&& ln -s /usr/lib/R/site-library/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
 	&& ln -s /usr/lib/R/site-library/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
 	&& install.r docopt \
 	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
 	&& rm -rf /var/lib/apt/lists/*

CMD ["bash"]

## Install rJava
RUN apt-get -y update && apt-get install -y \
   default-jdk  r-cran-rjava  r-cran-nloptr libssh2-1-dev

RUN apt-get install libcurl4-openssl-dev 
RUN yes | apt-get install libv8-dev
RUN yes | apt-get install libssl-dev
RUN yes | apt-get install libxml2-dev
RUN yes | apt-get install libudunits2-dev
RUN yes | apt install libgdal-dev
RUN yes | apt-get install libmagick++-dev

## Install extra R packages using requirements.R
RUN R -e "install.packages('JuliaCall',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('reticulate',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('tidyr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('XML', repos = 'http://www.omegahat.net/R')"
RUN R -e "install.packages('openssl',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('BiocManager'); library('BiocManager')"
RUN R -e "BiocManager::install('ShortRead')"
RUN R -e "BiocManager::install('dada2')"
RUN R -e "BiocManager::install('treeio')"
RUN R -e "BiocManager::install('ggtree')"
RUN R -e "BiocManager::install('phylobase')"
RUN R -e "BiocManager::install('randomcoloR')"
RUN R -e "BiocManager::install('tidyverse')"
RUN R -e "BiocManager::install('pegas')"
RUN R -e "BiocManager::install('foreach')"
RUN R -e "BiocManager::install('iterators')"
RUN R -e "BiocManager::install('doParallel')"
#RUN R -e "install.packages('devtools',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('optparse',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggplot2',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('reshape2',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('grid',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('nplr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('dplyr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('scales',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('gridExtra',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('RColorBrewer',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('stringr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('Biostrings',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('phylobase',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('phytools',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('phangorn',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('pegas',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('lubridate',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ape',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('data.table',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('vegan',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('cowplot',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('tibble',dependencies=TRUE, repos='http://cran.rstudio.com/')"



# pip install
RUN apt update && \
    apt install -y python3 ncbi-blast+ && \
    apt install -y python-biopython \
                   python3-pip \
                   python3-pysam \
                   wget \
                   unzip && \
    pip3 install biopython \
                 ete3 \
                 pycosat \
                 PyYAML \
                 requests \
                 numpy \
                 pandas \
                 bokeh \
				 matplotlib \
				 regex

# Install mafft
RUN cd /usr/local/ && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /usr/local/miniconda && \
    rm Miniconda3-latest-Linux-x86_64.sh && \
    ln -s /usr/local/miniconda/bin/conda /usr/local/bin/ && \
    conda init bash && \
    /bin/bash -c "source /root/.bashrc" && \
    conda install -c bioconda bowtie2 samtools=1.7 bedtools bwa mafft bcftools tabix fasttree && \
    conda clean -afy

RUN R -e "BiocManager::install('BiocGenerics')"
RUN R -e "install.packages('foreach', repos='http://R-Forge.R-project.org')"
RUN R -e "install.packages('iterators', repos='http://R-Forge.R-project.org')"
RUN R -e "BiocManager::install('doParallel')"

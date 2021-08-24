FROM rocker/tidyverse:latest

# Run the
RUN sudo apt-get update -y \
&& sudo apt-get install -y libpython2-dev libpython3-dev  dpkg-dev zlib1g-dev \
&& sudo apt-get install -y libssl-dev libffi-dev libglpk-dev curl nano less gawk \
&& sudo apt-get install -y openssh-client libzmq3-dev \
&& sudo apt-get install -y build-essential libcurl4-openssl-dev chrpath libssl-dev \
&& sudo apt-get install -y libxft-dev libfreetype6 libfreetype6-dev \
&& sudo apt-get install -y libfontconfig1 libfontconfig1-dev php php-common gcc git \
&& sudo apt-get install -y imagemagick php-imagick tabix bcftools awscli vsftpd \
&& sudo apt-get install -y vcftools coreutils r-cran-rjava software-properties-common \
&& sudo apt-get install -y g++ software-properties-common libpcre++-dev \
&& sudo apt-get install -y libcurl4-openssl-dev liblzma-dev libbz2-dev \
&& sudo apt-get install -y openjdk-11-jdk libxml2-dev glpk-utils libglpk-dev \
&& sudo apt-get install -y glpk-doc libxml2-dev

# Setup Phantom JS
RUN export PHANTOM_JS="phantomjs-1.9.8-linux-x86_64" \
&& wget https://bitbucket.org/ariya/phantomjs/downloads/$PHANTOM_JS.tar.bz2 \
&& tar xvjf $PHANTOM_JS.tar.bz2 \
&& mv $PHANTOM_JS /usr/local/share \
&& ln -sf /usr/local/share/$PHANTOM_JS/bin/phantomjs /usr/local/bin

# Install Python Dependencies And Code Repositories
RUN sudo apt-get install -y python \
&& sudo curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py \
&& sudo python3 get-pip.py \
&& pip install synapseclient \
&& pip install aws-mfa \
&& pip install awsmfa \
&& git clone https://github.com/jgockley62/igraphNetworkExpansion.git \
&& git clone https://github.com/jgockley62/networktracing.git

# Install R Dependencies
RUN R -e "install.packages('bcv')" \
&& R -e "install.packages('BiocManager')" \
&& R -e "BiocManager::install('biomaRt')" \
&& R -e "BiocManager::install('annotate')" \
&& R -e "BiocManager::install('paxtoolsr')" \
&& R -e "BiocManager::install('clusterProfiler')" \
&& R -e "BiocManager::install('ComplexHeatmap')" \
&& R -e "install.packages('corpcor')" \
&& R -e "BiocManager::install('cqn')" \
&& R -e "install.packages('data.table')" \
&& R -e "install.packages('doParallel')" \
&& R -e "BiocManager::install('DOSE', ask=F)" \
&& R -e "install.packages('dplyr')" \
&& R -e "BiocManager::install('edgeR')" \
&& R -e "install.packages('foreach')" \
&& R -e "install.packages('ggplot2')" \
&& R -e "install.packages('ggrepel', repos = 'http://cran.us.r-project.org')" \
&& R -e "install.packages('glmnet', repos = 'http://cran.us.r-project.org')" \
&& R -e "BiocManager::install('graphics')" \
&& R -e "BiocManager::install('grDevices')" \
&& R -e "install.packages('httr')" \
&& R -e "install.packages('igraph')" \
&& R -e "install.packages('kableExtra', repos = 'http://cran.us.r-project.org')" \
&& R -e "install.packages('knitr')" \
&& R -e "BiocManager::install('limma')" \
&& R -e "install.packages('optmatch', repos = 'http://cran.us.r-project.org')" \
&& R -e "BiocManager::install('org.Hs.eg.db')" \
&& R -e "install.packages('parallel')" \
&& R -e "install.packages('pheatmap', repos = 'http://cran.us.r-project.org')" \
&& R -e "install.packages('plotrix', repos = 'http://cran.us.r-project.org')" \
&& R -e "install.packages('plyr')" \
&& R -e "install.packages('RColorBrewer')" \
&& R -e "BiocManager::install('RCy3')" \
&& R -e "install.packages('reshape')" \
&& R -e "install.packages('reticulate', repos = 'http://cran.us.r-project.org')" \
&& R -e "install.packages('Rglpk')" \
&& R -e "install.packages('RJSONIO')" \
&& R -e "install.packages('rmarkdown')" \
&& R -e "BiocManager::install('scater')" \
&& R -e "BiocManager::install('stats')" \
&& R -e "BiocManager::install('sva')" \
&& R -e "install.packages('tidyr')" \
&& R -e "BiocManager::install('utils')" \
&& R -e "install.packages('vbsr', repos='http://cran.us.r-project.org')" \
&& R -e "BiocManager::install('WGCNA')" \
&& R -e "devtools::install_github('brian-bot/githubr')" \
&& R -e "devtools::install_github('blogsdon/spike/spike')" \
&& R -e "devtools::install_github('rlbarter/superheat')" \
&& R -e "devtools::install_github('jhrcook/HotNetvieweR')" \
&& R -e "devtools::install_github('jgockley62/igraphNetworkExpansion')" 

#RUN sudo apt-get update
#RUN sudo apt-add-repository -y ppa:marutter/rrutter
#RUN sudo apt-get update

RUN sudo R CMD javareconf

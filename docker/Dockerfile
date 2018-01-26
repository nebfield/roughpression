# R3.4.3_Bioc3.6
FROM bioconductor/release_base2@sha256:12bcb76ef632ea1fc89843bd7e970409895b892e6ab8b6f92d0c6f3150a53444
LABEL maintainer "wingfield-b@ulster.ac.uk"

RUN apt-get update && apt-get -y install fastx-toolkit parallel default-jdk ca-certificates-java 

RUN apt-get -y install python3 python3-dev python3-pip \
  && pip3 install cutadapt

RUN wget http://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz -P /tmp/ \
  && tar -xf /tmp/edirect.tar.gz -C /tmp/ \
  && ./tmp/edirect/setup.sh 

RUN wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.2-1/sratoolkit.2.8.2-1-ubuntu64.tar.gz -P /tmp/ \
  && tar -xf /tmp/sratoolkit.2.8.2-1-ubuntu64.tar.gz -C /tmp/

RUN R -e 'source("https://bioconductor.org/biocLite.R"); biocLite(c("dada2","phyloseq"))'

# sed fixes weird SSL cert bug 
RUN cd /tmp/ \
 && git clone https://github.com/wgromniak/mahout-extensions.git \
 && cd mahout-extensions/ \
 && sed -i -e 's/https/http/' gradle/wrapper/gradle-wrapper.properties \
 && ./gradlew jar4StandaloneReducts

ENV PATH $PATH:/tmp/edirect/
ENV PATH $PATH:/tmp/sratoolkit.2.8.2-1-ubuntu64/bin/

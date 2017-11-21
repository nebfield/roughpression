# Digest for reproducibility  
FROM bioconductor/release_base2@sha256:140fc150854b7c8c095f018f2b184e49b7771775714dce37a6debd8a1e9fa924

LABEL maintainer "wingfield-b@ulster.ac.uk"


RUN apt-get update && apt-get -y install fastx-toolkit parallel

RUN wget http://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz \
  && tar -xf edirect.tar.gz \
  && ./edirect/setup.sh 

ENV PATH $PATH:/edirect/

RUN wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.2-1/sratoolkit.2.8.2-1-ubuntu64.tar.gz \
  && tar -xf sratoolkit.2.8.2-1-ubuntu64.tar.gz 

ENV PATH $PATH:/sratoolkit.2.8.2-1-ubuntu64/bin/


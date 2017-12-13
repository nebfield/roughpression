#!/usr/bin/env nextflow

process get_jiang_dataset {
  storeDir "$baseDir/cache/download"

  output:
  file "*.fastq" into jiang_raw_fastq

  """
  esearch -db sra -query "srp042956[accn]" |\
    efetch --format runinfo |\
    cut -d ',' -f 1 |\
    grep SRR |\
    parallel fastq-dump --outdir \$PWD/ --skip-technical \
    --readids --read-filter pass --dumpbase --split-files --clip
  """
}

jiang_bc = file("$baseDir/bootstrap/jiang-bc.txt")
process prep_jiang_dataset {
  storeDir "$baseDir/cache/demuxed"

  input:
  file jiang_raw_fastq
  file bc from jiang_bc

  output:
  file "*.fastq" into jiang_split_fastq mode flatten

  """
  cat $jiang_raw_fastq | fastx_barcode_splitter.pl -bcfile $bc -bol --prefix jiang_ --suffix .fastq
  find . -name '*.fastq' -size 0 -print0 | xargs -0 rm # remove zero byte file (stats)
  """
}

process trim_jiang {
  cache true

  input:
  file fastq_file from jiang_split_fastq

  output:
  file "*.fastq.gz" into trimmed_jiang 

  """
  cutadapt -a barcode_linkedprimer=NNNNNNNNNNTTACCGCGGCTGCTGGCAC...CTGAGCCAGGATCAAACTCT ${fastq_file} > ${fastq_file.baseName}.fastq.gz
  """
}

process dada2_jiang {
  cache true

  input:
  file fastas from trimmed_jiang.collect()

  output:
  file 'dada2.Rdata' into gut_dada2

  """
  dada2_jiang.R $fastas
  """
}

process phyloseq_jiang {
  cache true
  echo false

  input:
  file gut_dada2
  
  output:
  file 'gut.rds' into gut_ps

  """
  phyloseq_jiang.R $gut_dada2
  """
}

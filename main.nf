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
  storeDir "$baseDir/cache/dada2_jiang"

  input:
  file fastas from trimmed_jiang.collect()

  output:
  file 'dada2.Rdata' into gut_dada2

  """
  dada2_jiang.R $fastas
  """
}

process phyloseq_jiang {
  storeDir "$baseDir/cache/phyloseq_jiang"

  input:
  file gut_dada2
  
  output:
  file 'gut.rds' into gut_ps

  """
  phyloseq_jiang.R $gut_dada2
  """
}

// TODO: change from hardcoded 
oral_dir = Channel.fromPath('/home/ben/data/oral_microbiome/combined-cohort/*.fastq.gz')

process dada2_oral {
  storeDir "$baseDir/cache/dada2_oral"
  stageInMode 'copy'
  
  input:
  file fastq from oral_dir.collect()

  output:
  file 'dada2.RData' into oral_dada2
  
  """
  dada2_oral.R 
  """
}

samp_data = file("$baseDir/bootstrap/sample_data.tsv")

process phyloseq_oral {
  storeDir "$baseDir/cache/phyloseq_oral"
  
  input:
  file oral_dada2
  file samp_data

  output:
  file 'oral.rds' into oral_ps

  """
  phyloseq_oral.R $oral_dada2 $samp_data
  """
}

process reduct_gut {
  input:
  file gut_ps
  
  output:
  file 'short_names.tsv' into gut_shortnames
  file 'gut.arff' into gut_arff
  
  """
  rr-gut.R $gut_ps
  # java -Xmx32g -jar /tmp/mahout-extensions/build/libs/mahout-extensions-standalone-reducts.jar -i gut.csv -numSub 5000 -subCard 40 -seed 0451 > gut_feature_ranks.txt
  # rr-select.R $gut_ps gut_feature_ranks.txt
  """
}

process rs_gut {
  publishDir "$baseDir/results/gut"
  
  input:
  file gut_arff
  
  output:
  file 'rules.txt' into gut_rules
  file 'rule-support.csv' into gut_rule_support 
  
  """
  java -cp /tmp/mbrs.jar uk.ac.ulster.rs.Microbiome $gut_arff
  """
}

process reduct_oral {
  input:
  file oral_ps

  """
  rr-oral.R $oral_ps
  # java -Xmx32g -jar /tmp/mahout-extensions/build/libs/mahout-extensions-standalone-reducts.jar -i oral.csv -numSub 5000 -subCard 40 -seed 0451 > oral_feature_ranks.txt
  """
}

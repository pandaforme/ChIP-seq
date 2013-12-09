#!/bin/bash

source $1

function getFasta
{
  path=$1
  name=$2

  if [ "${IS_INTERSECT_GTF}" = "true" ]; then
    intersectBed -wo -a ${path}/${name}.location -b ${ANNOTATION_FILE} | awk -v OFS="\t" '{print $1,$2,$3,"*","*",$10}' | uniq > ${path}/${name}.bed
    
    if [ ! -s ${path}/${name}.bed ]; then
      echo "  Error：Fetch sequences from the sense strand fail."
      exit -1
    fi
  fi  
  
  fastaFromBed -fi reference/ref.fasta -bed ${path}/${name}.location -fo ${path}/${name}.fa
}

#generate significant peak
function generateSignificantPeak
{
  echo "Info：Generating significant peak ..."
  
  grep -ve '^#' -ve '^$' -ve '^chr' MACS/${NAME}/${NAME}_peaks.xls | awk -v FDR_RATE=${FDR_RATE} '{if ($9 <= FDR_RATE) print}' > significantPeak/${NAME}_sig_peaks.xls
  
  if [ -e "MACS/${NAME}/${NAME}_negative_peaks.xls" ]; then
    grep -ve '^#' -ve '^$' -ve '^chr' MACS/${NAME}/${NAME}_negative_peaks.xls > significantPeak/${NAME}_negative_peaks.xls
  fi
}

#generate the input data for MEME
function generateMEMEData
{
  echo "Info：Generating the input data for MEME ..."
  
  sort -k8,8 -n -r significantPeak/${NAME}_sig_peaks.xls | head -${BEST_REGIONS} | awk '{summit = $2 - 1 + $5; start = summit - 51; if (start < 0) {start = 0;} print $1 "\t" start "\t" summit + 50}' > MEME/input/${NAME}_bestPeaks.location
  
  getFasta "MEME/input" "${NAME}_bestPeaks"

  if [ -e "significantPeak/${NAME}_negative_peaks.xls" ]; then
    awk '{summit = $2 - 1 + $5; start = summit - 51; if (start < 0) {start = 0;} print $1 "\t" start "\t" summit + 50}' significantPeak/${NAME}_negative_peaks.xls > MEME/input/${NAME}_negative_peaks.location
    
    getFasta "MEME/input" "${NAME}_negative_peaks"
    
    psp-gen -pos MEME/input/${NAME}_bestPeaks.fa -neg MEME/input/${NAME}_negative_peaks.fa -alpha dna -revcomp > MEME/input/${NAME}_priors.psp
  fi
}

#generate the input data for CentriMo
function generateCentriMoData
{
  echo "Info：Generating the input data for CentriMo ..."
  
  awk '{summit = $2 - 1 + $5; start = summit - 151; if (start < 0) {start = 0;} print $1 "\t" start "\t" summit + 150}' significantPeak/${NAME}_sig_peaks.xls > CentriMo/input/${NAME}_sig_peaks_summit.location
  
  getFasta "CentriMo/input" ${NAME}_sig_peaks_summit

  fasta-get-markov -nostatus -m 1 < CentriMo/input/${NAME}_sig_peaks_summit.fa 1>CentriMo/input/${NAME}_sig_peaks_summit.fa.bg

  if [ -e "significantPeak/${NAME}_negative_peaks.xls" ]; then
    awk '{summit = $2 - 1 + $5; start = summit - 151; if (start < 0) {start = 0;} print $1 "\t" start "\t" summit + 150}' significantPeak/${NAME}_negative_peaks.xls > CentriMo/input/${NAME}_negative_peaks.location
    
    getFasta "CentriMo/input" "${NAME}_negative_peaks"

    fasta-get-markov -nostatus -m 1 < CentriMo/input/${NAME}_negative_peaks.fa 1>CentriMo/input/${NAME}_negative_peaks.fa.bg
  fi   
}

#generate the input data for PeakAnnotator
function generatePeakAnnotatorData
{
  echo "Info：Generating the input data for PeakAnnotator ..."
  
  awk '{print $1 "\t" $2 "\t" $3}' significantPeak/${NAME}_sig_peaks.xls > PeakAnnotator/input/${NAME}_sig_peaks_PAinput.txt
}

#run MACS1.4
function runMACS14
{
  echo "Info：Running MACS 1.4 ..."
  
  mkdir -p MACS/${NAME}
  
  cd MACS/${NAME}
  if [ "${TREAT}" != "" ] &&  [ "${CONTROL}" != "" ]; then
    macs14 -t ${TREAT} -c ${CONTROL} -g ${GENOME_SIZE} -n ${NAME} -w 1>../../log/MACS_${NAME}.log 2>&1
  elif [ "${TREAT}" != "" ] &&  [ "${CONTROL}" = "" ]; then
    macs14 -t ${TREAT} -g ${GENOME_SIZE} -n ${NAME} -w 1> ../../log/MACS_${NAME}.log 2>&1
  fi

  if [ $? -ne 0 ]; then
    echo "  Error：Please look up log/MACS_${NAME}.log"
    exit -1
  fi

  Rscript ${NAME}_model.r 1>>../../log/MACS_${NAME}.log 2>&1
  if [ $? -ne 0 ]; then
    echo "  Error：Please look up log/MACS_${NAME}.log."
    exit -1
  fi

  cd ../../
}

#run MEME
function runMEME
{
  echo "Info：Running MEME ..."

  if [ -e "MEME/input/${NAME}_priors.psp" ]; then
    meme  MEME/input/${NAME}_bestPeaks.fa -p 24 -revcomp -dna -nmotifs ${NMOTIFS} -mod zoops -maxsize ${MAXSIZE} -psp MEME/input/${NAME}_priors.psp -oc MEME/output/${NAME}_bestPeaks_meme 1>log/MEME_${NAME}.log 2>&1 
  else  
    meme  MEME/input/${NAME}_bestPeaks.fa -p 24 -revcomp -dna -nmotifs ${NMOTIFS} -mod ${MOD} -maxsize ${MAXSIZE} -oc MEME/output/${NAME}_bestPeaks_meme 1>log/MEME_${NAME}.log 2>&1
  fi  

  if [ $? -ne 0 ]; then
    echo "  Error：Please look up log/MEME_${NAME}.log."
    exit -1
  fi
}

#run CentriMo
function runCentriMo
{
  echo "Info：Running CentriMo ..."

  centrimo --oc CentriMo/output/${NAME}_peaks_motif_centrimo --bgfile CentriMo/input/${NAME}_sig_peaks_summit.fa.bg CentriMo/input/${NAME}_sig_peaks_summit.fa MEME/output/${NAME}_bestPeaks_meme/meme.txt 1>log/CentriMo_${NAME}.log 2>&1

  if [ $? -ne 0 ]; then
    echo "  Error：Please look up log/CentriMo_${NAME}.log."
    exit -1
  fi

  if [ -e "CentriMo/input/${NAME}_negative_peaks.fa" ]; then    
    centrimo --oc CentriMo/output/${NAME}_negative_peaks_motif_centrimo --bgfile CentriMo/input/${NAME}_negative_peaks.fa.bg CentriMo/input/${NAME}_negative_peaks.fa MEME/output/${NAME}_bestPeaks_meme/meme.txt 1>>log/CentriMo_${NAME}.log 2>&1
  fi

  if [ $? -ne 0 ]; then
    echo "  Error：Please look up log/CentriMo_${NAME}.log."
    exit -1
  fi
}

#run PeakAnnotator
function runPeakAnnotator
{
  echo "Info：Running PeakAnnotator ..."

  java -jar ${PEAKANNOTATOR_FILE} -u ndg -p PeakAnnotator/input/${NAME}_sig_peaks_PAinput.txt -a ${ANNOTATION_FILE} -g all -o PeakAnnotator/output 1>log/PeakAnnotator_${NAME}.log 2>&1

  if [ $? -ne 0 ]; then
    echo "  Error：Please look up log/PeakAnnotator_${NAME}.log."
    exit -1
  fi

  awk '{if(NR!=1 && length($6) != 0){print $6}}' PeakAnnotator/output/${NAME}_sig_peaks_PAinput.ndg.txt | uniq > PeakAnnotator/output/${NAME}_FW_uniq_gene.txt 
  awk '{if(NR!=1 && length($9) != 0){print $9}}' PeakAnnotator/output/${NAME}_sig_peaks_PAinput.ndg.txt | uniq > PeakAnnotator/output/${NAME}_REV_uniq_gene.txt
}

#remove files
function removeFiles
{
  echo "Info：Cleaning all files ..."

  rm -rf significantPeak/*
  rm -rf MACS/*
  rm -rf MEME/input/*
  rm -rf MEME/output/*
  rm -rf CentriMo/input/*
  rm -rf CentriMo/output/*
  rm -rf PeakAnnotator/input/*
  rm -rf PeakAnnotator/output/*
  rm -rf log/*  
}

#make directory
function makeDirectory
{
  echo "Info：Making directories ..."

  mkdir -p log
  mkdir -p significantPeak
  mkdir -p MACS
  mkdir -p MEME/input
  mkdir -p MEME/output
  mkdir -p CentriMo/input
  mkdir -p CentriMo/output
  mkdir -p PeakAnnotator/input
  mkdir -p PeakAnnotator/output
}

#check
function check
{
  echo "Info：Checking annotation file ..."
  
  if [ ! -e ${ANNOTATION_FILE} ]; then
    echo "  Error：Please set a correct path of a annotation file."
    exit -1;  
  fi

  echo "Info：Checking reference files ..."
  if [ ! -e ${REFERENCE_FILE} ]; then
    echo "  Error：Please set a correct path of a reference file ."
    exit -1;  
  fi
  
  if [ ! -e ${REFERENCE_INDEX_FILE} ]; then  
    echo "  Error：Please set a correct path of a reference index file."
    exit -1;  
  fi 

  if [ ! -e ${PEAKANNOTATOR_FILE} ]; then  
    echo "  Error：Please set a correct path of a PeakAnnotator jar file."
    exit -1;  
  fi   
}

if [ "${IS_MAKE_DIRECTORY}" = "true" ]; then  
  makeDirectory
fi

if [ "${IS_REMOVE}" = "true" ]; then  
  removeFiles
fi

check

if [ "${IS_RUN_MACS}" = "true" ]; then
  runMACS14
fi

if [ "${IS_GENERATE_INPUT_FILES}" = "true" ]; then
  generateSignificantPeak

  generateMEMEData

  generateCentriMoData

  generatePeakAnnotatorData
fi

if [ "${IS_RUN_MEME}" = "true" ]; then
  runMEME
fi

if [ "${IS_RUN_CENTRIMO}" = "true" ]; then
  runCentriMo
fi

if [ "${IS_RUN_PEAKANNOTATOR}" = "true" ]; then
  runPeakAnnotator
fi
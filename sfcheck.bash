#!/bin/bash -f

if [[ -e $1 ]] 
then
  a="OK"
else
  echo "Input file $1 not found"
  echo " "
  echo "Please provide a name of a sequence file to be checked as an argument "
  echo "  "
  exit 1
fi

export LC_ALL=C
export PATH=/appl/molbio/emboss/6.5.7_gnu/bin:${PATH}


fformat=$(file $1 |awk '{print $2}')

#bam tiedostyt pitää lukea kokonaan
if [[ $fformat == "gzip" ]] 
then
  format=$(infoseq -nowarning -noerror -nodie -only -usa -nohead -filter $1 | awk -F "::" '{print $1}' | uniq | head -1)
else
  format=$(head -800 $1 | infoseq -nowarning -noerror -nodie -only -usa -nohead -filter | awk -F "::" '{print $1}' | uniq | head -1)
  if [[ $format == "" ]] 
  then
      format=$(head -99999 $1 | infoseq -nowarning -noerror -nodie -only -usa -nohead -filter | awk -F "::" '{print $1}' | uniq | head -1)
  fi
fi

  
if [[ $format == "fastq" ]] 
then
   q_count=$(head -800 $1 | grep "^+" -A 1 | grep -v "^+" | grep -v '\-\-' | wc -c)
   Sa_Il18=$(head -800 $1 | grep "^+" -A 1 | grep -v "^+" | grep -v '\-\-' | tr -d '!"#$%&()*+,-./0123456789:' | wc -c )
   So_1l13_Il15=$(head -800 $1 | grep "^+" -A 1 | grep -v "^+" | grep -v '\-\-' | tr -d 'KLMNOPQRSTUVWXYZ[\]^_abcdefgh' |wc -c )

   #echo Sanger or Illumina 1.8
   if [[ $q_count != $Sa_Il18 ]] 
   then
     Il_18=$(head -800 $1 | grep "^+" -A 1 | grep -v "^+" | grep -v '\-\-' | tr -d '!#$%&()*+,-./0123456789:;<=>?@ABCDEFGHI' | wc -c )
     if [[ $Il_18 -gt 0 ]] 
     then 
        #echo "Fastq format is: Illumina 1.8+"
        format=("fastq Phred+33 illumina18+")
     else
       #echo "Fastq format is: Sanger or Illumina 1.8+"
       format=("fastq Phred+33")
     fi
  fi

  #echo Solexa or Illumina 13+ or Illumina 15+ 
  if [[ $q_count != $So_1l13_Il15 ]]
  then
    So=$(head -800 $1 | grep "^+" -A 1 | grep -v "^+" | grep -v '\-\-' | tr -d ';<=>?' | wc -c ) 
    if [[ $q_count -ne $So ]] 
    then
      #echo "Fastq format is: Solexa" 
      format=("fastq Solexa+64")
    else
      So_Il13=$(head -800 $1 | grep "^+" -A 1 | grep -v "^+" | grep -v '\-\-' | tr -d '@A' | wc -c )
      if [[ $q_count -ne $So_Il13 ]] 
         then
         #echo "Fastq format is: Solexa or Illumina 13+" 
         format=("fastq Phred+64")
      else
         #echo "Fastq format is: Solexa or Illumina 13+ or Illumina 15+ " 
        format=("fastq Phred+64")
      fi
    fi 
  fi
fi

if [[ "$format" == "" ]] 
then
  echo "Not an EMBOSS compatible sequence file"
else
  echo $format
fi

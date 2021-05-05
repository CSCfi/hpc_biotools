#!/bin/bash -f
# A script to fetch sequences from a sequence file or database
# 21.1. 2011 KM

export LC_ALL=C
#setenv PLPLOT_LIB /v/linux26_x86_64/appl/molbio/emboss/6.5.7/share/EMBOSS
export PATH=/appl/molbio/emboss/6.5.7_gnu/bin:${PATH}
#setenv LD_LIBRARY_PATH /v/linux26_x86_64/appl/molbio/emboss/libharu/2.0.8/lib
osformat=("fasta")

if [[ "$1" == "" ]]
 then
   echo "Syntax:"
   echo "infoseq_summary sequence_file"
   exit 1
fi

if [[ -e $1 ]]
then
   infile=$1
   printf "%s\t%s\n" " File name:                           " $1
   size=$(du -sh "$1" | awk '{print $1}' )
   printf "%s\t%s\n" " File size:                          " $size
   sformat=$(sfcheck.bash "$1")
   printf "%s\t%s\n" " File format:                         " $sformat 
   if [[ "$sformat" == "Not an EMBOSS compatible sequence file" ]] 
   then
     exit 0
   fi
else 
#   seqret "$1" $TMPDIR/infoseq_summary_tmp_$$
#   set infile=$TMPDIR/infoseq_summary_tmp_$$
    echo ' USA definition:                      '"\t""$1" 
fi

infoseq -nowarning -nocolumn -delimiter "::"  -nohead "$1" -only -usa -name -type -length -filter | awk -F "::" 'BEGIN { s = 10000000000} { a = a +$NF} \
{ if ( $NF > l) l = $NF } { if ( $NF == l) ln = $3 }  { if ( $NF < s) s = $NF} { if ( $NF == s) sn = $3} {ka = a / NR} \
END { if ( $4 == "N")  print " Sequence type:                       \tNucleotide"} \
END { if ( $4 == "P")  print " Sequence type:                       \tProtein"} \
END { print " Number of sequences:                  \t" NR } \
END { print " Longest (or one of equally long):     \t" ln "\t" l  } \
END { print " Shortest (or one of equally short): \t" sn "\t"s } \
END { print " Average length:                      \t" ka } \
END {print  " Total amount of nucleotides/residues:\t" a } \
END { if ( NF > 5) print " Note: Sequence names seem to contain douple-douple point (::) characters!"}'



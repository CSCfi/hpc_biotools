#!/bin/bash
# pb blast- ohjelma puhti-palvelimelle
# K.M. 16.08. 2019
#
#
pblogfile=("/tmp/pb_blast_log")
gid=$(pwd | awk -F "/" '{print $3}' | awk -F "_" '{print $2}')
pname=$(pwd | awk -F "/" '{print $3}')


#echo "GID: $gid"
pbtmproot=("/scratch/$pname")
#pbtmproot=("/scratch/project_$gid")
start_time=$(date)
user=$(whoami)
#max_jobs=$(sacctmgr show -p assoc format=MaxSubmit where user=$user account=$pname partition=small | tail -1 | rev | cut -c2- |rev)
max_jobs=300
(( n_split = max_jobs / 2 - 1 ))

if [[ $# -le 1 ]]
then
    echo "--------------------------------------------------------------------"
    echo "pb program is designed to run effectively BLAST searches"
    echo "for large amounts of query sequences."
    echo "pb accepts blast search commands and all their arguments."
    echo " "
    echo "Usage: pb blast_command e.g.:"
    echo "pb blastn -db embl -query myseqs.fasta -out results"
    echo " "
    echo "pb specific arguments:"
    echo "  -dbprot [seqfile]  Use the given protein sequce file as search database "
    echo "  -ensembl_dna [species_name]  Retrieve the genomic sequence of the given species from ensembl database and use it as the database"
    echo "  -ensembl_cdna [species_name]  Retrieve the coding sequences of the given species from ensembl database and use it as the database"
    echo "  -ensembl_prot [species_name]  Retrieve the peptide sequences of the given species from ensembl database and use it as the database"   
    echo "  -dbnuc [seqfile]  Use the given nucleotide sequce file as search database "
    echo "  -batch submit the job and write out instuctions how to follow the job and collect the results"
    echo "  -no_slurm Run the job in local machine and not through the batch job system"
    echo "  -mem_limit [megabytes] Reserve at least the specified amount of memory for each subjob"
    echo "  -subjob_run_time [hours] Maximum run time for a subjob (default 48 h, max 72 h) "
    echo "  -post_process [post_process_script] defines name of a post processing script that is applied to the resuts of a sub-job"
    echo "--------------------------------------------------------------------"
    source /appl/soft/bio/blast/setup_blast.bash
    exit
fi

# Tarkista onko ohjelma sallittu Tarkista onko ohjelma sallittu
case "$1" in
       'blastp')
            blast_type=$1
        ;;
       'blastn')
            blast_type=$1
        ;;
       'blastx')
            blast_type=$1
        ;;
        'rpsblast')
            blast_type="blastp"
        ;;
       'psiblast')
            blast_type="blastp"
        ;;
        'rpsblastn')
            blast_type="blastx"
        ;;
       'tblastn')
            blast_type=$1
        ;;
       'tblastx')
            blast_type="blastn"
        ;;
       'deltablast')
            blast_type="blastp"
        ;;
        *)
            echo "   "
            echo "$argv[1] can not be used with pb"
            echo "pb can only be used to run commands:"
	    echo "  blastn"
	    echo "  blastp"
	    echo "  blastx"
            echo "  deltablast"
            echo "  psiblast"
            echo "  rpsblast"
            echo "  rpsblastn"
            echo "  tblastn"
            echo "  tblastx"
	    echo "  "
            echo "Usage: pb blast_command, e.g.:"
            echo "pb blastn -db embl -query myseqs.fasta -out results"
            exit 1
        ;;
esac



#module load emboss
command=" "


#Tarkista input, output ja tietokanta
inputflag=0
outputflag=0
databaseflag=0
ensemblflag=0
dbbuildflag=0
db_type="none"
pb_output_type=0
mode="interactive"
debug=0
parse_seqids="-parse_seqids"  
blastdb_orig=$BLASTDB
mem_limit=2048
run_time=48
post_process=0
use_gi_list=0
use_taxid_list=0

while [[ $# -ge 1 ]]
do
  case "$1" in
             '-query')
             # query file
                  inputseq="$2"
                  if [[ ! -e $inputseq ]] 
                  then
                     echo ""
                     echo " Query sequence file: $inputseq does not exist"
                     echo ""
                     echo "-----------------------------------------------------------"
                     exit 1
                  fi
                  inputflag=1
		  inp_type=$(head -40 $inputseq | infoseq -nohead -only -type -filter | tail -1 | tr -d " " ) 
                  shift
                  shift
                ;;
               '-out')                   
               # result file
                  blast_output="$2"
                  outputflag=1
                  if [[ -e $blast_output ]] 
                  then                        
		      echo " "
		      echo "Output file: $blast_output already exists."
		      echo ""
		      exit 1                      
    		  fi
                  shift
                  shift 
                ;; 
                '-db')
                 # database
                  blast_database="$2"
                  blast_database_name=($blast_database)
                  databaseflag=1
                  if [[ -e $BLASTDB/$blast_database.nal ]] 
                  then 
		    db_type="nuc"
		  fi
		  if [[ -e $BLASTDB/$blast_database.nsq ]]
                  then 
		    db_type="nuc"
		  fi
		  if [[ -e $BLASTDB/$blast_database.pal ]] 
                  then 
		    db_type="prot"
		  fi
		  if [[ -e $BLASTDB/$blast_database.psq ]] 
                  then 
		    db_type="prot"
		  fi
                  if [[ -e $blast_database.psq ]] 
                  then 
		    db_type="prot"
		  fi
                   if [[ -e $blast_database.nsq ]] 
                  then 
		    db_type="nuc"
		  fi
                  shift
                  shift
                ;;  
                '-dbprot')
                 # database
                  blast_database="$2"
                  dbbuildflag=1
                  databaseflag=1
                  db_type="prot"
                  db_letter=$(head $blast_database | grep -v ">" | infoseq -nohead -only -type -filter -sformat asis | tail -1 | tr -d " \t")
                  if [[ $db_letter != "P" ]]
                  then
		      echo "ERROR: Protein database expected. "
                      echo "Your database sequence does not look like a protein sequence." ; echo ""
                      echo "First rows of your database sequence file $blast_database looked like:" ; echo ""
                      head $blast_database 
                      exit -1
                  fi 
                  shift 
                  shift 
                ;;
		'-dbnuc')
                # database
                  blast_database="$2"
                  dbbuildflag=2
                  databaseflag=1
		  db_type=nuc
                  db_letter=$(head $blast_database | grep -v ">" | infoseq -nohead -only -type -filter -sformat asis | tail -1 | tr -d " ")
                  if [[ $db_letter != "N" ]] 
                  then
                      echo "ERROR: Nucleotide database expected "
                      echo "Your database sequence does not look like a nucleotide sequence."; echo ""
                      echo "First rows of your database sequence file $blast_database looked like:"; echo ""
                      head $blast_database 
                      exit 1
                  fi
                  shift
                  shift 
                ;;
                '-ensembl_dna')
                # database
                  blast_database=$(ensemblfetch -type dna $2 | tee | tail -1 )
                  if [[ -e $blast_database ]] 
                  then
                     echo "The genome data file $blast_database succesfully retrieved"
                  else
                     echo "The Ensembl data was not foud"
		     echo "Please check that species name $2 exists in the Ensemble or Ensembl genomes database"
                     exit 1
                  fi
                  dbbuildflag=2
                  databaseflag=1
                  ensemblflag=1
		  db_type="nuc"
                  db_letter=$(head $blast_database | grep -v ">" | infoseq -nohead -only -type -filter | tail -1 | tr -d " ")
                  if [[ $db_letter != "N" ]] 
                  then
                      echo "ERROR: Nucleotide database expected "
                      echo "Your database sequence does not look like a nucleotide sequence."; echo ""
                      echo "First rows of your database sequence file $blast_database looked like:"; echo ""
                      head $blast_database 
                      exit 1
                  fi 
                  shift 
                  shift
                ;; 
                '-ensembl_cdna')
                # database
                   blast_database=$(ensemblfetch -type cdna $2 | tail -1 )
                   if [[ -e $blast_database ]] 
                   then
                     echo "The genome data file $blast_database succesfully retrieved"
                   else
                     echo "The Ensembl data was not foud"
		     echo "Please check that species name $2 exists in the Ensemble or Ensembl genomes database"
                     exit 1
                   fi
                   dbbuildflag=2
                   databaseflag=1
                   ensemblflag=1
		   db_type="nuc"
                   db_letter=$(head $blast_database | grep -v ">" | infoseq -nohead -only -type -filter | tail -1 | tr -d " ")
                   if [[ $db_letter != "N" ]] 
                   then
                      echo "ERROR: Nucleotide database expected "
                      echo "Your database sequence does not look like a nucleotide sequence."; echo ""
                      echo "First rows of your database sequence file $blast_database looked like:"; echo ""
                      head $blast_database 
                      exit 1
                   fi 
                   shift 
                   shift
                 ;; 
   
                '-ensembl_prot')
                 # database
                   blast_database=$(ensemblfetch -type pep $2 | tail -1 )
                   if [[ -e $blast_database ]] 
                   then
                    echo "The protein data file $blast_database succesfully retrieved"
                   else
                     echo "The Ensembl data was not foud"
		     echo "Please check that species name $2 exists in the Ensemble or Ensembl genomes database"
                     exit 1
                   fi
                   dbbuildflag=1
                   databaseflag=1
                   ensemblflag=1
                   db_type="prot"
                   db_letter=$(head $blast_database | grep -v ">" | infoseq -nohead -only -type -filter | tail -1 |  tr -d " " )
                   if [[ $db_letter != "P" ]] 
                   then
		      echo "ERROR: Protein database expected. "
                      echo "Your database sequence does not look like a protein sequence." ; echo ""
                      echo "First rows of your database sequence file $blast_database looked like:" ; echo ""
                      head $blast_database 
                      exit
                   fi 
                   shift 
                   shift
                ;;  
		'-outfmt')
                 f_len=$(echo $2 | wc -c)
                 #If the argument is longer than 3 then we assume that this is format 6,7 or 10 with definitions
                 if [[ $f_len -le 4 ]]
                 then
                   if [[ $2 -eq 19 || $2 -eq 20 ]] 
                   then
                     fromatstring=('"6 qseqid sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore"')
                     command="$command $1 $fromatstring"
                     pb_output_type="$2"
                   fi
                   if [[ $2 -eq 21 ]] 
                   then
                     fromatstring=('"6 sacc sstart send sseq"')
                     command="$command $1 $fromatstring "
                     pb_output_type=($2)                    
                   fi
                   if [[ $2 -eq 22 ]] 
                   then
                      fromatstring="10"
                      command="$command $1 $fromatstring "
                      pb_output_type=($2)
                    fi  
                    if [[ $2 -eq 23 ]] 
                    then
                      fromatstring="6"
                      command="$command $1 $fromatstring  "
                      pb_output_type=($2)
                    fi 
                    if [[ $2 -eq 17 ]] 
                    then
                       echo " "
                       echo "pb blastn does not support output format 17 (SAM alignment)"
                       echo "If you wish to use this output format, you should execute your job as a plain blastn job"
                    #   exit 1
                    fi 
                    if [[ $2 -lt 19 ]] 
                    then
                       command="$command $1 \"$2\"" 
                       pb_output_type=($2) 
                       echo $command
                    fi 
                 else
                     #parsing format types 6,7,10                    
                     formatarray=($2)
                     pb_output_type=${formatarray[0]}
                     format_n_items=${#formatarray[*]}
                     command="$command $1 \"$pb_output_type"
                     for ((item=1; item<$format_n_items; item++))
                     do
                       command="$command ${formatarray[$item]}"
                     done
                     command="$command\" "
                 fi  
                 shift
                 shift  
                ;;
                '-num_threads')
                    echo " "
                    echo "Selection: -num_threads $2  ignored."
                    echo "gb defines this value automatically."
                    echo " "
                    shift 
                    shift 
                ;;
                '-help')
		  $blast_type -help
                  echo " *** pb spesific options:"
                  echo "  -dbprot  Use the given protein sequce file as search database "
                  echo "  -dbnuc  Use the given nucleotide sequce file as search database "
                  echo "  -ensembl_dna [species_name]  Retrieve the genomic sequence of the given species from ensembl database and use it as the database"
                  echo "  -ensembl_cdna [species_name]  Retrieve the coding sequences of the given species from ensembl database and use it as the database"
                  echo "  -ensembl_prot [species_name]  Retrieve the peptide sequences of the given species from ensembl database and use it as the database"
                  echo "  -batch submit the job and write out instuctions how to monitor the job and collect the results"
                  echo "  -no_seqid_parsing  With this option the database indexing for a sequence set defined with -dbnuc or -dbprot is done without -parse_seqids option" 
                  echo "  -no_slurm Run the job in local machine and not through the batch job system"
                  echo "  -mem_limit [megabytes] Reserve at least the specified amount of memory for each subjob"
                  echo "  -subjob_run_time [hours] Define maximum run time for a subjob (default 48 h, max 72 h)"
                  echo "  -project [project_number] Define the billing project used. By default the project name is set based on the name of the current scartch directory"
                  echo "  -post_process [post_process_script] defines name of a post processing script that is applied to the resuts of a sub-job"
                  echo "                The post process scriupt is extpected to be written so that it reads the name of the blast result file"
                  echo "                to be processed as the first argument and then wtites the results to standard output"             
                  echo "Extra modes for -outfmt option"
                  echo "  -outfmt 19  Print out a list of uniq hit sequece identifiers" 
                  echo "  -outfmt 20 print the hit sequences in fasta format"
                  echo "  -outfmt 21  Print out the matching regions of hit sequences in fasta format" 
                  echo " "
                  source /appl/bio/blast/setup_blast.bash
                  exit
		;;
                '-nsplit')
                  n_split=$2
                  shift
                  shift
                ;;	 
                '-batch')
                  mode="batch"
                  shift 
                ;;
                '-no_slurm')
                  mode="no_slurm"
                  shift 
                ;; 
                '-debug')
                  debug=1
                  shift
                ;;
                '-mem_limit')
                  mem_limit=$2
                  shift
                  shift 
                ;;
                '-subjob_run_time')
                  run_time=$2
                  if [[ $run_time -gt 72 ]]
                  then
                    echo "subjob_run_time must be less than 72 h"
                    exit 1
                  fi                    
                  shift
                  shift
                ;;
	        '-post_process')
                # query file
                  pp_script="$2"
                  if [[ ! -e $pp_script ]] 
                  then
                     echo ""
                     echo " Post processing script: $pp_script does not exist"
                     echo ""
                     echo "-----------------------------------------------------------"
                     exit 1
                  fi
                  post_process=1 
                  shift
                  shift
                ;; 

                '-gilist')
                  gi_list="$2"
                  if [[ ! -e $gi_list ]] 
                  then
                     echo ""
                     echo " GI list file: $gi_list does not exist"
                     echo ""
                     echo "-----------------------------------------------------------"
                     exit 1
                  fi
                  use_gi_list=1
                  shift
                  shift
                ;;     

                '-taxidlist')
                  taxid_list="$2"
                  if [[ ! -e $taxid_list ]] 
                  then
                     echo ""
                     echo " Taxonomy ID list file: $taxid_list does not exist"
                     echo ""
                     echo "-----------------------------------------------------------"
                     exit 1
                  fi
                  use_taxid_list=1
                  shift
                  shift
                ;;     

                '-no_seqid_parsing')
                  parse_seqids=" "
                shift
                ;;   
                '-archive')
                 inputseq="$2"
                 blast_database="$2"
                  inputflag=1 
                  inp_type=$(file $inputseq | awk -F : '{print $2}')                              db_type=(archive)
                shift
                shift
                ;;
                '-project')
                pname=$2
                gid=$(echo $2 | awk -F "_" '{print $NF}')
                if [[ $(groups | grep -c $gid) -ne 1 ]]; then
                    echo ""
                    echo "ERROR: Project $2 not found."
                    exit 1
                fi                    
                shift
                shift
                ;; 
                *)
            	command="$command $1"
                shift                       # No more switches
                ;;
    esac
done

if [[ $inputflag -eq 0 ]]
then
   echo "   "
   echo "Inputfile not defined"
   echo "You must define the file containing input sequences with option: -query file_name"
   echo "  "
   exit 1
fi

if [[ $databaseflag -eq 0 || $db_type == "none" ]] 
then
   echo "   "
   echo "Database not defined"
   echo "You must define the query database with option:" 
   echo "          -db database_name"
   echo "or feed in fasta formatted sequence file to be used as database with options:" 
   echo "          -dbprot protein_sequence_file"
   echo "or "
   echo "          -dbnuc nucleotide_sequece_file"
   echo "or you can use options -ensembl_dna, -ensembl_cdna, -ensembl_prot to retrieve "
   echo "a target genome from the Ensembl databases"
   exit 1
fi

if [[ $blast_type == "blastn" ]] 
then
     if [[ $db_type == "prot" ]]
     then
	echo "    "
	echo "Wrong database type."
	echo "blastn can only be used for searches against nucleotide databases."
	echo "    "
        exit 1
     fi
     if [[ $inp_type == "P" ]]
     then
	echo "    "
	echo "Wrong query sequence type."
	echo "blastn can only be used for searches with nucleotide sequences."
	echo "    "
        exit 1
     fi
fi

if [[ $blast_type == "blastp"  ]] 
then
     if [[ $db_type == "nuc" ]]
     then
	echo "    "
	echo "Wrong database type."
	echo "blastp can only be used for searches against protein databases."
	echo "    "
        exit 1
     fi
     if [[ $inp_type == "N" ]] 
     then
	echo "    "
	echo "Wrong query sequence type."
	echo "blastp can only be used for searches with protein sequences."
	echo "    "
        exit 1
     fi
fi

if [[ $blast_type == "blastx" ]] 
then
     if [[ $db_type == "nuc" ]]
     then
	echo "    "
	echo "Wrong database type."
	echo "blastx can only be used for searches against protein databases."
	echo "    "
        exit 1
     fi
     if [[ $inp_type == "P" ]] 
     then
	echo "    "
	echo "Wrong query sequence type."
	echo "blastx can only be used for searches with nucleotide sequences."
	echo "    "
        exit 1
     fi
fi
if [[ $blast_type == "tblastn" ]] 
then
     if [[ $db_type == "prot" ]] 
     then
	echo "    "
	echo "Wrong database type."
	echo "tblastn can only be used for searches against nucleotide databases."
	echo "    "
        exit 1
     fi
     if [[ $inp_type == "N" ]] 
     then
	echo "    "
	echo "Wrong query sequence type."
	echo "tblastn can only be used for searches with protein sequences."
	echo "    "
        exit
     fi
fi

if [[ $blast_type == "tblastx" ]] 
then
     if [[ $db_type == "prot" ]]
     then
	echo "    "
	echo "Wrong database type."
	echo "tblastx can only be used for searches against nucleotide databases."
	echo "    "
        exit 1
     fi
     if [[ $inp_type == "P" ]] 
     then
	echo "    "
	echo "Wrong query sequence type."
	echo "tblastx can only be used for searches with nucleotide sequences."
	echo "    "
        exit 1
     fi
fi

if [[ $outputflag == 0 ]] 
then
     blast_output=(blast_results)
fi

#Create a temporary directory
if [[ $gid == "" ]]; then
  if [[ -z $SLURM_JOB_ACCOUNT ]];then
     echo "Could not resolve project name"
     echo "Please give project name with option:"
     echo ""
     echo " -project project_name "
     exit 1
  else
     pname=$SLURM_JOB_ACCOUNT
  fi
fi

pbtmproot=("/scratch/$pname")

mkdir $pbtmproot/pb_$$_tmpdir



#mkdir $pbtmproot/pb_$$_tmpdir

#Move and check the inputseq
cp $inputseq $pbtmproot/pb_$$_tmpdir/
inputname=$(ls  $pbtmproot/pb_$$_tmpdir/)
inputseq=($pbtmproot/pb_$$_tmpdir/$inputname)

#Move te gi-list if used
if [[ $use_gi_list -eq 1 ]]
then
  cp $gi_list $pbtmproot/pb_$$_tmpdir/gi_list_file
  command="$command -gilist $pbtmproot/pb_$$_tmpdir/gi_list_file"
fi

#Move te gi-list if used
if [[ $use_taxid_list -eq 1 ]]
then
  cp $taxid_list $pbtmproot/pb_$$_tmpdir/taxid_list_file
  command="$command -taxidlist $pbtmproot/pb_$$_tmpdir/taxid_list_file"
fi

#Move post processing script if used
if [[ $post_process -eq 1  ]]
then
  cp $pp_script $pbtmproot/pb_$$_tmpdir/pp_script
  chmod u+x $pbtmproot/pb_$$_tmpdir/pp_script
fi

#In Batch mode, create a temporaray version of database to avoid
#problems caused by the database update
if [[ $mode == "batch" ]]
then
  if [[ $dbbuildflag -lt 1 ]] 
  then
      if [[ $outputflag -ne 0 ]] 
      then  
        echo "  "
        echo "Making temporary copy of $blast_database."
        echo " "
      fi
      blastdb_orig=($BLASTDB)
      
      # check if the database refres to a .pal or .nal file
      if [[ -e  "$BLASTDB"/"$blast_database".pal ]] 
      then 
       cp  "$BLASTDB"/"$blast_database".pal $pbtmproot/pb_$$_tmpdir/
          for subdb in $(grep DBLIST "$BLASTDB"/"$blast_database".pal | sed s/"DBLIST "/""/ )
          do
            location=$(pwd)
            cd $BLASTDB     
            cp "$subdb".*  $pbtmproot/pb_$$_tmpdir/
            cd $location
          done
      else
         if [[ -e  "$BLASTDB"/"$blast_database".nal ]]
         then 
           cp "$BLASTDB"/"$blast_database".nal $pbtmproot/pb_$$_tmpdir/
           for subdb in $(grep DBLIST "$BLASTDB"/"$blast_database".nal | sed s/"DBLIST "/""/)
           do
             location=$(pwd)
             cd $BLASTDB     
             cp "$subdb".*  $pbtmproot/pb_$$_tmpdir/
             cd $location
           done  
           else  
           cp "$BLASTDB"/"$blast_database".*  $pbtmproot/pb_$$_tmpdir/
        fi
      fi
      export BLASTDB=$pbtmproot/pb_$$_tmpdir 
      #ls $BLASTDB
  fi
fi

#Build BLAST-index if needed

#First check and fix possible duplicate database entries
if [[ $dbbuildflag -gt 0 ]] 
then
  grep ">" $blast_database | cut -f1 -d " " | sort | uniq -c | awk '{if ($1 > 1) print $0}'  >  $pbtmproot/pb_$$_tmpdir/non_unique_db_seq_names
  nonuniq=$(cat  $pbtmproot/pb_$$_tmpdir/non_unique_db_seq_names | wc -l )
  if [[ $nonuniq -gt 0 ]] 
  then
    echo "Following sequence names existed in the database sequence file more than once"
    cat $pbtmproot/pb_$$_tmpdir/non_unique_db_seq_names | awk '{print $2 " ocurrences: " $1}' 
    echo "Fixing names"
    cp $blast_database $blast_database"_tmp_"$$
    for non_uniq_name in $(awk '{print $2}' $pbtmproot/pb_$$_tmpdir/non_unique_db_seq_names)
    do
        echo "fixing $non_uniq_name"
        awk 'BEGIN{n=1}{ if ( $1 == "'"$non_uniq_name"'" ) print $1"_"n" "$2" "$3" "$4" "$5" "$6} {if ( $1 == "'"$non_uniq_name"'" ) n = n +1}{ if ( $1 != "'"$non_uniq_name"'" ) print $0} ' $blast_database"_tmp_"$$ >>!  "$blast_database"_"$$"_fixing
        rm -f $blast_database"_tmp_"$$
        mv "$blast_database"_"$$"_fixing $blast_database"_tmp_"$$
     done
     blast_database=("$blast_database"_tmp_"$$")
  fi
fi


if [[ $dbbuildflag -gt 0 ]] 
then

  if [[ $dbbuildflag -eq 1 ]]
  then
    dbtype=(prot)
  elif [[ $dbbuildflag -eq 2 ]] 
  then
    dbtype=(nucl)
  fi

  if [[ $outputflag -ne 0 ]] 
  then  
    echo "  "
    echo "Building BLAST indexes for sequences in $blast_database."
    echo " "
  fi
   
  num_db_seq=$(grep -c ">" $blast_database)
 

  #salloc -n 1 -t 8:00:00
  #srun seqret $blast_database $WRKDIR/pb_$$_tmpdir/blast_database.ncbi.fasta 
  #seqret $blast_database $WRKDIR/pb_$$_tmpdir/blast_database.ncbi.fasta
  cp $blast_database $pbtmproot/pb_$$_tmpdir/blast_database.ncbi.fasta 

  makeblastdb -dbtype $dbtype -in $pbtmproot/pb_$$_tmpdir/blast_database.ncbi.fasta $parse_seqids -out $pbtmproot/pb_$$_tmpdir/pb_$$_blast_database.tmp -blastdb_version 5 > $pbtmproot/pb_$$_tmpdir/makeblastdb_log

  blastdb_num=$(grep "Adding sequences from" $pbtmproot/pb_$$_tmpdir/makeblastdb_log | awk '{print $6}')
  if [[ $num_db_seq -ne $blastdb_num ]] 
  then
    cat $pbtmproot/pb_$$_tmpdir/makeblastdb_log
    exit 1
  elif [[ $outputflag -ne 0 ]] 
  then
     echo "$num_db_seq sequences indexed from $blast_database"
  fi  

  blastdb_orig=($BLASTDB)
  export BLASTDB=$pbtmproot/pb_"$$"_tmpdir  
  blast_database_name=($blast_database)                    
  blast_database=(pb_"$$"_blast_database.tmp)
fi

# Running the job
if [[ $outputflag -ne 0 ]]
then  
   echo "  "
   echo "Running:"
   if [[ $blast_type == "blast_formatter" ]] 
   then 
      echo $command -archive $inputseq -out $blast_output
   else 
      echo "$command -query $inputseq -db $blast_database -out $blast_output" 
   fi
   echo " "
fi

#check the database size

if [[ -e $BLASTDB/$blast_database.pal ]]
then
    db_koko=(0)
    for osa in $(grep DBLIST $BLASTDB/$blast_database.pal | cut -c7-900 | sed s/"\/qfs2\/biodb\/blast\/"//g | tr -d "\"" )
    do 
      osa=$( echo $osa | awk -F "/" '{ print $NF}' )    
      lisa=$( ls -l ${BLASTDB}/$osa* | awk '{ a = (a + $5)}{print a}' | tail -1 )
      (( db_koko = db_koko + lisa ))
    done
elif [[ -e $BLASTDB/$blast_database.nal ]] 
then
    db_koko=(0)
    for osa in $(grep DBLIST $BLASTDB/$blast_database.nal | cut -c7-900 | sed s/"\/qfs2\/biodb\/blast\/"//g | tr -d "\"" ) 
    do
      osa=$(echo $osa | awk -F "/" '{ print $NF}' ) 
      lisa=$(ls -l $BLASTDB/$osa* | awk '{ a = (a + $5)}{print a}' | tail -1 )
      (( db_koko = db_koko + lisa ))
    done
else
    db_koko=$(ls -l $BLASTDB/$blast_database* | awk '{ a = (a + $5)}{print a}' | tail -1)
fi

#Syötedatan käsittely tehdään eri tavalla jos kyse on archive-tiedostosta
if [[ $blast_type != "blast_formatter" ]] 
then
   #hatakorjaus 7.7. 2015 KM.
   inp_format=$(head -123  $inputseq | infoseq -nohead -only -USA -filter | tail -1 | awk -F : '{print $1}')
   #if [[ $inp_format != "fasta" ]]
   #then
   #   echo "Converting $inp_format formated query file to fasta format with seqret"     
   #   seqret $inputseq "$inputseq".fasta >& /dev/null
   #   inputseq=(""$inputseq".fasta")
   #fi 
   infoseq $inputseq -only -usa -nohead -outfile ${pbtmproot}/pb_"$$"_tmpdir/pb_"$$"_input_list >& /dev/null  
else
   # datan käsittely blast_formatterissa
   location=$(pwd)
   cd $pbtmproot/pb_"$$"_tmpdir/
   tar xf $inputseq
   num=$(ls *blast_database* | wc -l ) 
   if [[ $num > 1 ]] 
   then
      export BLASTDB=${pbtmproot}/pb_"$$"_tmpdir/  
   fi
   cd $location
   ls ${pbtmproot}/pb_"$$"_tmpdir/*result > ${pbtmproot}/pb_"$$"_tmpdir/pb_"$$"_input_list 
   cp ${pbtmproot}/pb_"$$"_tmpdir/pb_"$$"_input_list $pbtmproot/pb_"$$"_tmpdir/pb_"$$"_chunk_list
fi


r=$(cat $pbtmproot/pb_"$$"_tmpdir/pb_"$$"_input_list | wc -l)

if [[ $r -eq 0 ]] 
then
   echo "  "
   echo "Error: Input sequence $inputseq not found"
   exit 1
fi

if [[ $pb_output_type -eq 13 || $pb_output_type -eq 14 ]]
then
   if [[ $r -gt 5000 ]]
   then
     echo "Using output format $pb_output_type would produce in this case $r output files"
     echo "This is not permitted"
     echo "Please use some other output format"
     exit 1
   fi
fi 

# Less splitting if the database size is really small
if [[ $db_koko -lt 1000000 ]]; then
   n_split=10
fi

(( koko = r / n_split))
if [[ $blast_type == "blast_formatter" ]]
then 
     # Kaikki OK
   echo "blast_formatter"
else
   # No splitting in no_slurm mode
   if [[ $mode == "no_slurm" ]]
   then 
     cp $inputseq  ${pbtmproot}/pb_"$$"_tmpdir/pb_chunk_00001.fasta
   else
     echo "split"
     if [[ $koko -gt 2 ]] 
     then
         split_fasta.pl -p ${pbtmproot}/pb_"$$"_tmpdir/pb_chunk_ -s .fasta -o 1 $koko <  $inputseq  >& /dev/null
     else
         split_fasta.pl -p ${pbtmproot}/pb_"$$"_tmpdir/pb_chunk_ -s .fasta -o 1 2 <  $inputseq >& /dev/null
     fi
   fi
   ls ${pbtmproot}/pb_"$$"_tmpdir/pb_chunk_*.fasta > $pbtmproot/pb_$$_tmpdir/pb_$$_chunk_list
fi

n_split=$(cat ${pbtmproot}/pb_"$$"_tmpdir/pb_"$$"_chunk_list | wc -l )

if [[ $outputflag -ne 0 ]] 
then  
    echo "Total amount of query sequences: $r"
    echo "The job is split into $n_split pieces"
fi

location=$(pwd)
cd ${pbtmproot}/pb_"$$"_tmpdir/


#######
#Kirjoitetaan era-ajotiedosto
######


# 2 * 1048576 / 3 = 699051  
(( mem_request = db_koko / 699051 ))

#2.10.0 versiossa näyttäisi olevan pieni kokovaatimus
((  mem_request = mem_request / 2 ))

if [[ $mem_request -lt $mem_limit ]] 
then
  mem_request=$mem_limit
fi

if [[ $mem_request -gt 32000 ]] 
then 
  mem_request=32000
fi

if [[ $outputflag != 0 ]] 
then  
   echo "requesting memory $mem_request MB"
fi

#Set the number of cores:
if [[ $blast_type == "blast_formatter" ]]
then
  num_of_cores_to_use=(1)
else
  num_of_cores_to_use=(4)
fi

(( mem_request = mem_request / num_of_cores_to_use ))

#
#Create the batch job script
#
echo '#!/bin/bash -l' > $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo '#SBATCH -J pb_blast' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo '#SBATCH -o pb_blast_out_%A_%a.txt' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo '#SBATCH -e pb_blast_err_%A_%a.txt' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo '#SBATCH -t '$run_time':00:00' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo '#SBATCH --mem-per-cpu='"$mem_request" >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo '#SBATCH --array=1-'"$n_split" >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo "#SBATCH --account=${pname}" >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo '#SBATCH -p small' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo '#SBATCH -n 1' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo '#SBATCH --cpus-per-task='"$num_of_cores_to_use" >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo ""  >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo 'pbtmproot=("'$pbtmproot'")' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
#echo "cd $pbtmproot/pb_"$$"_tmpdir/jobs/job_\$SLURM_ARRAY_TASK_ID" >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo 'blasttest=$(which blastx 2>/dev/null | wc -c)' >>  $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo 'if [[ $blasttest -eq 0  ]]  ' >>  $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo 'then'  >>  $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo '  module load biokit 2> /dev/null ' >>  $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo 'fi' >>  $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
#echo "module load biokit 2> /dev/null " >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo "export BLASTDB=$BLASTDB" >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh

if [[ $mode == "no_slurm" ]]
then
echo "SLURM_ARRAY_TASK_ID=1" >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo "SLURM_CPUS_PER_TASK=3" >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
fi

echo 'nimi=(`sed -n "$SLURM_ARRAY_TASK_ID"p $pbtmproot/pb_'$$'_tmpdir/pb_'$$'_chunk_list `)' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo 'start_time=(`date +%s`)' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
if [[ $blast_type == "blast_formatter" ]]
then
    echo $command '-archive $nimi -out $nimi.result' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh

else
    echo $command '-num_threads $SLURM_CPUS_PER_TASK -query $nimi -db' $blast_database '-out $nimi.result.tmp' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
    echo 'mv $nimi.result.tmp $nimi.result'  >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
    echo 'rm -f $nimi' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
fi 

#  #add missing newline character to XML formatted outputfiles
#  if ( $pb_output_type == 5) then 
#      echo 'echo "" >> $nimi.result' >>  $WRKDIR/pb_$$_tmpdir/pb_blast_launch_tmp_$jobindex
#  endif

 # ID listan luonti tuloksista
if [[ $pb_output_type -eq 19 || $pb_output_type -eq 20 ]] 
then    
    echo 'numhits=$(cat $nimi.result | wc -l )' >>  $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
    echo 'if [[ $numhits -gt 0 ]] ; then ' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
    echo '  awk '\''{print $2}'\'' $nimi.result | sort | uniq  > $nimi.result.list' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
    echo '  rm -f $nimi.result' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
    echo '  mv $nimi.result.list $nimi.result' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
    echo 'fi' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh

fi 


#if ( $pb_output_type == 12) then
#    echo 'rm -f $nimi.result' >> $WRKDIR/pb_$$_tmpdir/pb_blast_launch_tmp_$jobindex
#    echo 'mv $nimi.result.list $nimi.result' >> $WRKDIR/pb_$$_tmpdir/pb_blast_launch_tmp_$jobindex
#endif
  
echo 'end_time=(`date +%s`)' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo '((duration = end_time - start_time))' >>  $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
echo 'echo '$jobindex' Duration: $duration'  >>  $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh


if [[ $pb_output_type -eq 21 ]] 
then 
  echo 'numhits=(`cat $nimi.result | wc -l `)' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
  echo 'if [[ $numhits -gt 0 ]]; then  ' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
  echo '  sort $nimi.result | uniq > $nimi.result.uniq' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
  echo '  rm -f $nimi.result' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh
  echo '  mv $nimi.result.uniq $nimi.result' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
  echo 'fi ' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
fi

if [[ $pb_output_type -eq 22 ]]
then
   echo 'numhits=(`cat $nimi.result | wc -l `)' >>  $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo 'if [[ $numhits -gt 0 ]] ; then  ' >>  $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo '  awk -F "," '\''{print $2}'\'' $nimi.result > $nimi.list' >>  $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo '  blastdbcmd -db '$blast_database' -entry_batch $nimi.list -outfmt %t -out $nimi.desc ' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo '  perl /appl/bio/blast/ncbi-blast-2.2.30+/bin/yhd.pl $nimi.result $nimi.desc COMMA > $nimi.ext ' >>  $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo '  rm -f $nimi.result ' >>  $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo '  mv $nimi.ext $nimi.result' >>  $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo 'fi ' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
fi

if [[ $pb_output_type -eq 23 ]]
then
   echo 'numhits=(`cat $nimi.result | wc -l `)' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo 'if [[ $numhits -gt 0 ]]; then  ' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo '  awk '\''{print $2}'\'' $nimi.result > $nimi.list' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo '  blastdbcmd -db '$blast_database' -entry_batch $nimi.list -outfmt %t -out $nimi.desc ' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo '  perl  /appl/bio/blast/ncbi-blast-2.2.30+/bin/yhd.pl $nimi.result $nimi.desc TAB > $nimi.ext ' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo '  rm -f $nimi.result ' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo '  mv $nimi.ext $nimi.result' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
   echo 'fi ' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
fi


if [[ $post_process -eq 1 ]]
then
  echo $pbtmproot/pb_"$$"_tmpdir/'pp_script $nimi.result >  $nimi.result.pp' >> $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
  echo 'mv -f $nimi.result.pp $nimi.result' >>  $pbtmproot/pb_"$$"_tmpdir/pb_batch.sh 
fi

#############
#Tyon lahetys
##############


cat <<EOF >  $pbtmproot/pb_"$$"_tmpdir/pb_blast.conf

mode=$mode
n_split=$n_split
location=$location
pb_output_type=$pb_output_type
blast_output=$blast_output
dbbuildflag=$dbbuildflag
pbtmproot=$pbtmproot
blastdb_orig=$blastdb_orig
ensemblflag=$ensemblflag
blast_database_name=$blast_database_name
blast_database=$blast_database
user=$user
start_time="$start_time"
pblogfile=$pblogfile
inputflag=$inputflag
outputflag=$outputflag
databaseflag=$databaseflag
dbtype=$dbtype
db_type=$db_type
debug=$debug
parse_seqids=$parse_seqids
blastdb_orig=$blastdb_orig
mem_limit=$mem_limit
run_time=$run_time
r=$r
EOF


cd  $pbtmproot/pb_"$$"_tmpdir/

if [[ $mode == "no_slurm" ]] 
then
  source pb_batch.sh > pb_blast_out 2> pb_blast_err
else 
 jobid=$(sbatch pb_batch.sh | grep "Submitted" | awk '{print $NF}')

 if [[ $jobid -gt 1  ]]
 then  
     echo "Batch job ID: $jobid"
 else
     echo "ERROR"
     echo "Automatic batch job submission failed!"
     echo "echo jobid value assigned:  $jobid"
     echo "Job submissin command provides following output:"
     echo ""
     echo "sbatch pb_batch.sh"
     sbatch pb_batch.sh
     exit 1
 fi

cat <<EOF >>  $pbtmproot/pb_"$$"_tmpdir/pb_blast.conf
jobid=$jobid
EOF
echo "starting blast_clusterrun with job ID: $$"
 blast_clusterrun -jobid $$ 

fi

exit 0


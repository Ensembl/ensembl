#!/bin/bash

#################################################
#                                               #
# Sanger Ensembl specific script to compute and #
# store Alternative splicing events in a core   #
# database.                                     #
# Author: Gautier Koscielny                     #
# e-mail: ensembl-dev mailing list              #
#                                               #
#################################################
 
port=3306
host=""
password=""
user=""
species=""
db=""
core=""
output_dir="/tmp"

print_help () { 
    echo "Usage:";
    echo "    as_event_computations.sh -h <dbhost> [-P <dbport>] -u <dbuser> [-p <dbpass>] -s <species> [-d <dbname>] [-o <output_dir>]";
    echo "";
    echo "If the species name is passed, the script will find the corresponding core database on <dbhost>.";
    echo "If the database name <dbname> is passed, the script will use this database as the core database.";
    echo "By default, all intermediate results will be written in the /tmp directory.";
    echo "Please use the -o parameter to pass a different existing writable directory.";
}

echo "Parsing parameters..."

while getopts ":h:u:p:P:s:d:o:" optname
  do
    #echo $optname
    case "$optname" in
      "h")
        echo "Hostname=$OPTARG";
				host=$OPTARG;
				echo $host
        ;;
      "o")
        echo "Output directory=$OPTARG";
				output_dir=$OPTARG;
        ;;
      "P")
        echo "Port=$OPTARG";
				port=$OPTARG;
        ;;
      "p")
        echo "Password=$OPTARG";
				password=$OPTARG;
        ;;
      "u")
        echo "User=$OPTARG";
				user=$OPTARG
        ;;
      "d")
        echo "Database=$OPTARG";
				db=$OPTARG
        ;;
      "s")
        echo "Species=$OPTARG";
				species=$OPTARG
        ;;
      "?")
        echo "Unknown option $OPTARG";
				print_help;
        ;;
      ":")
        echo "No argument value for option $OPTARG"
				print_help;
        ;;
      *)
      # Should not occur
        echo "Unknown error while processing options"
        ;;
    esac
    #echo "OPTIND is now $OPTIND"
  done

if [[ -n "$host" && -n "$user" ]]
then
  if [ -z "$db" ]
  then 
    if  [ -n "$species" ]
    then
      ## Find a core database for this species on the specified server.
	  if [[ -n "$password" && $password != "" ]]
	  then
				core=`mysql -h ${host} -P ${port} -u ${user} -p ${password} -s --skip-column-names -e "show databases like '${species}_core_%'"`
	  else
				core=`mysql -h${host} -P${port} -u${user} -s --skip-column-names -e "show databases like '${species}_core_%'"`
	  fi
    else
				echo "Species or/and database name are required." ;
				print_help;
				exit 1;
    fi
  else
			core=${db};
  fi
else
		echo "Hostname and username are required."
		print_help;
		exit 1;
fi

# count how many databases match this name

y=0;

for X in ${core}
do
		y=$[$y+1];
done

if [[ ${y} -eq 0 ]]
then
		echo "Check your parameters, there is no database matching species '${species}'";
		exit 1;
fi


if [[ ${y} -gt 1 ]]
then
		echo "Check your parameters, there are ${y} databases matching species '${species}':";
		echo "$core";
		exit 1;
fi


# Otherwise, start the pipeline

NOW=$(date +"%Y-%m-%d-%H-%M-%S")

echo "Starting pipeline with timestamp '${NOW}'";
if [ -n "$password" ]
then
		bsub -q normal -o ${output_dir}/as_gff_${species}_${NOW}.out -e ${output_dir}/as_gff_${species}_${NOW}.err -J as_gff_${species}_${NOW} "perl Fetch_gff.pl -dbname ${db} -dbhost ${host} -dbport ${port} -dbuser ${user} -dbpass ${password} -o ${output_dir}/${species}_${NOW}_variants.gff"
else
		bsub -q normal -o ${output_dir}/as_gff_${species}_${NOW}.out -e ${output_dir}/as_gff_${species}_${NOW}.err -J as_gff_${species}_${NOW} "perl Fetch_gff.pl -dbname ${db} -dbhost ${host} -dbport ${port} -dbuser ${user} -o ${output_dir}/${species}_${NOW}_variants.gff"
fi

bsub -q small -o ${output_dir}/as_compute_${species}_${NOW}.out -e ${output_dir}/as_compute_${species}_${NOW}.err -w "done(\"as_gff_${species}_${NOW}\")" -J as_compute_${species}_${NOW} "altSpliceFinder -i ${output_dir}/${species}_${NOW}_variants.gff -o ${output_dir}/${species}_${NOW}_events.gff --relax --statistics"

if [ -n "$password" ]
then
bsub -q normal -o ${output_dir}/as_populate_${species}_${NOW}.out -e ${output_dir}/as_populate_${species}_${NOW}.err -w "done(\"as_compute_${species}_${NOW}\")" -J as_populate_${species}_${NOW} "perl load_alt_splice_gff.pl -file ${output_dir}/${species}_${NOW}_events.gff -host ${host} -user ${user} -pass ${password} -dbname ${core}"
else 
		echo "Sorry, you did not provide any password. The script can't populate the database with the Alternative splicing information.";
		echo "However, the results are available in ${output_dir}/${species}_${NOW}_events.gff";
fi



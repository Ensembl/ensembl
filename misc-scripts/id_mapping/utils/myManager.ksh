#!/bin/ksh

function usage {
  cat <<EOT
Usage:
  $0 -h dbhost -P dbport -u dbuser -p dbpass -d dbname
EOT
}

if [[ ! -e ./manage_id_mapping_tables.pl ]]; then
  print -u2 "Expected to find the following executable file here:"
  print -u2 "\tmanage_id_mapping_tables.pl"
  exit
fi

dbport=3306

while getopts 'h:u:p:a:j:d:r:P:oci' opt; do
  case ${opt} in
    h)  dbhost=${OPTARG}    ;;
    P)  dbport=${OPTARG}    ;;
    u)  dbuser=${OPTARG}    ;;
    p)  dbpass=${OPTARG}    ;;
    d)  dbname=${OPTARG}    ;;
    *)  usage; exit         ;;
  esac
done

if [[
  -z ${dbhost} || -z ${dbport} ||
  -z ${dbuser} || -z ${dbpass} ||
  -z ${dbname}
]]; then
  usage
  exit
fi

./manage_id_mapping_tables.pl \
  -host ${dbhost} \
  -port ${dbport} \
  -user ${dbuser} \
  -pass ${dbpass} \
  -dbname ${dbname}


# $Id$

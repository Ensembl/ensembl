#!/bin/ksh
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


scriptname=$0

# THIS IS WHERE PERL WILL BE PICKED UP FROM:
export PATH=/software/perl-5.8.8/bin:${PATH}

function usage {
  cat <<EOT

Usage:

  ${scriptname} -h dbhost [-P dbport] -u dbuser -p dbpass -d dbname

EOT
}

if [[ ! -e ./manage_id_mapping_tables.pl ]]; then
  print -u2 "Expected to find the following executable file here:"
  print -u2 "\tmanage_id_mapping_tables.pl"
  exit
fi

dbport="3306"

while getopts 'h:P:u:p:d:' opt; do
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

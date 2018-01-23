#!/bin/bash
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

scan_and_replace() {
  search=$1
  replacement=$2
  echo "About to scan $(pwd) for files to replace '$search' with '$replacement'"

  for file in $(grep -R --files-with-matches "$search" .); do

    if [ "$(basename $file)" == "email_blatter.sh" ]; then
      echo "Skipping the email_blatter.sh script"
    else
      echo "Replacing email in $file"
      if [ "$(uname)" = "Darwin" ]; then
        LC_CTYPE=C LANG=C sed -i '' -e "s/$search/$replacement/g" $file
      else
        sed --in-place -e "s/$search/$replacement/g" $file
      fi
    fi
  done
}

if [[ "$#" -eq 0 ]]; then
  dirs=$(pwd)
else
  dirs=$@
fi

original_wd=$(pwd)

for var in $dirs; do

  if [ ! -d $var ] ; then
    echo "$var is not a directory. Skipping"
    continue
  fi

  cd $var

  scan_and_replace 'helpdesk@ensembl.org' 'http:\/\/www.ensembl.org\/Help\/Contact'
  scan_and_replace 'helpdesk&#64;ensembl.org' 'http:\/\/www.ensembl.org\/Help\/Contact'
  scan_and_replace 'dev@ensembl.org' 'http:\/\/lists.ensembl.org\/mailman\/listinfo\/dev'
  scan_and_replace 'dev&#64;ensembl.org' 'http:\/\/lists.ensembl.org\/mailman\/listinfo\/dev'
  scan_and_replace 'announce@ensembl.org' 'http:\/\\/lists.ensembl.org\/mailman\/listinfo\/announce'
  scan_and_replace 'announce&#64;ensembl.org' 'http:\/\\/lists.ensembl.org\/mailman\/listinfo\/announce'
  scan_and_replace 'vega@sanger.ac.uk' 'v@somewhere.com'
  
  cd $original_wd
done

#!/bin/bash
# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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
    echo "Replacing email in $file"
    sed -i -e "s/$search/$replacement/g" $file
  done
}

if [ -z "$@" ]; then
  dirs=$(pwd)
else
  dirs=$@
fi

original_wd=$(pwd)

for var in "$dirs"; do

  if [ ! -d $var ] ; then
    echo "$var is not a directory. Skipping"
    continue
  fi

  cd $var

  scan_and_replace 'helpdesk@ensembl.org' 'http:\/\/www.ensembl.org\/Help\/Contact'
  scan_and_replace 'dev@ensembl.org' 'http:\/\/lists.ensembl.org\/mailman\/listinfo\/dev'
  scan_and_replace 'announce@ensembl.org' 'http:\/\\/lists.ensembl.org\/mailman\/listinfo\/announce'
  
  cd $original_wd
done
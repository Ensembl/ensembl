#!/bin/bash
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2023] EMBL-European Bioinformatics Institute
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

  year=$(date "+%Y")
  last_year=$(($year - 1))

  search="^\(.*\)\\[\([0-9]*\)\(-*[0-9]*\)\\] EMBL-European Bioinformatics Institute"
  replacement="\1[\2-$year] EMBL-European Bioinformatics Institute"

  likely_copyright_line=".*[0-9].* EMBL-European Bioinformatics Institute"
  exception_line="Copyright \\[1999-2015\\] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute"

  echo "About to scan $(pwd) for files to replace '$search' with '$replacement'"

  for file in $(grep -r --files-with-matches "$likely_copyright_line" --exclude-dir=.git .); do

    if [[ $(grep "$search" $file) ]]; then
      echo "Replacing date in $file"
    fi

    if [ "$(uname)" = "Darwin" ]; then
      LC_CTYPE=C LANG=C sed -i '' -e "s/$search/$replacement/g" $file
    else
      sed --in-place -e "s/$search/$replacement/g" $file
    fi

    # If the line matches the more general $likely_copyright_line, but not the more specific $search regex that we
    # know how to update, and isn't the $exception_line, then we warn the user that the line is likely a copyright
    # line that we don't know how to update.
    grep "$likely_copyright_line" $file | while read -r line ; do
      if [[ $(echo $line | grep -v "$search") && $(echo $line | grep -v "$exception_line") ]]; then
        echo "Unable to replace date in $file.  The line is: $line"
      fi
    done

  done

  cd $original_wd
done

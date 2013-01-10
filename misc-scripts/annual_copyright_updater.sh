#!/bin/bash

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

  year=$(date "+%Y")
  last_year=$(($year - 1))

  search="Copyright (c) 1999-${last_year}"
  replacement="Copyright (c) 1999-${year}"

  echo "About to scan $(pwd) for files to replace '$search' with '$replacement'"

  for file in $(grep -R --files-with-matches "$search" .); do
    echo "Replacing date in $file"
    sed -i '' -e "s/$search/$replacement/g" $file
  done

  cd $original_wd
done
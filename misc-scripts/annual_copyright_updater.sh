#!/bin/bash

function confirm()
{
    echo -n "$@ "
    read -e answer
    for response in y Y yes YES Yes Sure sure SURE OK ok Ok
    do
        if [ "_$answer" == "_$response" ]
        then
            return 0
        fi
    done
	
    return 1
}

year=$(date "+%Y")
last_year=$(($year - 1))

search="Copyright (c) 1999-${last_year}"
replacement="Copyright (c) 1999-${year}"

confirm "About to scan $(pwd) for files to replace '$search' with '$replacement'. Ok?"

for file in $(grep -R --files-with-matches "$search" .); do
	echo "Replacing date in $file"
	sed -i "s/$search/$replacement/g" $file
done


#!/bin/ksh -ex

# $Id$
#
# Creates delta files between all consecutive revisions of all
# databases on ftp.ensembl.org using build.pl and apply.pl.
#
# Author:  Andreas Kahari <andreas.kahari@ebi.ac.uk>
#

export LANG=C

ftpsite='ftp.ensembl.org'
ftppass='${LOGNAME}@$(hostname).$(domainname)'

dbdir='./databases'
deltadir='./deltas'

build_cmd='./build.pl'
apply_cmd='./apply.pl'
time_cmd='/usr/bin/time'
xdelta_cmd='./xdelta.osf'
perl_cmd='/usr/local/ensembl/bin/perl -w'

trapsigs="INT HUP TERM"

#-------------------------------------------------------------

# Function: file_list
# Usage: file_list
#
# Downloads the ls-lR.Z file off the FTP site into the current
# working directory.  Extracts the names of the files that we are
# interested in from it and outputs them on standard output.  The
# format of the output is "path dbname version", where "path" is
# the path of the FTP directory, "dbname" is the name of the
# database, and "version" is the version of the database.  The
# output is sorted on "dbname", then on "version".

function file_list
{
    ftp -i -n -v <<-EOT >/dev/null
	open ${ftpsite}
	user anonymous ${ftppass}
	binary
	get ls-lR.Z
	EOT

    gunzip -c ls-lR.Z |
    grep 'data/mysql/.*[0-9][0-9]*_[0-9][0-9]*' | grep -v 'pub/NEW' |
    sed -n 's/^\(.*\)\/\([^\/]*\)_\([0-9][0-9]*_[0-9][0-9]*.*\):$/\1 \2 \3/p' |
    sort -k2,2 -k3,3
}

#-------------------------------------------------------------

# Function: cleanup
# Usage: cleanup dbname [version]
#
# Remove a downloaded database from ${dbdir}.  If "version" is
# ommited, remove all versions of the database.

function cleanup
{
    typeset dbname=$1
    typeset version=${2:-'*'}

    if [[ ! -d ${dbdir} ]]; then
	return
    fi

    rm -f -r ${dbdir}/${dbname}_${version}
}

#-------------------------------------------------------------

# Function: fetch_db
# Usage: fetch_db path dbname version
#
# Fetches version "version" of the database "dbname" at the
# path "path" off the ${ftpsite}.  The database will be stored
# in ${dbdir}.

function fetch_db
{
    typeset path=$1
    typeset dbname=$2
    typeset version=$3

    if [[ -d ${dbdir}/${dbname}_${version} ]]; then
	return
    fi
    mkdir -p ${dbdir}/${dbname}_${version}

    trap "rm -rf ${dbdir}/${dbname}_${version}; exit 1" ${trapsigs}
    (
	cd ${dbdir}
	ftp -i -n -v <<-EOT
		open ${ftpsite}
		user anonymous ${ftppass}
		binary
		cd ${path}
		mget ${dbname}_${version}/*
		EOT
    )
    trap - ${trapsigs}
}

#-------------------------------------------------------------

# Function: build_delta
# Usage: build_delta dbname opath oversion path version
#
# Records the changes between version "oversion" and "version"
# of database "dbname".  Also tests the generated delta files.

function build_delta
{
    typeset dbname=$1
    typeset opath=$2
    typeset oversion=$3
    typeset path=$4
    typeset version=$5

    typeset outdir=${deltadir}/to_${version%_*[0-9]*}
    mkdir -p ${outdir}

    typeset bout=${outdir}/${dbname}_${oversion}_delta_${version}_build.out
    typeset aout=${outdir}/${dbname}_${oversion}_delta_${version}_apply.out

    if [[ ! -f ${bout} ]]; then
	fetch_db ${opath} ${dbname} ${oversion}
	fetch_db ${path} ${dbname} ${version}

	trap "rm ${bout}; exit 1" ${trapsigs}
	${time_cmd} ${perl_cmd} ${build_cmd} -c ${xdelta_cmd} \
	    -s ${dbdir} -d ${outdir} \
	    ${dbname} ${oversion} ${version} | tee ${bout}
	trap - ${trapsigs}
    fi

    if [[ ! -f ${aout} ]]; then
	trap "rm ${aout}; exit 1" ${trapsigs}
	${time_cmd} ${perl_cmd} ${apply_cmd} -c ${xdelta_cmd} \
	    -d ${dbdir} -s ${outdir} \
	    ${dbname} ${oversion} ${version} | tee ${aout}
	trap - ${trapsigs}
    fi
}

#-------------------------------------------------------------

file_list |
while read path dbname version; do
    if [[ -n ${odbname} ]]; then
	if [[ ${odbname} != ${dbname} ]]; then
	    cleanup ${odbname}

	    opath=${path}
	    odbname=${dbname}
	    oversion=${version}

	    continue
	fi

	build_delta ${dbname} ${opath} ${oversion} ${path} ${version}
	cleanup ${dbname} ${oversion}
    fi

    opath=${path}
    odbname=${dbname}
    oversion=${version}
done

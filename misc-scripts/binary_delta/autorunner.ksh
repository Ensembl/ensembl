#/bin/ksh -x

# $Id$
#
# Creates delta files between all consecutive revisions of all
# databases on ftp.ensembl.org using build.pl and apply.pl.
#
# Author:  Andreas Kahari <andreas.kahari@ebi.ac.uk>
#

# Gets a database from the FTP site, puts it in ./databases/
function getdb
{
    typeset path=$1
    typeset db=$2
    typeset ver=$3

    typeset dbver=${db}_${ver}

    if [[ ! -d databases/${dbver} ]]; then
	trap "rm -rf databases/${dbver}; exit 1" INT
	mkdir -p databases/${dbver}
	( cd databases/${dbver}
	  ftp -i -n -v ftp.ensembl.org <<EOT
user anonymous ak@ebi.ac.uk
cd ${path}/${dbver}
bin
mget *
bye
EOT
	)
	trap - INT
    fi
}

# Builds delta files, and applies them for verification
function do_delta
{
    typeset db=$1
    typeset v1=$2
    typeset v2=$3

    typeset path1=$4
    typeset path2=$5

    if [[ ! -d deltas ]]; then
	mkdir deltas
    fi

    typeset build_out=deltas/${db}_${v1}_delta_${v2}_build.out
    typeset apply_out=deltas/${db}_${v1}_delta_${v2}_apply.out

    if [[ ! -f $build_out ]]; then
	getdb $path1 $db $v1
	getdb $path2 $db $v2

	trap "rm $build_out; exit 1" INT

	/usr/bin/time perl -w ./build.pl -c ./xdelta.osf \
	    -s databases -d deltas \
	    $db $v1 $v2 2>&1 | tee $build_out
	trap - INT
    fi
    if [[ ! -f $apply_out ]]; then
	trap "rm $apply_out; exit 1" INT
	/usr/bin/time perl -w ./apply.pl -c ./xdelta.osf \
	    -d databases -s deltas \
	    $db $v1 $v2 2>&1 | tee $apply_out
	trap - INT
    fi
}

# Removes a database that was fetched
function cleandb
{
    typeset db=$1
    typeset ver=$2

    if [[ -d databases && -n $db && -n $ver ]]; then
	rm -rf databases/${db}_${ver}
	rm -rf databases/${db}_${ver}.????
    fi
}

# ---------------  Main

# For debugging
typeset -ft getdb
typeset -ft do_delta
typeset -ft cleandb

# A regular expression that should be avoided
avoid_re='mart'

# A regular expression that should be required
require_re='1[0-5]_'

version_re='[0-9][0-9]*_[0-9][0-9]*'

# Use ftp://ftp.ensembl.org/ls-lR.Z to figure out what files are
# available
lynx -source ftp://ftp.ensembl.org/ls-lR.Z | \
    sed -n 's/^\.\(.*data\/mysql\)\/\(.*\)_\('"$version_re"'\):$/\1 \2 \3/p' | \
    grep -v $avoid_re | grep $require_re | sort -k2 >ls-lR

while read path db ver; do
    if [[ $db != $this_db ]]; then
	cleandb $this_db $ver
	cleandb $this_db $old_ver
	this_db=$db
	old_ver=$ver
	old_path=$path
	continue
    fi

    do_delta $this_db $old_ver $ver $old_path $path
    cleandb $this_db $old_ver

    old_ver=$ver
    old_path=$path
done <ls-lR

rm -rf databases/*

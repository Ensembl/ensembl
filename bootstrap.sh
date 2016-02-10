#!/usr/bin/env bash

join_array() { local d=$1; shift; echo -n "$1"; shift; printf "%s" "${@/#/$d}"; }

pathremove () {
        local IFS=':'
        local NEWPATH
        local DIR
        local PATHVARIABLE=${2:-PATH}
        for DIR in ${!PATHVARIABLE} ; do
                if [ "$DIR" != "$1" ] ; then
                  NEWPATH=${NEWPATH:+$NEWPATH:}$DIR
                fi
        done
        export $PATHVARIABLE="$NEWPATH"
}

pathprepend () {
        pathremove $1 $2
        local PATHVARIABLE=${2:-PATH}
        export $PATHVARIABLE="$1${!PATHVARIABLE:+:${!PATHVARIABLE}}"
}

pathappend () {
        pathremove $1 $2
        local PATHVARIABLE=${2:-PATH}
        export $PATHVARIABLE="${!PATHVARIABLE:+${!PATHVARIABLE}:}$1"
}

# Fetch the git tools so start the install process
echo "Fetching Ensembl install tools"
git clone --depth 1 https://github.com/Ensembl/ensembl-git-tools.git

# Put the git tools in the path
pathprepend $PWD/ensembl-git-tools/bin
pathprepend $PWD/ensembl-git-tools/advanced_bin

# Parameters for git
params=()
extra_repos=()

# Find what environment we're installing
GROUP=${GROUP:-"api"}

# What release are we targeting, if not default
if [ ! -z "$RELEASE" ]; then
    if [[ $RELEASE =~ ^-?[0-9]+$ ]]; then
	RELEASE="release/$RELEASE"
    fi
    params+=("--branch $RELEASE")
fi

if [ ! -z "$TEST_MODULE" ]; then
    params+=("--ignore_module $TEST_MODULE")
    ENS_TEST="true"
fi

# If we're setting up a test environment, we'll need
# the testing repo
if [ ! -z "$ENS_TEST" ] && [ "$ENS_TEST" = 'true' ]; then
    echo "yes"
    extra_repos+=('ensembl-test')
fi

# Prepare the extra parameters for the git call
param_str=$(join_array ' ' ${params[@]})
extra_repos_str=$(join_array ' ' ${extra_repos[@]})

echo "Fetching Ensembl repositories"
git ensembl --clone $param_str --depth 1 $GROUP $extra_repos_str

echo "Installing dependencies"
git ensembl --cmd install-dep $GROUP $extra_repos_str

echo "Installing BioPerl"
wget http://bioperl.org/DIST/BioPerl-1.6.1.tar.gz
mkdir bioperl-live && tar zxvf BioPerl-1.6.1.tar.gz -C bioperl-live --strip-components 1

echo "Installing common Ensembl dependencies"
git clone --branch master --depth 1 https://github.com/samtools/tabix.git
cd tabix
make
cd perl
perl Makefile.PL
make
cd ..
ln -sf perl/blib blib
cd ../
git clone --branch master --depth 1 https://github.com/samtools/htslib.git
cd htslib
make
export HTSLIB_DIR=$(pwd -P)
cd ../
git clone --branch master --depth 1 https://github.com/Ensembl/Bio-HTS
cd Bio-HTS
perl Build.PL $HTSLIB_DIR
./Build

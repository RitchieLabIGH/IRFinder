#!/bin/bash

if [[ "$1" != "local" ]] && [[ "$1" != "" ]] && [[ "$1" != "remove" ]]; then
    echo -e "Usage: \nGlobal installation:\tsudo install.sh\nLocal installation:\tinstall.sh local" >&2
    echo -e "Uninstall all:\tsudo install.sh remove\nUninstall local:\tinstall.sh remove local\nUninstall local:\tsudo install.sh remove global" >&2
    exit 1
fi

function checkDependencies(){
    local distro=$(cat /proc/version )
    local deps=0
    echo "Checking dependencies..."
    for pkg in $@ ; do
        if [[ "${distro}" =~ Ubuntu|Debian ]]  ; then
            if ! dpkg -s $pkg >/dev/null 2>/dev/null ; then 
                echo "Dependency $pkg not found." >&2
                deps=1
            fi
        else
            if ! rpm -q $pkg >/dev/null 2>/dev/null ; then
                echo "Dependency $pkg not found." >&2 
                deps=1
            fi
        fi
    done
    if [ $deps -eq 1 ]; then
        exit 1
    fi
}


if [[ $1 == "remove" ]]; then
    if [[ "$2" != "global" ]] && [[ "$2" != "local" ]] && [[ "$2" != "" ]] ; then
        echo "Error: $2 not recognized. Use 'local' or 'global' or leave empty" >&2
        exit 1
    fi
    if [[ "$2" == "" || "$2" == "global" ]]; then
        if [ "$EUID" -ne 0 ]; then 
            echo "Please run as root"
            exit
        fi
        if [ -d /usr/local/IRFinder ]; then
            rm -fr /usr/local/IRFinder /usr/bin/IRFinder
            echo "Removed system installation"
        else
            echo "Global installation of IRFinder not found"
        fi
    fi
    if [[ "$2" == "" || "$2" == "local" ]] ;then
        if [ -d ~/.local/IRFinder ] ; then
            rm -fr ~/.local/IRFinder ~/.local/bin/IRFinder
            echo "Removed local installation"
        else
            echo "Local installation of IRFinder not found."
        fi
    fi
    exit
fi


if [[ "${1}" != "local" ]]; then
    if [ "$EUID" -ne 0 ]; then 
        echo "Please run as root or to install IRFinder locally call"
        echo "./install.sh local"
        echo ""
        exit 1
    fi
fi

checkDependencies "make bedtools samtools gzip gawk libboost-iostreams-dev zlib1g"


ORIGINAL_FOLDER=$(realpath $PWD)
BASE_FOLDER=$(dirname "$(readlink -nf "$BASH_SOURCE")")

cd $BASE_FOLDER/src/trim/
make clean
make
cp ./trim $BASE_FOLDER/bin/util/trim
make clean
cd ../winflat
make clean
make
cp ./winflat $BASE_FOLDER/bin/util/winflat
cd ../irfinder/Release
make clean
make
cp ./irfinder $BASE_FOLDER/bin/util/irfinder
make clean
cd $BASE_FOLDER
chmod -R a+x ./bin
if [[ "${1}" == "local" ]];then
    if [ -d ~/.local/IRFinder ]; then
        rm -fr ~/.local/IRFinder ~/.local/bin/IRFinder
    fi
    cp -r $BASE_FOLDER ~/.local/IRFinder
    ln -s $(realpath ~/.local/IRFinder/bin/IRFinder) ~/.local/bin/IRFinder
else
    if [ -d /usr/local/IRFinder ]; then
        rm -fr /usr/local/IRFinder /usr/bin/IRFinder
    fi
    cp -r $BASE_FOLDER /usr/local/IRFinder 
    ln -s /usr/local/IRFinder/bin/IRFinder /usr/bin/IRFinder
fi

cd $ORIGINAL_FOLDER


if ! which suppa.py >/dev/null 2>/dev/null ; then
  echo "SUPPA2 not found. To use the RunMode Diff, install it: https://github.com/comprna/SUPPA " >&2
fi

if ! which STAR > /dev/null 2> /dev/null ; then
  echo "STAR not found. To use the RunMode FastQ and to produce your own mapability files during the reference build, install it: https://github.com/alexdobin/STAR " >&2
fi

if ! which minimap2 > /dev/null 2> /dev/null ; then
  echo "minimap2 not found. To use the RunMode Long, install it: https://github.com/lh3/minimap2 " >&2
fi







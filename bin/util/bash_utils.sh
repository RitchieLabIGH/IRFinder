#!/bin/bash
## Useful functions for IRFinder's utils

export IRFINDER_BASH_UTILS_IMPORTED=1
export VERSION=2.0.0
export LC_ALL=C
export LANG=C


function versionAlert(){
   echo "IRFinder version: $VERSION" 
   exit
}

function checkFile() {
    if [ ! -f "${1}" ]; then
        echo "Error: file $1 doesn't exists" >&2
        exit 1
    fi
}

function checkSamtools() {
    STVERSTR=`samtools --version`
	STVER=$(echo $STVERSTR|cut -d" " -f2)
	STVERMAIN=$(echo $STVER|cut -d"." -f1)
	STVERMINOR=$(echo $STVER|cut -d"." -f2)
	if [[ ! "$STVERMAIN" -ge 1 ]]; then
		echo "Error: Samtools $STVER: version too old (>=1.4 required)." >&2
		exit 1
	elif [[ ! "$STVERMINOR" -ge 4 ]]; then
		echo "Error: Samtools $STVER: version too old (>=1.4 required)." >&2
		exit 1
	fi
}

function getMem(){
    local MEMK=`awk '($1 ~ /^MemTotal:/) {print $2}' < /proc/meminfo`
    echo $(($MEMK/1000))
}

function getSortingMem() {
    MEMM=$(getMem)
    SORTMEM=$(echo "${MEMM}" | awk '{if ($1 > 10000 ) { print 10000 } else { if ($1 < 500) { print 500 } else { print $1 }  }   }' )
    echo $SORTMEM 
}

function checkStar(){
    MEMM=$(getMem)
    if [ "${MEMM}" -lt 32000 ]; then
        echo "System limitation: Minimum required RAM is 32GB. This software uses STAR for RNA mapping. RAM requirement is approximately 30GB for the human genome." >&2
        echo "  RunModes: BAM and BuildRefDownload, may be completed on servers with more RAM." >&2
        exit 2
    fi
    if [[ "$1" != "" ]]; then
        STAREXEC="$1"
    fi 
    if [[ "${STAREXEC}" == "" ]]; then
        STAREXEC="STAR"
    fi
    "$STAREXEC" --version &>/dev/null
    if [ ! $? -eq 0 ]; then
        echo "Error: STAR version is too old. --version parameter returns an error. Minimum version of 2.4.0 required." >&2
        exit 2
    fi
}


function checkMinimap(){
    if [[ "$1" != "" ]]; then
        MINIMAP_EXEC="$1"
    fi 
    if [[ "${MINIMAP_EXEC}" == "" ]]; then
        MINIMAP_EXEC="minimap2"
    fi
    if ! which $MINIMAP_EXEC > /dev/null 2> /dev/null ; then
      echo "minimap2 not found ( executable: $MINIMAP_EXEC ). To use the RunMode Long, install it: https://github.com/lh3/minimap2 " >&2
      exit 1
    fi
    MINIMAP_VERSION=$("$MINIMAP_EXEC" --version)
    if [[ $(echo ${MINIMAP_VERSION/-*/} | awk '{if ( $1 > 2.0 ) {print "ok" } else { print "no" }}') != "ok" ]]; then
        echo "Error: Minimap version is too old. Minimum version of 2.0.0 required. ${MINIMAP_VERSION} detected" >&2
        exit 2
    fi
}


function checkSuppa(){
    if ! which suppa.py >/dev/null 2>/dev/null ; then
        echo "SUPPA2 not found ( executable: suppa.py ). To use the RunMode Diff, install it: https://github.com/comprna/SUPPA " >&2
        exit 1
    fi
}

function checkDeseq(){
	if ! which Rscript > /dev/null 2>/dev/null; then
		echo "Rscript not found."
		exit 1
	fi
	DESeqVersion=$(Rscript -e 'installed.packages()' | awk  'BEGIN {v=0} $1=="Version" {v=1; } v==1 && $1 == "DESeq2" { gsub("\"", ""); print $2;v=0 } ' )
	
	if [[ "${DESeqVersion}" == "" ]]; then
		DESeqVersion=$(Rscript -e 'installed.packages()' | awk  'BEGIN {v=0} $NF=="Version" {v=1; } v==1 && $1 == "DESeq2" { gsub("\"", ""); print $NF;v=0 } ' )
		if [[ "${DESeqVersion}" == "" ]]; then
			echo "DESeq2 not installed. "
			exit 1
		fi
	fi
	logger "DESeq2 version $DESeqVersion"
}

function setThreads(){
    if [[ "${THREADS}" == "" ||  $THREADS == 0 ]]; then
	    THREADS=`grep -c ^processor /proc/cpuinfo`    
    	if [ ! -n $THREADS ] | [ $THREADS -eq 0 ]; then
        	THREADS=`awk 'BEGIN {FS=":"} ($0 ~ /^physical id/ ) { printf $2 " --"} ($0 ~ /^core id/) {print $2}' < /proc/cpuinfo | sort -u | wc -l`    	
    	    if [ ! -n $THREADS ] | [ $THREADS -eq 0 ]; then
        	    THREADS=1
        	fi
	    fi
    fi
}


function checkRef(){
    if [ ! "$1" ]; then
    	echo "Argument error: -r is required." >&2
    	exit 1	
    fi
    if [ ! -f "$1/IRFinder/ref-cover.bed" ]; then
		echo "Argument error: -r $1, Does not appear to be a valid IRFinder reference. Could not find $1/IRFinder/ref-cover.bed" >&2
		exit 1
	fi
}

function checkOutDir(){
    local OUTPUTDIR=$1
    if [ -d "$OUTPUTDIR" ]; then
		if [ -e "$OUTPUTDIR/IRFinder-IR-nondir.txt" ]; then
			echo "Argument error: -d $OUTPUTDIR, output directory contains files from a previous IRFinder run. Will not overwrite." >&2
			exit 1
		else
    		mkdir -p "$OUTPUTDIR/logs/"
		fi
	else
		mkdir -p "$OUTPUTDIR/logs/"
		if [ ! -d "$OUTPUTDIR" ]; then
			echo "Argument error: Output directory $OUTPUTDIR does not exist, and could not be created." >&2
			exit 1
		fi
	fi
}

function logger() {
    LOGOUT="./irfinder.stdout"
    if [[ "$OUTPUTDIR" != "" ]]; then
        if [ ! -d ${OUTPUTDIR}/logs/ ]; then 
            mkdir -p ${OUTPUTDIR}/logs/
        fi
        LOGOUT="$OUTPUTDIR/logs/irfinder.stdout" 
    fi
    if [[ "$1" == "init" ]] && [[ $# == 1 ]]; then
        > $LOGOUT
        LOG_MESSAGE="\n --------------------\n|  IRFinder v. $VERSION | \n --------------------\n"
    else
        LOG_MESSAGE="${@}"
    fi    
    if [[ "${VERBOSE}" == "1"  ]] ; then
        echo -e "${LOG_MESSAGE}" | tee -ai $LOGOUT 
    else
        echo -e "${LOG_MESSAGE}" >> $LOGOUT
    fi
}

function startMessage(){
    ## Check if the startMessage was called by the BAM mode after the FastQ or Long analysis
    if [[ "${IRF_RUNMODE}" == "" ]]; then
        logger "---" 
        logger "IRFinder version: $VERSION " 
        logger "IRFinder start: " `date` 
        logger "IRFinder runmode: $RUNMODE"
        logger "IRFinder user@host: $USER @ $HOSTNAME" 
        logger "IRFinder working dir: " `pwd` 
        logger "IRFinder reference: $REF" 
        n=1
        for f in $@; do
            logger "IRFinder file ${n}: $f" 
            n=$((n+1))
        done
        logger "---"  
        START_MESSAGE=1
    fi
}

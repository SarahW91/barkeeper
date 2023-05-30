#!/bin/bash

echo 'checking for bioconda'
bioconda_check=$(which conda)
if [ -z "$conda_check" ]
then
	echo $'conda not installed\nInstalling now.'

else
	echo $'conda already installed\nproceeding...'
fi

echo 'checking for ccs'
ccs_check=$(which ccs)
if [ -z "$ccs_check" ] 
then
	echo $'ccs not installed\nInstalling now.'
	conda install -c bioconda pbccs
else
	echo $'ccs already installed\nproceeding..'
fi

echo 'checking for lima'
lima_check=$(which lima)
if [ -z "$lima_check" ] 
then
	echo $'lima not installed\nInstalling now.'
	conda install -c bioconda lima
else
	echo $'lima already installed\nproceeding..'
fi

echo 'checking for samtools'
samtools_check=$(which samtools)
if [ -z "$samtools_check" ] 
then
	echo $'samtools not installed\nInstalling now.'
	wget 'https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2'
	tar -xvjf samtools-1.9.tar.bz2
	cd samtools-1.9/
	./configure
    make
    make install
    export PATH=/where/to/install/bin:$PATH 
else
	echo $'samtools already installed\nproceeding..'
fi

echo 'checking for seqtk'
seqtk_check=$(which seqtk)
if [ -z "$seqtk_check" ] 
then
	echo $'seqtk not installed\nInstalling now.'
else
	echo $'seqtk already installed\nproceeding..'
fi

echo 'checking for usearch'
usearch_check=$(which usearch)
if [ -z "$usearch_check" ] 
then
	echo $'usearch not installed\nInstalling now.'
else
	echo $'usearch already installed\nproceeding..'
fi

echo 'checking for blastn'
blastn_check=$(which blastn)
if [ -z "$blastn_check" ] 
then
	echo $'blastn not installed\nInstalling now.'
else
	echo $'blastn already installed\nproceeding..'
fi
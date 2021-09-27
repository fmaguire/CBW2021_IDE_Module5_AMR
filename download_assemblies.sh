#!/bin/env bash

# Sample IDs selected for demo purposes

sourcelink="http://ftp.ebi.ac.uk/pub/databases/ENA2018-bacteria-661k"
ref="http://ftp.ebi.ac.uk/pub/databases/ENA2018-bacteria-661k/sampleid_assembly_paths.txt"

neisseria=(SAMD00099400 SAMD00099401 SAMD00099402 SAMD00099403 SAMD00099404 SAMD00099405 SAMD00099406 SAMD00099407 SAMD00099408 SAMD00099409 SAMD00099410)
pseudomonas=(SAMD00019033 SAMD00019034 SAMD00019035 SAMD00019036 SAMD00019037 SAMD00019038 SAMD00019039 SAMD00019040 SAMD00019041 SAMD00019042 SAMD00019043)
salmonella=(SAMN06286049 SAMN06286050 SAMN06286051 SAMN06286052 SAMN06286053 SAMN06286054 SAMN06286055 SAMN06286056 SAMN06286057 SAMN06286058 SAMN06286059)

if [ -d "assemblies" ]; then
	echo "Directory 'assemblies' already exists! Remove or rename and try again!"
	exit 1
else
	mkdir -p assemblies/{neisseria,pseudomona,salmonella}

fi

if ! [ -f "assemblies/ref_ids.txt" ]; then
	wget $ref -O assemblies/ref_ids.txt
fi

# FIX ARRAY CALLING
for species in {'neisseria','pseudomonas','salmonella'}; do
	if [ $species = "neisseria" ]; then
		for id in ${neisseria[@]}; do
			path=$(grep "${id}" assemblies/ref_ids.txt | cut -f 2 | rev | cut -d/ -f-1,-2,-3 | rev)
			#echo $path
			wget "${sourcelink}/${path}" -P assemblies/neisseria
		done
	fi
	if [ $species = "pseudomonas" ]; then
		for id in ${pseudomonas[@]}; do
			path=$(grep "${id}" assemblies/ref_ids.txt | cut -f 2 | rev | cut -d/ -f-1,-2,-3 | rev)
                	#echo $path
                	wget "${sourcelink}/${path}" -P assemblies/pseudomonas
		done
	fi
	if [ $species = "salmonella" ]; then
		for id in ${salmonella[@]}; do
			path=$(grep "${id}" assemblies/ref_ids.txt | cut -f 2 | rev | cut -d/ -f-1,-2,-3 | rev)
                	#echo $path
                	wget "${sourcelink}/${path}" -P assemblies/salmonella
		done
	fi
done

rm assemblies/ref_ids.txt

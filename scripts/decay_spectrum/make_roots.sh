#!/bin/bash
#runs single decays and creates root-files from the resulting .tsv files

declare -a arr=("element1" "element2" "element3")

for i in "${arr[@]}"
do
    cd ../..
    ./decay_chain 12 10 50 0 $i
    mv $i scripts/decay_spectrum
    cd scripts/decay_spectrum
    ./create_root $i
done


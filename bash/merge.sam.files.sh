#!/bin/bash
clear
num_cores=4

ids=$((for f in $(find data-raw/sam.files/all -name '*.sam'); do
  basename ${f} | grep -o '^z[[:digit:]]*'
done) | sort | uniq)

echo --- Merge and Index SAM ---
  for i in ${ids}; do
  echo ${i}
  samtools merge -O SAM data-raw/sam.files/all/${i}.merged.sam \
    $(find data-raw/sam.files/all -name "${i}*.sam") \
    -f \
    -@ ${num_cores}
  samtools sort -O SAM $(find data-raw/sam.files/all -name "${i}*.merged.sam") \
  -o data-raw/sam.files/all/${i}.merged.sam \
    -@ ${num_cores}
done

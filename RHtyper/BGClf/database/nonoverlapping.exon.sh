#!/bin/bash

### generate non-overlapping exon positions

module load CASAVA/1.6

#script="nonoverlappingExonCoords.pl" ### from CASAVA
#$script <  hg19.refFlat.txt > hg19.refFlat.nonoverlapping.txt
#$script <  hg38.refFlat.txt > hg38.refFlat.nonoverlapping.txt


### my own code, CDS exon table downloaded from UCSC
#./nonoverlappingExonCoords.R hg38.exon.gz
./nonoverlappingExonCoords.R hg19.exon.gz


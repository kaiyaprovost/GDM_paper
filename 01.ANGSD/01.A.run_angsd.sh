#!/bin/bash
## STEP 0: INITIALIZE
spp="bellii" # #NOTE: replace species 
ref=Corv_chrom/pseudochromosomesSHORT.fasta ## NOTE: replace reference per species
bamlist=${spp}.bamlist ## NOTE: this will get overwritten during this pipeline if you run it
names=("AMN_245111_P002_WA06" "AMN_245111_P002_WA07" "AMN_245111_P002_WB07" "AMN_245111_P002_WB11" "AMN_245111_P002_WC01" "AMN_245111_P002_WD02" "AMN_245111_P002_WD07" "AMN_245111_P002_WD11" "AMN_245111_P002_WE08" "AMN_245111_P002_WF02" "AMN_245111_P002_WF03" "AMN_245111_P002_WF05" "AMN_245111_P002_WG05" "AMN_245111_P002_WH08" "AMN_245112_P001_WA06" "AMN_245112_P001_WB06" "AMN_245112_P001_WC06" "AMN_245112_P001_WE02")
## NOTE: replace names per species
## STEP 1: REHEADER
echo "Start Reheadering"
for i in *${spp}*.dedup.bam; do
echo $i;
findname=$(echo $i | cut -f 1 -d '.' | cut -f 3 -d "/")
replace=$(echo $i | cut -d "_" -f 2-5 | cut -f 3 -d "/"); 
if [ ! -f ${findname}*reheadered.bam ] ; then
echo "Reheadering"
~/nas3/samtools/samtools view -H $i | sed -e "s/SM:${findname}/SM:${replace}/" | ~/nas3/samtools/samtools reheader - $i > $i.dedup.reheadered.bam
else
echo "Bam already reheadered"
fi
done; 
## STEP 2: CLIP OVERLAPS TO REALIGN
for reheaderbam in *${spp}*.dedup.reheadered.bam; do
name=${reheaderbam%.dedup.reheadered.bam}
name=`echo $name | rev | cut -d'/' -f 1 | rev`
echo "#######################"
echo $name
echo "#######################"
echo "Indel Realigner"
bam clipOverlap --in $reheaderbam --noeof --out $name.clipped.bam --stats --params --unmapped --storeOrig
done
## STEP 3: MAKE REALIGNER TARGETS
echo "RealignerTargetCreator"
command="time java -jar GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${ref} "
for prebam in *${spp}*.dedup.reheadered.bam; do 
postbam=${prebam%.dedup.reheadered.bam}
name=`echo $postbam | rev | cut -d'/' -f 1 | rev`
inputbam=${name}.clipped.bam
inputbai=${inputbam%.bam}.bai
if [ ! -f $inputbai ]; then
time java -jar picard.jar BuildBamIndex I=$inputbam O=$inputbai
fi
command="${command}-I ${inputbam} ";
done
command="${command} -o ${spp}-allsamplesClipped-final.intervals > ${spp}-make-intervals.log 2>&1 "
echo
echo $command
eval $command
## STEP 4: REALIGN THE INDELS
for clippedbam in *${spp}*.clipped.bam; do
name=${clippedbam%.bam}
echo "#######################"
echo $name
echo "#######################"
echo Indel Realigner
time java -jar GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T IndelRealigner -R $ref -I $clippedbam -targetIntervals ${spp}-allsamplesClipped-final.intervals -o $name.realigned.bam
done
## STEP 5: COMBINE BAMS FROM MULTIPLE RUNS IF NEEDED
for name in ${names[@]}; do
echo "#######################"
echo $name
echo "#######################"
echo Merging Bam Files
command="time java -jar picard.jar MarkDuplicates TMP_DIR=tmp "
for namebam in *$name*.realigned.bam; 
do echo $namebam; 
command="${command}I=${namebam} "; 
done;
merged="${name}.merged.realigned.bam"
command="${command} O=${merged}"
command="${command} METRICS_FILE=$name.merged.realigned.metrics.txt REMOVE_DUPLICATES=false MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1024 TAGGING_POLICY=All"
echo "COMMAND IS:"
echo $command
eval $command
done
ls *merged.realigned.bam > $bamlist
## STEP 6: CALCULATE GENOTYPE LIKELIHOODS
date
time
time angsd -GL 2 -doMajorMinor 1 -doMaf 1 -doGlf 2 -minMaf 0.05 -SNP_pval 0.01 -minInd 4 -minMapQ 20 -minQ 20 -nThreads 16 -bam $bamlist -out ${spp}
## STEP 7: CALL SNPS
date
time
time angsd -GL 2 -bam $bamlist -doGeno 5 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf -minInd 4 -minMaf 0.05 -minMapQ 20 -minQ 20 -nThreads 16 -postCutoff 0.95 -SNP_pval 0.01 -out ${spp}-called 

spp="bellii" ## change as needed
ref=Corv_chrom/pseudochromosomesSHORT.fasta ## change as needed
bamlist=${spp}.bamlist
sonlist=SON_${spp}.bamlist
chilist=CHI_${spp}.bamlist
time realSFS ${spp}-taj1.saf.idx > ${spp}-taj.sfs
time realSFS SON_${spp}-taj1.saf.idx > SON_${spp}-taj.sfs
time realSFS CHI_${spp}-taj1.saf.idx > CHI_${spp}-taj.sfs
time angsd -GL 2 -bam $bamlist -doThetas 1 -doSaf 1 -pest ${spp}-taj.sfs -anc $ref -out ${spp}-taj2 
time angsd -GL 2 -bam $sonlist -doThetas 1 -doSaf 1 -pest SON_${spp}-taj.sfs -anc $ref -out ${spp}-taj2 
time angsd -GL 2 -bam $chilist -doThetas 1 -doSaf 1 -pest CHI_${spp}-taj.sfs -anc $ref -out ${spp}-taj2 
#Estimate for every Chromosome/scaffold
time thetaStat do_stat ${spp}-taj2.thetas.idx
time thetaStat do_stat SON_${spp}-taj2.thetas.idx
time thetaStat do_stat CHI_${spp}-taj2.thetas.idx
#Do a sliding window analysis based on the output from the make_bed command.
time thetaStat do_stat \
${spp}-taj2.thetas.idx \
-win 100000 \
-step 10000  \
-outnames ${spp}-taj2.thetasWindow.gz
time thetaStat do_stat \
SON_${spp}-taj2.thetas.idx \
-win 100000 \
-step 10000  \
-outnames SON_${spp}-taj2.thetasWindow.gz
time thetaStat do_stat \
CHI_${spp}-taj2.thetas.idx \
-win 100000 \
-step 10000  \
-outnames CHI_${spp}-taj2.thetasWindow.gz
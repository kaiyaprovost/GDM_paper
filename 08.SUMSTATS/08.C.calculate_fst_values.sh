spp="bellii" ## change as needed
ref=Corv_chrom/pseudochromosomesSHORT.fasta ## change as needed
bamlist=${spp}.bamlist
sonlist=SON_${spp}.bamlist
chilist=CHI_${spp}.bamlist
date
time
## CALCULATE AN SFS
time angsd -GL 2 -dosaf 1 -fold 1 -minInd 4 -minMapQ 20 -minQ 20 -nThreads 16 -ref $ref -anc $ref -bam $bamlist -out ${spp}-sfs1
## get the global sfs
time realSFS -nSites 50000000 ${spp}-sfs1.saf.idx
## get the by population sfs 
angsd -b $sonlist -anc $ref -out SON_${spp}-sfs1 -dosaf 1 -gl 1
angsd -b $chilist -anc $ref -out CHI_${spp}-sfs1 -dosaf 1 -gl 1
## now get the 2d sfs for each 
command="time realSFS -nSites 50000000 SON_${spp}-sfs1.saf.idx CHI_${spp}-sfs1.saf.idx > SON_CHI_${spp}.ml"
eval $command
## make sure the ML file is only a single row -- if not, you must sum the rows
echo "Prep for window analysis"
angsd sites index SON_${spp}-sfs1.saf.pos.gz
angsd sites index CHI_${spp}-sfs1.saf.pos.gz
# prepare the fst for easy window analysis etc
time realSFS fst index SON_${spp}-sfs1.saf.idx CHI_${spp}-sfs1.saf.idx SON_CHI_${spp}.ml -fstout SON_CHI_${spp}_FST 
echo "Get global estimate"
# get the global estimate
time realSFS fst stats SON_CHI_${spp}_FST.fst.idx 
echo "Sliding Windows" 
time realSFS fst stats2 SON_CHI_${spp}_FST.fst.idx -win 100000 -step 10000 > SON_CHI_${spp}_FST_slidingwindow.fst
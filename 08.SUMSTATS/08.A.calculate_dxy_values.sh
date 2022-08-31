spp="bellii" ## change as needed
bamlist=${spp}.bamlist
sonlist=SON_${spp}.bamlist
chilist=CHI_${spp}.bamlist
## GENERATE GENOTYPE LIKELIHOOD FOR ALL THREE POPULATIONS
time angsd -GL 2 -doMajorMinor 1 -doMaf 11 -doCounts 1 -doGlf 2 -minMaf 0.05 -SNP_pval 0.01 -skipTriallelic -minInd 4 -minMapQ 20 -minQ 20 -nThreads 16 -bam $bamlist -out ${spp}-DXY
time angsd -GL 2 -doMaf 11 -doCounts 1 -doMajorMinor 1 -skipTriallelic -ref $ref -sites CHROMS_FROM_DXY_RAWSIZES.txt -nThreads 16 -bam $sonlist -out ${spp}-SONDXY
time angsd -GL 2 -doMaf 11 -doCounts 1 -doMajorMinor 1 -skipTriallelic -ref $ref -sites CHROMS_FROM_DXY_RAWSIZES.txt -nThreads 16 -bam $chilist -out ${spp}-CHIDXY
## CALCULATE DXY ACROSS SLIDING WINDOWS
sonmaf=${spp}-SONDXY.mafs
chimaf=${spp}-CHIDXY.mafs
sonlen=`cat $sonmaf | wc -l`
chilen=`cat $chimaf | wc -l`
gunzip -f $sonmaf.gz
gunzip -f $chimaf.gz
## don't think you need totLen to generate the Dxy thing 
echo "son maf file is: ${sonmaf}"
echo "length: ${sonlen}"
echo "chi maf file is: ${chimaf}"
echo "length: ${chilen}"
Rscript calcDxy.R --popA $sonmaf --popB $chimaf --totLen $sonlen
mv Dxy_persite.txt ${spp}_Dxy_persite.txt

spp="bellii" ## change as needed
GL=$spp.beagle.gz
date
time
minMaf=0.05
maxiter=2000
misTol=0.2
minInd=6
tol=0.000001 
for k in {1..4}; do
	echo                          
	echo "#########################"
	echo $k
	echo "#########################"
	echo                          
	mkdir $k
	cd $k
	for seed in {120410..120419}; do
		echo $seed
		outfile="$k.$seed.$minMaf.$minInd.$misTol.$maxiter.$tol.${spp}_NGSadmix"
		if [ ! -f $outfile ]; then
			time NGSadmix \
			-likes $GL \
			-K $k \
			-P 16 \
			-o $outfile \
			-minMaf $minMaf \
			-seed $seed -maxiter $maxiter \
			-minInd $minInd -misTol $misTol -tol $tol
		fi
	done;
done

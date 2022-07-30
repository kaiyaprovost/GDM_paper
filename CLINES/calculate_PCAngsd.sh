spp="bellii"

for GL in *${spp}*.beagle.gz; do
	echo "###################"
	echo $GL
	outname=${GL%.beagle.gz}
	echo $outname
	for k in {1..4}; do
		echo $k
		# Estimate covariance matrix and individual admixture proportions
		echo "RUN"
		date
		time
		python pcangsd.py \
		-beagle $GL \
		-admix -admix_save -e $k \
		-inbreed 2 \
		-selection 1 -sites_save \
		-o ${outname}.${spp}_PCAngsd.$k \
		-threads 32
	done
done

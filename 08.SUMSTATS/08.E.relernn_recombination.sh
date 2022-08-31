spp="bellii" ## change as needed
SIMULATE="ReLERNN_SIMULATE"
TRAIN="ReLERNN_TRAIN"
PREDICT="ReLERNN_PREDICT"
BSCORRECT="ReLERNN_BSCORRECT"
CPU="1"
MU="2.21e-9"
RTR="1"
for vcffile in *{spp}*.vcf; do
	DIR="$PWD"
	VCF="$vcffile"
	GENOME="$vcffile.bed"
	MASK="$vcffile.mask"
	## check if file already exists 
	if [ -f "${vcffile%/vcf}.PREDICT.txt" ]; then
		echo "SKIPPING"
	else
		if [ ! -f "$GENOME" ]; then
			echo "##### MAKING BED FILES 1 #####" >> simulate.log;
			python3 generate_relernn_bed_files.py $vcffile
		fi
		if [ ! -f "$MASK" ]; then
			echo "##### MAKING BED FILES 2 #####" >> simulate.log;
			python3 generate_relernn_bed_files.py $vcffile
		fi
		## for this to work MUST have this header line:
		## ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
		## Simulate data
		echo "##### SIMULATING #####" >> simulate.log;
		${SIMULATE} \
		--vcf ${VCF} \
		--genome ${GENOME} \
		--unphased \
		--projectDir ${DIR} \
		--assumedMu ${MU} \
		--upperRhoThetaRatio ${RTR} \
		--forceDiploid \
		--nTrain 12800 \
		--nVali 2000 \
		--nTest 100 \
		--nCPU ${CPU} >> simulate.log 2>&1
		## Train network
		echo "##### TRAINING #####" >> simulate.log;
		${TRAIN} \
			--projectDir ${DIR} \
			--nEpochs 2 \
			--nValSteps 2 >> simulate.log 2>&1
		## Predict
		echo "##### PREDICTING #####" >> simulate.log;
		${PREDICT} \
			--vcf ${VCF} \
			--projectDir ${DIR} >> simulate.log 2>&1
		## Parametric Bootstrapping
		echo "##### BOOTSTRAPPING #####" >> simulate.log;
		${BSCORRECT} \
			--projectDir ${DIR} \
			--nCPU ${CPU} \
			--nSlice 2 \
			--nReps 2 >> simulate.log 2>&1
		echo "DONE!"
		gzip $vcffile
	fi
done;
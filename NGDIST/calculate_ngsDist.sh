spp="bellii"
GL=${spp}.beagle.gz
lines=`zcat $GL | wc -l`
sites="$(($lines-1))"

date
time

ngsDist \
--n_threads 16 \
--geno $GL \
--n_ind 18 \
--n_sites $sites \
--n_boot_rep 300 \
--boot_block_size 500 \
--out ${spp}_ngsdist
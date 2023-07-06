#! /usr/bin/env bash

set -e
# clone repos used for simulation
git clone git@github.com:bguo068/posseleff_simulations.git
git -C posseleff_simulations checkout c58f9d988c9827ce870f859fa1230a0679fc4e32

# create env and compile some c/c++ tools
mamba env create -f posseleff_simulations/env.yaml 
conda activate simulation
cd posseleff_simulations
python init.py
cd ..

# install hmmibd that allows to change recombination rate
git clone git@github.com:bguo068/hmmibd-rs.git
cd hmmibd-rs
x86_64-conda-linux-gnu-gcc -o hmmIBD_rec -O3 -lm -Wall hmmIBD.c
cargo build --release --bin hmmibd2
cd ..
cp hmmibd-rs/hmmIBD_rec  posseleff_simulations/bin/
cp hmmibd-rs/target/release/hmmibd2 posseleff_simulations/bin/
rm -rf hmmibd-rs

# modified some code to allow output heterzygous diploid vcf files
sed -i  's/write_peudo_homozygous_vcf(simulator/# write_peudo_homozygous_vcf(simulator/' posseleff_simulations/bin/sim_single_pop.py
echo '
    with open(f'\''{args.genome_set_id}_{args.chrno}.vcf'\'', "w") as f:
        simulator.mutated_trees.write_vcf(f, contig_id=f"chr{args.chrno}")
' >>  posseleff_simulations/bin/sim_single_pop.py

# download hapibd jar file
wget https://faculty.washington.edu/browning/hap-ibd.jar

# simulation and get tree/vcf files
for chrno in {1..3}; do 
    echo $chrno
    posseleff_simulations/bin/sim_single_pop.py \
        --chrno $chrno \
        --seqlen 1500000 \
        --selpos  500000 \
        --num_origins 1 \
        --N 10000 \
        --s 0.0 \
        --h 0.5 \
        --g_sel_start 80 \
        --r 6.6666667e-7 \
        --sim_relatedness 0 \
        --g_ne_change_start 200 \
        --N0 1000 \
        --u 1e-8 \
        --nsam 1000 \
        --genome_set_id 03
done

if [ ! -d ibd_tskibd ]; then mkdir ibd_tskibd; fi
if [ ! -d ibd_hapibd ]; then mkdir ibd_hapibd; fi
if [ ! -d ibd_hmmibd ]; then mkdir ibd_hmmibd; fi
if [ ! -d ibd_hmmibd2 ]; then mkdir ibd_hmmibd2; fi
if [ ! -d map ]; then mkdir map; fi
if [ ! -d bcf ]; then mkdir bcf; fi
if [ ! -d vcf_filt ]; then mkdir vcf_filt; fi

for chrno in {1..3}; do 
    mv 3_${chrno}.trees ibd_tskibd/

    # NOTE: it is important that the sample node ids and individual Ids are not in the same order
    # fix order
    tskit nodes ibd_tskibd/3_${chrno}.trees | awk '$2==1 && NR%2==1 {print "tsk_" $5}'  > sample_order.txt
    # make the sample consistent across chromosomes
    cat sample_order.txt  | awk -v OFS='\t' '{print $1, "tsk_" NR-1}' > sample_name_map.txt

    bcftools view -S sample_order.txt 3_${chrno}.vcf | bcftools reheader -s sample_name_map.txt | bcftools view -Ob -o bcf/sel_chr${chrno}.bcf
    bcftools index -f bcf/sel_chr${chrno}.bcf

    # ## fix ALT/REF due to slim modification
    bcftools view -m2 -M2 -q 0.01:minor bcf/sel_chr${chrno}.bcf |  awk -v OFS='\t' '/^#/{print}/^[^#]/{$4="A"; $5="C"; print}' | bcftools view -Oz -o vcf_filt/sel_chr${chrno}.vcf.gz
    bcftools index -f vcf_filt/sel_chr${chrno}.vcf.gz

    # ## call tskibd
    posseleff_simulations/bin/tskibd ${chrno} 15000 150 2 ibd_tskibd/3_${chrno}.trees
    mv ${chrno}.ibd ibd_tskibd/
    echo "chr${chrno} . 0.0 1" > ${chrno}.map
    echo "chr${chrno} . 100.0 1500000" >> ${chrno}.map
    mv ${chrno}.map map/chr${chrno}.map

    # ## call hapibd
    java -jar hap-ibd.jar gt=vcf_filt/sel_chr${chrno}.vcf.gz map=map/chr${chrno}.map min-output=2 out=ibd_hapibd/chr${chrno}

    # convert file
    echo | awk -v OFS='\t' 'END{$1="CHROM"; $2 = "POS"; for (i=0; i<1000; i++) $(3+i)=i; print}' > hmmibdinput.txt
    bcftools query -f "%CHROM\t%POS[\t%GT]\n" vcf_filt/sel_chr${chrno}.vcf.gz | tr '|' '\t' | sed 's/^chr//' >> hmmibdinput.txt
    ./posseleff_simulations/bin/hmmIBD_rec -i hmmibdinput.txt -n 100 -m 5 -r  6.6666667e-7 -o ./ibd_hmmibd/$chrno
    ./posseleff_simulations/bin/hmmibd2 -i hmmibdinput.txt -n 100 -m5 --rec-rate  6.6666667e-7 --max-all 2 --filt-ibd-only  --filt-min-seg-cm 2.0 -o ./ibd_hmmibd2/$chrno
done

# clear intermediate files
rm 1* 2* 3*

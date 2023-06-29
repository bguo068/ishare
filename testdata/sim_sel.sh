# clone an existing repo
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
cd ..
cp hmmibd-rs/hmmIBD_rec  posseleff_simulations/bin/


# modified some code to allow output heterzygous diploid vcf files
sed -i  's/write_peudo_homozygous_vcf(simulator/# write_peudo_homozygous_vcf(simulator/' posseleff_simulations/bin/sim_single_pop.py
echo '
    with open(f'\''{args.genome_set_id}_{args.chrno}.vcf'\'', "w") as f:
        simulator.mutated_trees.write_vcf(f, contig_id=f"chr{args.chrno}")
' >>  posseleff_simulations/bin/sim_single_pop.py

# simulation and get tree/vcf files
for chrno in {1..3}; do 
    echo $chrno
    posseleff_simulations/bin/sim_single_pop.py \
        --chrno $chrno \
        --seqlen 100000000 \
        --selpos 30000000 \
        --num_origins 1 \
        --N 10000 \
        --s 0.3 \
        --h 0.5 \
        --g_sel_start 80 \
        --r 1e-8 \
        --sim_relatedness 0 \
        --g_ne_change_start 200 \
        --N0 2000 \
        --u 1e-8 \
        --nsam 2000 \
        --genome_set_id 03
done

if [ ! -d ibd_tskibd ]; then mkdir ibd_tskibd; fi
if [ ! -d ibd_hapibd ]; then mkdir ibd_hapibd; fi
if [ ! -d ibd_hmmibd ]; then mkdir ibd_hmmibd; fi
if [ ! -d map ]; then mkdir map; fi
if [ ! -d bcf ]; then mkdir bcf; fi
if [ ! -d vcf_filt ]; then mkdir vcf_filt; fi
wget https://faculty.washington.edu/browning/hap-ibd.jar

for chrno in {1..3}; do 
    bcftools view 3_${chrno}.vcf -Ob -o bcf/sel_chr${chrno}.bcf
    bcftools index -f bcf/sel_chr${chrno}.bcf

    # ## fix ALT/REF due to slim modification
    bcftools view -m2 -M2 -q 0.01:minor bcf/sel_chr${chrno}.bcf |  awk -v OFS='\t' '/^#/{print}/^[^#]/{$4="A"; $5="C"; print}' | bcftools view -Oz -o vcf_filt/sel_chr${chrno}.vcf.gz
    bcftools index -f vcf_filt/sel_chr${chrno}.vcf.gz

    # ## call tskibd
    posseleff_simulations/bin/tskibd ${chrno} 1000000 10000 2 3_${chrno}.trees
    mv ${chrno}.ibd > ibd_tskibd/
    mv 3_${chrno}.trees ibd_tskibd/
    mv ${chrno}.map map/chr${chrno}.map

    # ## call hapibd
    java -jar hap-ibd.jar gt=vcf_filt/sel_chr${chrno}.vcf.gz map=map/chr${chrno}.map min-output=2 out=ibd_hapibd/chr${chrno}

    # convert file
    bcftools query -f "%CHROM\t%POS[\t%GT]\n" vcf_filt/sel_chr${chrno}.vcf.gz | tr '|' '\t' | sed 's/^chr//' > hmmibdinput.txt
    ./posseleff_simulations/bin/hmmIBD_rec -i hmmibdinput.txt -n 100 -r 1e-8 -o ./ibd_hmmibd/$chrno

done




rm 3_* tmp* *.log *.map
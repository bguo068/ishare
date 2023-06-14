import msprime
from subprocess import run
from pathlib import Path
import urllib.request

Path("bcf").mkdir(parents=True, exist_ok=True)
Path("vcf_filt").mkdir(parents=True, exist_ok=True)
Path("map").mkdir(parents=True, exist_ok=True)
Path("ibd1").mkdir(parents=True, exist_ok=True)
Path("ibd2").mkdir(parents=True, exist_ok=True)
if not Path("hap-ibd.jar").exists():
    urllib.request.urlretrieve(
        "https://faculty.washington.edu/browning/hap-ibd.jar", "hap-ibd.jar"
    )

for i in range(3):
    # Simulate ancestry and mutations, generate vcf file
    ts = msprime.sim_ancestry(
        samples=1000,
        ploidy=2,
        sequence_length=100_000_000,
        random_seed=2 + 10 * i,
        recombination_rate=1e-8,
        population_size=10_000,
    )
    ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=2)
    with open("test_sim.vcf", "w") as f:
        ts.write_vcf(f, contig_id=f"chr{i+1}")

    # convert to bcf files
    run(
        f"""
        bcftools view test_sim.vcf -Ob -o bcf/chr{i+1}.bcf
        rm -f test_sim.vcf
        bcftools index -f bcf/chr{i+1}.bcf
        """,
        shell=True,
        check=True,
    )

    # write genetic map files
    with open(f"map/chr{i+1}.map", "w") as f:
        f.write(f"chr{i+1} . 0.0 1\n")
        f.write(f"chr{i+1} . 100.0 100000000\n")

    # filter vcf and run hap-ibd
    cmd = f"""
    bcftools view -m2 -M2 bcf/chr{i+1}.bcf | bcftools view -q 0.01:minor -Oz -o vcf_filt/chr{i+1}.vcf.gz
    java -jar hap-ibd.jar gt=vcf_filt/chr{i+1}.vcf.gz map=map/chr{i+1}.map out=ibd1/chr{i+1} \
        min-seed=3 min-output=3
    java -jar hap-ibd.jar gt=vcf_filt/chr{i+1}.vcf.gz map=map/chr{i+1}.map out=ibd2/chr{i+1} \
        min-seed=2 min-output=2
    """
    run(cmd, shell=True, check=True)


cmd = f"""
    rm -f hap-ibd.jar
    bcftools concat bcf/chr1.bcf bcf/chr2.bcf bcf/chr3.bcf -Ob -o bcf/all.bcf
    bcftools index -f bcf/all.bcf
    bcftools view -r chr1:1-10000000 -m2 -M2 bcf/all.bcf -Ob -o bcf/small.bcf 
    bcftools index -f bcf/small.bcf
    bcftools query -l bcf/all.bcf > samples.txt
"""

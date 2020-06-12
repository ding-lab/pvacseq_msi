import os
import pandas as pd
from pathlib import Path

configfile: "config.yaml"

SAMPLE = pd.read_table(config["samples"])["sample"].drop_duplicates().tolist()
REF_DICT = pd.read_table(config["samples"])[["sample","reference"]].drop_duplicates().set_index("sample")["reference"].to_dict()
DIS_DICT = pd.read_table(config["samples"])[["sample","dispath"]].drop_duplicates().set_index("sample")["dispath"].to_dict()
FA_DICT = pd.read_table(config["samples"])[["sample","refpath"]].drop_duplicates().set_index("sample")["refpath"].to_dict()
HLA_DICT = pd.read_table(config["samples"])[["sample","hla_type"]].drop_duplicates().set_index("sample")["hla_type"].to_dict()

rule chi:
    input:
        lambda wildcards: DIS_DICT[wildcards.sample]
    output:
        temp("{sample}.chi_dis")
    log:
        "logs/{sample}.chi.log"
    shell:
        "zcat {input} | python scripts/chi.py - {output}"

rule dis2vcf:
    input:
        "{sample}.chi_dis"
    output:
        temp("{sample}.chi.vcf")
    log:
        "logs/{sample}.dis2vcf.log"
    shell:
        "perl scripts/MSISensor2Vepvcf.pl {input} {output}"

rule vep:
    input:
        vcf = "{sample}.chi.vcf",
        ref = lambda wildcards: REF_DICT[wildcards.sample],
        fa = lambda wildcards: FA_DICT[wildcards.sample]
    output:
        "{sample}.chi.vep.vcf"
    log:
        "logs/{sample}.vep.log"
    singularity: 
        "docker://ensemblorg/ensembl-vep:release_97.4"
    shell:
        "vep --input_file {input.vcf} \
             --output_file {output} \
             --assembly {input.ref} \
             --fasta {input.fa} \
             --force_overwrite --species homo_sapiens --cache --check_ref --dir /diskmnt/Datasets/VEP/ --vcf --offline --symbol --terms SO --tsl --hgvs --plugin Downstream --plugin Wildtype --format vcf"

rule pvacseq:
    input:
        vcf = "{sample}.chi.vep.vcf",
        hla = lambda wildcards: HLA_DICT[wildcards.sample],
        sample = "{sample}"
    output:
        "{sample}"
    log:
        "logs/{sample}.pvacseq.log"
    singularity:
        "docker://griffithlab/pvactools"
    shell:
        "pvacseq run {input.vcf} {input.sample} {input.hla} \
                      NetMHCpan {output} -e 8,9,10 --iedb-install-directory /opt/iedb"

rule all:
    input: 
        expand("{sample}.chi.vep.vcf", sample=SAMPLE)

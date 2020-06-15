import os
import pandas as pd
from pathlib import Path

configfile: "config.yaml"

OUTPUT = config["path2output"]
VEP = config["path2vep"]


SAMPLE = pd.read_table(config["samples"])["sample"].drop_duplicates().tolist()
REF_DICT = pd.read_table(config["samples"])[["sample","reference"]].drop_duplicates().set_index("sample")["reference"].to_dict()
DIS_DICT = pd.read_table(config["samples"])[["sample","dispath"]].drop_duplicates().set_index("sample")["dispath"].to_dict()
FA_DICT = pd.read_table(config["samples"])[["sample","refpath"]].drop_duplicates().set_index("sample")["refpath"].to_dict()
#HLA_DICT = pd.read_table(config["samples"])[["sample","hla_type"]].drop_duplicates().set_index("sample")["hla_type"].to_dict()

rule unzip:
    input:
        lambda wildcards: DIS_DICT[wildcards.sample]
    output:
        temp(f"{OUTPUT}{{sample}}_dis")
    shell:
        "zcat {input} > {output}"

rule chi:
    input:
        f"{OUTPUT}{{sample}}_dis"
    output:
        temp(f"{OUTPUT}{{sample}}.chi_dis")
    log:
        f"{OUTPUT}logs/{{sample}}.chi.log"
    shell:
        "python scripts/chi.py {input} {output}"

rule dis2vcf:
    input:
        f"{OUTPUT}{{sample}}.chi_dis"
    output:
        f"{OUTPUT}{{sample}}.chi.vcf"
    log:
        f"{OUTPUT}logs/{{sample}}.dis2vcf.log"
    shell:
        "perl scripts/MSISensor2Vepvcf.pl {input} {output}"

#rule vep:
#    input:
#        vcf = f"{OUTPUT}{{sample}}.chi.vcf",
#        fa = lambda wildcards: FA_DICT[wildcards.sample]
#    output:
#        f"{OUTPUT}{{sample}}.chi.vep.vcf"
#    log:
#        f"{OUTPUT}logs/{{sample}}.vep.log"
#    params:
#        i = "/data/{sample}.chi.vcf",
#        o = "/data/{sample}.chi.vcf.vcf",
#        veppath = f"{VEP}",
#        ref = lambda wildcards: REF_DICT[wildcards.sample]
#    container:
#        "docker://ensemblorg/ensembl-vep:release_97.4"
#    shell:
#        "vep --input_file {params.i} \
#             --output_file {params.o} \
#             --assembly {params.ref} \
#             --fasta {input.fa} \
#             --force_overwrite --species homo_sapiens --cache --check_ref \
#             --dir {params.veppath} \
#             --vcf --offline --symbol --terms SO --tsl --hgvs --plugin Downstream --plugin Wildtype --format vcf"

#rule pvacseq:
#    input:
#        vcf = f"{OUTPUT}{{sample}}.chi.vep.vcf",
#        hla = lambda wildcards: HLA_DICT[wildcards.sample],
#        sample = "{sample}"
#    output:
#        f"{OUTPUT}{{sample}}"
#    log:
#        f"{OUTPUT}logs/{{sample}}.pvacseq.log"
#    singularity:
#        "docker://griffithlab/pvactools"
#    shell:
#        "pvacseq run {input.vcf} {input.sample} {input.hla} \
#                      NetMHCpan {output} -e 8,9,10 --iedb-install-directory /opt/iedb"

rule all:
    input: 
        expand(f"{OUTPUT}{{sample}}.chi.vcf", sample=SAMPLE)

# pvacseq_msi
Identify immunogenic indels from MSIsensor output

## Authour
Wen-Wei Liang (liang.w@wustl.edu)

## Run the pipeline on katmai
### Preparation
1. Create the conda envrionment:
```
conda create -n pvacseq_msi python=3.7 snakemake=5.19.2 pandas=1.0.1 singularity=3.5.3 scipy
```

2. Activate the environment: 
```
conda activate pvacseq_msi
```

3. Modify `config.yaml` as needed.

### When singularity is supported by lab cluster
Run the following steps if the user has been add to fakeroot list by adim.

4. Set a TMPDIR for singularity: 
```
export SINGULARITY_TMPDIR=/diskmnt/Projects/Users/wliang/MSIsensor2/
```

5. Change the `singularity.conf` at `${Path2Conda/env/pvacseq_msi/etc/singularity/singularity.conf}` as following. 
```
mksquashfs path = ${Path2Cond}/env/pvacseq_msi/bin/
```

6. Run the pipeline by:
```
snakemake -p all --cores 1 --use-singularity \
          --singularity-args "--bind /diskmnt/Datasets/VEP/:/diskmnt/Datasets/VEP/" \
          --singularity-args "--bind /diskmnt/Datasets/Reference/:/diskmnt/Datasets/Reference/" \
          --singularity-args "--bind /diskmnt/Projects/Users/wliang/MSIsensor2/03_msi2vcf/output/:/data/"
```

### When singularity is not supported

4. Convert dis files to vcfs
```
snakemake -p all --cores 40
``` 

5. Manually run VEP. Need to change the access of output folder so VEP can write.
```
docker pull wenwiliang/vep2pvacseq
chmod 777 {path2output} {path2reference}
docker run -v /diskmnt/Datasets/VEP:/diskmnt/Datasets/VEP -v /diskmnt/Projects/Users/wliang/MSIsensor2/03_msi2vcf/output:/diskmnt/Projects/Users/wliang/MSIsensor2/03_msi2vcf/output -v /diskmnt/Datasets/Reference:/diskmnt/Datasets/Reference wenwiliang/vep2pvacseq vep --input_file /diskmnt/Projects/Users/wliang/MSIsensor2/03_msi2vcf/output/01BR001.chi.vcf --output_file /diskmnt/Projects/Users/wliang/MSIsensor2/03_msi2vcf/output/01BR001.chi.vep.vcf --force_overwrite --species homo_sapiens --cache --assembly GRCh37 --fasta /diskmnt/Datasets/Reference/GRCh37-lite/GRCh37-lite-chr_with_chrM.fa --check_ref --dir /diskmnt/Datasets/VEP/ --vcf --offline --symbol --terms SO --tsl --hgvs --plugin Downstream --plugin Wildtype --format vcf
```

6. Manually run pVACseq. Noted that the HLA types has been predicted by OptiType.
```
docker pull griffithlab/pvactools
docker run -v /diskmnt/Projects/Users/wliang/MSIsensor2/03_msi2vcf/output:/diskmnt/Projects/Users/wliang/MSIsensor2/03_msi2vcf/output griffithlab/pvactools pvacseq run /diskmnt/Projects/Users/wliang/MSIsensor2/03_msi2vcf/output/01BR001.chi.vep.vcf 01BR008 HLA-A*02:05,HLA-A*23:01,HLA-B*58:01,HLA-B*49:01,HLA-C*07:01,HLA-C*07:01 NetMHCpan /diskmnt/Projects/Users/wliang/MSIsensor2/03_msi2vcf/output/01BR001 -e 8,9,10 -t 30 --iedb-install-directory /opt/iedb
```

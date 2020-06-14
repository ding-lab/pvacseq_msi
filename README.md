# pvacseq_msi
Identify immunogenic indels from MSIsensor output

## Authour
Wen-Wei Liang (liang.w@wustl.edu)

## Run the pipeline on katmai

1. Create the conda envrionment:
```
conda create -n pvacseq_msi python=3.7 snakemake=5.19.2 pandas=1.0.1 singularity=3.5.3 scipy
```

2. Activate the environment: 
```
conda activate pvacseq_msi
```

3. Modify `config.yaml` as needed.

The following steps could run if the user has been add to fakeroot list by adim.

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

# pvacseq_msi
Identify immunogenic indels from MSIsensor output

## Authour
Wen-Wei Liang (liang.w@wustl.edu)

## Run the pipeline on katmai

1. Create the conda envrionment:
```
conda create -n pvacseq_msi python=3.6 snakemake-minimal=5.9.1  pandas=1.0.1 singularity scipy
```

2. Activate the environment: 
```
conda activate pvacseq_msi
```
3. Modify `config.yaml` as needed.
4. Run the pipeline by:
```
snakemake -p all --use-singularity --singularity-args "--bind ${path2data}:/data" --cores 1

```

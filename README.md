## A snakemake pipeline for fungal pathogen gene annotation

Under development for Pst. Should work for other fungal pathogens too but needs further testing.

Runs on SLURM clusters with mySQL (for PASA), e.g. dayhoff RSB ANU. 

### What it does:


Annotates genes with [funannotate](https://github.com/nextgenusfs/funannotate), plus some other extra things for capturing more pathogen genes.

See dag.svg.

Transcriptome assembly and genome repeat-masking are NOT included, might be added soon.

### To-dos


- conda env handling
- TM detection in SignalP unions faa can be better

bla bla bla
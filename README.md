# A BRCAness phenotype without genetic or epigenetic inactivation of homologous recombination genes dictates response to chemotherapy in metastatic colorectal cancer

Snakemake pipeline to reproduce bioinformatics main figures panels.

Keep in mind that final figures were assembled manually (e.g. for
legend positioning) and that some rules produce more than one plot
(with/without legend) to ease
those manual fixes.

A Dockerfile with the same R version, R packages and overall linux
environment used to produce the Figures is available here:
https://github.com/vodkatad/snakemake_docker/blob/master/Dockerfiles/godot/Dockerfile

The relevant rules are found in dataset/figures/Snakefile.

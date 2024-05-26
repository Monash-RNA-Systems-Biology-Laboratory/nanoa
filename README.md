# nanoa

Tools for working with Nanopore reads in the Beilharz Lab. Under early development, usage may change.

Requires the [FASTA](https://github.com/wrpearson/fasta36) aligner to be installed, specifically the program "ssearch36" (tested with version 36.3.8i).

Sample data in `inst/extdata` is a selection of reads from experiment "2023-12-14-LB1-167".

You can install this package from GitHub using:

```{r}
install.packages("BiocManager")

BiocManager::install("Monash-RNA-Systems-Biology-Laboratory/nanoa")
```

See the "Usage" article for usage, and the reference pages for further optional arguments.

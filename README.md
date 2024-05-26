# nanoa

Tools for working with nanopore reads in the Beilharz Lab.

Requires the Fasta aligner to be installed, specifically the program "ssearch36".

Sample data in `inst/extdata` is a selection of reads from experiment "2023-12-14-LB1-167".



## Development

```{r}
devtools::document()
devtools::load_all(export_all=FALSE)

```
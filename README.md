# HiVE v1.10 (20260210)

High-dimensional lipid Visualization and Enrichment analysis (HiVE)

## Description

This package provides tools for comprehensive visualization and analysis of multiomics data characterizing lipids.    HiVE integrates underlying knowledge of curated lipid metabolic pathways for quick mapping and analysis of lipidomics data onto a predefined network structure.    Flexible data formats supported by HiVE also enable simultaneous or individual mapping of complementary proteomics or transcriptomics datasets containing lipid metabolizing enzymes or lipid-related genes.

## Getting Started

### Dependencies
* Windows 10-11, WSL Ubuntu 22.04 or higher, Linux Ubuntu 22.04 or higher, or macOS 12.7.1 or higher
* R version 4.5.0 or higher (https://cran.r-project.org/)
* (Optional) RStudio version 2023.06.2 or higher (https://posit.co/download/rstudio-desktop/)
* R-packages (downloaded from CRAN unless otherwise specified):
    * Suggests: 
        * knitr,
        * rmarkdown,
        * reticulate
    * Imports:
        * dplyr,
        * ggrepel,
        * ggplot2,
        * viridis,
        * ggsci,
        * RColorBrewer,
        * rcdk,
        * igraph,
        * ggraph,
        * sf,
        * concaveman

### Installation
* Run the following in a new R session on the command line or within R-Studio:

```
devtools::install_github(
  "cschasestevens/HiVE", 
  ref = "main", 
  build_vignettes = TRUE
)
```

## Help
* Browse vignettes by running the following:

```
browseVignettes("HiVE")
```

* Access function documentation by running the following:

```
# Type function name after package name
?HiVE::hive_base()
```

## Authors

* Nathanial Chase Stevens, PhD, University of North Carolina at Chapel Hill
* Email: Nathanial_Stevens@med.unc.edu
* Alternate email: cschasestevens@gmail.com
* LinkedIn: https://www.linkedin.com/in/nathanial-chase-stevens-phd-08775180/

## Version History
* 1.10
    * Updated Network structure for HiVE base and oxylipin pathways.
    * Added support for COMb S4 objects as plot inputs, which are now the default option for all network plot types.
    * Expression data from single-cell transcriptomics datasets are automatically filtered based on lipid-related genes present in HiVE edge data.
    * Plot functions now support single-cell transcriptomics data provided from Seurat objects
    * Updated tutorial to reflect changes in version 1.10
    * Bioconductor packages now automatically installed
* 1.03
    * Added HiVE enrichment network
* 1.02
    * Added example dataset
    * Added HiVE annotation network
* 1.01
    * Added HiVE base network
    * Integrated HiVE node and edge metadata within package
* 1.00
    * Initial Release

## License

This project is licensed under the GNU General Public License Version 3 - see the LICENSE.md file for details

## Acknowledgments

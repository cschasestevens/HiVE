# HiVE v1.11 (20260327)

High-dimensional lipid Visualization and Enrichment analysis (HiVE)

## Description

This package provides tools for comprehensive visualization and analysis of multiomics data characterizing lipids.    HiVE integrates underlying knowledge of curated lipid metabolic pathways for quick mapping and analysis of lipidomics data onto a predefined network structure.    Flexible data formats supported by HiVE also enable simultaneous or individual mapping of complementary proteomics or transcriptomics datasets containing lipid metabolizing enzymes or lipid-related genes.

This package is intended as an extension COMb (https://github.com/cschasestevens/COMb), which is used for analyzing untargeted metabolomics/lipidomics datasets. However, all functions can be run without using COMb for prior data analysis.

## Getting Started

### Dependencies
* Windows 10-11, WSL Ubuntu 22.04 or higher, Linux Ubuntu 22.04 or higher, or macOS 12.7.1 or higher
* R version 4.5.0 or higher (https://cran.r-project.org/)
* (Optional) RStudio version 2023.06.2 or higher (https://posit.co/download/rstudio-desktop/)
* R-packages (downloaded from CRAN unless otherwise specified):
    * Suggests: 
        * knitr,
        * rmarkdown,
        * reticulate,
        * Seurat
    * Imports:
        *   ggsci,
        *   viridis,
        *   RColorBrewer,
        *   ggplot2,
        *   igraph,
        *   ggraph,
        *   dplyr,
        *   sf,
        *   concaveman,
        *   ggrepel,
        *   tidyr,
        *   ggnewscale,
        *   shadowtext,
        *   SeuratObject,
        *   gtools,
        *   magrittr,
        *   reshape2,
        *   ggforce,
        *   rcdk

### Installation
* Run the following in a new R session on the command line or within R-Studio:

```
devtools::install_github(
  "cschasestevens/HiVE", 
  ref = "main", 
  build_vignettes = TRUE
)
```

* When installing spatial packages (e.g. sf), GDAL system libraries must be installed first. If running HiVE in Linux, install using the following command within the conda environment:

```
$ sudo apt-get install libgdal-dev gdal-bin libproj-dev libgeos-dev
```

* To fix potential error installing spatial dependencies (e.g. units), run the following within the conda environment:

```
$ sudo apt-get install libudunits2-dev
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
* 1.11
    * Updated hive_anno, hive_enrch, and hive_flux parameters for additional input data flexibility.
* 1.10
    * Added functions for generating individual plot components. This option is useful for preparing publication quality network plots.
    * Updated Network structure for HiVE base and oxylipin pathways.
    * Added function for plotting individual lipid/gene-level data.
    * Plot functions now support single-cell transcriptomics data.
    * Added support for plotting lipid flux comparisons.
    * Updated tutorial to reflect changes in version 1.10.
    * Bioconductor packages now automatically installed.
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

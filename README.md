# gnomAD-SV codebase

***Please note that this repository is currently being updated, and some code is missing***

This repository contains custom code and scripts developed as described in the gnomAD-SV preprint, [Collins\*, Brand\*, et al., bioRxiv (2019)](https://broad.io/gnomad_sv).  

For more information, please refer to [this blog post](https://broad.io/gnomad_sv).  

---

# Structure of this repository

The contents of this codebase are subdivided into several sections, for clarity.  

These sections are subdivided as follows:

### SV pipeline scripts

**Directory:** `gnomad_sv_pipeline_scripts`  

This subdirectory contains the individual scripts and custom code to generate SV callsets using the gnomAD-SV pipeline.  

By and large, these scripts need not be called individually; they have all been placed in the appropriate order and with the correct arguments & options within their corresponding WDLs (see below).

### SV pipeline WDLs

**Directory:** `gnomad_sv_pipeline_wdls`  

This subdirectory contains the WDLs used to generate SV callsets using the gnomAD-SV pipeline.  

If you are unfamiliar with WDL, please see the `Dependencies` section below for more information.  

### gnomAD-SV analysis scripts

**Directory:** `gnomad_sv_analysis_scripts`  

This subdirectory contains the individual scripts and custom code used to filter, perform quality-control on, and analyze the gnomAD-SV callset on FireCloud.

### gnomAD-SV analysis WDLs

**Directory:** `gnomad_sv_analysis_wdls`  

This subdirectory contains the WDLs used to filter, perform quality-control on, and analyze the gnomAD-SV callset on FireCloud.

### gnomAD-SV manuscript code

**Directory:** `gnomad_sv_manuscript_code`  

This subdirectory contains the individual scripts and custom code used to conduct the final analyses presented in the gnomAD-SV preprint, and to generate the graphs and other plots presented in the manuscript.  

Please note that due to certain datasets appearing under restricted access (such as individual-level genotype data), not every analysis or plot will be able to be reproduced using the code provided herein.   

### Ancillary data and reference files

**Directory:** `data_and_refs`  

This subdirectory contains miscellaneous data and reference files used throughout the gnomAD-SV callset generation, analysis, and manuscript preparation.  

---

### Important notes regarding the design & implementation of the gnomAD-SV discovery pipeline 

1. As described in the supplementary information provided with [the gnomAD-SV preprint](https://broad.io/gnomad_sv), the gnomAD-SV discovery pipeline has been written for cloud-based implementation via [FireCloud](https://software.broadinstitute.org/firecloud/) on large (>500-sample) cohorts of standard Illumina short-read whole-genome sequencing data. As such, this pipeline was not designed to be executed on local networks or computing clusters, on smaller cohorts (<500 samples), or on sequencing data types other than standard Illumina short-read whole-genome sequencing (e.g., this pipeline will not work for exome sequencing, PacBio/Nanopore long-read sequencing, etc.).   

2. Please note that we are in the process of optimizing and refactoring the code involved in the gnomAD SV discovery pipeline to make it more streamlined and efficient for application to other datasets. We anticipate this optimized version will be available soon. Once available, it will be linked here, and will supercede all code herein.

---

### Dependencies

  * Most of the components of this codebase are assigned to a Docker image, and in many cases image hashes are specified in the relevant WDLs that invoke each script.  

  * All code has been designed & implemented for FireCloud, a user-friendly interface for large-scale parallel computation and genomics projects in Google Cloud. For more information on FireCloud, please refer to [the FireCloud website](https://software.broadinstitute.org/firecloud/).

  * FireCloud uses the [Workflow Description Language (WDL)](https://software.broadinstitute.org/wdl/documentation/spec) and the [Cromwell Engine](https://cromwell.readthedocs.io/en/stable/) to interface directly with Google Cloud.  

  * Many operations in the gnomAD-SV pipeline rely on `svtk`, a python package that handles operations with SVs in VCF and BED formats. Source code and installation instructions for `svtk` are [provided via gitHub](https://github.com/talkowski-lab/svtk). More details on `svtk` can be found in [Werling et al., Nat. Genet. (2018)](https://www.nature.com/articles/s41588-018-0107-y).

---

# Other information

### Data Availability
To browse the gnomAD-SV callset, please refer to [the gnomAD Browser](https://gnomad.broadinstitute.org).  

The gnomAD-SV callset is available for download in VCF and BED format, with no restrictions on reuse or reanalysis. To download these files, please visit [the gnomAD website downloads page](https://gnomad.broadinstitute.org/downloads).

Consistent with all prior ExAC and gnomAD studies, the majority of raw sequencing data are available to approved investigators through various repositories such as dbGaP, etc. Access to the raw sequencing data is handled by the investigators of each contributing study, not through gnomAD directly. For details on sequencing data access, please refer to [the gnomAD flagship manuscript](https://broad.io/gnomad_flagship)

### Terms of Use
All code in this repository is released under the MIT license (see `LICENSE`).  

If you reuse the code hosted in this repo, or the SV data hosted on the gnomAD Browser, please cite [the gnomAD-SV preprint](https://broad.io/gnomad_sv).

### Citation
[**An open resource of structural variation for medical and population genetics.**](https://broad.io/gnomad_sv) Ryan L. Collins, Harrison Brand, Konrad J. Karczewski, Xuefang Zhao, Jessica Alföldi, Amit V. Khera, Laurent C. Francioli, Laura D. Gauthier, Harold Wang, Nicholas A. Watts, Matthew Solomonson, Anne O’Donnell-Luria, Alexander Baumann, Ruchi Munshi, Chelsea Lowther, Mark Walker, Christopher Whelan, Yongqing Huang, Ted Brookings, Ted Sharpe, Matthew R. Stone, Elise Valkanas, Jack Fu, Grace Tiao, Kristen M. Laricchia, Christine Stevens, Namrata Gupta, Lauren Margolin, The Genome Aggregation Database (gnomAD) Production Team, The gnomAD Consortium, John A. Spertus, Kent D. Taylor, Henry J. Lin, Stephen S. Rich, Wendy Post, Yii,Der Ida Chen, Jerome I. Rotter, Chad Nusbaum,†, Anthony Philippakis, Eric Lander, Stacey Gabriel, Benjamin M. Neale, Sekar Kathiresan, Mark J. Daly, Eric Banks, Daniel G. MacArthur, Michael E. Talkowski. _bioRxiv_ (2019).

### Credits
Acknowledgements & credits for all members of the gnomAD consoritum and gnomAD production, analysis, and SV teams are [listed on the gnomAD website](https://gnomad.broadinstitute.org/about).   

### Contact
If you have questions about the code hosted in this repo, please [email Ryan Collins](mailto:rlcollins@g.harvard.edu) or file an issue on GitHub.

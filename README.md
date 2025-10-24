# **Coastal ecosystem degradation driven by decades of unregulated terrestrial mining**

**This is a workflow to reproduce analysis conducted in Meyneng et al. *(Submitted in 2025)*.**

*Note that every figures have been finalised using a graphic software (AffinityDesigner) for visualisation details.*

Date: 07/10/2025

Authors: Mathisse Meyneng and Arthur Monjot

Dependencies: R; ITSx; cutadapt; vsearch; seqkit

First, clone github repository: `git clone https://github.com/Matmey/Meyneng_NC.git`

## **Rawdata availability**

### Amplicon sequencing targeting 18S-V4 rDNA

Amplicon sequencing targeting 18S-V4 rDNA come from ECOMINE project (https://www.cresica.nc/projet/ecomine) and are available as raw sequence (https://doi.org/10.12770/62cde5b9-1888-4b2c-b8a5-a490a696b078).

For this study these sequences has been processed using SAMBA v4 available here: https://gitlab.ifremer.fr/bioinfo/workflows/samba 

### Amplicon sequencing targeting ITS rDNA

Amplicon sequencing targeting ITS rDNA were performed by the sequencing platform of Tartu University and processed by NextITS pipeline (https://github.com/vmikk/NextITS.git). Raw sequencing data are available on demand (vladimir.mikryukov@ut.ee).

### Directory organisation

rawdata/ - Contains all raw data required for the analyses.

Main_results.r - R script for all main analyses and figure generation.

LongRead_analysis.md - Bash script for analysis of long-read data.

Process_LR_analysis_figure.Rmd - R script for network analysis of long read and short read.


    





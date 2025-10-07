# **Coastal ecosystem degradation driven by decades of unregulated terrestrial mining**

**This is a workflow to reproduce analysis conduced in Meyneng et al. (2025).**

Date: 08/08/2025

Authors: Mathisse Meyneng and Arthur Monjot

Dependencies: R; ITSx; networkx; seqkit

First, clone github repository: `git clone https://github.com/Matmey/Meyneng_NC.git`

##













## I. Install dependencies

```
conda install bioconda::itsx
conda install bioconda::seqkit
conda install bioconda::vsearch
```

## I. LR data

### I.A. Prepare rawdata

```
gunzip rawdata/OBAMA01_Boxes-4_FullSequences_with-OTU-IDs.fasta.gz
```

### I.B. Extract 18S with ITSx

```
mkdir -p Outputs/LongReads_analysis/region/
ITSx -i rawdata/OBAMA01_Boxes-4_FullSequences_with-OTU-IDs.fasta -o Outputs/LongReads_analysis/region/Thio_LR_extracted --preserve T --cpu 12 --save_regions all
cat Outputs/LongReads_analysis/region/Thio_LR_extracted.SSU.fasta | grep "^>" | wc -l
```

7055 SSU derived from LR sequences remain.
The result appear convenable

### I.C. Extract V4 with cutadapt and the SR reverse primer

The SR reverse primer: TAReukREV3

5"-ACTTTCGTTCTTGATYRA-3"
reverse complement: TYRATCAAGAACGAAAGT
"T[CT][AG]ATCAAGAACGAAAGT"

```
nohup cutadapt -a TYRATCAAGAACGAAAGT";min_overlap=18" -o Outputs/LongReads_analysis/region/Thio_LR_cutadapt_trimmed.V4.fasta --action retain -e 0.2 --discard-untrimmed Outputs/LongReads_analysis/region/Thio_LR_extracted.SSU.fasta > Outputs/LongReads_analysis/region/Thio_LR_cutadapt_trimmed.V4.log
cat Outputs/LongReads_analysis/region/Thio_LR_cutadapt_trimmed.V4.fasta | grep "^>" | wc -l
```

6841 sequences remain after the selection with cutadapt. This is pretty better than the solution using regular expression (i.e. 5725).

I.D. Check the mean length of the V4 trimmed from LR sequences

We use the seqkit tool:

```
seqkit stats Outputs/LongReads_analysis/region/Thio_LR_cutadapt_trimmed.V4.fasta > Outputs/LongReads_analysis/region/Thio_LR_cutadapt_trimmed.V4.stat
```

file                                                                | format | type | num_seqs |   sum_len | min_len | avg_len | max_len
Outputs/LongReads_analysis/region/Thio_LR_cutadapt_trimmed.V4.fasta | FASTA  | DNA  |    6,841 | 2,642,458 |     236 |   386.3 |     555


## II. SR data

### II.A. Prepare rawdata

We have a csv table with in column 1 the header and in column 2 the sequence. We have to transform this file in fasta file.
```
cat rawdata/database_18SV4_Thio.csv | awk -F";" '{print ">"$1"\n"$2}' | tail -n+3 > Outputs/LongReads_analysis/region/Thio_SR_OTUs.V4.fasta
cat Outputs/LongReads_analysis/region/Thio_SR_OTUs.V4.fasta | grep ">" | wc -l
```

They are 25250 SR OTUs sequences.

II.B. Extract the same region (i.e. V4) in SR

The LR forward primer: Euk575Fngs

5"-ASCYGYGGTAAYWCCAGC-3"
"A[CG]C[CT]G[CT]GGTAA[CT][AT]CCAGC"

The SR forward primer: TAReukFWD1

5"-"CCAGCASCYGCGGTAAT-3"
"CCAGCA[CG]C[CT]GCGGTAAT""

Since bot primers are very closed in their position on the SSU, we choose to not trim the SR

```
#nohup cutadapt -g ASCYGYGGTAAYWCCAGC -o Outputs/LongReads_analysis/region/Thio_SR_cutadapt_trimmed.V4.fasta --action retain -e 0.2 Outputs/LongReads_analysis/region/Thio_SR_OTUs.V4.fasta > Outputs/LongReads_analysis/region/Thio_SR_OTUs_cutadapt_trimmed.V4.log
```

II.C. Check the mean length of the V4 trimmed from SR sequences

We use seqkit tool:

```
seqkit stats Outputs/LongReads_analysis/region/Thio_SR_OTUs.V4.fasta > Outputs/LongReads_analysis/region/Thio_SR_OTUs.V4.stat
```

file                                                    | format | type | num_seqs |   sum_len | min_len | avg_len | max_len
Outputs/LongReads_analysis/region/Thio_SR_OTUs.V4.fasta | FASTA  | DNA  |   25,250 | 8,573,229 |     232 |   339.5 |     451


## III. Clustering LR and SR derived V4 region and SSN-based analysis

### III.A. Clustering

We use Vsearch to align all LR and SR-derived region.
We check, in the same time, if the query and subject sequence length is up to 80% in coverage.

```
mkdir Outputs/LongReads_analysis/Clustering/
cat Outputs/LongReads_analysis/region/Thio_SR_OTUs.V4.fasta Outputs/LongReads_analysis/region/Thio_LR_cutadapt_trimmed.V4.fasta > Outputs/LongReads_analysis/Clustering/All_V4_sequences.fasta
vsearch --allpairs_global Outputs/LongReads_analysis/Clustering/All_V4_sequences.fasta --id 0.8 --query_cov 0.8 --target_cov 0.8 --threads 12 --blast6out Outputs/LongReads_analysis/Clustering/All_V4_allvsall_cover80.tsv
```

It takes a while!

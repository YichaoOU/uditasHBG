UDiTaS v1.1
===========

### version 1.1 (include some key modifications)

Specifically design for analyzing g34 5kb deletion in HBG1 and HBG2 regions

multi-mapped reads were randomly assigned to a location if score is the same, this creates an over-estimation of large indels

### Changes (not foundamental)

`uditas.py` and `uditas_helper.py` are modified to solve the following issues:

1. fastq file names were harcoded.
2. sample_into.csv contained empty rows.
3. bowtie2 index file incorrect

### Installation (tested at St. Jude)

1. Download hg38.2bit

`wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit`

2. Download bowtie2 index

```

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz

tar zxvf GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz

mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.1.bt2 hg38.1.bt2
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.2.bt2 hg38.2.bt2
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.3.bt2 hg38.3.bt2
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.4.bt2 hg38.4.bt2
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.rev.1.bt2 hg38.rev.1.bt2
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.rev.2.bt2 hg38.rev.2.bt2

chmod a+rx *

```

3. add ENV

Note that the environment variable is the folder that containing 2bit or index files.

```
export BOWTIE2_INDEXES=/home/yli11/Data/Human/hg38/bowtie2/
export GENOMES_2BIT=/home/yli11/Data/Human/hg38/

```

4. install udias dependencies

```
git clone https://github.com/YichaoOU/uditas.git

cd uditas

module load conda3

conda env create -f uditas_env.yml

source activate uditas_env

conda install -c anaconda python-dateutil

conda install -c conda-forge tbb

python -m pip install matplotlib

conda install -c conda-forge icu

conda install -c anaconda qt

conda install -c anaconda pyqt

python setup.py install

module load bowtie2/2.2.9

pytest


```

The old `uditas_env.yml` has quite some dependency errors. I have modified the yaml file. Bowtie2 is also not installed successfully through conda, luckly our HPC has bowtie2 installed, so I don't need to install it again.

You should be able to see this message if test is passed.

```
(uditas_env) [yli11@nodecn231 uditas]$ pytest
========================================================= test session starts =========================================================
platform linux2 -- Python 2.7.13, pytest-3.2.1, py-1.8.1, pluggy-0.4.0
rootdir: /research/rgs01/home/clusterHome/yli11/test/uditas, inifile:
collected 3 items                                                                                                                      

test/test_fig2a.py .
test/test_fig2b.py .
test/test_fig2c.py .

===================================================== 3 passed in 540.71 seconds =====================================================

```

Overview
--------

UDiTaS(TM) stands for UniDirectional Targeted Sequencing, a novel sequencing method useful for measuring small indels as well as
structural rearrangements, like translocations, in a single reaction.

See details of the method in Giannoukos et al. BMC Genomics (2018) 19:212, https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4561-9


Systems Requirements
--------------------

UDiTaS has been tested in python 2.7.13 and requires the python packages and tools listed in the file uditas_env.yml

The code requires setting up two environmental variables

`BOWTIE2_INDEXES` contains the location of the bowtie2 index files, typically hg38.1.bt2, hg38.2.bt2, etc.

`GENOMES_2BIT` contains the location of the 2bit files for the genomes used in the analysis, eg hg38.2bit

To test the code create a virtual python environment with

`conda env create -f uditas_env.yml`

then activate using

`source activate uditas_env`

To install uditas as an executable run

`python setup.py install`

To test the installation run

`pytest`

This will process the test data in

```
data/fig2a
data/fig2b
data/fig2c
```

These data are a subsample of the data displayed in Fig 2 of the paper.

Usage
-----
`uditas` is the main command to launch the UDiTaS analysis. The required argument is: `dir_sample`, the path of the directory with the data to be analyzed.

`dir_sample` should contain the fastq.gz files for R1, R2, I1 and I2. Used in the demultiplexing step.

`dir_sample` should also contain the analysis sheet `sample_info.csv` containing the description of all the samples, with their barcodes and guides used. See examples in the folder `data`

Once the setup has been run, the code can be run as

`uditas ./data/fig2a`

The full list of options are:

```
usage: uditas [-h] [-folder_genome_2bit FOLDER_GENOME_2BIT]
              [-skip_demultiplexing SKIP_DEMULTIPLEXING]
              [-skip_trimming SKIP_TRIMMING]
              [-skip_genome_local_alignment SKIP_GENOME_LOCAL_ALIGNMENT]
              [-skip_genome_global_alignment SKIP_GENOME_GLOBAL_ALIGNMENT]
              [-process_amplicon PROCESS_AMPLICON]
              [-skip_amplicon_global_alignment SKIP_AMPLICON_GLOBAL_ALIGNMENT]
              [-check_plasmid_insertions CHECK_PLASMID_INSERTIONS]
              [-skip_plasmid_alignment SKIP_PLASMID_ALIGNMENT] [-ncpu NCPU]
              [-window_size WINDOW_SIZE]
              [-default_amplicon_window_around_cut DEFAULT_AMPLICON_WINDOW_AROUND_CUT]
              [-min_MAPQ MIN_MAPQ] [-min_AS MIN_AS]
              [-process_AMP_seq_run PROCESS_AMP_SEQ_RUN]
              dir_sample

Process UDiTaS data

positional arguments:
  dir_sample            Directory with the sample to be processed

optional arguments:
  -h, --help            show this help message and exit
  -folder_genome_2bit FOLDER_GENOME_2BIT
                        Folder containing the 2bit file(s) with the reference
                        genome being used (default: GENOMES_2BIT)
  -skip_demultiplexing SKIP_DEMULTIPLEXING
                        Skip demultiplexing? Options: 0, 1 (skip) (default: 0)
  -skip_trimming SKIP_TRIMMING
                        Skip adapter trimming? Options: 0, 1 (skip) (default:
                        0)
  -skip_genome_local_alignment SKIP_GENOME_LOCAL_ALIGNMENT
                        Skip genome-wide local alignment? Options: 0 , 1
                        (skip) (default: 1)
  -skip_genome_global_alignment SKIP_GENOME_GLOBAL_ALIGNMENT
                        Skip genome-wide global alignment? Options: 0 , 1
                        (skip) (default: 0)
  -process_amplicon PROCESS_AMPLICON
                        Select row number (0-based) of amplicon to process,
                        set to all to process all amplicons (default: all)
  -skip_amplicon_global_alignment SKIP_AMPLICON_GLOBAL_ALIGNMENT
                        Skip amplicon global alignment? Options: 0, 1 (skip)
                        (default: 0)
  -check_plasmid_insertions CHECK_PLASMID_INSERTIONS
                        Check for plasmid insertions. Options: 0 (skip), 1
                        plamid_name and plasmid_sequence required in
                        sample_info.csv (default: 1)
  -skip_plasmid_alignment SKIP_PLASMID_ALIGNMENT
                        Skip plasmid alignment? Note, just alignment. Counts
                        still evaluated. Options: 0, 1 (skip) (default: 0)
  -ncpu NCPU            Number of CPUs to use (default: 4)
  -window_size WINDOW_SIZE
                        Window size around cut sites used to grab UDiTaS reads
                        (default: 15)
  -default_amplicon_window_around_cut DEFAULT_AMPLICON_WINDOW_AROUND_CUT
                        Window size around cut sites used to create amplicons
                        (default: 1000)
  -min_MAPQ MIN_MAPQ    Minimum mapping quality to include a read (default: 5)
  -min_AS MIN_AS        Minimum alignment score to include a read (default:
                        -180)
  -process_AMP_seq_run PROCESS_AMP_SEQ_RUN
                        Set to 1 to process an AMP-seq run using GUIDE-seq
                        adapters (default: 0)
```

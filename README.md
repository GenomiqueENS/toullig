Toullig 0.2-alpha
=======

Toullig is (since January 17 2017) a reader, parser of Fast5 files of ONT (basecalled and not basecalled) and a fastq files writer and the first of march Toullig can trim fastq of ONT based on the RT adaptors.


# Author


Aurélien Birer, [birer@biologie.ens.fr](birer@biologie.ens.fr)



# Table of contents


* [Requirements](#requirements)
* [Installation](#installation)
    * [To install Maven](#to-install-maven)
    * [To install Toullig](#to-install-toullig)
* [General options](#general-options)
* [How it works](#how-it-works)
    * [Classification MinION run](#classification-minion-run)
    * [Chemistry available](#chemistry-available)
* [Fast5tofastq](#fast5tofastq)
    * [Understand the type of sequence](#understand-the-type-of-sequence)
    * [Sequence available for each run configuration](#sequence-available-for-each-run-configuration)
    * [Options Fast5tofastq](#options-fast5tofastq)
    * [Example](#example)
* [TrimFastq](#trimFastq)
    * [Options TrimFastq](#options-trimFastq)
    * [Example trim](#example-trim)
* [Developpment environnement](#developpment-environnement)
* [Repository](#repository)
* [License](#license)


# Requirements


You just need to have java 8 and maven installed on your computer. This alpha version work on Ubuntu (Unix distribution).


# Installation


### To install Maven


    sudo apt-get install maven


### To install Toullig



    git clone https://github.com/GenomicParisCentre/toullig.git
    cd toullig
    mvn clean install


# General options


    #Information

    -help | -h      #display help
    -version        #display version of toullig
    -about          #display information of toullig
    -license        #display license of toullig
    


# How it works


Toullig have 2 tools :

- fast5tofastq : read the rootDirectory/downloads of your Fast5 run minION after the step of basecalling (metrichor/albacore).
- trim : trim the reads of a ONT fastq with a sam file, based on the RT adaptors.

### Classification MinION run



>**Not Barcoded Directories Tree**

├── downloads <br>
│   ├── fail <br>
│   └── pass <br>
└── uploaded <br>

>**Barcoded Directories Tree (with 6 barcode)**

├── downloads <br>
│   ├── fail <br>
│   │   └── unclassified <br>
│   └── pass <br>
│       ├── BC01 <br>
│       ├── BC02 <br>
│       ├── BC03 <br>
│       ├── BC04 <br>
│       ├── BC05 <br>
│       └── BC06 <br>
└── uploaded <br>


### Chemistry available


<p align="center">
  <img src="images/chemistry_tested.png"/>
</p>


# Fast5tofastq


In the execution of toullig Fast5tofastq, the programm step :

 + List the .fast5 files.
 + Read a .fast5 file.
 + Write the .fastq sequence(s).
 + Make few some log informations.

### Understand the type of sequence

Actually, we use in developpement Metrichor for the basecalling of our .fast5 files.

But it's important to understand clearly the type of the 4 fastq sequences give by metrichor (in our case in 2D).


The template sequence is the first sequence basecalled, this sequence correspond to the own read sequenced in 1D. This sequence contains section as follow :

<p align="center">
  <img src="images/template_sequence.png"/>
</p>

The complement sequence is the second sequence basecalled (in 2D), this sequence correspond to the reverse of the template sequence. This sequence contains section as follow :

<p align="center">
  <img src="images/complement_sequence.png"/>
</p>

/!\ The sequence template and complement can be reverse (it's depends of how the leader-adaptor is fix).

The consensus sequence is the sequence result of the alignement of the template and the complement sequence (in 2D). This sequence contains section as follow :

<p align="center">
  <img src="images/consensus_sequence.png"/>
</p>

/!\ The consensus sequence can be bases on the template sequence or complement sequence.

The transcript sequence is the sequence result of the consensus sequence (in 2D) with a trim of the barcode sequence (/!\ barcode can be still). This sequence contains section as follow :

<p align="center">
  <img src="images/transcript_sequence.png"/>
</p>


/!\ The transcript sequence is only available for barcoded ONT run.

### Sequence available for each run configuration

Here, for each sequecing configuration the sequence available for the fast5tofastq conversion.

In bold, the type of sequencing that it's mostly interessant.

<p align="center">
  <img src="images/type_of_sequencing_workflow.png"/>
</p>


### Options Fast5tofastq


    
    #Options
    
    -status pass|fail|unclassified (default: pass)                                  # The status of fast5 file
    -type template|complement|consensus|transcript (default: transcript)            # The type of sequence
    -mergeSequence true|false (default: false)                                      # If you want merge all type of sequence whatever the status
    -compress GZIP|BZIP2 (default: none)                                            # Set the type of compression for the output fastq
    
    #Arguments

    -rootDirectoryFast5run /home/user/yourRootDirectoryFast5run
    -outputDirectoryFastq /home/user/yourOutputDirectoryFastq
    
    
### Example


I have a directory of a minION run in 2D barcoded.
If i want just get the fastq sequence of the 'template', the 'complement' and the 'consensus' for the fast5 files in the status/repertory 'fail'.


    bash ./target/dist/toullig-0.2-alpha/toullig.sh fast5tofastq -status fail -type template,complement,consensus /home/user/myRootDirectoryFast5run /home/user/myOutputDirectoryFastq
    
    I have a directory of a minION run in 2D not barcoded.
If i want just get the fastq default sequence for the fast5 files in the status/repertory 'pass'.


    bash ./target/dist/toullig-0.2-alpha/toullig.sh fast5tofastq -status pass /home/user/myRootDirectoryFast5run /home/user/myOutputDirectoryFastq

WARNING : This last example is a trap ! The default type of sequence is transcript or the run is not barcoded and the transcript sequence is not available on nont barcoded run.

# TrimFastq

In the case of cDNA protocols like SQK-MAP006.
One of the problem of the minION reads in the format fastq is that the read is not the transcript as we expected. The read still have the RT adaptor and, in some case, the barcode with adaptors (leader/hairpin).

For enhance the mapping quality, it's important to trim the reads to delete these unwanted sequences.

In the execution of toullig Trim, the programm step :

 + Repere the index between the outliers and the mRNA (mode P or SW).
 
 For cutadapt:
 
 + Write Fasta File for each Outlier.
 + Trim with cutadapt
 + Write the cutadapt output merge.
 
  For trimmomatic:
 
 + Read a sequence fastq.
 + Trim with trimmomatic
 + Write the sequence trimmed.

### Options TrimFastq

    
    #Options
    
    -trimmer cutadapt|trimmomatic (default: cutadapt)       # The trimmer tool use for trimming
    -mode P | SW (default: P)                               # The type of trimming the transcripts reads
    -stats true|flase (default: false)                      # If you want somes stats on the trimming
    
    #Options Trimming by Side-window mode
    
    - addIndexOutlier (default: 15)                         # A addition of bases of the outlier to have a better catch of the adaptor
    
    #Options Trimming by Side-window mode
    
    -thresholdSideWindow (default: 0.8)                             # The threshold for the Side-Window algorithm
    -lengthWindowsSW (default: 15)                          # The length for the the Side-Window algorithm
    
    #Options Cutadapt
    
    -errorRateCutadapt (default: 0.5)                       # The error rate for Cutadapt (mismatch + deletion)
    
    #Options Trimmomatic
    
    -seedMismatchesTrimmomatic (default:17 )                # The number base of mismatchs maximun for Trimmomatic
    -palindromeClipThresholdTrimmomatic (default: 30)       # The threshold of palindrome clip for Trimmomatic
    -simpleClipThreshold (default: 7)                       # The threshold of simple clip for Trimmomatic
    
    #Options Post-Process
    
    -minlen (default : 100)                                 # The threshold of minimum length to write trimmed fastq
    
    #Arguments
    
    -samFile            /home/user/yourSamFile
    -fastqFile          /home/user/yourFastqFile
    -fastqOutputFile    /home/user/yourFastqTrimOutput
    -adaptorFile        /home/user/yourAdaptorFile
    -workDir            /home/user/yourTmpRepertoryOfWork


### Example trim


    bash ./target/dist/toullig-0.2-alpha/toullig.sh trim /home/user/samFile.sam /home/user/fastqONTFile.fastq /home/user/myFastqTrim.fastq ~/toullig/config_files/adaptor_RT_sequence_modify_for_nanopore.txt /home/user/yourTmpRepertoryOfWork


# Developpment environnement


Ubuntu version: '16.04.4'

Java version: '1.8.0_121'

Maven version: '3.2.3'


# Repository

Currently the Git reference repository is [https://github.com/GenomicParisCentre/toullig](https://github.com/GenomicParisCentre/toullig).


# License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
[CeCILL](http://www.cecill.info/licences/Licence_CeCILL_V1.1-US.html)



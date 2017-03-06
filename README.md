Toullig 0.1-alpha-2
=======

Toullig is (since January 17 2017) a reader, parser of Fast5 files of ONT (basecalled and not basecalled) and a fastq files writer.


AUTHOR
-

Aurélien Birer, [birer@biologie.ens.fr](birer@biologie.ens.fr)


REQUIREMENTS
-

You just need to have java 8 and maven installed on your computer. This alpha version work on Ubutun (Unix distribution).


###TO INSTALL MAVEN



    sudo apt-get install maven


###TO INSTALL Toullig



    git clone https://github.com/GenomicParisCentre/Toullig.git
    cd Toullig


HOW IT'S WORK
-

Nanoporetools read the rootDirectory/downloads of your Fast5 run minION after the step of basecalling.

>**Not Barcoded Directories Tree**

├── downloads <br>
│   ├── fail <br>
│   └── pass <br>
└── uploaded <br>

>**Not Barcoded Directories Tree**

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

In the execution of nanoporetools, the programm step :

 + List the fast5 files.
 + Read a fast5 file.
 + Write the fastq sequence(s).
 + Make few some log informations.

###UNDERSTAND THE TYPE OF SEQUENCE

Actually, we use in developpement Metrichor for the basecalling of our .fast5 file.

But it's important to understand clearly the type of the 4 fastq sequences give by metrichor (in our case in 2D).


The template sequence is the first sequence basecalled, this sequence correspond to the own read sequenced in 1D. This sequence contains section as follow :

leader-adaptor|barcodePOS0|Switch-PCR-RT-Primer|GGG|transcript|NVTTT{n}|PCR-RT-Primer|barcodePOS1|hairpin-adaptor

The complement sequence is the second sequence basecalled (in 2D), this sequence correspond to the reverse of the template sequence. This sequence contains section as follow :

leader-adaptor|barcodePOS0|Switch-PCR-RT-Primer|CCC|transcript|NBAAAAA{n}|PCR-RT-Primer|barcodePOS1|hairpin-adaptor

The consensus sequence is the sequence result of the alignement of the template and the complement sequence (in 2D). This sequence contains section as follow :

leader-adaptor|barcodePOS0|Switch-PCR-RT-Primer|GGG|transcript|NVTTT{n}|PCR-RT-Primer|barcodePOS1|hairpin-adaptor

The transcript sequence is the sequence result of the alignement of the template and the complement sequence (in 2D) trim of the barcode sequence. This sequence contains section as follow :

Switch-PCR-RT-Primer|GGG|transcript|NVTTT{n}|PCR-RT-Primer

OPTIONS
-

    #Information
    -help | -h      #display help
    -version        #display version of nanoporetools
    -about          #display information of nanoporetools
    -license        #display license of nanoporetools
    
    #Options
    -status pass|fail|unclassified|passbarcode (default : pass)                  #The status of fast5 file
    -type template|complement|consensus|transcript (default : transcript)   #The type of sequence
    -mergeSequence true|false (default : false)                                 #If you want merge all type of sequence whatever the status
    
    #Arguments
    -rootDirectoryFast5run /home/user/yourRootDirectoryFast5run
    -outputDirectoryFastq /home/user/yourOutputDirectoryFastq
    
    
###Example


I have a directory of a minION run in 2D with barcode.
If i want just get the fastq sequence of the 'template', the 'complement' and the 'consensus' for the fast5 files in the status/repertory 'fail'.


    bash ./target/dist/toullig-0.1-alpha-2/toullig.sh -status fail -type template,complement,consensus /home/user/myRootDirectoryFast5run /home/user/myOutputDirectoryFastq



DEVELOPPEMENT ENVIRONNEMENT
-

Ubutun : '16.04.4'

JAVA version : '1.8.0_121'

Maven version : '3.2.3'


REPOSITORY
-

Currently the Git reference repository is [https://github.com/GenomicParisCentre/Toullig](https://github.com/GenomicParisCentre/Toullig).

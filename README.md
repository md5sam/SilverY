# SilverY

SilverY is a tool for shortlisting regions that are shared across Y chromosomes of different species. Specifically, it can be used for isolating homologous regions for graph genome construction. Branching points in the graph genome can be added based on regions shared between most, and not all species. 

### Usage  

    python silverY.py
  	
Important parameters for the user to choose are : 


**kmer-size** : 
- the size of k used while iterating through every read 
- this must be the same as DSK's kmer-size
- usually optimal in the range [25, 31] but can be determined by running kmerGenie

**RMasking** :
- it is recommended that the Y chromosomes be RepeatMasked prior to computing set intersections, as this may artificially inflate the kmer proportions shared


### Installation 

	git clone https://github.com/md5sam/SilverY.git
	cd SilverY


### Dependencies     

Numpy and Biopython

    pip install numpy
    pip install biopython
    

### Input

The following input files are required in ./data folder
    	
	
	r1.fastq : Enriched raw reads (first in pair) 
	r2.fastq : Enriched raw reads (second in pair) 
	kmers_from_reads : kmer counts from DSK for r1.fastq
	trusted_kmers : kmer counts from DSK for human Y single copy genes

The input folder and filenames can be changed by the user within the program. 


### Output 

The ./op_dir folder contains kmers common to each set of comparisons :

 	species1_species2_species3
  species1_species2
  species1_species3
  species2_species3
  species1
  species2
  species3
  



### License
This program is released under the MIT License. Please see LICENSE.md for details


### Citation
Please cite this Github repository if you use this tool in your research. Thanks !
https://github.com/md5sam/SilverY

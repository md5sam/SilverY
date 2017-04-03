# SilverY

SilverY is a tool for shortlisting regions that are shared across Y chromosomes of different species. Specifically, it can be used for isolating homologous regions for graph genome construction. Branching points in the graph genome can be added based on regions shared between most, and not all species. Further, exact proportions of k-mer intersections can also be used to find a lower bound on estimated time of divergence between species or sub-species. 

### Usage  

    python silverY.py
  	
Important parameters for the user to choose are : 


**kmer-size** : 
- the size of k used while iterating through every read 
- this must be the same as DSK's kmer-size
- usually optimal in the range [25, 31] but can be determined by running kmerGenie

**RMasking** :
- it is recommended that the Y chromosomes be RepeatMasked prior to computing set intersections, as the presence of shared repeats may artificially inflate the kmer proportions shared


### Installation 

	git clone https://github.com/md5sam/SilverY.git
	cd SilverY


### Dependencies     

Numpy and Biopython

    pip install numpy
    pip install biopython
    

### Input

The following input files are required in the data folder (assuming a 3-way comparison) : 
    	
	uniq_kmers_species1_Y : kmer counts from DSK from Y chromosome of species1
	uniq_kmers_species2_Y : kmer counts from DSK from Y chromosome of species2
	uniq_kmers_species3_Y : kmer counts from DSK from Y chromosome of species3

The input folder and filenames can be changed by the user within the program. 


### Output 

The output folder contains kmers common to each pair or triplet of comparisons (assuming a 3-way comparison) :

    species1_species2_species3
    species1_species2
    species1_species3
    species2_species3
    



### License
This program is released under the MIT License. Please see LICENSE.md for details


### Citation
Please cite this Github repository if you use this tool in your research. Thanks !
https://github.com/md5sam/SilverY

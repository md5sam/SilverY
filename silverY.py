
# coding: utf-8

# In[1245]:

from __future__ import division
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import os
import itertools


# In[1246]:

def reverse_complement(seq):
    seq_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return "".join([seq_dict[base] for base in reversed(seq)])


# In[1247]:

def kmerize (ip_string, kmer_size) :
    return [ip_string[i:i+kmer_size] for i in range(0, len(ip_string)-kmer_size+1, 1)]


# In[1248]:

# function that uses BioPython to convert a fasta file to a dict
def parse_ctgs_fasta (ip_fasta_file) :
    with open(ip_fasta_file, "r") as handle :
        contigs_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    return contigs_dict


# In[1249]:

def make_set_from_kmer_abundance (ip_file, kmer_size) :
    ''' THIS FUNCTION HANDLES REVERSE COMPLEMENTS AS WELL'''
    with open(ip_file,'r') as file_handle1 :
        list1 = [line[:kmer_size] for line in file_handle1]
    with open(ip_file,'r') as file_handle2 :
        list2 = [reverse_complement(line[:kmer_size]) for line in file_handle2]   
    list1.extend(list2)
    #print set(list1)
    return set (list1)


# In[1250]:

possible_species = ['hum', 'chimp', 'gor', 'bono', 'orang']
def identify_species(kmer_file):
    for species in possible_species :
        if species in os.path.basename(kmer_file) :
            return species


# In[1251]:

species_kmers_dict = defaultdict(set)
def make_species_dict(source_folder, kmer_files_by_species, k_size) :
    for kmer_file in kmer_files_by_species : 
        global species_kmers_dict
        species = identify_species(kmer_file)
        # do you want to have a sentinel here for non-empty species only?
        # if species_kmers_dict[species] not None :
        full_kmer_file_path = source_folder+kmer_file
        species_kmers_dict[species] = make_set_from_kmer_abundance(full_kmer_file_path , k_size)
    return species_kmers_dict


# In[1252]:

def intersect_kmers(*species) :
    '''ONLY EVER GIVE SPECIES NAMES TO THIS, NOT THE KMER FILE NAME'''
    # First convert the ip_tuple into a list with correct # of elements
    species_list = list(itertools.chain(*species))
    
    # sort this species_list for consistency
    # i.e. human_chimp_gorilla should be same as chimp_gorilla_human
    sorted_species_list = sorted(species_list)
   
    # create a new name for the species you are intersecting
    name = '_'.join(sorted_species_list)
    
    to_intersect = []
    # to_intersect will soon become a list of sets 
    for specie in species_list : 
        for k, v in species_kmers_dict.iteritems() :
            if specie == k :
                to_intersect.append(v)
                
    # use list unpacking to find intersection of a list of sets 
    species_kmers_dict[name] = set.intersection(*to_intersect)
    return name, species_kmers_dict[name]
   


# In[1253]:

def main () :
    k_size = 31 
    kmer_files = []
    folder_to_look_in = "./data/"
    op_dir = "./op_dir"
    
    print "Started main" 

    # os.walk prints dirpath, dirnames, filenames
    # get the kmer_files you're interested in
    for _, _, files in os.walk(folder_to_look_in) :
        kmer_files.extend(files)
    
    # start using only the species short_name henceforth to avoid confusion
    # i.e. orang_kmers, orang_file, orang_kmerfile etc. becomes just 'orang'
    species_shortnames = []
    for kmer_file in kmer_files :
        species_shortnames.append(identify_species(kmer_file))
    
    # by default we do an all-v-all comparison 2,3,4,5 at a time
    # comparisons will look like : ['bono', 'chimp'], ['bono', 'chimp', 'gor'], etc. 
    comparisons = []
    for i in range (2, len(species_shortnames)+1) :
        curr_comparisons = itertools.combinations(species_shortnames,i)
        for group in curr_comparisons : 
            comparisons.append(group)
    
    print "Generated comparison names"
    for group in comparisons :
	print group

    # by default, for every species we have a kmer_file for, we make a dict
    # dict where key is species_shortname and value is all its kmers in kmerfile
    make_species_dict(folder_to_look_in, kmer_files, k_size)
    
    print "Made species dicts for all" 
    
    # make a new directory to write o/p files
    if not os.path.exists(op_dir):
        os.makedirs(op_dir)
    
    # now we compare every group in comparisons
    for group in comparisons :
        list_group = list(group)
        file_name, file_contents = intersect_kmers(list_group)
        op_name = op_dir+'/'+file_name
        with open(op_name,"w") as fp_out : 
            if file_contents : 
                for element in file_contents : 
                    fp_out.write(str(element))
                    fp_out.write('\n')
        print "Generated kmer-intersection for group : " , file_name
	print "Number of common kmers (including rev_comp) is : ", len(file_contents)

# In[1254]:

if __name__ == "__main__": main()


# In[ ]:




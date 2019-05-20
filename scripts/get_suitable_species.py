import subprocess
from Bio import SeqIO
import os

working_dir = os.environ['HOME'] + "/Desktop/projects/HB_reborn/"
eutils_dir = os.environ['HOME'] + "/edirect/"


##########################################################################
# CREATE DATABASE DICT WITH ALL 28S WITH SPECIFIC LENGTH
sequences_file = SeqIO.parse(working_dir + "Eukaryota28Ssequences.fasta", "fasta")
database = {}

min_28S_length = 2000
min_sra_size = 1000
max_sra_size = 4000

print("\nReading input file")
print("\n\nCreating database with all 28S sequences longer than " + str(min_28S_length) + "bp...")

counter = 0 
for entry in sequences_file:
    counter +=1
    #FIND SPECIES NAME FROM FASTA HEADER
    description_split = entry.description.split(";")
    species_container = description_split[-1].split(" ")
    if len(species_container) > 1:
        current_species = species_container[0] + " " + species_container[1]
        #ADD TO DATABASE IF SEQUENCE SATISFIES LENGTH CONDITION
        if len(entry.seq) >= min_28S_length:
            database[current_species] = [entry.id, str(entry.seq), description_split]

print("Sequences checked: " + str(counter))
print("Species with 28S longer than " + str(min_28S_length) + ": " + 
      str(len(database)) + "\n\n\n")
print("Will now check all " + str(len(database)) + " species for suitable SRA datasets\n")
print("Only SRA datasets smaller than " + str(max_sra_size) + "MB and larger than " + str(min_sra_size) + "MB will be kept.")



###########################################################################
# SEARCH NCBI FOR SRA DATASET FROM THE SPECIES FROM PREVIOUS STEP
# SAVE THE SPECIES THAT HAVE SUITABLE 28S AND SRA DATASET IN A NEW FILE CALLED SUITABLE_SPECIES.TSV

suitable_species = []

counter2 = 0
for species in database:
   counter2 += 1
   species_stripped = species.strip()
   
   result = str(subprocess.check_output(eutils_dir + 
                                        'esearch -db sra -query "' + 
                                        species_stripped + 
                                        '[ORGN] AND PAIRED[LAY] AND biomol rna[PROP]" | ' + eutils_dir + "efetch --format runinfo | cut -d ',' -f 1 | grep SRR | head -1",
                                        shell=True))

   print("\nCurrently checking: " + species_stripped + 
         " (" + str(counter2) + "/" + str(len(database)) + ")")
   
   if result[2] == "S":
       try:
           size = str(subprocess.check_output(eutils_dir + 'epost -db sra -id ' + 
                                              result[2:-3] + 
                                              ' -format acc | ' + eutils_dir + 'esummary -format runinfo -mode xml | ' + eutils_dir + 'xtract -pattern Row -element size_MB',
                                              shell=True))
           size_int = int(size[2:-3])
           
           print("Found SRA dataset for " + species_stripped + ": " 
                 + result[2:-3] + ", size: " + str(int(size[2:-3])) + "MB")
           
           if size_int > max_sra_size:
               print("SRA dataset too large.")
               
           elif size_int < min_sra_size:
               print("SRA dataset too small.")
               
           else:
               print("Saving data...")
               suitable_species.append([species_stripped, database[species_stripped][0], database[species_stripped][1], result[2:-3], int(size[2:-3]), database[species_stripped][2]])
       
       except ValueError:
           continue
       
   else:
       print("No SRA dataset found")
       
   result = ""
   
output_tsv = open(working_dir + "suitable_species.tsv", "w")

for entry in suitable_species:
    
    for i in range(len(entry)):
        
        output_tsv.write(str(entry[i]) + "\t")
    output_tsv.write("\n")
    

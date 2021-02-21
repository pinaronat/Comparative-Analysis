import math
#uses FASTA type alingment file
#gets the residues in given site_no for all sequences in the file
#WARNING: the function will not work if site_no is bigger than the alignment length
def get_site(align_File, site_no):
    alignment_file = open(align_File, "r")
    lines = alignment_file.readlines()
    residues = ""
    count = 0
    for line in lines:
        line = line.replace("\n","")
        #get which line the sequence is on in one sequence
        which_line = math.ceil(site_no/50)
        #if the first element is >, we start a new sequence. Reset the count
        if line[0] == ">":
            count = 0
            continue
        #for each line looked at, increase the count
        count = count + 1
        #if the count is equal to the line number the residue is on, get the residue for this sequence and add it to the residues string
        if which_line == count:
            #print(site_no)
            site_no_on_line = site_no - (50*(which_line-1))
            site = line[site_no_on_line-1]
            residues = residues + site
    alignment_file.close()
    return residues

#calculates the shannon entropy for given residues string
def shannon_entropy(residues):
    nr_of_residues = len(residues)
    unique_residues = list(set(residues))
    #if more than 30% of the residues are gaps, the site is very unlikely to be functionally important
    if residues.count("-") >= len(residues)*30/100:
        print("too many gaps")
        return 50
    #proportion of gaps, for gap penalty
    p_j = (residues.count("-"))/nr_of_residues
    total = 0
    for i in range(len(unique_residues)):
        if unique_residues[i] == "-":
            continue
        if "-" in unique_residues:
            p_i = (residues.count(unique_residues[i]))/nr_of_residues
            total = total + (p_i*math.log2(p_i)/(1-p_j))
        else:
            p_i = (residues.count(unique_residues[i]))/nr_of_residues
            total = total + p_i*math.log2(p_i)
    return -total

#get the conserved sites (shannon entropy < 1.3)
#gets as input a FASTA type alignment file
def get_conserved_sites(align_File):
    alignment_file = open(align_File, "r")
    lines = alignment_file.readlines()
    #get the alignment length to iterate through it and calculate shannon entropy for each site in alignment
    alignment = ""
    count = 0
    #get the alignment for the first sequence, and get its length for alignment bps
    for line in lines:
        if line[0] == ">":
            count = count+1
            continue
        if count == 2:
            break
        alignment = alignment + line.replace("\n","")
    align_length = len(alignment)
    conserved_dict = {}
    for i in range(1, align_length-1):
        residues = get_site(align_File, i)
        s_entropy = shannon_entropy(residues)
        if s_entropy <= 1.3:
            ###add it to dict (add the residues, site_no, what else??)
            conserved_dict[i] = residues + "\t" + str(s_entropy)
    alignment_file.close()
    return conserved_dict

#######main#######
conserved_sites = get_conserved_sites("seqdump.fasta.txt")
newFile = open("conserved_sites.txt","w")
for k,v in conserved_sites.items():
    newFile.write(str(k) + "-" + v + "\n")
newFile.close()

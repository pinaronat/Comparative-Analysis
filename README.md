# Comparative-Analysis
A comparative analysis between the proteome of Mycoplasma californicum and Mycoplasma agalactiae


### get.orthologs.py
contains two functions: *get_bidirectional_best_hits()* and *get_inparalogs()*.

*get_bidirectional_best_hits()* finds the BBH for two species by inputting their BLAST results for each other(BLAST species1 against species2 & BLAST species2 against species1).

*get_inparalogs()* removes possible inparalogs from the BBH list of the given species.


### entropy.py
contains three functions: *get_site()*, *shannon_entropy()* and *get_conserved_sites()*.

*get_site()* gets the residues in given site no for all sequences in the FASTA type alignment file.

*shannon_entropy()* gets as input a string of residues (output of *get_site()*), and calculates the Shannon Entropy.

*get_conserved_sites()* gives as output a list of conserved sites and their Shannon Entropy score for the given FASTA type alignment file. Conserved sites are defined as: Shannon Entropy < 1.3

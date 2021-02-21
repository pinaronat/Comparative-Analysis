#function get_bidirectional_best_hits: Finds the bidirestional best hits between two species' protein alignments (proteins that both have the highest match score for each other)

#inputs BLAST results of one species against the other one's database
def get_bidirectional_best_hits(crossfileName1, crossfileName2):
    file1 = open(crossfileName1,"r")
    lines1 = file1.readlines()
    best_hits_1 = {}
    #get the best hits fot first file
    for line in lines1:
        #if the first element in line != # AND the hit is not the protein itself
        if line[0] != "#" and line.split("\t")[2] != "100.000":
            line_elements = line.split("\t")
            queryID = line_elements[0].split("|")[1]
            subjectID = line_elements[1]
            bitScore = line_elements[11].replace("\n","")
            #if the queryID is not included yet in the dictionary, add it with its corresponding value(subjectID + bit score)
            if not(queryID in best_hits_1.keys()):
                value = subjectID + "\t" + bitScore
                best_hits_1[queryID] = value

    file2 = open(crossfileName2, "r")
    lines2 = file2.readlines()
    best_hits_2 = {}
    best_bidirectional = {}
    #get the best hits for second file
    for line in lines2:
        #if the first element in line != # AND the hit is not the protein itself
        if line[0] != "#" and line.split("\t")[2] != "100.000":
            line_elements = line.split("\t")
            queryID = line_elements[0].split("|")[1]
            subjectID = line_elements[1]
            bitScore = line_elements[11].replace("\n","")
            #if the queryID is not included yet in the dictionary, add it with its corresponding value(subjectID + bit score)
            if not(queryID in best_hits_2.keys()):
                value = subjectID + "\t" + bitScore
                best_hits_2[queryID] = value
                #if this query is also in best_hits_1, AND it has the same subjectID matched in best_hits_1 AND
                #if this query hasn't been added to bidirectional best hits before, add it to bidirectional hits dictionary
                if (subjectID in best_hits_1.keys()) and (queryID == best_hits_1[subjectID].split("\t")[0]) and not(queryID in best_bidirectional.keys()):
                    bidirect_value = subjectID + "\t" + bitScore
                    best_bidirectional[queryID] = bidirect_value
    file1.close()
    file2.close()
    return best_bidirectional

#function get_inparalogs: gets the bidirectional best hits for 2 species, and removes possible inparalogs. A protein is defined to be an inparalog if it has the highest match score
#with a protein in its own proteome rather than with a protein from the opposite species. An inparalog is most likely the result of a duplication of another gene, so it should be removed
#from the list, since the original protein is already matched with an ortholog.

#this function gets 2 inputs: crossFile, which is the BLAST output of one species with the other, and inFile, which is the BLAST output of one species with itself
#The bit-score gives the same value for hits in databases of different sizes and hence can be used for searching in a constantly increasing database
def get_inparalogs(inFileName1, inFileName2, crossFileName1, crossFileName2):
    #get the BBH dictionary from function. keys will be queryIDs of crossFileName1
    main_ortho1 = get_bidirectional_best_hits(crossFileName2,crossFileName1)
    #get the BBH dictionary from function (the other way). keys will be queryIDs of crossFileName2
    main_ortho2 = get_bidirectional_best_hits(crossFileName1,crossFileName2)
    file_in1 = open(inFileName1, "r")
    file_in2 = open(inFileName2, "r")
    lines_i1 = file_in1.readlines()
    lines_i2 = file_in2.readlines()

    BBH_no_inparalog1 = main_ortho1
    #get the best hits for in_file
    for line_i1 in lines_i1:
        #if the first element in line != # AND the hit is not the protein itself
        if line_i1[0] != "#" and line_i1.split("\t")[2] != "100.000":
            line_elements = line_i1.split("\t")
            queryID = line_elements[0].split("|")[1]
            subjectID = line_elements[1]
            in_bitScore = line_elements[11].replace("\n","")
            #look at overlap cutoff (threshold > 50%, to eliminate domains)
            seq_end = int(line_elements[9])
            seq_start = int(line_elements[8])
            q_end = int(line_elements[7])
            q_start = int(line_elements[6])
            align_length = int(line_elements[3])
            if (seq_end-seq_start) >= (q_end-q_start):
                cutoff = align_length / (seq_end-seq_start)
            if (seq_end-seq_start) < (q_end-q_start):
                cutoff = align_length / (q_end-q_start)
            #if threshold is below 0.5, remove the matches from BBH dictionary and move on to the next line/hit
            if cutoff < 0.5 and (queryID in BBH_no_inparalog1.keys()):
                del BBH_no_inparalog1[queryID]
                continue
            if queryID in BBH_no_inparalog1.keys():
                #get the bitscore for in main BBH list
                cross_bitScore = BBH_no_inparalog1[queryID].split("\t")[1]
                #if the found subjectID is an inparalog of the query
                if in_bitScore >= cross_bitScore and (subjectID in BBH_no_inparalog1.keys()):
                    #remove the inparalog of this query from the BBH dictionary
                    del BBH_no_inparalog1[subjectID]


    BBH_no_inparalog2 = main_ortho2
    #get the best hits for in_file
    for line_i2 in lines_i2:
        #if the first element in line != # AND the hit is not the protein itself
        if line_i2[0] != "#" and line_i2.split("\t")[2] != "100.000":
            line_elements = line_i2.split("\t")
            queryID = line_elements[0].split("|")[1]
            subjectID = line_elements[1]
            in_bitScore = line_elements[11].replace("\n","")
            #look at overlap cutoff (threshold > 50%, to eliminate domains)
            seq_end = int(line_elements[9])
            seq_start = int(line_elements[8])
            q_end = int(line_elements[7])
            q_start = int(line_elements[6])
            align_length = int(line_elements[3])
            if (seq_end-seq_start) >= (q_end-q_start):
                cutoff = align_length / (seq_end-seq_start)
            if (seq_end-seq_start) < (q_end-q_start):
                cutoff = align_length / (q_end-q_start)
            #if threshold is below 0.5, remove the matches from BBH dictionary and move on to the next line/hit
            if cutoff < 0.5 and (queryID in BBH_no_inparalog2.keys()):
                del BBH_no_inparalog2[queryID]
                continue
            if queryID in BBH_no_inparalog2.keys():
                #get the bitscore for in main BBH list
                cross_bitScore = BBH_no_inparalog2[queryID].split("\t")[1]
                #if the found subjectID is an inparalog of the query
                if in_bitScore >= cross_bitScore and (subjectID in BBH_no_inparalog2.keys()):
                    #remove the inparalog of this query from the BBH dictionary
                    del BBH_no_inparalog2[subjectID]

    #combine the two dictionaries, get only the BBHs that are found in both
    final_dict = {}
    for key1 in BBH_no_inparalog1.keys():
        for key2 in BBH_no_inparalog2.keys():
            match_seq = BBH_no_inparalog2[key2].split("\t")[0]
            if key1 == match_seq:
                final_dict[key1] = BBH_no_inparalog1[key1]
                break
    file_in1.close()
    file_in2.close()
    return final_dict

####main###
p = get_inparalogs("m_agal.out","m_cali.out","m_cali_db_aga_query.out","m_aga_db_cali_query.out")
writeFile = open("BBH_without_inparalogs.txt","w")
for k, v in p.items():
    writeFile.write(k + "-" + v + "\n")
writeFile.close()

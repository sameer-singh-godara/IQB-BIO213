INT_MAX = pow(10,6)
print("=====================================Assignment 2=====================================")
print("Sameer Singh Godara")
print("2022439")
# Chou and Fasman parameters for Alpha Helix detection
alpha_dictionary = {
    'E' : [1.53, 1],  # Glutamic Acid
    'A' : [1.45, 1],  # Alanine
    'L' : [1.34, 1],  # Leucine
    'H' : [1.24, 1],  # Histidine
    'M' : [1.20, 1],  # Methionine
    'Q' : [1.17, 1],  # Glutamine
    'W' : [1.14, 1],  # Tryptophan
    'V' : [1.14, 1],  # Valine
    'F' : [1.12, 1],  # Phenylalanine
    'K' : [1.07, 0.5],  # Lysine
    'I' : [1.00, 0.5],  # Isoleucine
    'D' : [0.98, 0],   # Aspartic Acid
    'T' : [0.82, 0],   # Threonine
    'S' : [0.79, 0],   # Serine
    'R' : [0.79, 0],   # Arginine
    'C' : [0.77, 0],   # Cysteine
    'N' : [0.73, -1],  # Asparagine
    'Y' : [0.61, -1],  # Tyrosine
    'P' : [0.59, -1],  # Proline
    'G' : [0.53, -1],  # Glycine
}

# Chou and Fasman parameters for Beta strand detection
beta_dictionary = {
    'M' : [1.67, 1],   # Methionine
    'V' : [1.65, 1],   # Valine
    'I' : [1.60, 1],   # Isoleucine
    'C' : [1.30, 1],   # Cysteine
    'Y' : [1.29, 1],   # Tyrosine
    'F' : [1.28, 1],   # Phenylalanine
    'Q' : [1.23, 1],   # Glutamine
    'L' : [1.22, 1],   # Leucine
    'T' : [1.20, 1],   # Threonine
    'W' : [1.19, 1],   # Tryptophan
    'A' : [0.97, 0],   # Alanine
    'R' : [0.90, 0],   # Arginine
    'G' : [0.81, 0],   # Glycine
    'D' : [0.80, 0],   # Aspartic Acid
    'K' : [0.74, -1],  # Lysine
    'S' : [0.72, -1],  # Serine
    'H' : [0.71, -1],  # Histidine
    'N' : [0.65, -1],  # Asparagine
    'P' : [0.62, -1],  # Proline
    'E' : [0.26, -1],  # Glutamic Acid
}


seq = input("Enter sequence (# Note: The sequence must be in single line): ")

seq_len = len(seq)

alpha_mapping = [] # H
beta_mapping = [] # S
ambigious_mapping = [] # will be used for H or S, but will be decided later
final_seq = [] # final sequence

# identify alpha helices;
# for each 6 residue window; 
# check if score is 4 or more than 4;
# if yes, then check the left and right residues;
# if left and right residues are also 4 or more than 4;
# then extend the window size of the helix;
# if not, then check the next window;
def alpha_position(seq, dictionary):
    window_size = 6
    minimum_score = 4
    array = []
    # taking 6 residue window
    for i in range (seq_len-window_size+1):
        window = seq[i:i+window_size]
        score = 0
        count_helix_breaker = 0
        count = 0
        for x in window:
            score += dictionary.get(x)[1]
            if dictionary.get(x)[1] == -1:
                count_helix_breaker += 1
            if dictionary.get(x)[1] >= 0.5:
                count += 1
        if count >= 4:
            for j in range(window_size):
                if not(i+j in array):
                    array.append(i+j)
            
            # extend the window size of the helix right side
            score_for_extension = pow(10,6)
            index_of_extension = i+window_size
            while score_for_extension >= minimum_score:
                if index_of_extension < seq_len:
                    extended_seq = seq[index_of_extension-3:index_of_extension+1]
                    score_for_extension = 0
                    for y in range(len(extended_seq)):
                        score_for_extension += dictionary.get(extended_seq[y])[0]
                    if score_for_extension >= minimum_score:
                        if not(index_of_extension in array):
                            array.append(index_of_extension)
                else:
                    break
                index_of_extension += 1

            # extend the window size of the helix left side
            score_for_extension = INT_MAX
            index_of_extension = i-1
            while score_for_extension >= minimum_score:
                if index_of_extension >= 0:
                    extended_seq = seq[index_of_extension:index_of_extension+4]
                    score_for_extension = 0
                    for y in range(len(extended_seq)):
                        score_for_extension += dictionary.get(extended_seq[y])[0]
                    if score_for_extension >= minimum_score:
                        if not(index_of_extension in array):
                            array.append(index_of_extension)
                else:
                    break
                index_of_extension -= 1
    return array         

# identify beta strands;
# for each 5 residue window;
# check if score is 3 or more than 3;
# if yes, then check the left and right residues;
# if left and right residues are also 3 or more than 3;
# then extend the window size of the strand;
# if not, then check the next window;
def beta_position(seq, dictionary):
    window_size = 5
    minimum_score = 3
    array = []
    # taking 5 residue window
    for i in range (seq_len-window_size+1):
        window = seq[i:i+window_size]
        score = 0
        count_strand_breaker = 0
        count = 0
        for x in window:
            score += dictionary.get(x)[1]
            if dictionary.get(x)[1] == -1:
                count_strand_breaker += 1
            if dictionary.get(x)[1] >= 0.5:
                count += 1
        if count >= 3:
            for j in range(window_size):
                if not(i+j in array):
                    array.append(i+j)
            
            # extend the window size of the strand right side
            score_for_extension = pow(10,6)
            index_of_extension = i+window_size
            while score_for_extension >= minimum_score:
                if index_of_extension < seq_len:
                    extended_seq = seq[index_of_extension-2:index_of_extension+1]
                    score_for_extension = 0
                    for y in range(len(extended_seq)):
                        score_for_extension += dictionary.get(extended_seq[y])[0]
                    if score_for_extension >= minimum_score:
                        if not(index_of_extension in array):
                            array.append(index_of_extension)
                else:
                    break
                index_of_extension += 1

            # extend the window size of the strand left side
            score_for_extension = INT_MAX
            index_of_extension = i-1
            while score_for_extension >= minimum_score:
                if index_of_extension >= 0:
                    extended_seq = seq[index_of_extension:index_of_extension+3]
                    score_for_extension = 0
                    for y in range(len(extended_seq)):
                        score_for_extension += dictionary.get(extended_seq[y])[0]
                    if score_for_extension >= minimum_score:
                        if not(index_of_extension in array):
                            array.append(index_of_extension)
                else:
                    break
                index_of_extension -= 1
    return array

# Set all elements of final sequence to 0
for i in range(seq_len):
    final_seq.append("*")

alpha_mapping = alpha_position(seq, alpha_dictionary)
beta_mapping = beta_position(seq, beta_dictionary)

#print(alpha_mapping)
#print(beta_mapping)

# Update final sequence based on alpha_mapping and beta_mapping
for i in range(seq_len):
    if (i in alpha_mapping) and (i in beta_mapping):
        ambigious_mapping.append(i)
    elif not (i in alpha_mapping) and (i in beta_mapping):
        final_seq[i] = "S"
    elif (i in alpha_mapping) and not (i in beta_mapping):
        final_seq[i] = "H"
    else:
        final_seq[i] = "*"


# initialise the variable
index = 0
counter = 0
alpha_score = 0
beta_score = 0
# print(ambigious_mapping)
# update the final sequence based on ambigious_mapping
while index < len(ambigious_mapping):
    storing_value_of_i = index
    value_of_ambigious = ambigious_mapping[index]
    index += 1
    if index >= len(ambigious_mapping):
        break

    while ambigious_mapping[index] == value_of_ambigious + 1:
        value_of_ambigious += 1
        index += 1
        if index >= len(ambigious_mapping):
            break
    counter += 1
    # print(ambigious_mapping)
    extended_seq = seq[ambigious_mapping[storing_value_of_i]:ambigious_mapping[index-1]+1]

    # calculating the score of the conflicting sites, and one with greater score will be selected
    for x in range(len(extended_seq)):
        alpha_score += alpha_dictionary.get(extended_seq[x])[0]
        beta_score += beta_dictionary.get(extended_seq[x])[0]
    
    if alpha_score > beta_score:
        for y in range (index-storing_value_of_i):
            final_seq[ambigious_mapping[storing_value_of_i]+y] = 'H'
    else:
        for y in range (index-storing_value_of_i):
            final_seq[ambigious_mapping[storing_value_of_i]+y] = 'S'


# printing the final output
# part a - display the sequence regions that are helical in nature
print("\nPart 1: Regions which have tendency of alpha-helix")
H_string = ""
for i in range(seq_len):
    if i in alpha_mapping:
        H_string += "H"
    else:
        H_string += "*"

# print seq and H_stirng in proper alignment

i = 0
while i < seq_len:
    print(f"{i+1}-{min(i+60, seq_len)}")
    print("Sequence:                   ", seq[i:min(i+60, seq_len)])
    print("Regions containing helices: ", H_string[i:min(i+60, seq_len)])
    print("\n")
    i += 60

# part b - display the sequence regions that are beta-strand in nature
print("\nPart 2 : Regions which have tendency of beta-strand")
S_string = ""
for i in range(seq_len):
    if i in beta_mapping:
        S_string += "S"
    else:
        S_string += "*"

i = 0
while i < seq_len:
    print(f"{i+1}-{min(i+60, seq_len)}")
    print("Sequence:                       ", seq[i:min(i+60, seq_len)])
    print("Regions containing beta-strand: ", S_string[i:min(i+60, seq_len)])
    print("\n")
    i += 60


# part c - display the show the conflicting regions and then decide whether the final structure
print("\nPart 3: Final Structure")

print("Code for printing conflicting regions is commented in the code\n")
# code for writing conflicting regions
# 
# flag = 0
# initial_index = 0
# ending_index = 0
# for i in range(seq_len):
#     if H_string[i] == "H" and S_string[i] == "S" and flag == 0:
#         initial_index = i
#         flag = 1
#     if i+1 < seq_len:
#         if (H_string[i+1] != "H" or S_string[i+1] != "S") and flag == 1:
#             ending_index = i
#             flag = 0
#             print(f"Conflicting region: {initial_index+1}-{ending_index+1}")
# 


Final_string = ""
for i in range(seq_len):
    Final_string += final_seq[i]

# print("The sequence is as follows:             ", seq)
# print("The helical regions are as follows:     ", H_string)
# print("The beta-strand regions are as follows: ", S_string)
# print("The final structure is as follows:      ", Final_string)

i = 0
while i < seq_len:
    print(f"{i+1}-{min(i+60, seq_len)}")
    print("Sequence:                       ", seq[i:min(i+60, seq_len)])
    print("Regions containing helices:     ", H_string[i:min(i+60, seq_len)])
    print("Regions containing beta-strand: ", S_string[i:min(i+60, seq_len)])
    print("Final Structure:                ", Final_string[i:min(i+60, seq_len)])
    print("\n")
    i += 60

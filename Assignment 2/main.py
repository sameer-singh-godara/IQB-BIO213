print("=====================================Assignment 2=====================================")
print("Sameer Singh Godara")
print("2022439")
# Chou and Fasman parameters for Alpha Helix detection
alpha_dictionary = {
    'E' : 1.53,  # Glutamic Acid
    'A' : 1.45,  # Alanine
    'L' : 1.34,  # Leucine
    'H' : 1.24,  # Histidine
    'M' : 1.20,  # Methionine
    'Q' : 1.17,  # Glutamine
    'W' : 1.14,  # Tryptophan
    'V' : 1.14,  # Valine
    'F' : 1.12,  # Phenylalanine
    'K' : 1.07,  # Lysine
    'I' : 1.00,  # Isoleucine
    'D' : 0.98,  # Aspartic Acid
    'T' : 0.82,  # Threonine
    'S' : 0.79,  # Serine
    'R' : 0.79,  # Arginine
    'C' : 0.77,  # Cysteine
    'N' : 0.73,  # Asparagine
    'Y' : 0.61,  # Tyrosine
    'P' : 0.59,  # Proline
    'G' : 0.53,  # Glycine
}

# Chou and Fasman parameters for Beta strand detection
beta_dictionary = {
    'M' : 1.67,  # Methionine
    'V' : 1.65,  # Valine
    'I' : 1.60,  # Isoleucine
    'C' : 1.30,  # Cysteine
    'Y' : 1.29,  # Tyrosine
    'F' : 1.28,  # Phenylalanine
    'Q' : 1.23,  # Glutamine
    'L' : 1.22,  # Leucine
    'T' : 1.20,  # Threonine
    'W' : 1.19,  # Tryptophan
    'A' : 0.97,  # Alanine
    'R' : 0.90,  # Arginine
    'G' : 0.81,  # Glycine
    'D' : 0.80,  # Aspartic Acid
    'K' : 0.74,  # Lysine
    'S' : 0.72,  # Serine
    'H' : 0.71,  # Histidine
    'N' : 0.65,  # Asparagine
    'P' : 0.62,  # Proline
    'E' : 0.26,  # Glutamic Acid
}
sequence = "MNASSEGESFAGSVQIPGGTTVLVELTPDIHICGICKQQFNNLDAFVAHKQSGCQLTGTSAAAPSTVQFVSEETVPATQTQTTTRTITSETQTITVSAPEFVFEHGYQTYLPTESNENQTATVISLPAKSRTKKPTTPPAQKRLNCCYPGCQFKTAYGMKDMERHLKIHTGDKPHKCEVCGKCFSRKDKLKTHMRCHTGVKPYKCKTCDYAAADSSSLNKHLRIHSDERPFKCQICPYASRNSSQLTVHLRSHTASELDDDVPKANCLSTESTDTPKAPVITLPSEAREQMATLGERTFNCCYPGCHFKTVHGMKDLDRHLRIHTGDKPHKCEFCDKCFSRKDNLTMHMRCHTSVKPHKCHLCDYAAVDSSSLKKHLRIHSDERPYKCQLCPYASRNSSQLTVHLRSHTGDTPFQCWLCSAKFKISSDLKRHMIVHSGEKPFKCEFCDVRCTMKANLKSHIRIKHTFKCLHCAFQGRDRADLLEHSRLHQADHPEKCPECSYSCSSAAALRVHSRVHCKDRPFKCDFCSFDTKRPSSLAKHVDKVHRDEAKTENRAPLGKEGLREGSSQHVAKIVTQRAFRCETCGASFVRDDSLRCHKKQHSDQSENKNSDLVTFPPESGASGQLSTLVSVGQLEAPLEPSQDL"

seq = sequence

seq_len = len(seq)


def alpha_position(seq, dictionary):
    H_array = ["-"]*seq_len
    window_size = 6
    minimum_count = 4

    for i in range(seq_len-window_size+1):
        count = 0
        for j in range(i, i+window_size):
            if dictionary[seq[j]] >= 1:
                count += 1
        
        if count >= minimum_count:
            for j in range(i, i+window_size):
                H_array[j] = "H"
            
            # extend the window size of the helix right side
            sites = []
            score_for_extension = 0
            for j in range(i+3, seq_len-3):
                sites.append(j)
                score_for_extension = 0
                for k in range(j, j+4):
                    score_for_extension += dictionary[seq[k]]
                if score_for_extension >= 4:
                    H_array[j:j+4] = ["H"] * 4
                else:
                    break
                    
            # extend the window size of the helix left side
            sites = []
            for j in reversed(range(i+3)):
                sites.append(j)
                score_for_extension = 0
                for k in range(j-3, j+1):
                    score_for_extension += dictionary[seq[k]]
                if score_for_extension >= 4:
                    H_array[j-3:j+1] = ["H"] * 4
                else:
                    break
    return H_array         


def beta_position(seq, dictionary):
    window_size = 5
    minimum_count = 3
    S_array = ["-"]*seq_len

    for i in range (seq_len - window_size + 1):
        count = 0
        for j in range(i, i+window_size):
            if dictionary[seq[j]] >= 1:
                count += 1
        
        if count >= minimum_count:
            for j in range(i, i+window_size):
                S_array[j] = "S"

            # extend the window size of the beta right side
            sites = []
            score_for_extension = 0
            for j in range(i+3, seq_len-3):
                sites.append(j)
                score_for_extension = 0
                for k in range(j, j+4):
                    score_for_extension += dictionary[seq[k]]
                if score_for_extension > 4:
                    S_array[j:j+4] = ["S"] * 4
                else:
                    break
            
            # extend the window size of the beta left side
            sites = []
            score_for_extension = 0
            for j in reversed(range(i+3)):
                sites.append(j)
                score_for_extension = 0
                for k in range(j-3, j+1):
                    score_for_extension += dictionary[seq[k]]
                if score_for_extension >= 4:
                    S_array[j-3:j+1] = ["S"] * 4
                else:
                    break
            

    return S_array

final_mapping = []
for i in range(seq_len):
    final_mapping.append("-")

alpha_mapping = alpha_position(seq, alpha_dictionary)
beta_mapping = beta_position(seq, beta_dictionary)

conflicting_area = ["-"]*seq_len
for i in range(seq_len):
    if alpha_mapping[i] == "H" and beta_mapping[i] == "S":
        conflicting_area[i] = "^"

# print(alpha_mapping)
# print(beta_mapping)
# print(conflicting_area)

alpha_strand_score = 0
beta_strand_score = 0
flag = 0
start = 0

index_of_conflict = []
index_of_non_conflict = []
for i in range(seq_len):
    if conflicting_area[i] == "^": # when there is a conflict
        index_of_conflict.append(i)
        alpha_strand_score += alpha_dictionary[seq[i]]
        beta_strand_score += beta_dictionary[seq[i]]
        
        if flag == 0:
            start = i
            flag = 1
    
    else: # only come here when there is no conflict
        index_of_non_conflict.append(i)
        if alpha_mapping[i] == "H" and beta_mapping[i] == "-":
            final_mapping[i] = "H"
        elif alpha_mapping[i] == "-" and beta_mapping[i] == "S":    
            final_mapping[i] = "S"
        
        if flag == 1:
            if alpha_strand_score < beta_strand_score:
                for j in range(start, i):
                    final_mapping[j] = "S"
            else:
                for j in range(start, i):
                    final_mapping[j] = "H"
            flag = 0
            alpha_strand_score = 0
            beta_strand_score = 0

    if i == seq_len-1 and flag == 1:
        if alpha_strand_score < beta_strand_score:
            for j in range(start, i+1):
                final_mapping[j] = "S"
        else:
            for j in range(start, i+1):
                final_mapping[j] = "H"


# print H, S, * in the final_mapping

alpha_mapping_string = ("".join(alpha_mapping)) 
beta_mapping_string = ("".join(beta_mapping)) 
conflicting_area_string = ("".join(conflicting_area)) 
final_mapping_string = ("".join(final_mapping)) 
print("Part 1: Alpha Helix Prediction")
# write in fragments of 60 and also write the indexes minimum of length and 60 
for i in range(0, seq_len, 60):
    print("Index: ", i+1, " to ", min(i+60, seq_len))
    print("Original Sequence:      ", seq[i:i+60])
    print("Alpha Helix Prediction: ", alpha_mapping_string[i:i+60])
    print("\n")

print("Part 2: Beta Strand Prediction")

for i in range(0, seq_len, 60):
    print("Index: ", i+1, " to ", min(i+60, seq_len))
    print("Original Sequence:      ", seq[i:i+60])
    print("Beta Strand Prediction: ", beta_mapping_string[i:i+60])
    print("\n")


print("Part 3: Conflicting Area")
for i in range(0, seq_len, 60):
    print("Index: ", i+1, " to ", min(i+60, seq_len))
    print("Original Sequence:      ", seq[i:i+60])
    print("Conflicting Area:       ", conflicting_area_string[i:i+60])
    print("\n")


print("Part 4: Final Prediction")

for i in range(0, seq_len, 60):
    print("Index: ", i+1, " to ", min(i+60, seq_len))
    print("Original Sequence:      ", seq[i:i+60])
    print("Alpha Helix Prediction: ", alpha_mapping_string[i:i+60])
    print("Beta Strand Prediction: ", beta_mapping_string[i:i+60])
    print("Final Prediction:       ", final_mapping_string[i:i+60])
    print("\n")

print("=========================================================================================")
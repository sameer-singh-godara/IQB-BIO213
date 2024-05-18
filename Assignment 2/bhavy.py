#calculate scores for each amino acid in the sequence based on predefined alpha-helix and beta-sheet propensity values
def calculate_alpha_score(sequence):  
    alpha_values = {
        'A': 1.45, 'C': 0.77, 'D': 0.98, 'E': 1.53, 'F': 1.12,
        'G': 0.53, 'H': 1.24, 'I': 1.0, 'K': 1.07, 'L': 1.34,
        'M': 1.2, 'N': 0.73, 'P': 0.59, 'Q': 1.17, 'R': 0.79,
        'S': 0.79, 'T': 0.82, 'V': 1.14, 'W': 1.14, 'Y': 0.61
    }
    return sum(alpha_values[i] for i in sequence)

def calculate_beta_score(sequence):  
    beta_values = {
        'A': 0.97, 'C': 1.30, 'D': 0.80, 'E': 0.26, 'F': 1.28,
        'G': 0.81, 'H': 0.71, 'I': 1.60, 'K': 0.74, 'L': 1.22,
        'M': 1.67, 'N': 0.65, 'P': 0.62, 'Q': 1.23, 'R': 0.90,
        'S': 0.72, 'T': 1.20, 'V': 1.65, 'W': 1.19, 'Y': 1.29
    }
    return sum(beta_values[i] for i in sequence)

#identify regions in the sequence likely to form alpha-helices and beta-sheets, respectively, based on the calculated scores along with extensions
def find_alpha_regions(sequence):
    helix_structure = ['_' for _ in range(len(sequence))]
    i = 0
    while i < len(sequence)-5:
        n = sum(1 for j in range(i, i+6) if calculate_alpha_score(sequence[j]) >= 1)
        if n >= 4:
            for j in range(i, i+6):
                if helix_structure[j] != 'H':
                    helix_structure[j] = 'H'
            p1, p2 = i+6, i-1
            while p1 < len(sequence) and calculate_alpha_score(sequence[p1-3:p1+1]) >= 4:
                helix_structure[p1] = 'H'
                p1 += 1
            while p2 >= 0 and calculate_alpha_score(sequence[p2:p2+4]) >= 4:
                helix_structure[p2] = 'H'
                p2 -= 1
        i += 1

    return "".join(helix_structure)

def find_beta_regions(sequence): 
    beta_sheet_structure = ["_" for _ in range(len(sequence))]
    i = 0
    while i < len(sequence) - 4:
        if sum(1 for j in range(i, i + 5) if calculate_beta_score(sequence[j]) >= 1) >= 3:
            for j in range(i, i+5):
                if beta_sheet_structure[j] != 'S':
                    beta_sheet_structure[j] = 'S'
            p1, p2 = i+5, i-1
            while p1 < len(sequence) and calculate_beta_score(sequence[p1-3:p1+1]) >= 4:
                beta_sheet_structure[p1] = 'S'
                p1 += 1
            while p2 >= 0 and calculate_beta_score(sequence[p2:p2+4]) >= 4:
                beta_sheet_structure[p2] = 'S'
                p2 -= 1
        i += 1
    return "".join(beta_sheet_structure)

#resolve conflicts between predicted alpha-helix and beta-sheet regions by assigning secondary structure based on the relative scores of the conflicting regions
def resolve_conflicts(seq1, seq2, sequence):
    resolved_structure = ["" for i in range(len(sequence))]
    sequence_length = len(sequence)
    i = 0
    while i < sequence_length:
        if (seq1[i] == 'H' and seq2[i] == '') or (seq1[i] == '' and seq2[i] == 'H'):
            resolved_structure[i] = 'H'
            i += 1
        elif (seq1[i] == '' and seq2[i] == 'S') or (seq1[i] == 'S' and seq2[i] == '') :
            resolved_structure[i] = 'S'
            i += 1
        elif seq1[i] == '' and seq2[i] == '':
            resolved_structure[i] = 'X'
            i += 1
        else:
            n = 0
            while (seq1[i] == 'H' and seq2[i] == 'S') or (seq1[i] == 'S' and seq2[i] == 'H'):
                n += 1
                if i < sequence_length-1:
                    i += 1
                else:
                    break
            if i == sequence_length-1:
                i += 1 
            alpha_score = calculate_alpha_score(sequence[i-n:i])
            beta_score = calculate_beta_score(sequence[i-n:i])
            if alpha_score > beta_score:
                for k in range(i-n,i):
                    resolved_structure[k] = 'H'
            else:
                for k in range(i-n,i):
                    resolved_structure[k] = 'S'
    return "".join(resolved_structure)              

#apply conflict resolution to determine the final secondary structure of the sequence.
def resolve_conflict_case(seq1, seq2, sequence):
    resolved_structure = resolve_conflicts(seq1, seq2, sequence)
    return resolved_structure

# visualize the predicted secondary structure of the sequence
def plot_secondary_structure(sequence, structure):
    colors = {'H': 'red', 'S': 'blue', '_': 'white', 'X': 'green'}
    fig, ax = plt.subplots(figsize=(10, 0.5))
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    plt.xlabel("Secondary Structure Visualization")
    ax.set_title('Secondary Structure Prediction')
    ax.set_xlim(-0.5, len(sequence) - 0.5)
    i = 0
    for aa, ss in zip(sequence, structure):
        ax.fill_betweenx([0, 1], i - 0.5, i + 0.5, color=colors[ss], edgecolor='black', linewidth=0.5)
        i += 1
    plt.show()

#define the sequence here 
sequence = 'MNASSEGESFAGSVQIPGGTTVLVELTPDIHICGICKQQFNNLDAFVAHKQSGCQLTGTSAAAPSTVQFVSEETVPATQTQTTTRTITSETQTITVSAPEFVFEHGYQTYLPTESNENQTATVISLPAKSRTKKPTTPPAQKRLNCCYPGCQFKTAYGMKDMERHLKIHTGDKPHKCEVCGKCFSRKDKLKTHMRCHTGVKPYKCKTCDYAAADSSSLNKHLRIHSDERPFKCQICPYASRNSSQLTVHLRSHTASELDDDVPKANCLSTESTDTPKAPVITLPSEAREQMATLGERTFNCCYPGCHFKTVHGMKDLDRHLRIHTGDKPHKCEFCDKCFSRKDNLTMHMRCHTSVKPHKCHLCDYAAVDSSSLKKHLRIHSDERPYKCQLCPYASRNSSQLTVHLRSHTGDTPFQCWLCSAKFKISSDLKRHMIVHSGEKPFKCEFCDVRCTMKANLKSHIRIKHTFKCLHCAFQGRDRADLLEHSRLHQADHPEKCPECSYSCSSAAALRVHSRVHCKDRPFKCDFCSFDTKRPSSLAKHVDKVHRDEAKTENRAPLGKEGLREGSSQHVAKIVTQRAFRCETCGASFVRDDSLRCHKKQHSDQSENKNSDLVTFPPESGASGQLSTLVSVGQLEAPLEPSQDL'

#call the functions
alpha_helix_structure = find_alpha_regions(sequence)  
beta_sheet_structure = find_beta_regions(sequence)  

#print in a versitile manner
print("alpha helix sites:",alpha_helix_structure)
print()
print("beta helix sites:",beta_sheet_structure)

final_structure = resolve_conflict_case(alpha_helix_structure, beta_sheet_structure, sequence)  
print("----------------------------------------------")
print("Final structure is: " + final_structure)

#call to plot the final graph

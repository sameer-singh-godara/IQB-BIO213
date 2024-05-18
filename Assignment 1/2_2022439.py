import numpy as np
import pandas as pd

# Smith-Waterman algorithm for local sequence alignment. 
def sw(x, y, match=2, mismatch=-1, gap=-3):
    nx = len(x)
    ny = len(y)

    # Initialization of the matrix.
    F = np.zeros((nx + 1, ny + 1))

    # Pointers to trace through an optimal alignment.
    P = np.zeros((nx + 1, ny + 1), dtype=int)
    P[:, 0] = 3
    P[0, :] = 4

    max_score = 0
    max_i = 0
    max_j = 0

    # Matrix filling.
    max_positions = []
    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            match_score = match if x[i - 1] == y[j - 1] else mismatch
            a1 = F[i - 1, j - 1] + match_score
            a2 = F[i - 1, j] + gap
            a3 = F[i, j - 1] + gap
            F[i, j] = max(a1, a2, a3, 0)
            if (a1 == F[i, j]): # Match or mismatch
                P[i, j] += 2
            if (a2 == F[i, j]): # Gap in one of the sequence
                P[i, j] += 3
            if (a3 == F[i, j]): # Gap in the other sequence
                P[i, j] += 4

            if F[i, j] > max_score:
                max_score = F[i, j]
                max_positions = [(i, j)]
            elif F[i, j] == max_score:
                max_positions.append((i, j))

    # Print scoring matrix using pandas DataFrame
    print("Scoring Matrix:")
    df = pd.DataFrame(F, index=['-'] + list(x), columns=['-'] + list(y))
    print(df)

    # Iterative traceback
    alignments = []
    for max_i, max_j in max_positions:
        stack = [(max_i, max_j, '', '', 0)]
        # Using stack to avoid recursion depth limit.
        # Each element in the stack is a tuple (i, j, rx, ry, score).
        # i, j: current position in the matrix.
        # rx, ry: partial alignments.
        # score: current alignment score.
        # The stack is initialized with one of the maximum scores in the matrix.
        # The traceback is performed iteratively until the stack is empty.
        # The traceback is performed in a depth-first manner.
        # The traceback is finished when the current score is 0.
        while stack:
            i, j, rx, ry, score = stack.pop()
            if F[i, j] == 0:
                alignments.append((rx[::-1], ry[::-1], score))
            else:
                if P[i, j] in [2, 5, 6, 9]:
                    if x[i - 1] == y[j - 1]:
                        stack.append((i - 1, j - 1, rx + x[i - 1], ry + y[j - 1], score + match))
                    else:
                        stack.append((i - 1, j - 1, rx + x[i - 1], ry + y[j - 1], score + mismatch))
                if P[i, j] in [3, 5, 7, 9]:
                    stack.append((i - 1, j, rx + x[i - 1], ry + '-', score + gap))
                if P[i, j] in [4, 6, 7, 9]:
                    stack.append((i, j - 1, rx + '-', ry + y[j - 1], score + gap))

    # Print all optimal alignments.
    print("\nAll Optimal Local Alignments:\n")
    for alignment in alignments: # this will print the motif of the alignment which is common in both
        print('\n'.join(alignment[:2]))
        print("Score:", alignment[2])
        print()

    return

seq1 = "GATGCGCAG"
seq2 = "GGCAGTA"

sw(seq1, seq2)
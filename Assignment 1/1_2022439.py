import numpy as np
import pandas as pd

# Needleman-Wunsch algorithm for global sequence alignment.
def nw(x, y, match=2, mismatch=-3, gap=-1):
    nx = len(x)
    ny = len(y)

    # Initialization of the matrix.
    F = np.zeros((nx + 1, ny + 1))
    F[:, 0] = np.linspace(0, gap * nx, nx + 1)
    F[0, :] = np.linspace(0, gap * ny, ny + 1)

    # Pointers to trace through an optimal alignment.
    P = np.zeros((nx + 1, ny + 1), dtype=int)
    P[:, 0] = 3
    P[0, :] = 4

    # Matrix filling.
    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            match_score = match if x[i - 1] == y[j - 1] else mismatch
            a1 = F[i - 1, j - 1] + match_score
            a2 = F[i - 1, j] + gap 
            a3 = F[i, j - 1] + gap
            F[i, j] = max(a1, a2, a3)
            if (a1 == F[i, j]): # Match or mismatch
                P[i, j] += 2
            if (a2 == F[i, j]): # Gap in one of the sequence 
                P[i, j] += 3
            if (a3 == F[i, j]): # Gap in the other sequence
                P[i, j] += 4

    # Print scoring matrix using pandas DataFrame
    print("Scoring Matrix:")
    df = pd.DataFrame(F, index=['-'] + list(x), columns=['-'] + list(y))
    print(df)

    # Iterative traceback
    alignments = []
    stack = [(nx, ny, '', '', 0)]
    # Using stack to avoid recursion depth limit.
    # Each element in the stack is a tuple (i, j, rx, ry, score).
    # i, j: current position in the matrix.
    # rx, ry: partial alignments.
    # score: current alignment score.
    # The stack is initialized with the last element of the matrix. Last element is the bottom right corner of the matrix.
    # The traceback is performed iteratively until the stack is empty.
    # The traceback is performed in a depth-first manner.
    # The traceback is finished when the current position is (0, 0) and current score is 0.
    while stack:
        i, j, rx, ry, score = stack.pop()
        if i == 0 and j == 0:
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
    print("\nAll Optimal Global Alignments:\n")
    for alignment in alignments: # this will print the complete alignment end to end
        print('\n'.join(alignment[:2]))
        print("Score:", alignment[2])
        print()

    return

seq1 = "GATGCGCAG"
seq2 = "GGCAGTA"

nw(seq1, seq2)
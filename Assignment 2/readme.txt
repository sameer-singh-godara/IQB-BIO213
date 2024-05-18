Protein Secondary Structure Prediction using Chou and Fasman Algorithm

This Python script implements the Chou and Fasman algorithm for predicting alpha helix and beta strand secondary structures in protein sequences. The algorithm utilizes propensity parameters for amino acids to predict regions of secondary structure based on local sequence characteristics.

1. Files:
    - main.py: Python script for secondary structure prediction using the Chou and Fasman algorithm.
    - readme.txt: This file.

2. Usage:
    - Input the protein sequence in the sequence variable within the script.
    - Run the script using Python.
    - The script will output the predicted secondary structure, including alpha helix regions, beta strand regions, areas of conflict and the final predicted structure in fragments of 60 characters.

3. Functionality:
    - The script implements the Chou and Fasman algorithm, which uses propensity parameters for amino acids to predict alpha helices and beta strands.
    - Propensity parameters are used as given in the question.
    - The algorithm scans the protein sequence and evaluates each window of residues and based on the concept the structure of local sequence is predicted. 
    - If particular window size is accepted, then the window of residues is extended in both the direction.

4. References:
    - Theory taught in class, algorithm is mentioned in the slides

5. Disclaimer:
    - This script provides predictions based on the Chou and Fasman algorithm and may not always accurately reflect the actual secondary structure of a protein.
    - Users are encouraged to verify the predictions using experimental or more advanced computational methods.

For any questions or issues, please contact me at sameer22439@iiitd.ac.in

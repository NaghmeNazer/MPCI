# MPCI
Quantify methylation haplotypes in a region
________________________________________
# MPCI Calculation Script
This script calculates the Methylation Pattern Consistency Index (MPCI), a novel metric for quantifying consistent methylation patterns across sequencing reads in differentially methylated regions (DMRs). The script also includes helper functions for calculating signed Manhattan similarity and assigning weights based on methylation status.
________________________________________
# Requirements
•	R (version 4.0 or higher)
________________________________________
# Functions
1. give_sign_for_weights(x_1, x_2)

•	Purpose: Assigns weights based on the average methylation status of two rows.

•	Input: Two binary vectors (x_1, x_2) representing methylation status of CpG sites (1 = methylated, 0 = unmethylated).

•	Output:

o	1 if the average methylation is > 0.5 (more methylated).

o	-1 if the average methylation is < 0.5 (more unmethylated).

o	Randomly 1 or -1 if the average methylation is exactly 0.5.

2. signed_manhattan_sim(binary_dmr)

•	Purpose: Calculates the signed Manhattan similarity for a binary DMR matrix.

•	Input: A binary matrix (binary_dmr) where rows represent sequencing reads and columns represent CpG sites.

•	Output: A weighted average of pairwise similarities between rows.

3. MPCI(binary_dmr)

•	Purpose: Calculates the Methylation Pattern Consistency Index (MPCI) for a binary DMR matrix.

•	Input: A binary matrix (binary_dmr) where rows represent sequencing reads and columns represent CpG sites.

•	Output: The MPCI value, ranging from -1 to +1, where:

o	Positive values indicate consistent methylation.

o	Negative values indicate consistent unmethylation.

o	Values near zero indicate random methylation patterns.
________________________________________
# Usage
1. Synthetic Data
The script includes examples of synthetic binary DMR matrices to test the MPCI function:

•	fully_methylated: A fully methylated DMR.

•	mixed_methylation_1: A DMR with mixed methylation patterns.

•	mostly_unmethylated: A mostly unmethylated DMR with some random methylation.

•	fully_unmethylated: A fully unmethylated DMR.

•	mixed_methylation_2: A DMR with mixed methylation patterns.

•	mostly_methylated: A mostly methylated DMR with some random unmethylation.

To test MPCI on these matrices, simply run the corresponding lines in the script.

2. Real Data
To calculate MPCI for a real binary DMR matrix:

1.	Load the binary DMR data from a CSV file using read.csv().

2.	Remove unnecessary columns (e.g., real_binary_dmr$X <- NULL).

3.	Pass the cleaned matrix to the MPCI() function.
________________________________________
# Creating a Binary DMR Matrix from Bisulfite Sequencing Reads
To create a binary DMR matrix from bisulfite sequencing reads:

1.	Extract Methylation Calls

2.	Convert to Binary Format: For each read, create a binary vector where:

o	1 represents a methylated CpG.

o	0 represents an unmethylated CpG.

o	Missing or ambiguous calls can be represented as NA.

3.	Build the Matrix: Combine these binary vectors into a matrix where:

o	Rows represent sequencing reads.

o	Columns represent CpG sites.

4.	Save as CSV: Save the binary matrix as a CSV file for input into the MPCI script.
________________________________________
# Output
The script outputs the MPCI value for each binary DMR matrix.
________________________________________
# Notes

•	Ensure that the input binary DMR matrix is properly formatted (rows = reads, columns = CpG sites).

•	The script handles edge cases and returns NA if similarity cannot be calculated.
________________________________________
# Contact
naghme93@gmail.com


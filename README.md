DoubletFinder is an R package that:

Simulates artificial doublets by averaging expression profiles of random real cells.

Projects real and artificial cells into PCA space.

Uses k-nearest neighbors to determine how close each real cell is to the artificial doublets.

Calculates a doublet score (pANN) for each cell.

Classifies cells as “Singlet” or “Doublet” based on expected doublet rates (e.g., ~0.8% per 1,000 cells for 10X).

# Phase type theory for the coalescent with structure

In this repo, we mainly show the implementation of formula (5.22) in Coalescent Theory: An introduction by John Wakeley (2009), which generates the subintensity matrix for the coalescent with structure. Using a phase-type distribution, this can be used to model e.g. the time to the most recent common ancestor ($T_{MRCA}$) and provides an alternative to the first-step analysis formulas in Wakeleys book.
This implementation is given in `generate_subint_matrix.R`

Using the PhaseTypeR package (https://rivasiker.github.io/PhaseTypeR/), we simulate the total branch length of the trees and sprinkle on mutations, giving us simulated realizations of the number of mutations on the genealogy. Using method of moments and MLE, we attempt to recover the mutation rate and migration rate used to simulate the data. 
This is done in `Simulations.Rmd`.

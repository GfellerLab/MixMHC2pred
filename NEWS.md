# Release notes from MixMHC2pred

## Version 2.0.2 (08.12.2022)

* Definition files of the alleles are not tracked on github anymore but can
  be downloaded from <http://mixmhc2pred.gfellerlab.org/PWMdef>, including
  alleles from multiple species.
* Added the option `-e` to output intermediate scores in addition to the final
  *%Rank*.

## Version 2.0 (03.10.2022)

* Release of MixMHC2pred-v2.0 (as full version, not RC version anymore).
* No major change with respect to the 2.0.RC (see news of 10.06.2022). It was
  only retrained after minor changes on the training data used.

## Version 2.0.RC (10.06.2022)

* Major update of MixMHC2pred, allowing for pan-allele predictions for human and
  non-human MHC-II alleles.
* Implementation vastly changed; it is now based on multiple neural networks and
  it was trained on a much larger MHC-II peptidomics dataset.
* Please note that the format of the input and output files have changed, as
  well as the available input arguments.

## Version 1.2 (03.02.2020)

* Output contains now also the individual predictions from each allele in
  addition to the score of the best allele per peptide.
* Outputted core positions correspond now to the P1 of the binding core, counted
  from the N-terminus of the peptide (first amino acid has a value of 1).
* Added a "FAQ" document.
* There is an option to allow having unspecified amino acids in the peptide
  sequence (these should be written as "X" or "-").
* Added the option "--best_s" (see README).
* Slightly updated the %Rank computations, avoiding in general having scores of
  exactly 0, and allowing the predictions for peptides up to a size of 25 amino
  acids.
* Corrected a bug when input file was given in fasta format with multiple lines
  per peptide sequence.

## Version 1.1 (14.10.2019)

* Public release of MixMHC2pred accompanying the Nature Biotechnology paper.

## Version 1.0 (08.02.2019)

* Public release of MixMHC2pred accompanying our bioRxiv preprint.

# Frequently asked questions for MixMHC2pred

## What kind of predictions can I make with MixMHC2pred

MixMHC2pred is a predictor of MHC-II ligands displayed at the cell surface.
It combines predictions of affinity to MHC-II molecules together with other
features linked to antigen processing and presentation. It can work for any
human MHC-II allele (HLA-DR/-DP/-DQ) as well as for many alleles from
other species.

Peptides with lower %Ranks in MixMHC2pred predictions are more likely to
elicit CD4+ T cell recognition

## How should I interpret the output of MixMHC2pred

* *Peptide* and *Context* give the sequence of the peptides (respectively
  context), in the same order as provided in input.
* *BestAllele* gives the name of the best allele (based
  on the allele with best predicted score for the given peptide).
* *%Rank_...* give the predicted score for the best allele and for each allele
  asked. The score is given as a percentile rank, best score is about 0, worst
  score is 100. This value corresponds to the percent of random peptides that
  would have a better score than the peptide provided in input among peptides of
  sizes 12-21 amino acids.
* *Core_best* indicates the best predicted core binding sequence for each
  peptide towards its best allele.
* *CoreP1_...* give the most likely binding core position for the given peptide
  towards the given allele (this tells the position of the first amino acid from
  the binding core (which has a size of 9 aa in the predictions), starting at a
  value of 1 (i.e. if binding core corresponds to the 9 first amino acids
  from the peptide, this *CoreP1 = 1*)).
* *subSpec_...* tell in which sub-specificity the given peptide
  is likely bound toward the given allele. The value *1* corresponds to the
  main sub-specificity (the only one for multiple alleles). But for example
  for *DRB1\*08:01* allele a 2nd sub-specificity exists and is indicated by the
  value *2*. For alleles accomodating reverse binding, a value of *-1* indicates
  that the given peptide is bound in the reverse orientation.

Peptides that are too short (less than 12 amino acids), too long (more than
21 amino acids) or that contain non standard amino acids have *NA* values
instead of their scores.

## How can I rank my peptides based on MixMHC2pred predictions

The best way to rank your peptides is to use the global score with the best
allele (*%Rank_best*), the best predicted peptides have the lowest
scores.

## Can I use MixMHC2pred for commercial purposes

If you plan to use MixMHC2pred for commercial purposes, you are required to
obtain a separate license. To do so, please contact <nbulgin@lcr.org>
at the Ludwig Institute for Cancer Research Ltd.

## Who should I contact in case of a technical or other issue

Julien Racle ([julien.racle@unil.ch](mailto:julien.racle@unil.ch)). Please
provide as much details as possible and ideally send also an example input file
that is causing the issue.

## How should I cite MixMHC2pred

If you are using MixMHC2pred, please refer to

Racle, J., et al., Machine learning predictions of MHC-II specificities reveal
alternative binding mode of class II epitopes. *bioRxiv* (2022)
(available [here](https://doi.org/10.1101/2022.06.26.497561)).

and

Racle, J., et al. Robust prediction of HLA class II epitopes by deep motif
deconvolution of immunopeptidomes. *Nat. Biotechnol.* 37, 1283–1286 (2019)
(available [here](https://www.nature.com/articles/s41587-019-0289-6)).

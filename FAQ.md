# Frequently asked questions for MixMHC2pred

## What kind of predictions can I make with MixMHC2pred

MixMHC2pred is a predictor of HLA-II ligand displayed at the cell surface.
It combines predictions of affinity to HLA-II molecules including
binding core offset preferences, as well as peptide N-/C-terminal motifs and
peptide length preference (which are likely emerging from the peptide
processing steps).

Peptides ranking high in MixMHC2pred predictions are more likely to elicit
CD4+ T cell recognition.

## How should I interpret the output of MixMHC2pred

* *Peptide* gives the sequence of the peptides, in the same order as
  provided in input.
* *BestAllele* gives the name of the best allele (based
  on the allele with best predicted score for the given peptide).
* *%Rank_...* give the predicted score for the best allele and separately for
  each allele asked. The score is given as a percentile rank (i.e., the
  percent of random peptides that  would have a score higher than the peptide
  provided in input among peptides of sizes 12-25 amino acids; best score is
  about 0, worst score is 100).
* *%Rank_best_perL* is similar to *%Rank_best* but computed based only on peptides
  having the same length than the given peptide. This score thus doesn't follow
  the length distribution observed in naturally presented ligands.
* *Core_best* indicates the best predicted core binding sequence for each
  peptide towards its best allele.
* *CoreP1_...* give the most likely binding core position for the given peptide
  towards the allele (this tells the position of the first amino acid from the
  binding core (which has a size of 9 aa in the predictions), starting at a
  value of 1 (i.e. if binding core corresponds to the 9 first amino acids
  from the peptide, this *CoreP1 = 1*)).

Peptides that are too short (less than 12 amino acids), too long (more than
25 amino acids) or that contain non standard amino acids have *NA* values
instead of their scores.

## How can I rank my peptides based on MixMHC2pred predictions

The best way to rank your peptides is to use the global score with the best
allele (*%Rank_best*), the best predicted peptides have the lowest
scores.

## Can I use MixMHC2pred for commercial purposes

If you plan to use MixMHC2pred for commercial purposes, you are required to
obtain a separate license. To do so, please contact <eauffarth@licr.org>
at the Ludwig Institute for Cancer Research Ltd.

## Who should I contact in case of a technical or other issue

Julien Racle ([julien.racle@unil.ch](mailto:julien.racle@unil.ch)). Please
provide as much details as possible and ideally send also an example input file
that is causing the issue.

## How should I cite MixMHC2pred

If you are using MixMHC2pred, please refer to

Racle, J., et al. Robust prediction of HLA class II epitopes by deep motif
deconvolution of immunopeptidomes. *Nat. Biotechnol.* 37, 1283–1286 (2019).

It is available [here](https://www.nature.com/articles/s41587-019-0289-6).

# CRISPR CIRCLE-seq — research notes

## Assay target

The sequenced material is circularized human genomic DNA linearized in vitro at Cas9 cleavage sites. This protocol does not target, amplify or reconstruct a pathogen genome: every oligo it uses is a host-genome or library oligo (adapter duplexes and library PCR primers). The note is here so that the protocol's CRISPR, off-target and double-strand-break vocabulary is not read as pathogen work.

## Sources read

- 🟢 Tsai et al. 2017, DOI [10.1038/nmeth.4278](https://doi.org/10.1038/nmeth.4278).
  Round 2026-10-17 re-read: PMC5924695 JATS Methods lines 115-116 (shearing to 300 bp, oSQT1288
  with its phosphorothioate, Lambda Exonuclease/Exonuclease I, USER + T4 PNK, 5 ng/µl T4 ligase
  circularization, Plasmid-Safe ATP-dependent DNase, 90 nM SpCas9 and 90 nM gRNA on 250 ng DNA,
  KAPA HiFi PCR, 150-bp paired-end MiSeq) and Figure 1 legend line 450-451 (4-nt palindromic
  overhangs); Supplementary Protocol (Springer MOESM84) text lines 141-200 (exonuclease, USER/T4
  PNK, 16 °C 16 h circularization, Plasmid-Safe) and 296-320 (USER, NEBNext i5/i7 primers,
  KAPA HiFi HotStart, 22 cycles).

## Molecular order

🟢 Genomic DNA is sheared to about 300 bp, repaired, dA-tailed and ligated to the
uracil-containing oSQT1288 stem-loop `5'-P-CGGTGGACCGATGATCUATCGGTCCACCG*T-3'`. Lambda
Exonuclease plus E. coli Exonuclease I destroy molecules that kept a free end; USER plus
T4 PNK then excise the uracil and phosphorylate, exposing the 4-nt palindromic overhangs
(Fig. 1 legend), and T4 DNA ligase at 5 ng/µl (16 °C, 16 h) closes one-fragment circles.
Plasmid-Safe ATP-dependent DNase removes residual linear DNA. 90 nM SpCas9 with 90 nM
in-vitro-transcribed gRNA linearizes 250 ng of target-bearing circles. The fresh ends are
dA-tailed, ligated to the USER-openable hairpin adapter, opened with USER, amplified with
KAPA HiFi HotStart for 22 cycles using NEBNext i5/i7 primers, and sequenced 150-bp
paired-end on a MiSeq (index read lengths not stated).

**Model corrections (round 2026-10-17-breadth-first-08):** the graph previously went
straight from stem-loop ligation to USER opening; the Lambda Exonuclease/Exonuclease I
selection step and its digested waste are now represented, the named enzymes, buffers and
concentrations above are declared, the oSQT1288 bases are an exact input state checked by
the source-transcription guard, and the sequencing plan records the 150-bp paired-end
reads with the index lengths left unstated.

🟢 The post-cleavage adapter is named in the supplementary protocol as NEBNext Adaptor for
Illumina (#E7601A) and printed there in full. 🟡 The paper prints 64 nt, with
`ACACTCTTTCCTACACGACG`; the later NEB E7600 Set 1 manual prints 65 nt, with
`ACACTCTTTCCCTACACGACG`. The model uses the later NEB sequence in `lib/illumina.py`
because the paper names that product. NEB's manual revision history records a correction
to the adapter sequence in November 2017, after this paper, so the physical oligo used
in the study remains unverified. `SOURCE_OLIGO_EXCEPTIONS` in `tools/crispr_circle.py`
records why the model cannot match the paper's printed oligo. The differing base falls
in the loop rather than the 12-base stem, so duplex-pairing checks cannot resolve it.

🟢 The paper names the NEBNext Multiplex Oligos dual-index kit by catalogue number only,
so the i5/i7 primer bases come from the NEB kit rather than from this paper.

## Disambiguation

This method creates circles from ordinary genomic fragments to select Cas9 cleavage.
eccDNA Circle-Seq starts with naturally circular molecules and uses exonuclease selection
plus phi29 amplification. The shared word “Circle” does not make them one protocol.

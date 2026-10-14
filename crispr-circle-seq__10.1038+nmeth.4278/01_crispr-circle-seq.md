# CRISPR CIRCLE-seq — research notes

## Assay target

The sequenced material is circularized human genomic DNA linearized in vitro at Cas9 cleavage sites. This protocol does not target, amplify or reconstruct a pathogen genome: every oligo it uses is a host-genome or library oligo (adapter duplexes and library PCR primers). The note is here so that the protocol's CRISPR, off-target and double-strand-break vocabulary is not read as pathogen work.

## Sources read

- 🟢 Tsai et al. 2017, DOI [10.1038/nmeth.4278](https://doi.org/10.1038/nmeth.4278), PMC5924695, main text, Methods and supplementary protocol.

## Molecular order

🟢 Genomic DNA is sheared to about 300 bp, repaired, dA-tailed and ligated to the
uracil-containing oSQT1288 stem-loop. Exonucleases remove open molecules; USER plus PNK
opens the loop and dilute T4 ligation circularizes individual fragments. Plasmid-Safe
DNase removes residual linear DNA. Cas9–guide complexes linearize target-bearing circles.
The fresh ends are dA-tailed, ligated to a USER-openable hairpin adapter, opened, PCR
amplified and sequenced 2×150.

🟢 The post-cleavage adapter is named in the supplementary protocol as NEBNext Adaptor for
Illumina (#E7601A) and printed there in full. 🟡 Its printed transcription is 64 nt where
the ordered oligo is 65: `ACACTCTTTCC`**`C`**`TACACGACG` loses one C at position 43. The
page therefore draws the correct 65-base form, transcribed in `lib/illumina.py` from the
NEB E7780 manual, and `SOURCE_OLIGO_EXCEPTIONS` in `tools/crispr_circle.py` records why a
verbatim match against this paper fails for the right sequence. The dropped base falls in
the loop rather than the 12-base stem, so duplex-pairing checks cannot catch it.

🟢 The paper names the NEBNext Multiplex Oligos dual-index kit by catalogue number only,
so the i5/i7 primer bases come from the NEB kit rather than from this paper.

## Disambiguation

This method creates circles from ordinary genomic fragments to select Cas9 cleavage.
eccDNA Circle-Seq starts with naturally circular molecules and uses exonuclease selection
plus phi29 amplification. The shared word “Circle” does not make them one protocol.

# TrAC-looping — source notes

Defining paper: Lai et al., *Nature Methods* (2018), [doi:10.1038/s41592-018-0107-y](https://doi.org/10.1038/s41592-018-0107-y).

## Sources read

- 🟢 PMC full text, Methods and Fig. 1 workflow.
- 🟢 Supplementary Fig. 1 and Supplementary Table 2 oligonucleotide list.

## Chemistry boundary

- 🟢 Two 19-bp mosaic ends separated by a short spacer favor a tetrameric Tn5 complex; insertion in trans bridges two interacting accessible regions.
- 🟢 Bio67F and 67bpR form the 67-bp biotinylated bivalent linker. The table prints both strands and their modifications.
- 🟢 Supplementary Table 2 also prints `ME3ddCTPtop` (`/5phos/…/3ddCTP/`) and `5NH2MEbottom` (`/5AmC6/…`) for transposome assembly; they were absent from the ordering list. The same table gives named N5xx/N7xx indexing primers. The page now shows exact N501 and N701 examples; its other listed index members can be selected from the table.
- 🟢 A four-base cutter fragments DNA; biotinylated products are captured, circularized with T7 DNA ligase, amplified by TempliPhi RCA, and amplified with N501/N701 primers.
- 🟡 The page uses the canonical Nextera sequencing geometry because the printed N501/N701 primers contain those exact canonical components. The genomic order inside an RCA product is shown schematically rather than as one asserted circle traversal.

## Bounded source pass — 2026-10-17

Read published PMC7212307 JATS Methods131–163: transposase/transposome
context, assembly, chromatography and the complete TrAC-looping preparation.
Read published MOESM1 Supplementary Table2 text570–607: modified half-linker,
Bio67F/67bpR and indexed N5/N7 oligo rows. Source manifest identifies the
archived article and publisher supplement. Ran the full bundle oligo scraper
and verified used hits at those original passages. Gene-targeting capture
probes after that table region were not followed in this chemistry pass.

🟢 The printed Bio67F and 67bpR are reverse complements. The old shared
constructor repeated ME_RC at both ends, so its assembled sequence differed
from Bio67F despite passing a 67-base length check. The local declaration
now assembles ME_RC–source-derived spacer–ME and enforces equality to Bio67F.
The input duplex preserves both 5′ phosphates and locates internal biotin-dT
at nucleotide 35. This is a repository correction, not a publication error.

🟢 Post-transposition DNA is recovered after quench/protease/crosslink reversal,
the transposition gaps are repaired with T4 DNA polymerase, and AMPure-size
recovery precedes a split into separate NlaIII and MluCI digestions. The graph
now traces the NlaIII aliquot; the parallel MluCI path is explicit but not
reconstructed. Capture-eluted DNA is recovered before dilute T7 ligation,
then TempliPhi amplification precedes Phusion HF PCR with N501/N701 for 11
cycles and 250–650 bp E-gel recovery. The graph now includes these boundaries.
Release/recovery transitions remain warned summaries rather than exact
molecular mechanisms.

🟢 The endpoint and PCR inputs now derive from source N501/N701 indexes
and canonical arms, rather than a generic Nextera index shell. Source-listed
primer bases do not prove the heterogeneous RCA-to-PCR molecular relation.
Paired-end HiSeq 2500 is stated; no acquired cycle/index recipe is assigned
from canonical primer presence alone. Candidate non-contact Tn5 products have
unresolved fate rather than being presumed physically discarded at insertion.

Deferred: selected detailed source sections outside the read scope, linked
predecessor/vendor originals, complete non-transfer half-linker/Tn5 assembly
and insertion geometries, MluCI path, exact restriction/circle/RCA traversal,
individual recoveries and acquired run mask. Existing generic circle flanks
are schematic. This capped pass does not validate full experimental chemistry
or accept any missing-information claim.

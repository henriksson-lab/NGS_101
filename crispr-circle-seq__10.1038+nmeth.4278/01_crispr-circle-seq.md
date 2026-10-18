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
`ACACTCTTTCCCTACACGACG`. The migrated model now uses the defining supplementary protocol’s printed 64-nt
reference and keeps the actual historical reagent lot unresolved. NEB's manual revision history records a correction
to the adapter sequence in November 2017, after this paper, so the physical oligo used
in the study remains unverified. Both source values and the current source-printed reference choice are displayed
on the public page. The differing base falls
in the loop rather than the 12-base stem, so duplex-pairing checks cannot resolve it.

🟢 The paper names the NEBNext Multiplex Oligos dual-index kit by catalogue number only,
so historical physical i5/i7 primer sequences and the selected index pair remain
unassigned in this model; no generic eight-base dual-index endpoint is inserted.

## Disambiguation

This method creates circles from ordinary genomic fragments to select Cas9 cleavage.
eccDNA Circle-Seq starts with naturally circular molecules and uses exonuclease selection
plus phi29 amplification. The shared word “Circle” does not make them one protocol.

## Source-owned validation migration (2026-10-18)

The generic `batch_ngs.truseq_library` endpoint is removed. The final library now
records source-named NEBNext i5/i7-derived ends without assigning historical
primer sequences, index lengths or exact primer sites. The insert path follows
one Cas9-created end through the original circle-closure junction to the other
Cas9-created end; the former invented six-position internal cut-site payload is
removed. MiSeq Read1/Read2 acquisition remains 150 bp each; acquired index cycles
are unassigned.

The native graph distinguishes un-opened stem-loop-capped molecules, USER/PNK-opened
duplex substrates, two complementary closed circle strands, Cas9-opened circles,
un-opened post-cleavage hairpin caps, and the actual USER-opened PCR substrate.
It no longer calls a single-strand circularization constructor on duplex DNA, uses
an invented eight-position stem-loop scar, or uses a final amplified library as
the pre-PCR substrate. The final 0.7× AMPure cleanup is explicit. Uncut circles have
an unresolved physical fate; losing adapter-compatible ends does not physically
remove those molecules.

Both modified adapters are transcribed directly from the supplementary protocol,
including /5Phos/, /ideoxyU/ and the terminal C*T linkage. Their folded stem
relations are checked locally. The E7601A reference is the defining protocol’s
64-nt printed form; the previous later 65-nt vendor choice is no longer treated as
the physical reagent. The known published discrepancy stays visible. Complete
historical library-adapter ligation and PCR geometry remain unavailable. Duplex
states show a representative strand with their base-level complement geometry
unassigned, rather than manufacturing literal base pairs for role labels.

This capped pass read two originals: PMC5924695 JATS Methods lines 115–116 and the
published Supplementary Protocol MOESM84, pp.1–6 (modified oligos, end preparation,
exonucleases, USER/PNK, circularization, Cas9, adapter ligation, USER before PCR,
22-cycle PCR and final cleanup). Mandatory bundle scraper pointers were followed
into those actual source passages. The later NEB manual was not re-read; its
existing discrepancy evidence is retained. Other intermediate SPRI cleanups,
optimization supplements, historical kit lot, selected indexes and run metadata
remain deferred. This is not a full chemistry review or an absence conclusion.

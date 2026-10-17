# RAMPAGE

## 1. What it is

RAMPAGE selects 5′-complete capped transcripts and uses paired-end sequencing to connect a base-resolution promoter start to downstream transcript structure.

## 2. Sources and why this CAGE-class representative

- 🟢 Batut P et al. *High-fidelity promoter profiling reveals widespread alternative promoter usage and transposon-driven developmental gene expression.* Genome Research (2013). [doi:10.1101/gr.139618.112](https://doi.org/10.1101/gr.139618.112). Read: PMC3530677 JATS Methods "Library preparation and sequencing", "Sequencing reads alignment" and duplicate removal (`.xml.txt` lines 141–153).
- 🟢 Batut P & Gingeras TR. *RAMPAGE: promoter activity profiling by paired-end sequencing of 5′-complete cDNAs.* Current Protocols in Molecular Biology 104:25B.11 (2013). [doi:10.1002/0471142727.mb25b11s104](https://doi.org/10.1002/0471142727.mb25b11s104). Read: PMC4372803 author-manuscript JATS, Synthetic oligonucleotides, Sequence barcodes and Basic Protocol steps 1–42 (lines 186–735).
- 🟢 The typeset Current Protocols article (Supplement 104, pp. 25B.11.1–25B.11.16) is reproduced as Appendix 1 of [Philippe Batut's CSHL thesis](https://repository.cshl.edu/id/eprint/31200/1/Batut%20Thesis.pdf) (fetched 2026-10-17). Read: PDF pp. 122–124 (article pp. 25B.11.3–25B.11.5): the materials list and **Table 25B.11.1**, which prints all 40 TSOs rampage_TS_01–40 with their bold library barcodes. The rendered page was inspected, because the text extraction scrambles the name/sequence order. This is an author reproduction of the published article; the publisher's own Supplementary File 1 (PMC `NIHMS665960-supplement-Supp.pdf`) is still behind PMC's download gate.
- 🟡 RAMPAGE is represented rather than a “CAGE family” page because it has a protocol-specific defining paper, a complete open molecular protocol, and a distinct paired-end construct.
- 🟢 A later [Gingeras Lab ENCODE library-production document](https://www.encodeproject.org/documents/db4fd5dc-145b-44a7-b3ea-0cddb5b458cb/@@download/attachment/296007.pdf) (pp. 9–10) prints four sets of five indexes labelled sets 1, 2, 3 and 5 while its reaction text calls for set 4. It is a later implementation, not the authority for this page; the ENCODE sample's barcode remains unresolved there.

## 3. Oligos (verbatim, Current Protocols; identical in the typeset materials list)

- 🟢 rampage_TS_** (TSO): `5′-TAGTCGAACTGAAGGTCTCCAGCANNNNNNrGrGrG` — "The 6 N's … correspond to the library identification barcodes"; "The last 3 residues must be riboguanosines". Ordered as 250 nmol RNA oligo batches.
- 🟢 rampage_RT: `5′-TAGTCGAACGAAGGTCTCCGAACCGCTCTTCCGATCT(N)15`
- 🟢 rampage_F: `5′-AATGATACGGCGACCACCGAGATCTACACTAGTCGAACTGAAGG` (the typeset materials list names the same qPCR primer CAGEscan-erF)
- 🟢 rampage_R: `5′-CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT` (typeset: CAGEscan-erR)
- 🟢 Custom Read 1 primer rampage_r1: `5′-TAGTCGAACTGAAGGTCTCCAGCA`; Read 2: Illumina SBS8 `5′-CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT`.
- 🟢 Barcodes: {{= len(rampage.TSO_BARCODES) }} six-mers, designed with GC 20–80% and pairwise Hamming distance ≥3 (computed minimum for the transcribed set: {{= rampage.min_pairwise_hamming() }}). Table 25B.11.1, in order 01–40: {{= ", ".join(rampage.TSO_BARCODES) }}.

### How they interlock (computed)

- 🟡 rampage_F is P5 followed by the first {{= len(rampage.PCR_F_TAIL_3P) }} bases of the TSO handle, so it primes on the TSO complement at the 3′ end of first-strand cDNA.
- 🟡 rampage_R is P7 followed by the complete SBS8 site, but only its 3′ {{= rampage.RT_R_OVERLAP }} bases repeat the rampage_RT handle. The RT handle's first {{= len(rampage.RT_HANDLE) - rampage.RT_R_OVERLAP }} bases are therefore replaced by the PCR primer tail, and the SBS8 Read 2 site exists only after PCR. (Model correction 2026-10-17: the earlier model drew the SBS8 site itself as the RT-primer handle.)

## 4. Step by step (Current Protocols Basic Protocol)

1. 🟢 Start from 5 µg DNase I-treated total RNA (2–3 µg per sample when pooling after RT). Terminator (TEX) digest of 5′-monophosphate RNAs, 30 °C 90 min; RNAClean XP. The protocol's troubleshooting ties residual rRNA to an inefficient TEX digest. (Earlier model correction: the input is not separately ribosome-depleted, and the exonuclease acts before RT, not after cap trapping.)
2. 🟢 RT and template switching in one reaction: rampage_RT (400 µM, 1 µl) and rampage_TS_** (4 mM, 1 µl), 65 °C denaturation, SuperScript III with sorbitol/trehalose and betaine; 22 °C 1 min, 42 °C 30 min, 75 °C 15 min. RNAClean XP.
3. 🟢 Optional: qPCR quantification with rampage_F/R and equimolar pooling of barcoded samples, then RNAClean XP. The defining paper pooled all samples of a library right after RT.
4. 🟢 Cap oxidation: NaIO4 in NaOAc pH 4.5, 45 min on ice in the dark (oxidizes cap riboses and 3′-terminal riboses); glycerol and Tris stop; RNAClean XP.
5. 🟢 Biotinylation with biotin hydrazide (long arm), sodium citrate pH 6.0, 14–15 h at room temperature; no cleanup.
6. 🟢 RNase I, 37 °C 30 min (up to 60 min when pooling), 65 °C 5 min; RNAClean XP. The typeset figure legend: RNase I removes single-stranded RNA, leaving biotin only at the 5′ ends of capped transcripts hybridized to 5′-complete cDNAs.
7. 🟢 Streptavidin pulldown on tRNA-blocked MPG beads; washes with buffers 1–4; supernatants discarded.
8. 🟢 Elution with 50 mM NaOH, 10 min, neutralized into Tris-HCl pH 7.0; AMPure XP.
9. 🟢 PCR with rampage_F and rampage_R, Ex Taq HS, 16 cycles.
10. 🟢 Double-sided AMPure XP size selection: 0.52× (discard beads with large inserts), then 0.70× with enriched beads (keep). Expected library ~300–1,000 bp.

## 5. Final library and reads

- 🟢 Paired-end run; Read 1 uses rampage_r1, Read 2 uses SBS8. The defining study sequenced 2×101 on HiSeq (embryo, larva/pupa libraries) and 2×76 on GAIIx (adults). The model draws the 2×101 case.
- 🟢 The first 6 bases of Read 1 are the library barcode; a stretch of 3–4 G follows (troubleshooting note on low pass-filter rates). The barcode and RT-primer sequences were trimmed before mapping.
- 🟢 The random 15-mer at the start of Read 2 often primes with mismatches. It was used together with full alignment coordinates as a pseudo-random barcode for duplicate collapsing. It is modelled as a degenerate UMI-like feature, not a designed UMI.

## 6. Open questions

- 🔴 Which of the 40 barcodes each sample of the defining study used. Genome Research Supplemental Methods/Table S1 (genome.cshlp.org) were not obtained; a supplementary URL guess returned 404. The Europe PMC supplement zip contains only figure images.
- 🔴 The publisher's Supplementary File 1 was not read directly. The barcode list comes from the author's thesis reproduction of the typeset Table 25B.11.1, and its 40 strings match the TSO barcodes the thesis prints.
- 🟡 PCR duplex geometry is not reconstructed: the page shows the PCR product as an endpoint.

# PicoPLEX DNA-seq

## Defining source

🟢 Rubicon Genomics, *PicoPLEX™ DNA-seq Instruction Manual — Single Cell Library Preparation for Illumina® NGS Platforms*, QAM-106-002 (Cat# R300381, 48 reactions). The PDF is hosted by Takara Bio, which acquired Rubicon (file name `PicoPLEX-DNA-seq-Manual-QAM-106-002-updated(RB3381).pdf`), but every page of the manual is footed QAM-106-002 and carries Rubicon branding. SurePlex is a related commercial WGA implementation and is not merged into this page.

Sections read (2026-10-17, text twin plus rendered pp. 3, 7 and 17): Product Description, Kit Contents, Principle and Figure 1, Workflow, Getting Started, Starting Material, Dual Index Plate and Figure 4, Cell Lysis, Pre-Amplification and Library Amplification protocols, Library Processing (pooling, AMPure XP), Library Quantification, Sequencing Recommendations, Analysis of the Sequencing Reads, and Appendix 2 with Figure 7. The oligo scraper found no hits because the only printed sequences are 8-nt indexes.

## Molecular workflow

🟢 **Lysis.** One to ten cells (no more than 2.5 µL of carry-over PBS), or 6–60 pg of purified DNA, are made up to 5 µL with Cell Extraction Buffer. 5 µL of Extraction Enzyme Dilution Buffer/Cell Extraction Enzyme master mix is added, and the reaction is incubated at 75 °C for 10 min, 95 °C for 4 min, then held at 22 °C. The positive control is 15 pg of reference DNA; the negative control is 2.5 µL of PBS or TE.

🟢 **Pre-amplification**, in the same tube. 5 µL of Pre-Amp Buffer/Enzyme master mix is added (15 µL total). Cycling: 95 °C for 2 min, then 12 cycles of 95 °C 15 s, 15 °C 50 s, 25 °C 40 s, 35 °C 30 s, 65 °C 40 s and 75 °C 40 s, then a 4 °C hold. The manual says that "proprietary quasi-random primers bind to selective sites on the gDNA and, through a linear amplification, create a highly reproducible library". Figure 1 draws primers with a common tail, and an "accumulating hairpin library".

🟡 The manual does not name the Pre-Amp Enzyme or say that it displaces strands. The catalogue description "displacement DOP-PCR" therefore does not come from this manual, and earlier versions of this note were wrong to mark it 🟢.

🟢 **Library amplification**, in the same tube. 30 µL of Amplification Buffer/Enzyme master mix is added, with EvaGreen/fluorescein if amplification is monitored, followed by 5 µL from one well of the Dual Index Plate (50 µL total). Cycling: 95 °C for 4 min; 4 cycles of 95 °C 20 s, 63 °C 25 s and 72 °C 40 s; 7 cycles of 95 °C 20 s and 72 °C 55 s (fluorescence is read during these cycles); then a 4 °C hold.

🟢 **Dual Index Plate** (Figure 4 on p. 7, repeated as Figure 7 on p. 17; checked against the rendered figure). Rows A–H hold i5 D501–D508 and columns 1–6 hold i7 D701–D706, giving 48 pairs; wells A7–H12 are blank. The sequences below are transcribed in `picoplex_dna_seq.I5_INDEXES` / `I7_INDEXES`, and the self-test checks them against the source text:

- i5: {{= ", ".join(f"{k} {v}" for k, v in picoplex_dna_seq.I5_INDEXES.items()) }}
- i7: {{= ", ".join(f"{k} {v}" for k, v in picoplex_dna_seq.I7_INDEXES.items()) }}

🔴 The manual does not say whether these sequences are written as the oligo strand or as the base calls. It credits the sequences to Illumina ("Oligonucleotide sequences © 2007-2012 Illumina"), but this note does not take orientation from any Illumina document.

🟢 **Pooling and cleanup.** 10 µL of each library is pooled, and 100 µL of the pool is purified once with 100 µL of AMPure XP (1:1). The beads are washed twice with 300 µL of 80% ethanol, dried at 37 °C and eluted in 50 µL of TE. The manual warns against QIAquick and other silica columns because they leave primers behind. Libraries are quantified by KAPA qPCR, calculated at 300 bp, or on a Bioanalyzer high-sensitivity chip.

## Final library and reads

🟢 The library spans "~300 to 1000 bp total size (~200 to 900 bp insert size)". Sequencing uses standard Illumina reagents for multiplexed libraries; a suggested starting load is 16 pM on MiSeq v3 with at least 5% PhiX. The manual sets no read lengths. Indexes are 8 nt.

🟢 "The first 11 cycles of each read will contain quasi-random bases … either trim the initial 14 bases from each read or begin … at base position 15." The page shows an {{= picoplex_dna_seq.QUASI_RANDOM_CYCLES }}-nt quasi-random primer-derived region at each end of the insert. The 11 versus 14 difference appears on the public page as a discrepancy within the manual. A second wording discrepancy is also shown there: purified DNA input is "less than 6 pg to 60 pg" on p. 6 but "6-60 pg" on p. 8.

🔴 The quasi-random primers, their common tail, and the complete amplification/index oligos and P5/P7 ends are proprietary. The page therefore shows read roles and index lengths but cannot compute where the sequencing primers land.

## Open questions

- Pre-Amp Enzyme identity and mechanism, including whether it displaces strands (not stated in the manual).
- Common-tail sequence and how it relates to the read-one and read-two handles. Figure 1 suggests the same tail at both ends.
- Orientation of the printed i5/i7 sequences.
- Possible follow-ups (not read in this pass): US Patent 8,206,913, cited on the back page; the Illumina TruSeq Sample Preparation Pooling Guide 15042173 Rev A, cited for index compatibility; and the Takara-era revisions of this manual.

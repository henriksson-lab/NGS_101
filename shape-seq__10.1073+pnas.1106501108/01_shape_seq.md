# SHAPE-Seq

**Defining source:** Lucks *et al.*, *PNAS* (2011),
doi:10.1073/pnas.1106501108.

## Evidence used

- 🟢 Each designed RNA carries a DNA barcode at its 3′ end before pooling.
- 🟢 SHAPE chemistry acylates flexible 2′ hydroxyls; primer-extension stops encode the
  modification sites relative to a minus-reagent control.
- 🟢 A second adapter is ligated to cDNA and paired-end sequencing connects each stop
  with the barcode identifying the RNA template.
- 🟢 The original supplementary methods were recovered from the CaltechAUTHORS mirror.
  They print the v1 A_adapter_t, A_adapter_b, treated/control handle patterns and PCR
  construction; the shared model now uses those sequences.

## Published adapter and primer architecture

🟢 `A_adapter_t` / `PE_forward` is
`AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT`.
It is the P5 arm followed by the Read 1 site and forms the 5′ tail of each RT primer.

🟢 The RT primers then carry a four-base treatment handle—`RRRY` for the 1M7-treated
branch or `YYYR` for the solvent control—followed by a target-specific RT-primer sequence.
There is therefore no single universal terminal RT-primer sequence to assert.

🟢 The single-stranded ligation oligo `A_adapter_b` is printed as
`AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG`
with a 3′ C3 spacer, which prevents concatemerization. 🟡 A 5′ phosphate is chemically
required for ligation to the cDNA, but the SHAPE-Seq SI does not print that modification
on this named oligo; it should not be presented as source-verbatim.

🟢 The cited Kozarewa *et al.* adapter supplement prints its own `A_adapter_b`
beginning `GATC`, while the SHAPE-Seq SI prints `AGATC`. The SHAPE-Seq model follows
the defining SI's complete oligo, not the predecessor's adapter. Kozarewa's
supplement also marks a terminal phosphorothioate on its `A_adapter_t`; the
SHAPE-Seq SI does not mark one on its RT-primer tail.

🟢 The defining experiment used 2×50 cycles on an Illumina Genome Analyzer IIx.
Read 1 begins with the treatment handle and RNA-identity end; the paired read reports the
SHAPE-dependent primer-extension stop from the opposite end.

🟢 Illumina's [11 January 2011 adapter sequence letter](https://www.seqanswers.com/filedata/fetch?id=303482),
p. 6, prints the paired-end Read 1 primer (the `A_adapter_t` Read 1 site) and Read 2
primer `CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT`. The latter reverse-complements
exactly to the first 37 bases of `A_adapter_b`. These are source-supported binding
sites for the contemporary paired-end primer family; the paper does not state the
specific reagent lot used for its sequencing run.

## RNA bar-code position and read layout (2026-10-17 bounded pass)

🟢 SI "Synthesis of … RNAs" (SI p.1) states that every RNA except the non-bar-coded
RNase P RNA carries "a four nucleotide degenerate bar code between the reverse
transcriptase primer binding site and the 3′ linker region within the 3′ flanking
structure cassette". Fig. S1 draws the bar code between the two cassette hairpins, and
Table S1 shows each RNA ending `…GAUCCGCUUCGGCGGAUCC[NNNN]AAAUCGGGCUUCGGUCCGGUUC-3′`.
The RT primer therefore extends directly into the bar code.

🟡 Read 1 consequently reads the treatment handle, the target-specific RT-primer bases,
then the cDNA copy of the RNA bar code; the SHAPE stop lies at the A_adapter_b junction
and is reported by Read 2. This agrees with the main text ("the SHAPE modification
position and the RNA identity (bar code)—are on opposite ends of the cDNA molecules").
**Model correction:** the earlier schematic drew the bar code beside A_adapter_b.

🔴 The target-specific RT-primer bases are not printed in the SHAPE-Seq SI; the primer
extension procedure is cited to Mortimer & Weeks 2009 (*Nat Protoc* 4:1413), not read
in this pass. The PCR primers are described only as binding "the Illumina adapter
sequences" (SI p.1, citing Ingolia 2009); their bases are not printed there.

🟢 Other SI p.1 steps now in the model: 1M7 at 6.5 mM with a solvent (−) control; RT
with 9 pmol primer after a 95 °C step; RNA hydrolysis with 1 µL 4 M NaOH, 95 °C, 5 min;
CircLigase ligation of A_adapter_b (5 µM, 68 °C, 6 h, separate (+) and (−) reactions,
80 °C 10 min stop); AMPure XP removal of excess adapter; 9 or 12 PCR cycles (the main
text reports nine) and AMPure cleanup; no size selection.

**Model correction:** `A_adapter_t` was assembled as `P5 + Read 1 site`, duplicating the
shared `ACAC` (62 nt). It is now the 58-nt SI sequence.

## Construct logic

The original protocol is a multiplexed in-vitro accessibility assay for defined RNA
pools. It is distinct from SHAPE-MaP, where modifications are read as mutations rather
than termination sites.

# Digital Telomere Measurement (DTM)

## 1. What it is

DTM ligates barcoded capture/tether duplexes directly to native human telomeres, fills the terminal gap and measures individual telomeres by Oxford Nanopore sequencing.

## 2. Sources read

- 🟢 Sanchez et al. (2024), [doi:10.1038/s41467-024-49007-4](https://doi.org/10.1038/s41467-024-49007-4), full text and methods.
- 🟢 Supplementary Data 3 (MOESM6 XLSX; capture oligo and seqTether sequences) and Supplementary-file descriptions (MOESM3).
- 🟡 Makarova et al. (2025), [doi:10.3389/fmolb.2025.1725112](https://doi.org/10.3389/fmolb.2025.1725112), workflow comparison.

## 3. Construct and reaction order

Methods, "Telomere capture sequencing nanopore library preparation" (`s41467-024-49007-4_PMC11189511.xml.txt:158`):

1. 🟢 Equimolar barcoded capture oligo and seqTether are duplexed in low TE (95 °C 2 min, then 1 h at room temperature).
2. 🟢 About 3 µg HMW genomic DNA is ligated to freshly duplexed capture oligos (rCutSmart, 25 pmol duplex, 4000 U T4 DNA ligase, 0.1 mM ATP added) overnight at 37 °C, then heat-inactivated at 65 °C for 10 min.
3. 🟢 In the same tube, Sulfolobus DNA polymerase IV (4 U, ThermoPol, dNTPs, ATP; 56 °C 2 min then 72 °C 15 min) fills "potential gaps between the capture oligo and the double-strand/single-strand junction".
4. 🟢 ProNex beads (1.6×), two 80% ethanol washes, elution in 60 µL ONT elution buffer give the "post-capture sequencing library".
5. 🟢 The whole post-capture library is ligated to the ONT sequencing adapter (AMII for R9, NA for R10) with Quick T4 DNA ligase for 1 h at room temperature; the rest follows the whole-genome LSK-110 workflow. Multiplexed runs pool up to eight differently barcoded post-capture libraries (1000 ng total) before adapter ligation (not drawn).

🟡 Which nick T4 ligase seals is not stated. By polarity, the capture oligo annealed to the G-rich overhang has its 3′ end (CCTAACC) facing the C-rich 5′ terminus and its 5′-phosphorylated end distal, paired with seqTether. Fig. 1a draws the CCAATCC (3′→5′) arm abutting the C-rich strand and seqTether opposite the G-strand terminus, but is schematic. Before the 2026-10-17 bounded pass the page claimed that ligase "joins the capture 5′ phosphate at the native C-rich terminus"; that is geometrically impossible and was a model error, now removed.

🟢 Supplementary Data 3 gives `seqTether` as `AACCTTGGAGATGCACGGAGCAAGCAAT` and the representative 5′-phosphorylated p-t3-nb01 capture oligo as `TGCTCCGTGCATCTCCAAGGTTCACAAAGACACCGACAACTTTCTTCCTAACC`. 🟡 The table does not partition the oligos. Comparing its rows gives a 22-nt tether-complementary arm (shared prefix of all eight canonical oligos; its reverse complement is the first 22 nt of seqTether), a 24-nt barcode and a 7-nt `CCTAACC` telomere arm (shared suffix of the canonical oligos; the seven TSQ variants share the first 46 nt with p-t3-nb01 and differ only in the last 7). The model previously split this as 25 + 6; that was a repository error, corrected 2026-10-17, and `tools/selftest.py` now derives the partition from the table. Eight canonical-telomere barcode variants and seven telomere-sequence-variant capture oligos are provided; the page uses p-t3-nb01 as the representative molecule.

🟢 The core Methods used high-molecular-weight DNA directly. The paper discusses later/optional restriction-digest refinements, but the page no longer invents a mechanical-fragmentation step.

## 4. Final library and sequencing

A motor-loaded nanopore adapter feeds an intact tagged strand through the pore; no sequencing primer binds.

🔴 The defining Methods identify the historical adapters as AMII for R9.4.1 and NA for R10.4, but neither the paper nor its supplementary files print their complete motor-bearing assemblies. ONT's current chemistry document now prints current LA/NA oligonucleotide strands; it does not establish that those current sequences are the same reagents used with LSK-110 in this study. The page therefore does not backfill a later adapter sequence into the historical library.

🔴 The Methods also link the authors' [Telometer repository](https://github.com/santiago-es/Telometer). Its README points to `TelometerLibraryPrep.docx` for a current bench protocol. That linked DOCX could not be downloaded during the source review, so its contents and relationship to the 2024 library remain unchecked.

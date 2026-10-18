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

🟢 The Methods also link the authors' [Telometer repository](https://github.com/santiago-es/Telometer). The now-cached `gh_TelometerLibraryPrep.docx.txt` identifies itself at line 2 as “04/2025”. This is a later author protocol, rather than the defining 2024 preparation. This migration inspected its identity only; its digestion, polymerase and kit choices are excluded from the 2024 graph. The earlier failed-download statements in historical review passes describe those passes, not its present availability.

## 5. Bounded native graph migration (2026-10-18, batch 12)

Inspected three originals: published Methods `s41467-024-49007-4_PMC11189511.xml.txt:150–160` (whole-genome workflow incorporated by reference and capture-specific overrides), Supplementary Data 3 `s41467-024-49007-4_supp_41467_2024_49007_MOESM6_ESM.xlsx.txt:1–23` (oligos and modifications), and the identity header of the April 2025 author DOCX `gh_TelometerLibraryPrep.docx.txt:2`. Ran the mandatory bundle oligo scraper and checked its published-table hits directly. Superseded preprint chemistry and later author-protocol chemistry did not supply this graph.

The protocol now owns its native graph. Capture and seqTether enter as independent source-exact polymers, then form a computed antiparallel duplex. The capture/tether material is allocated once to the chromosome-end hybrid; the following ligation consumes that complex, rather than consuming the same capture aliquot again. The capture 5′ phosphate is an input modification, not a base copied onto a new strand. The sample barcode has structured sample-index identity. The drawn telomere lengths and repeat register are schematic, not measured properties of a native chromosome end.

The inspected Methods do not resolve the sealed capture nick sufficiently for this model to construct its covalent product. Post-ligation and post-gap-fill states therefore carry explicit unresolved assembly tokens instead of claiming the unchanged pre-ligation chromosome is the product or choosing a covalent connection. The gap-fill tract and historical AMII/NA motor geometry remain unassigned. These are model limits of this bounded reading, not accepted source-absence claims. No PCR or amplification is introduced.

Both source-described ProNex cleanups are explicit: 1.6× after gap filling and 1.1× after motor-adapter ligation, the latter following the whole-genome workflow incorporated by reference. Each keeps the eluate and physically discards supernatant/washes. Cleanup does not claim complete removal of untagged chromosomes or any other specific molecular class. No computational read rejection is modelled as physical removal. The complete 60-µL post-capture eluate enters adapter ligation. No aliquot split or sample pool is invented.

Deferred: the up-to-eight-sample 1000-ng multiplex pool; native end heterogeneity and sealed capture connectivity; exact Pol-IV tract; historical AMII/NA assembly; and full reading of the later April 2025 author protocol as a separate variant. Existing blocked missing-information status is unchanged. This capped migration is not a complete protocol source review.

### Code review corrections (2026-10-18; no additional source reading)

The reviewed source inventory is now an explicit curator-written list rather than automatic enumeration of graph nodes. Post-ligation, gap-filled and motor-ligated populations retain the source barcode as an **unplaced structured segment** with `capture_sample` sample-index identity. This metadata assigns neither a covalent attachment nor a sequence coordinate within the unresolved assembly; it is not drawn as a separate free oligo. Cleanup preserves it.

The same populations retain independent **possible residual classes**: unreacted genomic DNA, partially capture-adapted DNA, free capture and free/incompletely incorporated seqTether; after motor ligation, captured DNA without motor and free/incompletely attached motor adapter are also retained. These classes are inferred possibilities, not observed products or proportions. They remain possible in bead-cleanup eluates because cleanup has not established their complete removal. Their structured candidate status does not convert a computational filter into physical removal. `NICK_SEALING` classifies the source-stated ligation within the already assembled complex only; its structured geometry explicitly declines to assign any particular sealed nick.

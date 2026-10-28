# ARM-seq: AlkB-sensitive RNA recovery

Cozen et al., *Nature Methods* (2015), [doi:10.1038/nmeth.3508](https://doi.org/10.1038/nmeth.3508).

This capped source pass establishes the pretreatment and reported sequencing readout. The local graph ends at the selected small-RNA input to the historical commercial library kit. It intentionally fails sequencing publishability. Adapter identity and kit revision remain unreviewed, not confirmed absent.

## Source reading

- Published author-manuscript article `nmeth.3508_PMC4553111.pdf.txt`: abstract 28–35; Introduction 85–108; Discussion 272–297; Methods 315–402 (AlkB purification, biological RNA preparation, AlkB treatment, separate primer-extension assay, library preparation); read processing 403–426 and differential analysis 458–465.
- Published Springer figure supplement `nmeth.3508_41592_2015_BFnmeth3508_MOESM407_ESM.pdf.txt`: Supplementary Figure 1, lines 1–13, AlkB activity and fresh cofactors; Figures 2–5 captions, lines 15–65, fragment/modification interpretation. The author-manuscript DOC counterpart starts with the same figure captions; it is not an independent chemistry method. Remaining figure panels/captions were not fully inspected.
- Published workbook `nmeth.3508_41592_2015_BFnmeth3508_MOESM408_ESM.xlsx.txt`: first sheet heading and opening rows only; this is an expression/count result table, not an inspected library oligo table. Other sheets and MOESM409/410 workbooks are deferred. Supplementary Table 5, explicitly cited for primer-extension oligos, remains uninspected.

The full bundle oligo scraper returned zero hits. This is a retrieval aid result, not evidence that the original tables disclose no oligos. No image table was declared inspected based on empty extracted text.

## Chemical distinction

🟢 E. coli AlkB removes several RT-blocking methylations, including m1A, m3C and m1G, before conventional small-RNA sequencing. Differential abundance between treated and untreated samples exposes modification-sensitive recovery. It does not identify every modification chemically or quantify a modification from an invented per-base mismatch signal. The article explicitly distinguishes ARM-seq from methods relying on RT misincorporation (Discussion 292–297).

🟢 Supplementary Figure 1 also shows that stops at modifications outside this treatment's scope can become more visible after dominant blocks are removed. No blanket claim of complete demethylation or restoration of all RNA molecules belongs in the model.

## Source-defined preparation

1. 🟢 The paper prepares bulk yeast or human RNA. Independent cultures or technical replicates remain distinct biological designs; neither creates a different protocol chemistry.
2. 🟢 Treat 50 µg bulk RNA with 50 µg AlkB in 200 µl containing 50 mM HEPES-KOH pH 8, 75 µM ferrous ammonium sulfate pH 5, 1 mM α-ketoglutarate, 2 mM sodium ascorbate and 50 µg/ml BSA, for 100 min at 37 °C. Prepare the reaction buffer fresh. Matched untreated controls receive AlkB storage buffer.
3. 🟢 Stop with EDTA/ammonium acetate, extract with phenol, ethanol-precipitate and resuspend in water.
4. 🟢 MirVana selects RNA below 200 nt; concentration and DNase I precede column cleanup. The paper gives 1 µg RNA input to NEBNext Small RNA Library Prep Kit for Illumina.
5. 🟢 E-Gel size selection retains libraries corresponding to 18–120 nt RNA inserts, followed by BioAnalyzer QC and Illumina MiSeq paired-end sequencing (Methods 376–402). Reads are adapter/index trimmed and merged before mapping; analyzed merged inserts are at least 15 nt (407–410).

The partial graph retains physically rejected extraction and small-RNA selection material as waste and ends at RNA entering the kit. Its RNA states describe a variable population, not fabricated tRNA sequences or residue positions. The graph does not assume treatment repairs native RNA ends. It does not import DM-tRNA-seq's mutant enzyme, TGIRT reaction, template-switch starter or circularization.

## Inventory and read layout

The source names a commercial kit but the applicable historical catalogue number, edition and its exact oligos have not been established in this pass. No current NEBNext Low-bias or modern small-RNA manual is substituted. Named primer-extension oligos belong to the separate radiolabeled validation assay; they do not establish library primers.

🟢 MiSeq paired-end reads are reported. 🟡 Their analyzed content is overlapping RNA-insert sequence after adapter and barcode removal. Exact cycle counts, index read(s), barcode width, native end requirements and run-primer identities are deferred. No exact molecular library or cut-and-paste sample-sheet index is exposed by this partial model.

## Ancestry candidates and follow-up

No direct ancestry relation to DM-tRNA-seq is asserted: shared AlkB use and related publication dates are insufficient, and their downstream chemistries differ.

Candidate `arm-seq__10.1038+nmeth.3508` → `pandora-seq__10.1038+s41556-021-00652-7` is plausible only for demethylation-enabled small-RNA recovery. Existing PANDORA notes support the comparison but the defining PANDORA original's citation and inheritance statement were not re-read here, so this is an inferred candidate awaiting original-source review, with no exact buffer/enzyme/kit inheritance established.

Deferred: identify and retrieve the applicable historical kit original; inspect Supplementary Table 5 and remaining linked materials; model source-specific kit reactions and molecular states; model the matched control branch; declare actual read priming when justified. No confirmed-absence claim follows from this capped pass.

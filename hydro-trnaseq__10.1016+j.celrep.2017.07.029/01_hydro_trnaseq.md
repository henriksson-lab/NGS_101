# hydro-tRNAseq: size-enriched tRNA fragments after partial hydrolysis

Gogakos et al., *Cell Reports* (2017), [doi:10.1016/j.celrep.2017.07.029](https://doi.org/10.1016/j.celrep.2017.07.029).

This bounded pass scopes hydro-tRNAseq itself. The paper's SSB/La PAR-CLIP experiments are a complementary assay, not silently incorporated into the hydrolysis library. The current source-qualified graph continues from CIP/PNK through the explicitly inherited small-RNA ligation/RT/PCR DNA-library route. Whole occupied strands and exact study-specific variants remain unassigned; no complete molecular/readout or source-absence claim is asserted.

## Sources inspected and their versions

- `j.celrep.2017.07.029_PMC5564215.xml.txt` is an **accepted author manuscript**, not the final publisher version. Inspected Introduction 46–48; Results 52–54; Discussion 110–114; Experimental Procedures 118–128; source identity of the library predecessor 244–250; Fig. 1 legend 368. Detailed facts below are qualified to this manuscript pending final-publication comparison.
- The defining supplement's six XLSX files and PDF are actual linked originals, but all seven were gated in standard retrieval. Their individual URLs and actual failure reasons are in `MANIFEST.tsv`; their original contents were not read or judged irrelevant.
- The article explicitly cites Hafner et al. 2012, [doi:10.1016/j.ymeth.2012.07.030](https://doi.org/10.1016/j.ymeth.2012.07.030), for barcoded library preparation. Its `doi_10.1016+j.ymeth.2012.07.030/j.ymeth.2012.07.030_PMC3508525.xml.txt` is also an **accepted author manuscript**. Inspected overview 34–40; adapter and primer inventory 90–148; 3′ ligation/pooling 302–326; gel selection 338–392; RT/hydrolysis/PCR 398–442. This original supports component inventory, but every downstream detail is not automatically assigned to the 2017 experiments.

Full-bundle scrapers were run for both source bundles. The incomplete hydro-tRNAseq bundle returned zero hits; the predecessor returned 15, including a collapsed twenty-adapter family. The directly relevant bases and modifications were checked at their original table locations. No negative inference follows from the defining bundle's zero hits.

## What changes relative to other tRNA methods

🟢 The manuscript enriches mature and precursor tRNAs by size, then partially hydrolyzes them so individual fragments contain less stable structure and fewer RT-blocking modifications. Basic treatment also deacylates charged tRNAs, making the native 3′ end available for ligation. This does **not** enzymatically demethylate every tRNA, erase the modifications remaining inside a fragment, or establish complete nucleotide-level coverage for each molecule.

The two sequential adapter ligations select fragments that can be copied between both handles rather than recovering every RT-truncated product (Discussion 110). Aggregate overlapping fragments inform the full tRNA and precursor annotation; a short fragment library is not an intact full-length tRNA readout. ARM-seq and DM-tRNA-seq are explicitly discussed as comparisons, not inherited molecular constructors.

## Source-qualified preparation

1. 🟢 Isolate HEK293 total RNA with TRIzol. Resolve 20 µg per sample on 12% urea-PAGE and collect 60–100 nt RNA. The local graph exposes the physically rejected size fractions.
2. 🟢 Partially hydrolyze in 15 µl containing 10 mM Na₂CO₃ and 10 mM NaHCO₃, pH 9.7, at 65 °C: 10 min for replicate 1, 1 h for replicates 2–4. These are different source instances; the local `graph_for_replicate(...)` records the corresponding duration rather than averaging them.
3. 🟢 Treat with 10 U calf intestinal phosphatase for 1 h at 37 °C in the source's NaCl/Tris/MgCl₂/DTT/carbonate buffer.
4. 🟢 Rephosphorylate with 10 U T4 PNK and 1 mM ATP for 1 h at 37 °C. This ATP-containing operation is distinct from HydraPsiSeq's ATP-free 3′ phosphatase treatment.
5. 🟢 The manuscript says 19–35 nt fragments were converted into barcoded small-RNA cDNA libraries as described by Hafner 2012 and sequenced on HiSeq 2500. It does not establish a physical 19–35 nt selection immediately before the first ligation. The cited procedure performs gel selection on adapter-bearing products; that placement remains to be reconciled to the actual hydro-tRNAseq implementation.

Hydrolysis products and any deacylation-derived molecular material are retained in the compact graph until a declared physical separation. The model does not turn a cleaved fragment or a released aminoacyl-derived species into automatic waste. It also does not assume that CIP alone proves every possible 2′,3′-cyclic/2′/3′ phosphate species has become an adapter-compatible 3′ hydroxyl: the exact post-treatment end population remains unassigned pending source-qualified modelling.

## Direct predecessor inventory, with unassigned member

🟢 The 2012 original prints a family of twenty preadenylated 3′ DNA adapters with distinct pentamer barcodes. Example **Ad01 in that source**, not an asserted choice in hydro-tRNAseq: `rAppTCACTTCGTATGCCGTCTTCTGCTTG-L`. `rApp` is the activating 5′ adenylate; `L` is the terminal aminohexyl blocking group. The actual hydro-tRNAseq adapter member remains unassigned. The local inventory uses a structured `sample_index` feature for that unspecified pentamer, with the source family named as whitelist; it supplies no invented sample-sheet index value.

🟢 The same predecessor prints RNA 5′ adapter `GUUCAGAGUUCUACAGUCCGACGAUC`, forward PCR primer `AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA`, and reverse primer `CAAGCAGAAGACGGCATACGA`. These historical sequences are not extended with contemporary P5/P7 tails. Its JATS underline marks `CCACCGACAGGTTCAGAGTTCTACAGTCCGA` within the forward primer as the sequencing primer; this is supporting-source inventory, not yet an established hydro-tRNAseq acquired-primer assignment.

The predecessor describes Rnl2(1–249)K227Q without ATP for 3′ ligation, pooling after barcode attachment, 5′ adapter ligation, SuperScript III RT, KOH RNA hydrolysis, PCR and size-marker removal. Exact calibrator use, selected barcode members, pooling design and all marker/cleanup steps in the 2017 instance remain to be checked. No generic small-RNA recipe is imported into the local graph.

## Readout, lineage and remaining work

🟢 The manuscript reports HiSeq 2500 sequencing and adapter trimming; actual acquired cycles, read configuration and instrument-primer use have not been assigned in this pass. Biological tRNA modifications may contribute mismatch signatures; that is separate from PAR-CLIP's labeling-dependent T→C signature.

Explicit component ancestry: Hafner 2012 → hydro-tRNAseq, limited to barcoded small-RNA cDNA preparation, cited in Methods 124 and duplicated references 2012a/2012b at 244–250. The two reference labels resolve to the same DOI and do not identify different protocols. The predecessor currently has no registered directory. Comparisons to ARM-seq/DM-tRNA-seq do not establish direct inheritance of their AlkB/TGIRT chemistry. Earlier yeast hydrolysis-based methods are mentioned at Introduction 48 but were not followed in this capped pass.

Proposed metadata: RNA; bulk HEK293 RNA; size enrichment; partial alkaline fragmentation/deacylation; CIP plus ATP/PNK end treatment; dual RNA-adapter/RT/PCR library with inline sample barcode; HiSeq 2500; published instance. Keep `notes` and `pending_source_review` until the source/version and downstream-library boundaries are resolved.

Deferred: final publisher article/supplement comparison, gated supplements, exact treatment-dependent RNA-end species, stage of 19–35 nt selection, actual barcode assignment, complete ligation/RT/PCR/pooling graph and acquired read plan. A second independent reader and linked-source review are still required before any accepted missing-information claim.

## Source-qualified DNA-library completion, 30 October 2026

The normal public graph now has an actual DNA-producing path: inherited RNA-adapter ligations → reverse-primer carrier retained from product recovery → reverse transcription → RNA-template hydrolysis/neutralization → DNA PCR → source-positive finishing → selected DNA sequencing contribution. Sequencing is never attached directly to the end-treated RNA. The actual RNA→cDNA/PCR-library route is explicitly stated by 2017 Methods 124 and its cited 2012 source.

**Reference applicability is inferred by constructor.** The twenty-adapter family, RNA5 adapter and historical PCR primers are exact source-verified predecessor references, not established individual 2017 orders. No Ad01, fixed pool of twenty samples, calibrator/PmeI use, exact QC/no-template control schedule, acquired cycle count or modern P5/P7 extension is invented. The source guides prescribe these control/titration/marker practices; actual study implementation beyond the common core remains unassigned and visible. No complete VOR/manuscript equivalence is asserted.

PNK ATP is a molecular phosphate donor, not merely a condition. The source 50 µl CIP→20 µl PNK interface is explicitly unassembled; no whole-stock transfer, new cleanup or aliquot quantity is fabricated. Residual end classes and material remain unresolved. The source 19–35 nt statement does not create a fictitious separate preligation PAGE step. Guide adapter-product recovery is qualified to the common inherited route; exact selection/pooling schedule and accessory-buffer additives/carryover are not claimed complete.

The full graph, uncertain reference polymers and unresolved molecular checks remain visible under an expandable section. Existing legacy selftest entry points remain; a single wrapper-free page declaration was added to the existing protocol module. No source acquisition or source-obligation exclusion occurred.

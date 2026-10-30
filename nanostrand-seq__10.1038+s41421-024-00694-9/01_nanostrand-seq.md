# NanoStrand-seq — selected GM12878 preparation

Defining paper: [Simultaneous de novo calling and phasing of genetic variants at chromosome-scale using NanoStrand-seq](https://doi.org/10.1038/s41421-024-00694-9), 2024.

Private source-qualified partial model, not a finished schematic or complete molecular/absence review. Published main76/172–194, publishedSI FigS1caption and the whole TableS1 workbook context were inspected; three-original cap closed. Historical vendor kit originals and full tagged-product geometry remain unassembled.

🟢 The selected GM12878 route uses10µM RO-3306 for4h, followed by40µM BrdU24h; the80µM mouse routes are separate. After coldPBS washing, nuclei are isolated on ice10min and stained with Hoechst33258/PI10µg/ml each. Half-Hoechst G1 nuclei are sorted into2.5µl lysis buffer. Protease1mg/ml,0.3%Triton,20mMKCl/2mMTE act50°C3h→70°C30min; optional−80°C storage follows. Unselected nuclei remain outside the selected preparation.

🟢 Source main176 specifies the loaded Tn5 strands:

```
5′-phos-CTGTCTCTTATACACATCT-NH2-3′
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
```

Canonical mosaic-end/S5 bases come from `nextera`, with the source modifications retained locally. Stocks100µM are annealed to50µM adaptor. Vazyme complex assembly uses10µl enzyme(2µg/µl),7µl adaptor and33µl coupling buffer55°C10min; selected stock dilution is1:5000. The graph models independent oligo annealing/loading,1:5000 dilution and the1µl complex-stock allocation; remaining stock has an unresolved fate. The chosen loaded-stock aliquot supplies its enzyme and adapter components through a non-transforming same-stock representation, not new free ME/S5 additions to the nucleus. Exact loading stoichiometry/companion species remain unassembled. Tagmentation55°C10min precedes SDS55°C10min, Tween20 quenching5min at room temperature and Bst2.0 gap filling65°C30min.0.8×AMPure cleanup yields9µl.

🟢 Hoechst33258 plus365nm UVA3000J/m² nicks BrdU-associated daughter strands. The model retains labeled strand/remnant classes rather than declaring quantitative physical daughter-strand destruction. This is followed by two separate single-primer extension rounds, each94°C2min→62°C5min→68°C10min, with a distinct source24nt barcode in each primer. Each round has its own fresh Gflex addition, then2µlExoI37°C30min→80°C20min and0.8×AMPure cleanup(elute10µl first,11µl second).

The published whole ordering templates are:

```
ACACTCTTTCCCTACACGACGCTC[Bar1-24]TCGTCGGCAGCGTCAGATGTG
ACACTCTTTCCCTACACGACGCTC[Bar49-64]TCGTCGGCAGCGTCAGATGTG
```

The bracket labels denote sets of barcode names, not barcode lengths. The24nt widths are established by main188 and distinct first/second codes by main190. Their identifiers are cell barcodes; no UMI meaning is assigned.

**Apparent conflict in inspected source TableS1:** sheetTableS1 has only two columns, `Primer` and `Sequence`, dimensionA1:B47; all40 sequence entries forBar1–24 andBar49–64 literally repeat `AATGATACGGCGACCACCGAGATCT`. Original XLSX cells, headers, neighbouring ordering templates and all relevant barcode rows were checked; there are no hidden sheets, associated drawings, comments or sheet external links. Other sheets are indexed as benchmarking or SV-validation tables. This does not establish that the intended codes are absent from all linked material. The model leaves actual codes unassigned and does not use the repeated string as40 distinct codes. Main176 additionally specifies3′NH2 onME while TableS1 only writes its5′phosphate; the model retains the explicit main modification.

🟢 A single amplification-primer species `ACACTCTTTCCCTACACGACGCTC` is added for25µlPCR:94°C3min,21cycles98°C15s/62°C30s/68°C8min,68°C10min. Whole symmetric primer landing geometry remains unassembled.0.8×AMPure/Qubit precede equal-amount pooling of cells, then two separate0.7×AMPure long-fragment selections. Analytical sizeQC and the~1µg ONT input allocation remain distinct from unselected physical material. SQK-LSK109 finishing and PromethION48/R9.4.1 sequencing are source-known, but historical kit end/adapter geometry is not imported from another protocol.

🟢 Main76 explicitly develops the parental-template principle from Strand-seq. This supports a scoped ancestry link for BrdU/Hoechst/UVA strand selection. NanoStrand's single-adapter Tn5, sequential tagging and long-amplification preparation differ from the original MNase/Y-adapter route; no whole Strand-seq recipe inheritance is asserted.

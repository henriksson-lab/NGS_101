# LASER-Seq

Zinshteyn et al., *Assaying RNA structure with LASER-Seq*. DOI https://doi.org/10.1093/nar/gky1172 . Published online2018-11-22; issue2019. LASER-Seq RT-stop and LASER-MaP are branches of this one method.

## Sources and chemistry

🟢 Published VOR PMC6326810 main introduction/Fig1 and Methods52–80, original published TableS1 workbook and SI FigureS1–S8 captions were inspected. Mandatory bundle scraper reviewed all three text sources and original modification strings. No preprint substitute or missing-information conclusion.

NAz photochemistry produces covalent A/G C8 adducts. The model requires explicit actual purine sites; it does not predict accessibility, site occupancy, structure, mutation frequency or a deterministic copying response. Physical RNA bases/adducts remain distinct from source-selected DNA substitutions and RT-stop classes. Conditional stop copying returns the complementary3′-side RNA interval without declaring the modified base itself converted.

The in-vitro E.coli and yeast ribosome routes and K562 in-cell route retain their distinct conditions; no DMS/SHAPE control chemistry is borrowed. E.coli NAz reaction contains2µl300mMNAz in25µl, incubates5min37°C then3minUV. Yeast uses2.5µl3MNAz/25µl and3minUV after25°C incubation. K562 uses300mMNAz and1or3minUV, then Trizol/DirectZol extraction with DNase omitted. In-vitro reactions undergo physical isopropanol/PCA/chloroform/ethanol recovery. Opaque historical kit internals do not establish a perfect molecular partition.

## Library operations and physical constructs

Main76: ZnCl2 fragmentation10mM95°C5min, ice plus immediateEDTA20mM; ethanol recovery, denature95°C3min and10%ureaPAGE size selection (100–120bp relative to non-denatured DNA ladder), gel-elution/filter/isopropanol. T4PNK37°C1h has no ATP addition in the printed mixture; selected3′OH competence is reaction-derived, not a complete terminal distribution. Adenylated linker ligation uses truncatedRnl2 37°C3h, followed by isopropanol recovery.

🟢 TableS1 oBZ191: `/preA/CACTCGGGCACCAAGGA/ddC/`, purchased with5′phosphate then preadenylated in-house. **Terminal ddC is a cytidine nucleotide plus3′block**. The local molecule contains that C; processing80's `CACTCGGGCACCAAGGAC` is consistent with this nucleotide inventory. No discrepancy is asserted.

🟢 oBZ192: `/5Phos/AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGC/iSP18/CACTCA/iSp18/TTCAGACGTGTGCTCTTCCGATCTGTCCTTGGTGCCCGAGTG`. Both C18 spacers are distinct typed non-nucleotide backbone nodes. RT adds an explicitly unassigned RNA-derived DNA payload to the same ordered primer. The pre-hydrolysis RT mixture retains independent original-RNA and possible unextended-primer/enzyme/reagent companion rows; circularization retains all nodes and removes free-end topology. Calling exact_sequence on a spacer-containing backbone raises; no spacer-erased ordinary DNA circle is supplied to PCR. Circle-to-PCR traversal/polymerase behavior remains unassembled rather than inheriting an irCLIP mechanism from primer similarity.

E.coli RNA is split to TGIRTIII60°C1h or SSII/Mn2+42°C1h; other samples use SSII only. NaOH hydrolysis1µl5M95°C3min generates RNA remnants, followed by cDNA precipitation/denaturation/ureaPAGE. MaP isolates full-length cDNA; stop libraries isolate~20nt above primer to~20nt below full-length, including untreated controls. CircLigaseI60°C2h then80°C20min precedes8–12cyclePCR, gel purification and HiSeq2500SE50.

🟢 oBZ287 forward `AATGATACGGCGACCACCGAGATCTACAC`; source reverse primer family `CAAGCAGAAGACGGCATACGAGAT[index6]GTGACTGGAGTTCAGACGTGTGCTCTTCCG`, TableS1 rows6–29 specify actual sample-index candidates. A study member is not arbitrarily selected here. Source order/acquisition are retained through the endpoint, but exact product/primer landings, index layout and circle-to-PCR spacer handling remain unresolved. No finished-page marker.

## Lineage

Methods76 explicitly adapts library preparation from references22,23,31, not a whole imported family recipe. These predecessors and NAz-synthesis reference8 remain directly linked pending originals. Similarities to other two-spacer RT primers alone establish no ancestry. Proposed mRNA rRNA-depletion/polyA-selection extension in Discussion is not inserted into actual ribosome/K562 source routes. No confirmed absence or two-reader absence audit.

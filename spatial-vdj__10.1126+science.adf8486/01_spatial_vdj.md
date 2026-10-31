# Spatial VDJ: published SR-TRB preparation

Engblom, Thrane, Lin et al., *Spatial transcriptomics of B cell and T cell receptors reveals lymphocyte clonal dynamics*, [Science382,eadf8486](https://doi.org/10.1126/science.adf8486) (2023).

## Source identity and resumed scope

🟢 The primary source is now the **final published article**, [Science, DOI 10.1126/science.adf8486](https://www.science.org/doi/full/10.1126/science.adf8486). An authorized visible browser displayed full institutional access through Umeå University; its full article HTML was saved, source-manifested, validated and converted in the external store on 30 October 2026. Inspected: Results1580 (constant-region capture); spatial GEX Methods1623–1624, capture1625–1642, LR boundary1643–1644, selected SR preparation1645–1667 and supplementary index1767–1782 of `science_adf8486_final_browser.html.txt`. The mandatory bundle oligo scraper was rerun; sequence-like analysis barcodes and bulk-validation oligos are not reassigned to the spatial preparation.

The final published Methods directly support the selected short-read TRB workflow below. This establishes the inspected preparation against the publication; it does not assert equivalence of every accepted-manuscript or supplementary detail. The earlier [Zenodo13148082 accepted author version](https://zenodo.org/records/13148082) and its independent review remain historical evidence. No preprint or biological data archive supplies chemistry.

The user supplied both published downloads on 30 October 2026. The supplementary PDF and intact Data S1–S5 archive were moved from `/tmp/camilla` to the external source store, checksummed, manifested and revalidated. The ZIP contains five named XLSX workbooks. Only **Data S1** was selected for chemistry inspection; Data S2–S5 were extracted and converted but not read as chemistry evidence.

Inspected supplementary PDF: published title/DOI; Methods1–72 distinguishing the barcode-losing validation PCR and capture-enrichment evaluation from Spatial VDJ. The selected spatial oligos come from the actual Data S1 worksheets, not the separate bulk SS3 or validation recipes. The graph still begins with prepared full-length spatial cDNA, without inheriting a Visium recipe.

## Selected physical preparation

🟢 Choose the **short-read TRB** route, with the second experiment's single14-cycle postcapture PCR. TRA's14-cycle outer-PCR alternative, LR SMRTbell preparation, bulkSS3 validation and the FigS1F V/C screen are separate scopes. In particular, the FigS1F screen loses the spatial barcode; it is not substituted for SR-SpatialVDJ.

1. Start with full-length spatially barcoded Visium cDNA **before fragmentation**. Select10µl from the library, divide into10ng PCR subsamples whose number depends on concentration, then amplify with25µl KAPAHiFiHotStartReadyMix2X/7.5µl10xcDNA primers/17.5µl sample and water **per reaction**. The graph preserves a variable array of separate reactions through this opaque interface; it does not clone a single PCR or pool templates before amplification. Pool completed products, clean with0.6xSPRIselect and elute40µlQiagenEB.
2. The source delegates concentration to xGenv4AppendixA, then PlateProtocolstep13 onward. Coenrich IG and TRA/TRB constant-region transcripts with custom adapter blockers and IDTxGenHybridizationandWashKit1080584,16h overnight. ChooseIG:TR1:3. This is achieved by **triplicating TR probe sequences in the order and adding equal volumes of the produced IG/TR pools**, not by adding three volumes of TR pool. The graph does not invent bait biotin or a streptavidin formulation from the kit name.
3. Use the selected14-cycle postcapture amplification and purification. Four actual TRB outer-PCR replicates each receive2.5–5ng captured cDNA. OuterTRBV primers and partRead1 amplify for12cycles; source primer stock-dilution wording is100µM→40x→2.5µM, with final-reaction concentration convention unassigned. PreserveqPCR calibration/Qubit/Bioanalyzer analytical allocations separately. Clean each outer product with0.6xAMPure/two80%EtOH washes and eluteEB after15min37°C.
4. Select5ng from each outer stock for its own innerTRBV/partRead1 PCR. Run8cycles; retain the printed72°C30s followed by72°C90s in the cycle programme. Purify and quantify each replicate separately. Select10ng from each of the four inner stocks; **pool at the indexing stage**, assigning the same index to the four replicates.
5. Index-PCR8cycles, bead-purify/quantify and select sequencing allocations. The source locally names KAPAReadyMix for common/outer/inner reactions; its exact indexing formulation/packaging is not assigned from that earlier name. Aggregate index-PCR enzyme/nucleotide/buffer supply classes remain unassembled. Samples are pooled for NovaSeq6000SP withR1=35nt, one6-cycle index read andR2=259nt. Primer landing, index orientation and exact amplicons remain unresolved.

Known ReadyMix components are exposed from each actual named ReadyMix allocation; they are not fresh extra enzyme/nucleotide additions. Stock remainders and analytical aliquots are kept apart. Amounts not specified by the source remain unassigned.

## Sequence and validation boundaries

The earlier accepted text gives two short computational search patterns (`CGACGCTCTTCCGATCT` and `TCTGCGTTGATACCACT`) in **LR analysis**, followed by16spatial-barcode and12UMI bases. These are source analysis locators, not complete spatial primers or permission to import Visium/TruSeq constructs. No exact source oligo state is manufactured from a platform name or matching anchor. Spatial barcodes identify tissue spots, not single cells or physically paired receptor chains; receptor pairing is computational.

## Published Data S1 oligo inventories

🟢 Data S1 now supplies the complete source sequence fields for the selected preparation:

- `SR-Spatial VDJ PartRead1 primer`: the whole22nt PartRead1 order, source line4.
- `SR-Spatial VDJ Inner TRBV`:62 named whole orders and their explicit target/handle columns, lines55–119. All whole orders equal the published34nt Read2 handle followed by that row’s target sequence. The constructor checks this relation rather than guessing a handle from a kit name.
- `SR-Spatial VDJ outer TRBV`:62 named rows, lines171–235, with58 distinct sequence bodies. The worksheet’s interior caption says “Outer TRAV” while the worksheet name and individual IDs say TRBV; the model follows the worksheet/row identities and retains this label discrepancy. It does not silently relabel the original.
- `Hyb. capture - Probes`:174 Ig rows and20 TCR rows, all120nt, lines236–431. Ig has169 distinct bodies. Every row is preserved, including repeated sequences. Source row multiplicity does not independently establish pool molarity or the experiment’s separate triplicate TR-order instruction.
- `Hyb. capture - Blocking Oligos`:both whole patterns, lines432–435. Read1+polyT is82nt:PartRead1(22)+N28+T30+VN. The N/V positions are a degenerate blocking-reagent pattern, not assigned sample barcodes or UMIs. The second blocker, TSO, is30nt.

The source sequence cells do not assign terminal chemistry or occupied binding sites; no5′phosphate, blocked3′end, probe biotin or concrete N/V value is inferred from plain sequence text. These are oligo-pool inventory views, not annealed multimolecular complexes. The four actual replicate allocations retain their respective whole outer/inner pools and PartRead1.

Preamplification/postcapture cDNA primer preparations, the selected whole index-primer orders, historical kit geometry, captured transcript sequences and complete molecular products remain unassembled. This is an open modelling/source question, not a conclusion that those details are unpublished. The28N blocker window alone does not prove the earlier accepted-analysis16/12 barcode/UMI partition. No bulk SS3 or barcode-losing validation primer is substituted.

The final published preamplification programme (line1632) literally says “Repeat steps 2 to 5 six times for a total of five cycles.” Both printed values are retained, and executed cycle count is unassigned. The same ambiguity was previously observed in the accepted manuscript. The model does not choose a corrected cycle count.

The public source-qualified workflow records the final article’s positive source order and stock/fraction boundaries. It remains an incomplete library model: selected cDNA/index primer preparations, occupied transcript geometry, molecular products and primer landing remain unresolved. No complete source audit or accepted missing-information conclusion is claimed.

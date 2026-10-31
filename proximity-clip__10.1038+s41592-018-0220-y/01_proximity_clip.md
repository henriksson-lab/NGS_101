# Proximity-CLIP

Benhalevy et al.(2018), *Proximity-CLIP provides a snapshot of protein-occupied RNA elements in subcellular compartments*. DOI https://doi.org/10.1038/s41592-018-0220-y .

🟢 Three-original capped pass: associated accepted-manuscript PMC6289640 (provenance22; Results32–48 and Methods/analysis184), preliminary associated-manuscript protocol `NIHMS1508063-supplement-2.docx`, then the actual published publisher protocol MOESM3 PDF. Active operational authority is the published protocol, steps33–109/source499–880 plus labeling steps; earlier manuscript reading remains history, not a claim of complete version equivalence. Main VOR Methods and published oligo Table1 are not independently reconciled/read in this cap. Mandatory bundle scraper run; the31MB manuscript supplement was skipped by the default scan. No absence claim.

## Source-owned chemistry and routes

4SU metabolic labeling in APEX2-expressing cells precedes local protein biotinylation; UV312–365nm crosslinking occurs during quencher handling. This differs from APEX-RIP formaldehyde capture. A4SU–protein adduct is not a physical U-to-C mutation of RNA: diagnostic T-to-C is an RNA-oriented sequencing signature of selected copying events. The local relation preserves full RNA and distinguishes complementary firstDNA from reference-oriented report; it predicts no universal stall, mutation frequency or protein specificity.

Published protocol preparative RIPA lysate is clarified; streptavidin capture flowthrough is saved for analysis, not declared disposed. Washed beads split into MS, intact-RNAseq and footprint aliquots. Proteomics stays outside NGS scope. Intact RNA avoids RNaseT1, while footprints receive RNaseT1(1U/µl15min22°C), washes, CIP10min37°C, washes, T4PNK/γ32PATP30min37°C and coldATP100µM5min37°C, then washes. Published MOESM3 step46/source533–534 saves50µl from the first radioactive bead wash as a marker for subsequent Urea-PAGEs. The graph preserves this QC aliquot separately from discarded remaining washes; its exact composition and later disposal are unassembled. Radioactivity is a tracer, not a sequence barcode.

Both RNA branches undergo three sequential ProteinaseK additions/incubations50°C30min each, then source magnetic bead removal. Digestion alone neither physically removes peptides nor transfers protein biotin onto every RNA. Source phenol/chloroform, chloroform and ethanol recovery follow. Intact-RNAseq uses historical NEBNextUltraRNA kit with10min fragmentation; source rRNA depletion applies only to total-cell RNA controls, not automatically to captured RNA.

Footprints are gel-selected into source20–40/41–~80nt classes, extracted, ligated to an adenylated3′adapter by Rnl2(1–249)K227Q overnight on ice, gel-recovered, ligated to5′adapter using Rnl1/ATP1h37°C, and gel-recovered again. Exact oligo bodies/3′blocks/identifier geometry remain unassembled from the unreviewed Table1, not assumed absent. SSIII1h50°C precedes diagnostic cycle selection and preparative PlatinumTaqPCR, Zymo concentration and Pippin size selection. The diagnostic aliquot is a separate physical preparation; it is not reused as the preparative PCR template. Actual acquired read cycles/primer geometry remain unassembled; AM main reports HiSeq2500/3000 without assigning a library-specific run plan here.

## Lineage and open questions

Published protocol explicitly cites PAR-CLIP and its streamlined small-RNA cDNA predecessor plus APEX2 proteomic labeling protocol. This is component-strategy evidence; no complete predecessor recipe is imported. Exact historical NEBNext/primer/library geometry, published main applicability and directly cited detailed predecessors remain pending. No finished schematic or two-reader absence audit.


## Bounded publication follow-up (2026-10-30)

Bounded private correction: gamma-32P ATP and fresh cold ATP are molecular substrates, with the same retained PNK. RNA, adapter and actual ligation mixture pass through heating before fresh ligase. Source ATP-free 3′ versus ATP-containing 5′ supplies remain separate. GlycoBlue and RT dNTP contributions are explicit. Diagnostic PCR is six actual 10µl tubes (9/11/13/15/17/19 cycles), consumed analytically; preparative PCR is three actual 100µl tubes using a separately allocated cDNA contribution. Each tube has copresent primer/enzyme/nucleotide components from its actual aggregate supplies, without asserting commercial premixing. The downstream three-product workup is compound: pooling versus cleanup scheduling remains unassembled. Actual30 of70µl Pippin selection leaves a separate library remainder. Published main, original Table1 and historical kit/read geometry remain genuinely required/unreviewed; no public promotion or absence finding.

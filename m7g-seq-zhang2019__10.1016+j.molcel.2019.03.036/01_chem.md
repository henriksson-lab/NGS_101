# m7G-seq: Zhang et al. 2019

Defining article: [Transcriptome-wide Mapping of Internal N7-methylguanosine Methylome in Mammalian Messenger RNA](https://doi.org/10.1016/j.molcel.2019.03.036).

🟢 Read the actual associated author manuscript PMC6588483.1/NIHMS1526211, main211–215, Fig3B/main97 and the referenced proteinase recipe207. One original. This capped pass does not establish equivalence to the complete final published article. The mandatory primer scraper returned four candidates in that original. Cited MethodS1 (supplement1 PDF) is currently PMC-gated; other listed supplements were not read. No confirmed absence claim or full source audit.

The source-positive boundary is chemical reduction of internal m7G, hydrazide-biotin capture of abasic RNA, streptavidin enrichment and HIV reverse transcription. m7G-MeRIP, model-probe assays and the later m7G-MaP method are separate scopes. The private model selects the enriched mRNA branch. Source input and before-pulldown companion libraries remain unassembled here, rather than being silently folded into that branch.

## Source orderings and local ends

🟢 RNA 3′ linker: `5′rApp-AGATCGGAAGAGCGTCGTG-3SpC3`.

🟢 RT primer: `5′-ACACGACGCTCTTCCGATCT-3′`.

🟢 cDNA 3′ linker: `5′Phos-NNNNNNNNNNAGATCGGAAGAGCACACGTCTG-3SpC3`.

The source's spacer C3 is a terminal chemical block, not an extra cytidine. Canonical Illumina fixed bodies are imported and the whole source orderings are checked locally. The cDNA linker's ten random bases are retained as a placeholder. The inspected Methods give literal N10 degeneracy without an assigned purpose or read-cycle mapping here; no molecule-UMI meaning is imposed. Local Scenes show the source-known 5′ rApp or phosphate donor against RNA/cDNA 3′OH; the full sealed junctions and products remain unassembled.

## Selected preparation

🟢 Starting with4–6µg human mRNA after two poly(A) enrichments, fragment to50–100nt withAM8740. Decap with TobaccoDecappingPlus2, concentrate, perform PNK3′ dephosphorylation37°C1h and rSAP5′/3′ dephosphorylation37°C1.5h. The graph preserves a20µl RNA aliquot and untransferred material from the preceding50µl preparation. It does not import the separate MeRIP5′ phosphorylation step.

🟢 Add1.6µl100µM preadenylated RNA linker,70°C2min/ice, then ligate in50µl with truncatedKQ Rnl2, PEG and buffer25°C2h/16°C8h or longer. Dilute97µl, add3µl5′Deadenylase30°C1h then1.5µlRecJf37°C1h; column recover RNA. Retained primers/adapters are same-tube components, not repeat stock additions.

🟢 Reduce with freshly prepared1MNaBH4 used as a5× stock25°C45–60min; quench/recover with7volumes column binding buffer. Add5µl1MMESpH4.5 and10µl50mMhydrazide-biotin inDMSO to35µl RNA,55°C1h. Reduction and labelling are explicitly incomplete; the model does not turn observed mutation frequencies into stoichiometric conversion. Hydrazide-biotin is a molecular contribution, not just a condition.

🟢 Column recover; bind20µl streptavidinC1 resin4°C20min, wash five times, proteinaseK55°C30min and save the source flow-through/RNA column recovery. The referenced MeRIP paragraph includes another digestion, but this selected paragraph specifies one30min digestion: no extra second treatment is imported.

🟢 Anneal1µl2µM RT primer65°C2min/ice. Add HIVRT,2µl10mM dNTP mix, AMV buffer and RNaseOUT37°C1h/75°C5min; oligo column recover. dNTP is a molecular structural substrate, separate from enzyme/conditions, without guessing whether the concentration label is per nucleotide or total.

🟢 Denature with1µl50µM cDNA linker75°C2min; add Rnl1, ATP, PEG and buffer25°C12h. Amplify using historical NEBNext universal/indexed primers and sequence NextSeq500 single-end80bp. Exact PCR orderings/cycles and full read landing remain unassembled; no current kit recipe is borrowed.

## Remaining scope

The actual cited [MethodS1](https://pmc.ncbi.nlm.nih.gov/articles/instance/6588483/bin/NIHMS1526211-supplement-1.pdf) is needed for a later direct-source pass; its retrieval is gated, not evidence that its details are absent. Full final-version equivalence, complete junctions/read masks, control-library allocation and quantitative recovery remain open. Private partial graphs and checks do not publish a complete schematic or establish reaction yield.

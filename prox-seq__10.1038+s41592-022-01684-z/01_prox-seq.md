# Prox-seq — research notes

Defining source: Vistain *et al.*, *Nature Methods* (2022), [10.1038/s41592-022-01684-z](https://doi.org/10.1038/s41592-022-01684-z) ([open author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC11289786/)). Detailed protocol: Vistain *et al.*, *Nature Protocols* (2024), [10.1038/s41596-024-01030-x](https://doi.org/10.1038/s41596-024-01030-x).

## 1. What it is

Proximity sequencing (Prox-seq) measures extracellular proteins, pairwise protein proximity and mRNA in the same single cells. Two antibody-conjugated DNA probes produce a sequenceable product only when a universal connector can bridge and ligate them.

## Sources read

- 🟢 The defining paper's probe design, reporter architecture, ligation and single-cell library sections were read in the open manuscript.
- 🟢 The later authoritative protocol confirms probe-conjugation and branch-amplification order.
- 🟢 The paper's Supplementary Methods and Supplementary Tables 7–10 (MOESM1 PDF and MOESM3 XLSX) print the connector, branch-specific library and sequencing primers, index choices and bead oligos. The downloaded `.txt` twins were checked against the named table rows.
- 2026-10-17 bounded pass (round breadth-first-15): published JATS PMC11289786 lines 53, 131, 145–173; Supplementary Methods MOESM1 extracted lines 233–345 and 393–415; MOESM3 Tables 7–10 (lines 138–296). The Nature Protocols 2024 detailed protocol and the 10x CG000315/CG000317 Rev B guides were not re-read in this pass.
- The source fetch command is `python3 tools/get_sources.py --doi 10.1038/s41592-022-01684-z`.

## Reporter architecture

- 🟢 Every protein target has probe A and probe B, each an antibody conjugated to a single-stranded oligo. Antibodies are activated with DBCO-PEG4-NHS and reacted overnight at 4 °C with an azide-functionalized PLA oligomer (published Methods, antibody–DNA conjugation, JATS line 131).
- 🟢 Supplementary Table 10 (MOESM3 XLSX, sheet "Supplementary Table 10", extracted lines 193–296) prints one A and one B oligo per 8-nt target barcode. The fixed parts are identical in every pair:
  - Probe A: `/5AzideN/CGCATTGCATCGTCTCGTGGGCTCGGCHHHHACHHHHACHHHNGCAG[8-nt barcode]GATCGCTAAATCGTG-3'`
  - Probe B: `/5Phos/TCGTGTCGTGTCTAAAGTCC[8-nt barcode]ACAAAAAAAAAAAAAAAAAAAAAAAAAAAAA/3AzideN/`
  - The rows are named `PLA.seq.Bead.N.A/B`, but the azide ends, the Methods conjugation passage and the connector pairing below identify them as the antibody PLA oligomers (🟡 identification). An earlier version of this note called them bead oligos; that was a repository misreading.
- 🟢 The ligation connector is `5′-TTTCACGACACGACACGATTTAGGTC-3′` (published plate-based Methods ligase solution, JATS line 149; Supplementary Methods). 🟡 Computed with the model: its central 20 bases pair with probe A's 3′ arm `CTAAATCGTG` and probe B's 5′-phosphorylated arm `TCGTGTCGTG`, leaving 3-nt single-stranded connector ends. This agrees with the article's "20 of those bases hybridized to a connector" (JATS line 53) and with the expected connector-region read sequence `TCGTGTCGTGTCGTGTCTAAAG` used for read validation (Supplementary Methods, extracted lines 407–409). The connector is a splint, not inserted sequence.
- 🟢 The product contains protein barcodes A and B, a UMI, a PCR-primer site and a 3′ poly(A) tail (JATS line 53). In the printed oligos the PCR site is probe A's 5′ region containing U_fwd `GCATCGTCTCGTGGGCTC` (Tables 7–8), and the poly(A) is probe B's 3′ end. 🟡 The Supplementary Methods (line 401–402) say the UMI is extracted from probe A; the `HHHHACHHHHACHHHN` window is taken as that UMI. The schematic previously placed the UMI on probe B; that model error is corrected.
- ⚠️ Published discrepancy, shown on the page: the article says the complete PLA product spans 119 bases (JATS line 53), but the printed Table 10 probe A and probe B are 70 nt and 59 nt, i.e. a 129-nt ligated product (computed by the page model from the transcribed sequences). The schematic draws the printed sequences; how the 119 was counted is not stated.
- 🟢 The poly(A) tail makes the reporter compatible with the same oligo(dT) capture used for cellular mRNA. Unligated probes are said to be "automatically discarded during the library preparation step" (JATS line 53); the page keeps them as an unresolved off-product because no physical separation is described.

## Construct-changing steps

1. Bind antibody–oligo probe pairs to cell-surface proteins.
2. 🟢 Wash, then hybridize connector and ligate proximal A/B probes in 50 mM HEPES pH 7.5, 10 mM MgCl2, 1 mM rATP, 9.5 nM connector, 130 U/ml T4 ligase (NEB), 3 h at 37 °C for the plate experiment and 30 min for PBMCs (JATS lines 149, 153). This is the defining covalent junction and is explicitly marked on the schematic.
3. 🟢 Process cells by Drop-seq, Smart-seq2 or 10x. The single-cell barcode joins the PLA and transcriptome measurements computationally.
4. 🟢 PCR the PLA and cDNA branches and add platform indexes.

## Final library and uncertainty

- 🟢 Each PLA read identifies both targeted proteins and its UMI; reporter counts encode protein abundance and proximity.
- 🟢 Supplementary Table 7 gives droplet-branch custom sequencing primers: Read 1 `GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC`, Read 2 `AGATGCATCGTCTCGTGGGCTCGG`, and i7 `CCGAGCCCACGAGACGATGCATCT`. The paper's Methods state 20/85/8 cycles for Read 1/Read 2/i7.
- 🟢 Supplementary Table 8 gives plate-branch custom primers: Read 1 `TATTTGGCATCGTCTCGTGGGCTCGG`, i5 `CCGAGCCCACGAGACGATGCCAAATA`, and i7 `ACTCTGCGTTGATACCACTGCTTCGGC`. The Methods state 75/8/8 cycles for Read 1/i5/i7.
- 🟢 The distinct 10x branch in Supplementary Methods (MOESM1, extracted lines 317–344) uses CG000315 Rev B for cDNA and CG000317 Rev B for PLA. Its cDNA run is Read 1 28, i7 10, i5 10, Read 2 90 cycles; its PLA run is Read 1 28, i7 8, Read 2 75 cycles, with custom Read 2 and i7 primers spiked into the Illumina primers. These are not the Drop-seq 20/8/85 recipe.
- 🟢 Table 7 also prints `TSO_RNAhybrid` as `AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG` and the 10x PCR primer with terminal phosphorothioates as `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGC*T*C`; Table 10 prints the azide-modified probe A/B oligos described above. The Table 7–9 primers are specific to the published back ends and should not be merged into a single synthetic library.
- 🟢 The 10x branch's supplementary Methods cite CG000315 Rev B for 3′ gene expression and CG000317 Rev B for PLA/feature library construction. Both [CG000315 Rev B](https://assets.ctfassets.net/an68im79xiti/1C16trEdzy1Folq5xbOijE/7e6fb1f504e130bd561d898384da99d9/CG000315_ChromiumNextGEMSingleCell3-_GeneExpression_v3.1_DualIndex__RevB.pdf) and [CG000317 Rev B](https://assets.ctfassets.net/an68im79xiti/1W6c5yj7FAiBQdOKabc6sc/06e962ecdff018fb20b6979859cf1f5e/CG000317_ChromiumNextGEMSingleCell3-v3.1_CellSurfaceProtein_RevB.pdf) are available as original guides. Their appendices print the fixed gel-bead oligo and library shells; the Prox-seq-specific reporter remains defined by its own supplementary Tables 7–10.
- 🟡 The method supports several single-cell back ends. The current schematic represents invariant reporter roles and does not select one of the printed readout branches; its generic sequencing-primer limitation concerns that unspecialized schematic, not an absence of published branch-specific primer sequences.

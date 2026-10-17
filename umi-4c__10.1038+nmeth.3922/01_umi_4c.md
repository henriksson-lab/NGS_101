# UMI-4C — research notes

## Citation

Schwartzman O et al. *UMI-4C for quantitative and targeted chromosomal contact profiling.* Nature Methods 13, 685–691 (2016). [doi:10.1038/nmeth.3922](https://doi.org/10.1038/nmeth.3922)

## Sources read

- 🟢 The defining paper's publisher record and abstract. The full article and Online Methods are not openly available (Europe PMC: PMID 27376768, no PMC copy; `tools/get_sources.py` leaves it as a manual row). They remain unread.
- 🟢 Supplementary Information PDF (publisher MOESM330). Identity: every page footer reads `Nature Methods: doi:10.1038/nmeth.3922`. Read: Supplementary Figure 1 legend (`.txt` lines 79–85), Supplementary Table 1 legend (2193–2195), Supplementary Figure 10 legend (1529–1530) and the complete Supplementary Table 4 (2529–2641). The oligo scraper was run on the bundle and its hits were checked against Table 4.
- 🟢 The author package usage document from the publisher ZIP (MOESM331, `umi4c-usage.html`). Read: Background (lines 74–83) and the baits-table description (190–209).
- 🟡 GEO [GSM2054403](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2054403) belongs to series GSE77553. That series is the C-walks study (Olivares-Chauvet et al., Nature 2016, PMID 27919068), not the UMI-4C paper, and this sample is a C-walk library. Its "3way 4C" extraction protocol applies UMI-4C library construction within that study. Earlier notes called it the authors' deposited UMI-4C description, which was a repository attribution error. Its reaction details are used here only as labelled related-study context.

## 1. What it is

🟢 UMI-4C makes quantitative, viewpoint-targeted profiles from sonicated 3C material. The random sonication endpoint counts molecules (usage Background; Supplementary Figure 1a shows the DNA before and after sonication). The "UMI" is a physical genomic coordinate, not a synthesized N-base barcode.

## 2. Construct-changing steps

1. 🟡 Crosslinking, DpnII digestion, proximity ligation and crosslink reversal. DpnII comes from the usage document's split on the restriction site and the pad sequences ending in `gatc` in Supplementary Table 4. The defining conditions are in the unread Online Methods.
2. 🟢 Sonication of 3C DNA (Supplementary Figure 1a). 🟡 End repair and dA tailing are inferred: the printed adapter strand 1 ends in a single unpaired 3′ T (computed below), and the related-study description includes both steps.
3. 🟡 Dephosphorylation of 5′ ends after dA tailing comes from the related-study description only. It explains why only the phosphorylated strand 2 can ligate. The page labels this step as related-study context.
4. 🟢 Ligation of the forked adapter (Supplementary Table 4). 🟡 Heat release of the non-ligated strand 1 at 95 °C, followed by 1x SPRI, comes from the related-study description only.
5. 🟢 Two rounds of PCR (Supplementary Figure 1c: "two rounds of PCR amplifications"). Universal Primer 2 is the reverse primer in every reaction (Table 4). Every downstream (DS) bait primer carries the printed Illumina prefix, and reads begin with the DS primer followed by the pad (Supplementary Table 1 legend). 🟡 The upstream (US) primer is therefore taken to belong to the first round and DS to the second. That assignment is inferred from the nested design.
6. 🟢 Supplementary Figure 1c quantifies linear products from PCR with only the reverse universal primer. The legend says that linear products are not compatible with Illumina sequencing. The graph shows them as a candidate off-product with an unresolved fate, because no physical separation is described.

## 3. Oligos (Supplementary Table 4, verbatim)

- Forked adapter strand 1: `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT` (no modification printed). This is identical to the DS Illumina prefix.
- Forked adapter strand 2: `(PO4)-GATCGGAAGAGCACACGTCTGAACTCCAGTCAC(NNNNNN)ATCTCGTATGCCGTCTTCTGCTTG`.
- Universal Primer 2: `CAAGCAGAAGACGGCATACGA`.
- G1p1 example bait: US `CATCCCCTCAACTCCACCTTAC`; DS `CTCAACTCCACCTTACCCACCAT` after the Illumina prefix; pad `gatc`. Twelve further baits plus separate 4C-seq primers are listed in the table.

🟡 Computed in `tools/chem.py`, which raises an error if the transcription drifts: strand 1 is P5 (less its last four bases) followed by TruSeq Read 1. Strand 2 is the TruSeq Index 1 primer, a six-base window and the P7 reverse complement. Universal Primer 2 is a prefix of P7. The 5′ end of strand 2 pairs with the 3′ end of strand 1, which leaves the final T of strand 1 unpaired as a single-base cohesive end. The page shows the paired core and its length.

## 4. Final library and reads

🟡 P5/Read 1 (from the DS primer) – DS bait – pad – GATC – reporter – sonication endpoint – A – Index 1/Read 2 arm – six-base window – P7. The page computes the primer landing sites.

🟢 Reads are paired-end: the usage document requires R1 and R2 FASTQs. Read 1 starts with the DS primer and the pad (Supplementary Table 1). 🟡 Read 2 starts at the adapter junction and the sonication endpoint. Read lengths are not given in the inspected sources. Nor do they say whether the six-base window was read as a sample index. The model declares the Index 1 site but does not list it as acquired.

## 5. Open questions (bounded pass 2026-10-17)

- The defining Online Methods remain unread (gated). Still to be checked there: enzymes, end repair and dA tailing, any dephosphorylation, strand release, PCR polymerase and cycles, cleanups, and acquired read lengths.
- Whether the related-study 3way 4C details (C-walks, GSE77553) match the defining UMI-4C experiments.
- Whether the six-base adapter window was used as a sample index in the defining runs.

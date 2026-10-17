# scRibo-seq

## Purpose

scRibo-seq couples ribosome-footprint generation to one-pot small-RNA library construction in individual wells.

## Source read

🟢 VanInsberghe et al., *Single-cell Ribo-seq reveals cell cycle-dependent translational pausing*, Nature 597, 561–565 (2021), [doi:10.1038/s41586-021-03887-4](https://doi.org/10.1038/s41586-021-03887-4). The publisher's version is paywalled at the journal; it was read (2026-10-17, deferred follow-up) from the open copy deposited in the Utrecht University repository ([hdl:1874/445587](https://dspace.library.uu.nl/handle/1874/445587), publisher version, Taverne open access), saved as `uu_dspace_1874-445587_s41586-021-03887-4.pdf` in the source bundle. Sections read in its `.txt`: Fig. 1a schematic text (lines 90–200), main-text method summary (lines 32–41), Methods "Library construction", "Cell lysis and footprint digestion", "Small-RNA library preparation", "Pooling and purification", "Sequencing" (lines 1082–1158), "Read processing" and "Read parsing" (lines 1171–1195), Data/Code availability (lines 1258–1267).

🟢 Supplementary Table 1 (`41586_2021_3887_MOESM3_ESM.xlsx`, Sheet1) was re-read directly for the oligo-ordering strings (2026-10-17). The Reporting Summary (MOESM1) and Peer Review File (MOESM2) were checked; neither prints adapter, primer or read-structure chemistry.

## Oligos from Supplementary Table 1

🟢 `OMV630_miRNA4_3App`: `/5rApp/TGGAATTCTCGGGTGCCAAGG/3ddC/`; the 5′ adenylation and 3′ dideoxy-C are part of the ordered adapter.

🟢 `OMV572_miRNA4_RT`: `GCCTTGGCACCCGAGAATTCCA`.

🟢 `OMV632_miRNA5_5A_10U`: `rCrUrArCrArCrGrArCrGrCrUrCrUrUrCrCrGrArUrCrUr(N1:25252525)r(N1)r(N1)r(N1)r(N1)r(N1)r(N1)r(N1)r(N1)r(N1)`. This retains the source's modified RNA and random-base notation; the ten random positions follow the fixed 5′ RNA sequence.

🟢 The well-indexed forward primer family is printed as individual `miRv6-PCR_F-cbc` rows (for example A1: `AATGATACGGCGACCACCGAGATCTACACGATTCACGACACACTCTTTCCCTACACGACGCTCTTCCG`). The plate-indexed reverse family is printed as `RPI` rows (for example RPI24: `CAAGCAGAAGACGGCATACGAGATGCTACCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA`). The XLSX is the authority for each chosen well/plate member; these are examples, not universal index choices.

## Construct-changing path

🟢 (Methods lines 1088–1099) Single cells are sorted into 384-well plates holding 5 µl mineral oil and 50 nl lysis buffer with cycloheximide. 50 nl MNase is added (37 °C, 30 min); a stop mix of thermolabile Proteinase K, EGTA, EDTA and guanidinium thiocyanate ends digestion (37 °C 30 min, 55 °C 10 min).

🟢 (lines 1100–1113) One-pot small-RNA library preparation: T4 PNK end repair (37 °C, 1 h; the mix as printed contains 4.1 mM uridine triphosphate), then 3′ ligation of the pre-adenylated `OMV630_miRNA4_3App` adapter with T4 RNA ligase 2 truncated KQ in 35.5% PEG-8000 (4 °C, 18 h).

🟢 (lines 1114–1129) `OMV572_miRNA4_RT` is pre-annealed to the 3′ ligation products (65 °C 1 min, 37 °C 2 min, 25 °C 2 min) **before** `OMV632_miRNA5_5A_10U` is joined with T4 RNA ligase 1 (Ambion; 37 °C, 2 h). Maxima H Minus RT then generates cDNA (50 °C 1 h, 85 °C 5 min). The graph now carries this as a hybridization step before the 5′ ligation, with RT extending the pre-annealed primer.

🟢 (lines 1130–1136) Each well receives one of the 384 well-specific `miRv6-PCR_F-cbc` forward primers; the PCR brew carries an `RPI-` reverse primer; Q5 Hot Start, 10 cycles (98 °C 15 s, 65 °C 30 s, 72 °C 30 s).

🟢 (lines 1138–1155) Plates are pooled, the oil removed, the aqueous phase concentrated with n-butanol/diethyl ether and cleaned with diluted AMPure XP beads. Half the pool is run on a 7% polyacrylamide gel and the 175–185 region, "corresponding to an insert size of ~30–40 nt", is excised. Counted reads (lines 1191–1193) are those whose trimmed length is 30–45 nt.

## How the oligos interlock (computed in `scribo_seq._assemble_sites`)

🟡 The well-indexed forward primer is P5 + 10-nt cell index + `ACACTCTTTCCCTACACGACGCTCTTCCG`; its tail overlaps the DNA equivalent of the 5′ RNA adapter's fixed part (`CUACACGACGCUCUUCCGAUCU`) so that the amplified library carries the complete {{= len(scribo_seq.READ1_SITE) }}-nt TruSeq Read 1 site `{{= scribo_seq.READ1_SITE }}` immediately before the UMI. The construction raises if this is not so.

🟡 The RPI 3′ region (after P7 and the 6-nt index) is the reverse complement of `{{= scribo_seq.SMALL_RNA_I1 }}`, which begins with the 3′ adapter bases; the RT primer `OMV572_miRNA4_RT` is the reverse complement of the 3′ adapter bases plus one C, i.e. it also pairs with the adapter's terminal `/3ddC/`.

**Model correction (2026-10-17).** Earlier versions of the model and page drew the TruSeq small-RNA RA5 sequence (`GTTCAGAGTTCTACAGTCCGACGATC`) as the Read 1 site and its reverse complement as the i5 primer. That sequence occurs nowhere in Supplementary Table 1; the published 5′ adapter and forward primers carry TruSeq Read 1. This was a repository modelling mistake, not a publication error.

## Read layout

🟢 (lines 1156–1158, 1171–1177) NextSeq 500, v2.5 chemistry: Read 1 75 cycles, i7 six cycles (plate index), i5 ten cycles (cell index). bcl2fastq is run with `--use-bases-mask Y*,I*,Y*`, so the i5 cell index is written as a read; the first ten bases of Read 1 are the UMI and are concatenated to the cell barcode. `TGGAATTCTCGGGT` is trimmed from Read 1 with cutadapt.

🟡 The sequencing primers are not printed in Supplementary Table 1, and the Methods Sequencing paragraph names no custom primer. Read 1 and i5 (reverse-complement workflow on NextSeq 500) are modelled as the standard TruSeq primers that land on the computed Read 1 site; the i7 read primer is modelled as the small-RNA index primer `{{= scribo_seq.SMALL_RNA_I1 }}` derived from the RPI 3′ region, because the standard TruSeq Index 1 primer would pair over only its last 12 nt.

## Open questions

- Sequencing primers: neither Table 1 nor the Methods print them; the i7 read primer remains inferred from the RPI 3′ region.
- The Methods cite existing protocols (refs 5–7) and a code repository (`github.com/mvanins/scRiboSeq_manuscript`); these were not followed in the 2026-10-17 capped pass.

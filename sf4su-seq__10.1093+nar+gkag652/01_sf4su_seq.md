# sf4sU-Seq — research notes

Gerlevik *et al.*, *Nucleic Acids Research* (2026),
[10.1093/nar/gkag652](https://doi.org/10.1093/nar/gkag652).

## Source scope

🟢 Original article PMC13309788, Methods 118–125, and original Supplementary Methods
297–455 establish the shared no-sonication capture front end and the distinct sf4sU branch.
The complete source bundle was oligo-scraped before reading (17 hits in two files).
The NEB manual is **Version 9.0_07/25**; it is a positive dated kit reference, not an
assertion of the exact study revision or lot. Full linked-source review remains open.

## Source-owned route

1. Living-cell 4sU incorporation is followed by QIAzol/chloroform RNA extraction and
   precipitation. The SNU-derived front end omits TT sonication. Label-containing RNA
   can include sequence synthesized before the pulse; molecular age is not inferred
   from whole-transcript capture.
2. HPDP-biotinylation is followed by Phase-Lock/chloroform cleanup and precipitation,
   µMACS streptavidin binding/washes, **DTT reductive elution**, and miRNeasy recovery
   with on-column DNase I. RNA without the pulse-installed label is a separate class;
   affinity capture is not asserted to be quantitatively perfect.
3. Denaturing TBE-urea gel selection is followed by dialysis and isopropanol recovery.
   The 2026 paragraph says the small-RNA region; a precise universal excision interval
   is not assigned here from an earlier protocol version.
4. RppH treatment is followed by the actual EDTA/heat stop and a second isopropanol
   recovery. The resulting pool is not assumed to consist exclusively of ligation-compatible
   5′-phosphate/3′-OH substrates.
5. The conditional dated E7300 reference ligates the 3′ adapter, hybridizes the RT primer
   to excess free 3′ adapter **before 5′ ligation**, then performs RT and PCR. The graph
   shows the blocked-adapter duplex and its unassigned later fate.

The sf4sU branch does not receive SNU bacterial polyadenylation/QuantSeq or TT directional
library choices. Local declarations own all stages and fates; no shared 4sU recipe or
fixed eight-base index is selected.

## Dated kit oligos

The official E7300 manual p.16 prints:

| Role | As ordered, 5′→3′ |
|---|---|
| 3′ SR adapter | `rAppAGATCGGAAGAGCACACGTCT-NH2` |
| 5′ RNA adapter | `rGrUrUrCrArGrArGrUrUrCrUrArCrArGrUrCrCrGrArCrGrArUrC` |
| RT primer | `AGACGTGTGCTCTTCCGATCT` |
| SR PCR primer | `AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCG*A` |
| Index-primer family | `CAAGCAGAAGACGGCATACGAGAT[6-base index]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC*T` |

The manual's bases and modifications are source-positive. Their mapping to an unverified
study revision is conditional: drawn kit-reference bases are marked inferred, with a
preamble and explicit captions. The exact study index member and complete pre-PCR/final
assembly remain open questions. The kit substrate requirement does not prove that every
source RNA is compatible or recovered.

## Adapter dimers and readout

🟢 Original Supplementary Figure 4C explicitly identifies the desired top bands and a
second-highest adapter-dimer band after ligation/PCR. The graph retains the dimer product
with an unresolved study-specific fate. The kit offers QC-dependent cleanup choices;
this model does not invent complete dimer removal or a universal sf4sU excision window.

🟢 The sf4sU method names NextSeq 500 and computes alignment endpoints using bedtools
`-5` and `-3`. Those operations do not by themselves establish paired-end acquisition or
cycle counts. Deposited sample GSM5452295 now supplies a positive SINGLE, 75-base
genomic-read declaration, detailed below; technical/index cycles remain unassigned.
Computational mapping/filtering remains separate from molecular purification.

## Published wash-buffer difference

Both original PDFs were inspected independently twice:

- The 2026 Supplementary Methods p.18 (text 338) visibly prints **1 mM NaCl**.
- Dölken *et al.* 2008, [10.1261/rna.1136108](https://doi.org/10.1261/rna.1136108),
  original Methods p.1970/PDF p.13 (text 966–967), visibly prints **1 M NaCl**.

The 2026 method explicitly includes modifications. This is a published source difference,
not a proven typo. The page shows both locations/values; the model leaves wash salt
numerically unassigned and does not import the ancestor's other conditions or recovery kit.

## Remaining obligations

The author analysis repository and remaining study-linked metadata scope, the original
TT/capture predecessor chain, remaining supplement/figure scope and relevant kit/QC instructions
still require bounded original-source dispositions. The currently fetched bioRxiv v2 is
posted **24 October 2025**, despite the DOI containing 2021; version identity must be
preserved. No missing-information claim is accepted from this bounded checkpoint.

## Deposited size-selected library and version differences

🟢 Original GSE179306 sample metadata identifies **GSM5452295**,
`PL_size_selected_SNU_WT`, as the HEK293 size-selected library. Its treatment/extraction
and linked **SRX11439177** experiment confirm the size-selection/RppH/NEBNext route.
The source window varies:

| Original record | Stated window |
|---|---|
| 2026 article Fig. 2D legend (JATS line 168) | size-fractionated (≈ 100 nt) RNA |
| GSE179306 series Overall Design, line 11 | 50–100 bp |
| GSM5452295 extraction, line 70; SRA construction protocol | 70–100 bp |

The 2026 published Methods paragraph describes the small-RNA region without repeating an
interval. Round 18 removed the superseded bioRxiv v1/v2 values from this comparison: under
the published-version-first rule an ordinary preprint revision is not a competing source.
The public page now shows this discrepancy and leaves the exact gel window unassigned.
A protocol revision, specific-library choice or transcription difference remains possible;
none is silently selected as the explanation.

🟢 Both original ENA experiment XML and NCBI full SRA metadata explicitly declare
**SINGLE** for SRX11439177. All four archived runs, **SRR15130590–SRR15130593**, have
75-base populated reads: original filenames are lane L001–L004 `_R1_001.fastq.gz`, and
read-statistics index 0 has average 75, standard deviation 0, and count equal to spots;
the other statistics slot has count 0. The page now declares this sample's **Read 1 = 75**.
This is the deposited genomic read length, not a recovered complete instrument/index mask.
No raw reads were downloaded. Alignment `-5`/`-3` operations remain separate from acquisition.

🟡 The 2026 Results (JATS line 178) say the short transcripts were "sequenced from both
ends (Supplementary Fig. S4C)". The deposited layout above is SINGLE with an empty second
read slot. The page shows both statements; the model keeps one 75-nt Read 1 for
GSM5452295 and treats 5′/3′ ends as computed from alignments. Whether the phrase means
paired-end acquisition of other, undeposited runs or end recovery after adapter trimming
remains open.

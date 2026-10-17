# SMART-seq3xpress — oligo source review

Hagemann-Jensen M, Ziegenhain C, Sandberg R. "Scalable single-cell RNA sequencing from
full transcripts with Smart-seq3xpress." *Nat Biotechnol* 2022;**40**(10):1452–1457.
doi:[10.1038/s41587-022-01311-4](https://doi.org/10.1038/s41587-022-01311-4)

> **Evidence marking.** 🟢 verbatim from a source read · 🟡 derived or inferred ·
> 🔴 not published in any source read here.

## Sources read

- 🟢 **The defining paper's own Methods**, open access under CC BY as PMC9546772. It
  prints the two TSOs side by side: "Original Smartseq3 TSO
  (5′-Biotin-AGAGACAGATTGCGCAATGNNNNNNNNrGrGrG-3′; IDT). Improved TSO
  (5′-Biotin-AGAGACAGATTGCGCAATGNNNNNNNNWWrGrGrG-3′; IDT)", and both PCR primers,
  "Smart-seq3 forward primer (5′-TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGATTGCGCAATG-3′)" and
  "Smart-seq3 reverse primer (5′-ACGAGCATCAGCAGCATACGA-3′)". It also gives the UMI-read
  pattern, `ATTGCGCAATG` with up to two mismatches.
- 🟢 **The authors' protocols.io protocol** `smart-seq3xpress-bwh4pb8w` — the ordered oligo
  table, with supplier, purification and concentration:
  `/5BiosG/ACGAGCATCAGCAGCATACGA(T)30VN`, `/5BiosG/AGAGACAGATTGCGCAATGNNNNNNNNWWrGrGrG`,
  `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGATTGCGCAA*T*G` and `ACGAGCATCAGCAGCATAC*G*A`, with a
  footnote that `*` marks phosphorothioate bonds.
- 🟢 **Springer supplementary files** `41587_2022_1311_MOESM1/2_ESM.pdf` — supplementary
  figures; no additional oligo table.
- 🟢 **The original protocol-linked `384set_Nextera_BC_order.xls`**, sheets
  `UD_i5_indices` and `UD_i7_indices`: complete ten-base index primers; retrieval
  and the explicit illustrative pair are detailed below.
- 🟢 **The local core oligos and index workbook** establish the represented handles.
  Canonical strings in `lib/nextera.py` and `lib/illumina.py` are used for computation;
  their reuse alone does not establish an xpress reagent assignment.

The RT/TSO bases are transcribed from the sources above. The indexed endpoint is
now declared in this protocol, using the author-linked order workbook below. No secondary redrawing of the
SMART-seq family is used as evidence here, and none is needed: the paper and the author
protocol supply the represented core oligos directly. The article prints the TSOs
and PCR primers; the author protocol supplies the ordered oligo-dT handle and
the full modification notation.

## Oligos

```
Smartseq3_OligodT30VN  5'- /5BiosG/ ACGAGCATCAGCAGCATACGA TTTTT...TTTTT VN -3'   (30 T)  🟢
Smartseq3xpress_TSO    5'- /5BiosG/ AGAGACAGATTGCGCAATG [8 N positions] WW rGrGrG -3'     🟢
Fwd_PCR_primer         5'- TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGATTGCGCAA*T*G -3'            🟢
Rev_PCR_primer         5'- ACGAGCATCAGCAGCATAC*G*A -3'                                   🟢
```

🟢 The final paper distinguishes the original SMART-seq3 TSO from the improved
xpress TSO by the **`WW` spacer** between the UMI and the G tail.
`AGAGACAG` is the 3'-terminal 8 nt of the Nextera
mosaic end, so the stock Read 1 primer (`s5 + ME`) ends exactly where the 11-bp tag
`ATTGCGCAATG` begins; the forward PCR primer rebuilds the whole s5 arm as it amplifies.
`W` is A or T, so the spacer is a two-base degenerate window, not a fixed sequence.

The local declaration selects the final improved TSO, Maxima H Minus RT,
the published SeqAmp pre-amplification option and commercial Nextera tagmentation
(ATM from the XT kit, or TDE1) with Phusion index PCR. The paper separately evaluates other pre-amplification
polymerases and in-house Tn5. These are explicit source variants; matching adapter
patterns do not make their conditions or reagent assemblies interchangeable.

## What the page asserts, by role

| Role | Status | Where it comes from |
|---|---|---|
| biotinylated oligo-dT with its own handle | 🟢 | protocols.io ordered table; shared with SMART-seq3 |
| UMI-bearing TSO with `WW` spacer | 🟢 | defining paper's Methods, "Improved TSO"; protocols.io |
| 10-position UMI, comprising 8N + 2W | 🟢 | ordered TSO and deposited annotation `UMIsequence`/`UMIlength`; WW also acts as a spacer |
| 11-bp 5' tag `ATTGCGCAATG` | 🟢 | both sources, and the read-parsing pattern |
| forward / reverse PCR primers, incl. terminal phosphorothioates | 🟢 | protocols.io ordered table |
| Nextera ME / s5 / s7, read and index primer entry | 🟢/🟡 | own forward primer and complete workbook index oligos; sequencing sites computed on the schematic endpoint |
| P5 / P7 flowcell ends | 🟢 | complete workbook index oligos match the canonical strings |
| i5 / i7 index *members* | 🟢 illustrative | original order workbook, A1 on each sheet; no experimental cell assignment |

## Open roles

- 🟢 **Index oligos recovered.** The original author protocol (DOI
  `10.17504/protocols.io.yxmvmk1yng3p/v1`, integer ID 51484) names
  `384set_Nextera_BC_order.xls`. Its anonymous public API stores the attachment URL
  in `guidelines.entityMap[2].data.source`:
  `https://content.protocols.io/files/gh7u2zw6.xls`. The PDF names the file but has
  no clickable attachment annotation. The recovered original has sheets
  `UD_i5_indices` and `UD_i7_indices`, each with 384 full oligos. All index windows
  are ten bases. The defining paper Methods, "Sequence library preparation for
  Smart-seq3xpress", independently specifies ten-base dual indexes.
- 🟢 **Explicit illustrative pair.** The page now uses i5 A1 `TTGGTACGCG` and i7
  A1 `ATAGGCGCTC`. Their complete ordered primers are `P5 + i5 + s5` and
  `P7 + i7 + s7`; the top-strand endpoint therefore carries the reverse complement
  of the ordered i7 index. These are source-known available oligos, **not an
  assignment to any published cell or well**. The former eight-base endpoint was
  a repository model choice inherited from the shared recipe; it was not an error
  in the publication. The author's eight-base alternative remains a protocol option,
  while the represented published xpress route uses ten.
- 🟡 **Sequencing scope.** This page represents the Illumina-compatible library.
  The defining Methods also reports MGI conversion/circularization and custom
  sequencing primers; that distinct branch is not represented here. A library
  architecture does not by itself establish the run cycle mask of every experiment.
- 🟡 Exact tagmented-end and extension intermediates remain unassigned molecular
  states. The declaration represents the commercial Nextera option and the reported
  initial extension hold; it does not inherit a complete loaded-transposome assembly
  or a fixed gap geometry from a SMART-family helper. Other tagmentation products
  remain visible without an absolute claim that they cannot amplify.

## Processing-example discrepancy

🟢 The original author protocol v1 p33 (local text2252) gives the bcl2fastq
example `--use-bases-mask Y150N,I8,I8,Y150N`. Its adjacent p34 configuration
(local text2286–2290) requests `BC(1-10)` from both index FASTQ files. The final
article specifies ten-base dual indexes, and the original order workbook contains
ten-base indexes. Eight retained index cycles do not supply ten bases to that
configuration. The protocol also permits eight-base index alternatives, so this
is an inconsistency between the neighboring examples rather than evidence that
all eight-base use is wrong. The page displays both source locations and uses the
ten-base molecular endpoint; it selects no actual run mask and does not silently
repair the command. The experimental application of these examples remains uncertain.

🟢 The same p34 example extracts `UMI(12-21)` after the eleven-base tag.
Recovered original annotation tables go further: the v2/v5/v6/v7 TSO-test
tables label `NNNNNNNNWW` with `UMIlength=10`, and give complexity 262144
for that design. The represented identifier therefore spans all ten positions,
while preserving the distinct eight-N region and two-W spacer. There are eight
fully degenerate N positions, not ten. The MGI configuration repeats the same
extraction window. The earlier graph labelled only eight positions as UMI;
the deposited source metadata prompted this repository model correction.

## Recovered study annotations

🟢 The article-linked ArrayExpress records have original IDF/SDRF metadata
at E-MTAB-11488, E-MTAB-11452 and E-MTAB-11467. The first is explicitly the
10x v3.1 comparator. The other two name sample-level barcode annotations,
which were recovered through the current BioStudies file inventories after
the legacy processed-ZIP URLs failed. They are available originals, not
confirmed missing information.

The PBMC annotation contains 43,671 records with explicit `cellID`, `run`
and `barcode` fields. Fifteen optimization annotation files retain distinct
conditions and TSO designs. The two low-volume files have 16-base read-group
identifiers, while the other thirteen optimization files have 20-base
identifiers. These fields must not be converted into a universal ten-plus-ten
physical primer assignment without orientation and source-version evidence.
The schematic continues to use its explicit illustrative workbook pair.

## Local declaration boundary

The RT primer, improved TSO, pre-amplification product and selected reactions are
declared locally from the xpress originals. The amplified product includes the
source-known s5/ME overhang introduced by the forward primer. The former shared
helper left that overhang out of its amplified-state drawing; this was a repository
model error, not a publication error. Final 5′ and internal fragment classes remain
separate. The endpoint and exact primer-site checks do not prove the complete
RT, tagmentation or amplification transitions or an experiment-specific cycle recipe.


## Bounded pass 2026-10-17 (round breadth-first-18)

🟢 Final paper Methods (PMC9546772 JATS text lines 70–100): "Low-volume Smart-seq3
and Smart-seq3xpress library preparation", "After pre-amplification workflow",
"Sequence library preparation for Smart-seq3xpress" and "Sequencing". Author
protocol v1 PDF (May 30 2022) RT/pre-amplification tables (local text 690–960),
tagmentation/SDS/index PCR (1130–1290) and pooling/clean-up (1990–2025).

| # | Step (represented route) | Final Methods conditions |
|---|---|---|
| 1 | Lysis under 3 µL overlay, 72 °C 10 min | 0.1% Triton X-100, 5% PEG8000, 0.5 µM oligo(dT), 0.5 mM dNTP each, RNase inhibitor (adjusted to RT volume) |
| 2 | RT + template switch, Maxima H Minus 2 U/µL | 25 mM Tris-HCl pH 8.3, 30 mM NaCl, 0.5 mM GTP, 2.5 mM MgCl2, 8 mM DTT; 42 °C 90 min, 10×(50 °C 2 min, 42 °C 2 min); SeqAmp route 0.125/0.5 µM oligo-dT and 0.75/2 µM TSO |
| 3 | SeqAmp pre-amplification, 0.6 µL mix, 0.025 U/µL | 95 °C 1 min; 6–18×(98 °C 10 s, 65 °C 30 s, 68 °C 2–6 min); 72 °C 10 min |
| 4 | Add 9 µL water; transfer 1 µL to tagmentation | no bead cleanup on this route |
| 5 | Tagment in 2 µL: 1× buffer (10 mM Tris pH 7.5, 5 mM MgCl2, 5% DMF) with 0.025–0.5 µL ATM or 0.0002–0.01 µL TDE1 | 55 °C 10 min |
| 6 | Stop/strip with 0.5 µL 0.2% SDS | |
| 7 | Custom Nextera index primers (10-bp dual indexes) + Phusion 0.01 U/µL | 72 °C 3 min; 95 °C 30 s; 12–14×(95 °C 10 s, 55 °C 30 s, 72 °C 30–60 s); 72 °C 5 min |
| 8 | Spin-out pooling into a reservoir; 22% PEG beads 0.7:1 | |

The graph now declares the lysis and RT conditions, the dilution/transfer
bottleneck (remainder kept outside the library), the SDS strip and the pooled
0.7:1 bead cleanup. The protocol adds Tween-20 to the index-PCR mix to quench
leftover SDS; the article does not mention it, and the graph does not select it.

🟢 **Primer-concentration discrepancy shown on the page.** The Methods SeqAmp paragraph
gives "0.5 µM/1 µM Smartseq3 forward and reverse primer"; the protocol v1 step 9
table lists 0.5 µM each. The page names the primer pair but selects no
concentration. Other article-versus-protocol differences (cycle and elongation
ranges, index-primer volumes) are option ranges versus a default, not conflicts.

Open: experiment-specific pre-amplification and index-PCR cycles; which tagmentation
enzyme and amount each data set used; exact tagmented-end geometry.

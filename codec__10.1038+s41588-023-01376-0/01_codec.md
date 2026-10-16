# CODEC — research notes

Defining source: Bae *et al.*, [10.1038/s41588-023-01376-0](https://doi.org/10.1038/s41588-023-01376-0).

## Evidence used

- 🟢 The article defines CODEC as “Concatenating Original Duplex for Error Correction” and states that an adapter quadruplex plus strand-displacing extension physically links information from both original DNA strands in one Illumina read pair.
- 🟢 The read design recovers UMI, sample index and insert while retaining both complementary copies for a single-duplex consensus.
- 🟢 Supplementary Table 5 prints all four CODEC adapter oligos for eight matched
  sample-index sets. It defines `*` as a phosphorothioate linkage and `P-` as a 5′
  phosphate. The schematic uses set 1 as a representative molecule and preserves both
  modifications in its ordering list.
- 🟢 `LD4-brid5` contains the reverse complement of the standard TruSeq Read 1
  primer site; `LD4-brid7` contains the reverse complement of the standard TruSeq Read 2
  primer site. Each is preceded by a 3-nt UMI and an 18-nt inline sample index, consistent
  with the outward Read 1/Read 2 geometry in Fig. 1b and the CODECSuite description.
- 🟡 The paper draws the post-extension linked product by region rather than
  printing it as one contiguous strand. That state is therefore role-scaled; it is not a
  claim that the adapter oligos themselves are unavailable.

## Protocol boundary

This page covers CODEC library construction. Conventional Duplex Sequencing and NanoSeq remain separate protocols because they recover the two strands as independently amplified read families rather than concatenating them physically.

## Fresh original-source audit

The standard source fetcher now retrieves the final PMC10181940 article, all four
publisher supplements, and the linked 2021 preprint. Supplementary Table 5 is
available; the old ledger's inaccessible-table claim was stale. Full-bundle oligo
scrapers preceded reading (24 hits in three files, then 26 in four after the named
NEB manual was added).

Original locators in the CODEC source bundle:

| Original | Sections checked | What it establishes |
|---|---|---|
| `s41588-023-01376-0_PMC10181940.xml.txt` | 104–114; 198–220; 294–302; 378–384 | Internal outward reads; commercial E7103L baseline with explicit CODEC modifications; optional repair methods; expected byproducts; data/software links |
| `s41588-023-01376-0_supp_41588_2023_1376_MOESM1_ESM.pdf.txt` | 14–40; 484–545 | Demultiplexing/read-start definitions; all eight matched adapter sets and modified capture blockers |
| `s41588-023-01376-0_supp_41588_2023_1376_MOESM4_ESM.xlsx.txt` | All four sheet headers and unique repair/library/instrument/read-length fields | Experiments use multiple repair and CODEC versions; comparator workflows are separately labeled |
| `manuale7103-e7645.pdf.txt` | Title/version; 70–230 | Official E7103/E7645 v7.0 September 2022 end-preparation/ligation baseline; USER treatment is conditional on standard non-indexed NEB adapter, which CODEC replaces |

### Repository correction

A literal whole-oligo comparison found that the old shared `CODEC_BRID7_1`
constructor duplicated `GC` at its spacer/core boundary: it made 100 nt while
Table 5 prints 98 nt. The local constructor now uses the 12-nt
`CAATCTATAAGTT` spacer plus the complete 30-nt reverse-complement linker core.
The original bridge-5 is also 98 nt. All four complete modified set-1 ordering
sequences are now checked against Table 5. This was our model transcription
error, not a publication discrepancy. A Scene checks pairing of the two bridge
cores during construction.

### Implementation and read boundaries

The local model selects the paper's commercial end-repair/dA-tail baseline,
CODEC ligation with added 5′-deadenylase, phi29 extension, and KAPA amplification
with the named xGen 1077677 primer mix. It does not inherit NanoSeq restrictive
repair or the optional Duplex-Repair route. Ordered xGen primer strands and a
contiguous amplified product remain unassigned model data; this does not assert
that a relevant source omits them.

The displayed set 1 has three `N` positions followed by the 18-base inline sample
index and a terminal T. Supplementary Text and Extended Data Figure 1f also show
`NNNA`/`NNNT` read-start variants, so a uniform three-position identifier should
not be generalized to every matched set. Both sample identifiers are sequenced
within Read 1/Read 2. The graph selects the original paired 166-cycle NovaSeq
option and labels the paired 260-cycle alternative separately.

Extended Data Figure 2 describes putative blank-ligation, double-ligation and
intermolecular-extension mechanisms. The graph includes these off-products by
role and retains the source's conditional mechanism; it does not claim exact
byproduct nucleotide sequences or proven reaction routes.

### Remaining source scope

The original image-based reporting summary has been visually inspected in full
(three pages): statistical reporting, software/data, study design and human-sample
ethics; an independent second source-scope read remains pending. The peer-review file and linked preprint/supplements are retrieved
but not fully read in this pass. The named KAPA amplification manual, exact xGen
1077677 commercial-primer identity, and original CODECsuite/Zenodo software scope
need source-specific checks. Main-article references to optional Duplex-Repair,
restrictive ddBTP repair and targeted enrichment are distinct branches; their
recipes are not assigned to the selected baseline. No missing-information claim
has been accepted, and the audit remains in progress.

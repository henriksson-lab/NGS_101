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

### Bounded pass 2026-10-17-breadth-first-08

🟢 Re-read Methods lines 202–220 of the final JATS. Added to the graph the
source-stated steps the model had omitted: Covaris shearing to a mean of
150 bp (line 202; not applied to cfDNA), the 1.8× AMPure XP cleanup before
ligation with 40 µl low-TE elution (exception 1), 0.75× AMPure XP after the
phi29 extension (exception 4), and two 0.65× AMPure XP cleanups after PCR
(exception 6; for cfDNA the second is a 0.5×/0.65× double-sided selection).
Reagent labels now carry the stated adapter annealing (85 °C 3 min, −1 °C/min
to 20 °C, 12 h RT; 500 nM in ligation), 1 h ligation, phi29 M0269L at 30 °C
for 20 min with 0.2 mM dNTP, and KAPA HiFi with 2 min extension. Input is
2.5–20 ng. The KAPA cycle number is not stated in these lines and stays
unassigned.

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

## Native ownership migration, 2026-10-18

Sources read: Published PMC10181940 JATS .xml.txt: Results Adapter quadruplex and workflow/Fig1, lines100–118; Methods DNA samples/ERAT/CODEC/Sequencing197–220. Published MOESM1 supplementary PDF .txt Table5 lines481–548, allfour modified set1 ordering sequences and eightmatched sets. Preprint and current vendor manual not used.

Removed protocol_specs.get/Definition/workflow and rna_library_operations.molecule. Native graph owns actual genomic duplex and dA-tailed substrate, four-strand representative adapter (two bridge strands plus two phosphorothioate terminal-T arms), separate ligated Watson/Crick intermediate chains and actual linked single-chain PCR precursor. Sample/UMI feature identities preserved from adapter through local endpoint. Original Methods quantities/variant conditions retained, including pre-ligation cleanup absent from v0, cfDNA alternate final cleanup and166/260 paired-cycle alternatives. Candidate blank/double/intermolecular byproducts are unresolved material fates rather than unsupported physical discards.

Model boundaries: Native complete ligation branch geometry/contiguous extension product; xGen1077677 ordered physical primers/first-copy PCR proof; optional capture and ERAT alternative graphs remain outside selected baseline route. No source-absence claim is accepted by this capped pass. Bundle oligo scraper run and source passages checked.

### Local SID7 orientation correction, 2026-10-18

The final region-level reference now retains the source bridge7 orientation: reverse complement of SID7 → three UMI positions → reverse complement of the Read 2 primer site. Its previous SID7 segment had the sense ordering-oligo value adjacent to a Read2rc site, disagreeing with the exact bridge7 local region already verified against Supplementary Table 5. This was a repository orientation error, not a publication discrepancy. Read 2 extends toward that region's 5′ end and reports the source SID7 value; the sample-index feature identity is unchanged, and no alternate/transformed source index was invented to preserve the old computed readout.

The complete concatenated post-extension architecture and xGen PCR geometry remain explicitly unassigned. This implementation-only correction relies on the existing Table 5 source scope; no new original was investigated and no absence conclusion is made.

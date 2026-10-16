# NicE-seq — research notes

Ponnaluri et al., *Genome Biology* (2017),
[doi:10.1186/s13059-017-1247-6](https://doi.org/10.1186/s13059-017-1247-6).

## Sources read

- 🟢 The defining open-access article and methods (PMCID: PMC5488340).
- 🟢 Original GSE97889 family metadata: all 14 sample headers and the five unique
  extraction/seven processing paragraphs were inspected. Nine NicE-seq records share
  the deposited labeling/library paragraph; four WGBS and one ATAC record describe
  separate branches. The large platform annotation table was not audited.
- 🟢 New England Biolabs, *NEBNext Ultra DNA Library Prep Kit for Illumina*, E7370
  manual, for the discontinued commercial end-prep and ligation workflow named by the
  paper. The fetched manual is **Version 9.1_1/23**, a later companion; it requires
  separate oligo kits and does not establish which historical oligos the study used.
- 🟢 All S1–S17 captions and a visual survey of their 17 embedded original figures:
  peak-caller comparisons, genomic tracks, Venn diagrams, heatmaps and biological
  results. This was a source-scope inspection, not a statistical reanalysis.

## What carries the accessibility signal

1. 🟢 Cells are either used natively or crosslinked in 1% formaldehyde for 10 min and
   quenched with 125 mM glycine. NP-40 cytosolic buffer isolates nuclei.
2. 🟢 Nt.CviPII recognizes `CCD` (`D` = A, G or T) and nicks one strand where chromatin
   permits access. The new 3′-OH primes *E. coli* DNA polymerase I. Polymerase extension
   plus 5′→3′ exonuclease displacement translates the nick through a tract of DNA.
3. 🟢 The 200-µl labeling reaction contains 2.5 U Nt.CviPII, DNA polymerase I and 30 µM
   each dNTP, including 6 µM biotin-14-dATP and 6 µM biotin-16-dCTP. It runs for 2 h at
   37 °C with 800-r.p.m. mixing. Biotin is therefore internal to accessible DNA, not an
   adapter-end label.
4. 🟢 EDTA stops the reaction while RNase A removes RNA. Proteinase K and SDS reverse
   crosslinks at 65 °C, and phenol/chloroform extraction recovers covalently labelled
   genomic DNA.
5. 🟢 Covaris sonication makes approximately 150-bp fragments. The E7370 workflow
   end-repairs, dA-tails and ligates an Illumina-compatible adapter *before* enrichment.
6. 🟢 Adapter-ligated DNA is incubated with blocked streptavidin beads for 2 h in 2 M
   NaCl. Four B&W/Triton washes and a TE/Triton wash remove unlabelled molecules. PCR
   uses **4 µL of the 40-µL bead suspension**, usually ten cycles. The graph feeds
   only this aliquot into PCR and marks the remaining 36 µL as unresolved; it does not
   invent disposal or storage. The 25-cell branch instead sonicates all recovered DNA
   and uses 10 µL beads; that distinct variant is not merged into the standard aliquot.

This order matters: read endpoints are mechanical shear positions and do not identify
the original Nt.CviPII nick. Enrichment across the fragment reports proximity to an
internally biotinylated accessible-chromatin tract.

## Article/deposit differences

🟢 The article Methods and the nine deposited NicE-seq extraction paragraphs differ:

| Item | Article Methods | GSE97889 NicE-seq metadata, including GSM2579632 |
|---|---|---|
| Standard DNA polymerase I | 10 U | 50 U |
| SDS/proteinase K treatment at 65 °C | Overnight | 1 h |
| Triton X-100 in B&W washes | 0.05% | 0.005% |

The public page displays all three differences. It retains the shared operations without
selecting universal values. The separately described 25-cell branch uses 5 U polymerase I
in both sources; it is not evidence that either standard-reaction value is a typo. Whether
the differences represent run variants or transcription errors remains unresolved.

## Commercial-oligo and model boundaries

🟢 The article names E7370S for end repair, dA tailing, ligation and PCR. Its source-owned
order is ligation → **capture without further purification** → bead-bound PCR. The model
keeps that order and does not import the later vendor manual's intervening cleanup.
The 2023 vendor manual makes USER treatment conditional on use of NEBNext adapters;
its separate WGBS loop-adapter/USER recipe in GEO does not assign that adapter to NicE-seq.

The current model uses **role-level library ends**. The previous compatible six-base
single-i7 shell was a repository-selected example, not a recovered study oligo identity;
it has been removed. No index width, index-read count, exact adapter arms or base-level
sequencing primer sites are assigned. Adapter ligation and PCR geometry remain explicitly
unavailable pending reconciliation of the study-specific oligos. These are open model
questions, not accepted claims that the authors never supplied the information.

The graph serializes nick creation and nick translation to explain molecular roles;
the original performs them in the same labeling mixture. RNase digestion, protein removal
and DNA extraction are distinct stages; actual extraction and capture/wash fractions are
shown as discarded material. The graph does not equate every unbound molecule with an
experimentally guaranteed absence of biotin.

## Sequencing

🟢 The NicE-seq GEO records name NextSeq 500 and paired-end processing. The article's
explicit 75-bp run statement describes the separate ATAC-seq comparison. A study-specific
NicE-seq run mask remains unassigned in this model; linked experiment metadata have not
yet been read. No read-cycle length or index architecture is inferred from the instrument.

## Remaining source obligations

The review remains in progress. The actual linked NEB oligo hub and relevant historical
manual identity, study experiment metadata, and bounded processing links require their
own original-source dispositions. The registered E7370 product page returned HTTP 403;
that is an availability result, not a missing-sequence finding. No blanket absence claim
is accepted from the scraper's empty article/supplement output.

## Bounded curation round 2 (2026-10-16)

Inspected: Published original PMC5488340 JATS160–170.

No correction identified in inspected scope:150bp sonication, NEB E7370S adapters before immediate streptavidin capture,4µL aliquot of40µL beads and usually10 PCR cycles agree with current graph.25-cell branch remains separately qualified.

Deferred: Remaining original Methods/supplements, reagent-specific kinetics and vendor adapter identities; existing article/deposit discrepancy evidence was not independently re-reviewed. This capped pass does not establish a complete source review or publication-wide absence.

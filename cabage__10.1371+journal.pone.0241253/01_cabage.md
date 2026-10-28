# CaBagE — Cas9 protection followed by background exonuclease digestion

Wallace et al. (2021), *CaBagE: A Cas9-based Background Elimination strategy for targeted, long-read DNA sequencing*. [Published article](https://doi.org/10.1371/journal.pone.0241253).

## Bounded source pass

The published article, Methods and the two published review/response attachments
were inspected on 2026-10-28. This is a capped pass, not an exhaustive source review
or an accepted missing-information audit. Originals remain in the external source
store under `sources/cabage__10.1371+journal.pone.0241253/`.

- 🟢 `journal.pone.0241253_PMC8031414.xml.txt`, Introduction paragraphs 73–77,
  Results overview 81–90, and Methods 198–226 establish the molecular workflow.
- 🟢 `journal.pone.0241253_supp_pone.0241253.s007.docx.txt` is the published
  referee report, not a bench protocol. Its comparison of nCATS and CaBagE was
  inspected; it is not the authority for experimental reagent choices.
- 🟢 `journal.pone.0241253_supp_pone.0241253.s008.docx.txt`, author response to
  reviewer 1, acknowledges the earlier Negative Enrichment strategy and explains
  the scope of the paper's amended novelty claim.
- The full-bundle oligo scraper was run. Hits from the downloaded superseded
  preprint were excluded from active chemistry. The article identifies S1 Table
  as the guide/target-sequence table; its original XLSX was inspected in the follow-up model pass (S1 Table, all
  guide rows and gBlock entry). No claim of absent guide sequences is made.

## What distinguishes the chemistry

🟢 Two Cas9–guide complexes flank the target in **PAM-in** orientation. Their
stable association with the cut DNA protects the target from an exonuclease mixture.
Exonucleases I, III and lambda deplete susceptible background substrates before
Cas9 is released by heating. This is physical depletion of background DNA within
the reaction, followed by bead recovery; it is not nCATS's selective adapter
ligation after dephosphorylation. Degraded background products remain in the tube
until the cleanup. The method is PCR-free.

## Published native-DNA branch

1. 🟢 The genomic DNA is pretreated with NEB FFPE Repair Mix M6630S, using the
   manufacturer protocol for other user-supplied library reagents. This reduces
   nicks that could expose desired fragments to exonuclease. Exact repair chemistry
   is delegated to that named reagent and was not independently resolved here.
2. 🟢 Cas9 and guides are preassembled in CutSmart buffer at 25°C for 10 minutes.
   The source reports a roughly 10:10:1 Cas9:guide:target molar ratio and cleavage
   at 37°C for 15 minutes. S1 Table gives locus-specific guide sequences. The model records
   the C9orf72 examples without inventing a target/reference sequence.
3. 🟢 Add Exonuclease I (40 U), Exonuclease III (200 U), and lambda exonuclease
   (20 U) immediately after cleavage; incubate at 37°C for two hours, then 80°C
   for 20 minutes. The Results describes the heat step as both exonuclease
   inactivation and Cas9 dissociation.
4. 🟢 Add dATP and Taq DNA polymerase M0267S; incubate at 72°C for five minutes
   to A-tail the exposed target ends. The article calls its dATP stock “10mm”;
   the current model does not rely on interpreting this typography as a numerical
   concentration specification.
5. 🟢 Ligate adapters from **SQK-LSK109**, using its ligation buffer and adapter
   mix plus Quick T4 ligase E6057, for ten minutes at room temperature.
6. 🟢 Recover with 0.3× AMPure XP beads and two **Short Fragment Buffer** washes,
   then elute. Do not silently substitute nCATS's Long Fragment Buffer. The paper
   compares SFB and LFB in S1 Figure, which remains a separate uninspected original.
7. 🟢 Sequence on MinION **R9.4.1**, with the LSK109 sequencing-buffer, loading-bead
   and sequencing-tether mixture, for 48 hours. Reads are native-DNA pore reads,
   without an Illumina sequencing-primer or sample-index recipe in this branch.

## Molecular representation and limits

The source supports a two-ended Cas9-protected target, background degradation,
heat release, dA-tailing and adapter ligation. These states can be represented
without inventing a guide/target sequence or a specific target length. Protein
protection is a source-declared expectation, not a computed prediction of Cas9
residence time, exonuclease processivity or enrichment efficiency. The article
explicitly reports residual off-target molecules; expected target survival must
not be interpreted as perfect purification.

The complete guide-to-reference cut mapping, guide modifications, complete
LSK109 adapter geometry and complete FFPE-repair substrate effects remain outside
this bounded reading. Their omission from a model does not establish non-disclosure.

## Lineage evidence

🟢 The Results overview explicitly takes downstream library construction from the
adapter-ligation stage of the ONT Cas-mediated PCR-free protocol used with nCATS;
Methods also cites that predecessor for flow-cell preparation. This supports a
**component-scoped link** to nCATS, not inheritance of its dephosphorylation-based
enrichment chemistry. The author response acknowledges Negative Enrichment as
conceptually the same protection/depletion strategy, while distinguishing input,
exonucleases, digestion duration and downstream preparation. Its defining original
must be identified and read before curating a more precise relationship.

## Deferred source targets

- Published S1 Figure (`pone.0241253.s001.tif`): paired SFB/LFB comparison.
- The cited Negative Enrichment original and the specific historical LSK109
  protocol/vendor guide, with their own identities and revisions.
- FFPE Repair Mix M6630S protocol cited in Methods; no substitute current revision
  is assumed equivalent.
- Remaining experimental branches, linked material and independent review.

## Follow-up molecular pass (2026-10-28)

The published article overview 81–94 and Methods 198–226 were independently read
for this model pass; the source-boundary reviewer/author attachments above are
retained as prior-pass evidence, not reread as laboratory protocols. S1 Table was
read at `journal.pone.0241253_supp_pone.0241253.s004.xlsx.txt`, complete
Supplementary Table 1 sheet, guide rows and gBlock entry. The full-bundle scraper
was run again; only published source hits were used. The superseded preprint was
not read for active chemistry. Two originals were read in this follow-up.

🟢 C9orf72 guide entries, verbatim as written in S1 Table:

| Name | Guide sequence | Source target position | Strand |
|---|---|---|---|
| c9orf72_for | TAACGTAGAATAGAACCCGA | chr9:27572027 | + |
| c9orf72_rev | GTTCCCTGGGTGTGAGTCAA | chr9:27576071 | - |

These 20-base sequence entries are not complete manufactured sgRNA molecules.
No scaffold, chemical modification, PAM sequence, exact cut coordinate or target
length is inferred from them. The gBlock is an experimental validation substrate;
this model represents the native genomic branch, not that branch.

The local `ProtectedTarget` state records a duplex occupied at both flanking ends
by guide-directed Cas9 in PAM-in geometry. Its constructor enforces that declared
local scope; it does not compute Cas9 affinity, nuclease accessibility, product
yield or an exact footprint. The graph explicitly carries digested background
products to the later bead cleanup. Residual off-target material may remain in
the recovered library; the waste node does not imply perfect specificity.

The native branch graph passes the structural publication contract for its
**represented scope**, while catalogue status remains notes. In particular, a
passing graph does not establish full source review, exact adapter-end chemistry
or quantitative survival. The local preview is a review artifact. Its LSK109
ligation uses `AdapterLigationGeometry.unavailable`, not an inferred exact T:A
junction. No predecessor constructor supplies nCATS dephosphorylation, LFB washes
or other unsupported defaults.

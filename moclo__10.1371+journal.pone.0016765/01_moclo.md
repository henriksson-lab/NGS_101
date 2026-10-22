# MoClo: modular assembly of transcription units

Weber E, Engler C, Gruetzner R, Werner S, Marillonnet S (2011).
[A Modular Cloning System for Standardized Assembly of Multigene Constructs](https://doi.org/10.1371/journal.pone.0016765).
PLoS ONE 6(2): e16765.

## Source scope

Bounded source readings, 2026-10-21. A protocol-owned graph now renders assembly
from selected, sequenced level 0 modules through level 1 transcription units to a
level 2 construct. This is not a completed missing-information audit.

- Defining article: PMC3041749, verified against its own DOI. Read Results from
  “Technological background and general considerations” through “Design of multigene
  constructs”; Discussion on module libraries; Materials and Methods “MoClo cloning
  protocol”. Local JATS text lines 46–166, 236–254.
- Inspected original PDF page 3, Figure 1 panels A–C, including the fusion-site
  strings and the separate cytosolic and secreted-protein architectures. These are
  graphical source values, not a reconstruction from a primer scraper.
- Follow-up: inspected original PDF page 5, Figure 3, and page 7, Figures 5–6.
  Verified BsaI/BpiI site orientation, nested level 1 cassette boundaries, forward
  positions 1–2, and the matching basic end-linker. The extended end-linker branches
  are identified in the source but are not represented by the current graph.
- Supplement Table S1 (`PMC3041749_supp_pone.0016765.s001.doc`): inspected the module
  inventory and its source/accession column. Full source reconciliation remains open.
- Ran `tools/scrape_primers.py` on the retrieved source directory and explicitly on
  the article PDF text. Both returned zero candidate sequences. This does not
  establish absence: Figure 1 contains short fusion-site strings, and the article
  Methods contain spaced lowercase primer sequences outside the present model scope.

Originals and text twins are under
`$CHEM_DATA/sources/moclo__10.1371+journal.pone.0016765/`.
See [retrieval manifest](ref/MANIFEST.md).

## What is new for this catalogue

🟢 MoClo assembles standardized parts into transcription units, then assembles those
units into multigene constructs. It is a construct-library preparation method rather
than an NGS adapter installation method. The Discussion explicitly describes replacing
one selected part with a library of alternatives.

The model should distinguish a valid declared junction from successful experimental
assembly, clone identity, representation of every intended combination, and expression.
None of the latter follows from a compatible sequence boundary alone.

## Published assembly boundaries

🟢 Figure 1B defines the following module boundaries, written on the same reference
strand. Matching labels identify compatible junctions; each junction belongs once in
the assembled top-strand sequence.

| Module | Upstream fusion site | Downstream fusion site |
|---|---|---|
| Promoter | GGAG | TACT |
| 5′ untranslated region | TACT | AATG |
| Signal peptide | AATG | AGGT |
| Coding sequence following a signal peptide | AGGT | GCTT |
| Terminator | GCTT | CGCT |
| Cytosolic coding sequence, replacing signal peptide plus coding sequence | AATG | GCTT |

The model must keep the two coding-sequence cases separate. A generic illustration of
part roles does not identify a particular biological insert or its sequence length.

## Source-supported workflow

1. 🟢 Prepare level 0 modules. The article describes PCR-added fusion sites and BpiI
   sites, followed by BpiI cloning into the appropriate level 0 destination vector.
   Removal of internal interfering sites is a source-specific preparation choice;
   it must not be assumed to occur merely because a module is declared.
2. 🟢 Assemble compatible level 0 modules into a level 1 transcription unit using
   BsaI and T4 DNA ligase. The assembled unit retains a different, position-dependent
   pair of BpiI fusion sites for its next assembly level.
3. 🟢 Assemble level 1 units with a level 2 destination vector and the appropriate
   end-linker using BpiI. Destination-vector choice and unit position remain explicit.
4. 🟢 A basic end-linker finishes that assembly route. Alternative end-linkers retain
   a selectable cassette and BsaI or Esp3I sites for further assembly. These are
   distinct products; no new assembly round follows from a basic end-linker alone.
5. 🟢 Transformation, clone selection and plasmid examination follow assembly.
   The molecular graph should account for the discarded cassette/backbones without
   claiming experimentally measured proportions of correct and incorrect products.

## Output and readout

🟢 The assembly output is a cloned DNA construct. Article Results report analytical
restriction digestion of selected plasmids and GFP expression as demonstrations.
The Discussion also proposes sequencing an individual colony for verification.

No Illumina adapter, index read, or paired-end sequencing recipe is assigned by this
note. Such a downstream readout would need its own source and declared scope.

## Implemented model and validation

The graph in `tools/moclo.py` owns its part types, selected transcription-unit
architectures, position-specific vector boundaries, enzymes, clone-screening
choices and source inventory. Shared geometry receives those explicit declarations;
it supplies neither a MoClo part vocabulary nor an inferred assembly order.

The example uses the cytosolic and secreted architectures from Figure 1. These are
illustrative source-defined part classes, not particular experimental inserts or a
claim that two source-demonstration clones have these sequences. Payloads, end-linker
body and most vector sequence are role placeholders without invented lengths.

🟢 Original Figure 5 gives forward position 1 as `TGCC` → `GCAA`, position 2 as
`GCAA` → `ACTA`, and pELE-2 as `ACTA` → `GGGA`. The graph uses those boundaries
with pL2-1. It does not import another toolkit's fusion vocabulary.

The independent chemistry checker compares each ligation's molecular inputs to the
explicit geometry and checks its reference-strand product. Adjacent cohesive ends
must pair, fragment order must be retained and a shared fusion appears once in the
product. A different circular drawing origin does not change the molecule. The
product's closed-circle drawing and its opened duplex derive from that same model;
ligation junctions are marked on the strands.

These relation checks are conditional when cores contain unspecified sequence.
They do not establish full-vector restriction digests, end phosphorylation, enzyme
activity, absence of competing products, colony genotype or experimental coverage.
The input-to-fragment digestion transitions remain unresolved in the independent
checker rather than being certified from an enzyme name. Colony colour enriches
candidate clones; it does not prove the full expected sequence.

Graph accounting retains removed donor backbones/dropout material through the
reaction mixture until clone selection. No separate physical purification is
invented between the cleavage and sealing views of a one-pot reaction.

## Implementation and review queue

- The current construction scope begins with sequenced level 0 modules. PCR-based
  part entry/domestication and complete vector-map digestion are not yet modelled.
- Further source-specific vector records and inserted payloads need their own
  identities before the diagram can assert their exact bases or internal cut sites.
- Extend the graph to the published alternating pELB/pELR routes only after their
  complete retained cassettes and next-round cleavage geometry have been assigned.
- Keep a constructed-DNA endpoint. The graph does not designate its output as an
  NGS library or claim a sequencing acquisition.
- Follow directly relevant predecessor cloning sources in a later bounded pass;
  introductory references do not establish inheritance of their construct designs.
- A second independent reader and a complete linked-source checklist remain pending.
  There is no confirmed-missing sequence claim.

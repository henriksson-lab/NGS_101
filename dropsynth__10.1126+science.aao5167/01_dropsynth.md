# DropSynth: grouping oligos for pooled assembly

Plesa C et al. (2018).
[Multiplexed gene synthesis in emulsions for exploring protein functional landscapes](https://doi.org/10.1126/science.aao5167).
Science 359:343–347.

## Source scope

Initial bounded reading, 2026-10-21; a limited construction graph is now implemented below.
Read the defining article's published-version author manuscript, PMC6261299,
own DOI verified: abstract and main-text assembly description, local JATS text
lines 36–50; source/code statement at line 72 and supplementary inventory.
This entry concerns the original 2018 method, not DropSynth 2.0.

Ran `tools/scrape_primers.py` across the retrieved bundle. The bundle contains
automatically fetched preprint files. Their sequences are excluded from active
evidence: the published supplement must be inspected rather than silently
replaced by the preprint supplement.

The published Supplemental Information PDF is identified in the article, but
its PMC download route is gated and the PMC Cloud licence record supplies only
article XML/text, not that attachment. This is an unresolved access route, not
evidence that the original source omits a sequence. Other legitimate publisher
and author-copy routes remain to be checked before calling the source inaccessible.

## Source-supported architecture

🟢 A group of oligos belonging to one intended assembly shares a bead-association
barcode. Complementary bead-displayed strands gather those oligos from the pool.
Emulsion compartmentalization separates assembly groups before release and
polymerase cycling assembly. The emulsion is then broken to recover products.

🟢 The main text describes twelve-base bead barcodes. A bead barcode records an
assembly group; it must not be represented as a cell barcode, a sequencing sample
index, or proof that every compartment contains exactly the intended input set.

The model must represent capture as molecular association. Distinct fragments
remain separate molecules until a source-supported assembly step joins them.
Drawing all captured fragments on one covalent strand would falsely claim a
product before assembly.

## Remaining work

- Retrieve and read the **published** Supplemental Information, keeping its
  identity separate from the downloaded preprint materials.
- Resolve exact nick/release/ligation geometry beyond the already inspected original Figure1 graphic.
- Identify the source-specific fragment, bead-oligo and release-site definitions
  without inferring them from a matching sequence in another synthesis method.
- Distinguish bead association, compartment membership, assembled product identity
  and any post-assembly barcode assignment.
- Resolve source-owned cloning and sequencing/QC branches and their exact scope.
  No sequencing-primer or run-length assignment is made in this first-pass note.
- Check only the explicit fragment/junction and barcode-association relations that
  the supplied molecular states support. Do not infer specificity, representation,
  error rate, protein function or experimental yield from an assembled sequence.
- The source-linked repository is `https://github.com/kosurilab/DropSynth`; its
  relevant original files remain a named follow-up target, not reviewed evidence.
- A second independent reader and complete linked-source inventory remain pending.
  No confirmed-missing-information conclusion is accepted.

## Original figure and repository follow-up, 2026-10-21

Retrieved the university-hosted published PDF `Plesa_Kosuri_2018.pdf`; verified
title, authors, Science359:343–347, year2018 and own DOI. It is the same defining
article as the cached published-version author manuscript, not a newer DropSynth
method. Read the initial assembly description and inspected original Figure1A–D
from PDF page2 using a temporary `pdftoppm` rendering.

🟢 Figure1A separates array oligo amplification/nicking, selective bead pulldown,
emulsion assembly, recovery, and later assembly-barcode addition. Its lower
schematic shows nicking followed by loss of a short strand fragment to expose
the bead-association region. Terminal polarity and exact nick/release positions
are not assigned from that unlabelled schematic alone. The article identifies
Taq ligase for capture and streptavidin bead chemistry in the chosen route.
Its distinct twelve-base capture barcode and later random twenty-base clone
barcode must remain separate features; neither is a sequencing sample index.

The original repository landing page was retrieved and its file inventory/readme
inspected. It identifies oligo-design code and separate amplification-primer,
assembly-primer and barcode directories, not a published wet-lab supplement.
No executable was run, no target-gene sequence was selected and no design or
optimization script was incorporated into the model. Current repository content
is not automatically assigned to the2018 publication without version checks.

Reran `tools/scrape_primers.py` after adding the published PDF. Hits in fetched
preprint files remain excluded; no primer was promoted from a filename or scraper
result alone.

The publisher's old supporting-material page did not return a readable original.
The published SI remains unresolved; the preprint supplement stays excluded.
This pass covers the defining article and linked repository, not complete source
review. Actual exposed capture ends, bead oligos, nick/release chemistry,
fragment overlaps and cloning/sequencing branches remain deferred. A future graph
must keep captured fragments as separate molecules, and must not treat equal
unknown `N` strings as evidence of matching barcode sequences or occupancy.

## Draft implementation and remaining chemistry

`tools/dropsynth.py` now declares a protocol-owned construction graph and renders
native SVG with the normal protocol CLI. The catalogue entry remains `notes`;
its page and unplaced features are excluded from public navigation, site assembly
and the sequence-tool catalogue. This is not the completed sequencing protocol.

The schematic uses two explicitly illustrative capture identities, A/B, within
the source's four-fragment model class. They are not actual barcode values,
target genes or an assertion of exactly two experimental groups. Fragment
classes remain separate through amplification, capture and the compartment
model. The group constructors reject duplicate classes and foreign feature
identities, even though both unknown barcode strings consist of N placeholders.
Those are checks of declarations, not evidence of base pairing or specificity.

The twelve-base capture feature is unplaced: full oligo layout, orientation and
capture duplex geometry are unassigned. It is not drawn covalently between
fragments. Polymerase assembly produces a source-expected product class with
unknown overlap junctions; fragment strings and lengths are not concatenated.
The later random twenty-base clone feature remains independent and is carried
from the source barcode-bearing recipient into an inferred construct arrangement.
No target gene sequence, design algorithm or functional optimization is selected.

Displaced short-strand and bead-remnant fates are explicit unresolved terminals.
They are not silently discarded or assigned to a particular droplet. Source
post-assembly size selection is represented at class level; its precise placement
relative to unexpanded bulk processing still requires the published SI. Cloning
is an explicit unexpanded boundary: the word cloning does not establish a ligase
recipe or a complete vector map. Intended product classes and drawing arrangement
are marked inferred. Sequencing/QC remains an unfinished branch of the goal.

Seven local checks cover grouping invariants, separate identifier meanings and
the two lengths against their separate published main-text clauses. The source
oligo guard has no exact sequences to verify and reports three inputs without
exact molecular states. Chemistry claims remain unresolved; no fully specified
reaction, occupancy, correctness, coverage, abundance balance or yield is proved.

Additional source attempts: the public Europe PMC supplementary archive route
returned a response that failed ZIP validation; the candidate publisher SI PDF
returned HTTP403. Neither response entered ordinary source reading. These are
route failures, not an absence claim. The exact published SI remains the named
next source, alongside applicable bead/primer/overlap and sequencing definitions.

## Supplied published supplement, 2026-10-22 — current generic processing scope

The user supplied the published supplementary PDF `aao5167-plesa-sm.pdf`, its separately
named Tables S8–S14 workbook and the article PDF. DOI/title and supplement inventory
were checked; provenance and hashes are in the external MANIFEST.tsv. The workbook is
received, not an inspected collection of target sequences. Earlier acquisition failures
above no longer describe the availability of this supplied supplement.

🟢 Read the generic DropSynth protocol (printed p4), bead-preparation description p3,
and p5 statement that presence/absence of post-assembly size selection was tested.
Ran `scrape_primers.py` on the supplied supplement; target-specific assay/design and
functional-selection sequences are outside this generic construction pass.

🟢 The protocol specifies initial and bulk pool amplification with KAPA HiFi and15-mer
primers; the bulk forward is biotinylated. Nt.BspQI processing exposes the capture
region. M270 streptavidin beads remove the short biotinylated product, with column
cleanup before association with separately prepared barcoded capture beads and Taq
ligase. Generic emulsion assembly uses KAPA2G Robust, BtsI and Bio-Rad Droplet
Generation Oil, followed by chloroform emulsion break, purification and re-amplification.
The graph now records these reagent identities and the short-strand-removal branch.
Size selection is labelled as a tested option, not a universal protocol step.

The four-fragment A/B classes remain source-model schematic examples, not a declaration
of every experimental oligo design. Equal unknown barcode strings do not establish
physical binding. Bead anchor/ligation/barcode-oligo geometry, both pre-capture PCRs,
post-recovery PCR, exact selected-oligo cut/overlap states and downstream cloning/QC
remain unexpanded. Relevant original figures/tables are supplied and await direct review;
these are modelling/review tasks, not confirmed absent disclosures. Target-specific
sequence design, rescue/functional-selection experiments and activity optimization are
outside this generic graph. The protocol remains notes-only.

## Generic bead geometry and amplification stages, 2026-10-22

🟢 Read the generic SI p4 description and the separate initial/bulk/re-amplification
stages in the detailed emulsion protocol (Steps2,3,10, printed pp12–17). The latter
is explicitly the source's PPAT preparation: no condition is generalized to every
DropSynth experiment. The graph now separates initial amplification, the biotinylated
bulk forward-primer PCR and post-recovery re-amplification. Selected primer values,
cut coordinates and assembled products remain unassigned. No source target sequence,
functional-selection or optimization branch is imported.

🟢 The processed fragment class now shows the source5′12-nt capture overhang beside
a duplex body. Unknown barcode letters are unpaired; no equal-N comparison proves
physical capture. Source classes remain separate throughout recovery and re-amplification.

🟢 Inspected rendered original SI TableS6 p52, including the wrapped anchor and the
separate modification column, and read the bead-preparation description p3 and detailed
bead protocol pp10–11. The anchor is5′ dual-biotinylated; the ligation oligo has5′
phosphate and3′ biotin. The barcode oligo is the12-mer capture-barcode reverse
complement followed by the printed common arm. `Scene` verifies both20-nt common
arms against the source anchor. It shows the splinted nick, exact ligated junction,
subsequent T4-PNK barcode5′ phosphorylation, and streptavidin immobilization/washes.
A/B remain schematic unknown barcode classes; bead occupancy and physical pairing to
a particular fragment are not inferred from those placeholders. Biomolecule/bead
attachment coordinates and stoichiometry are unassigned.

The generic oligo text guard cannot reconstruct the anchor wrapped through unrelated
modification/amount columns or the split12-mer shorthand. Two explicit, source-located
exceptions document this extraction limitation. Independent protocol source checks
reconstruct the anchor from its two source tokens, compare the ligation/common arms,
and verify the12-mer reverse-complement notation and modifications. This is not an
exemption from source transcription verification or a missing-information claim.
Local suite now has12 checks; target-specific assembly and sequencing/QC remain deferred.

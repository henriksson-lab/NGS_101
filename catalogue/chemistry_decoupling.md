# Chemistry ownership migration queue

The [queue](chemistry_decoupling.tsv) records shared code that chooses protocol
chemistry or source evidence. It was inspected on 2026-10-16 after the PARIS,
scCC and Chromium source audits. The initial inventory recorded 47 entries covering 254 protocol
directories with direct imports of the listed symbols. These are static import
relationships, including local build wrappers and tests; they do not establish
that every imported recipe executes. Counts can change as callers migrate. Resolved rows retain their historical
migration scope; SEC-seq is included in the `tenx_rna` row because its dependency
was indirect.

The queue excludes neutral base operations, validators, rendering, canonical
Illumina/Nextera sequences, and simple dataclass/import wrappers. Some listed
modules contain both valid neutral operations and source-specific choices. Migrate
the named recipe or default, rather than deleting the whole module.

Priority 0 addresses shared source-absence claims, complete graph recipes and
implicit final-library selection. Priority 1 moves named oligos, library
assemblies and inherited enzyme/condition choices to their source owners.
Priority 2 reviews smaller operations that still select biological ordering,
lengths, fates or evidence without explicit caller declarations. Chromium recipes
are being decoupled; their queue entry records that work in progress.

The `tenx_rna` migration is resolved for Direct-seq, SEC-seq, TEA-seq and ZipSeq.
Each now owns its endpoint and sequencing declarations, checked against its own
article, supplement or cited vendor protocol. The unused helper was removed after
a repository-wide import check. Direct-seq uses its printed guide-PCR primers;
SEC-seq and TEA-seq distinguish their source-specific gene-expression and feature
architectures. ZipSeq keeps its primer-derived small-RNA endpoint inferred and
leaves the unresolved loaded sequencing primers explicit. This resolves chemistry
ownership; it does not mark unfinished source reviews complete.

The source audits give concrete reasons for this order. Original PARIS used two
library routes; its adapter route needs cDNA circularization, and its original
HeLa primer differs from the later protocol. scCC needs an enrichment PCR before
circularization and end repair/A-tailing/adapter ligation after capture. Their
former generic finish helpers omitted those operations. The lentiCRISPR v2
pooled-library front end uses Gibson assembly rather than the individual-guide
ligation route previously shown. A common endpoint or similar diagram did not
establish equal reaction sequences.

For each migration, the protocol owns its source-specific inputs, modifications,
intermediates, reaction order, side products, selection decisions, sequencing
acquisition and evidence statements. Source-known oligos should enter the graph
as molecular states. A model gap must say what remains unrepresented; it cannot
turn into a statement that the source withholds the information. Direct inheritance
from a predecessor is explicit and source/version checked, with molecular
constructors imported from the predecessor's owning model.

Shared code may still calculate complements, exact primer sites, end compatibility,
strand transformations and graph invariants from supplied declarations. It should
not choose a named protocol's enzyme, adapter, index length, reaction sequence or
source-review conclusion. Similarities should be detected from completed protocol
graphs and their molecular relations, retaining source/version boundaries.

The queue is a migration inventory, not a scientific error verdict or completed
source audit. Each row names the implementation choice to inspect, its direct
importers and a concrete next action. Recheck actual call sites and original
sources before changing a protocol. This inventory does not authorize copying one
recipe's inferred states into other callers.

Seven modules with a sole protocol owner now live under that protocol's `tools/`:
`bliss_chem`, `mario_chem`, `chil_chem`, `rage_seq_chem`, `seqscope_chem`,
`stlfr_chem` and `nano3p_chem`. The move preserves their declarations exactly and
removes those recipes from the shared library. The queue marks this bounded change
`ownership_resolved`; it does not resolve their other uses of generic endpoint or
workflow recipes, nor complete their source audits. Similar-looking declarations
can be compared without making either protocol depend on the other's assumptions.

A subsequent source check recovered the original Smart-seq3xpress index-order
workbook and found ten-base dual indexes where its inherited endpoint had eight.
The xpress endpoint and PCR inputs now use local declarations from that workbook;
its earlier RT/Tn5 reuse remains separate migration work. Tang 2009 likewise uses
fragmentation and selection ranges that differ from its later companion protocol.
These are additional reasons to keep source/version choices at their owners.

The GEM-X 3-prime v4 and 5-prime v3 gene-expression protocols now assemble their
own graphs and use local source-qualified cDNA primers. They no longer import the
shared Chromium reaction recipe or its assumed physical dT adapter end. The
manual Ligation Mix identity remains unresolved; an automated Kit C guide uses a
different part number and does not establish the manual reagent's strands. Other
Chromium callers remain migration work. The scSPRITE audit similarly removed a
shared generic TruSeq endpoint: source tag order, parser orientation and acquired
read lengths are separate claims, and none supplies absent physical adapter bases.
This exposed an additional queue entry, `sprite.final_library`, beyond the initial
47-entry inventory. Both callers now declare their endpoints locally and the
unused shared recipe has been removed. Original SPRITE uses its printed PCR
primer example, one i7 and source read minima; unassembled tag junctions remain
explicit. Its full source audit remains open independently of this ownership fix.

ONT Rapid V14 likewise declares its published RA strand and sequencing plan
locally. The known strand is a graph input; the complete motor/tag attachment is
a separate component-level model boundary. An identical sequence in another kit
does not establish matching bottom strands or physical assembly.

The two TruePrime input workflows now own their TthPrimPol/phi29 reaction
declarations, amplification durations and reported Covaris/NEBNext continuation.
The purified-DNA branch retains its three-hour schedule, while the intact-cell
branch retains three or six hours and its separate L2 lysis/version boundary.
The shared `primpol_mda` recipe was removed after checking all actual callers.
Unresolved downstream kit and primer geometry remain model limits; this
ownership change does not complete either source audit.

GRO-cap now owns the inherited RT, RNA-template digestion, Phusion PCR and PAGE
route, with local source-derived conditions and an unassembled endpoint. The
former shared shell assigned an eight-base index without evidence; even the
six-base RPI family printed in the adjacent GRO-seq method does not assign an
index to GRO-cap. The original mapping method supplies an RA3 trimming target,
which is distinct from a complete modified ordering oligo. Source transcription
checks read those originals directly. The other nascent-RNA callers remain open.

Kinnex no longer imports a generic template-switch illustration. The independently
checked SMRT Link instructions identify its printed sequences as primer-removal
references, including reverse-complemented 3-prime orientation without poly(T).
They are no longer labelled physical PCR primers. The endpoint also marks its
unknown ligation junction explicitly instead of implying blunt ends. Research
MAS-Seq and commercial Kinnex remain separate source declarations.

The inventory now also tracks the SMART predecessor's source-selecting
`page_parts`, not only modules under `lib/`. FLASH-seq and Smart-seq3xpress have
local declarations after fresh original-source checks. FLASH removes inherited
Smart-seq2 ISPCR/LNA choices: the original figure visibly crosses out ISPCR, a
detail lost in extracted text. Its standard branch declares SSIV/KAPA together,
commercial XT/NT/NPM, an illustrative source-table index pair, and 75/8/8/0.
Commercial ME assignments remain explicitly inferred. Xpress keeps its own
improved TSO, forward-primer overhang, ten-base indexes and source-annotated8N+2W
identifier. Remaining SMART callers still need review.

CODEC now owns its published quadruplex, products and amplification route. A
whole-oligo source check found an extra GC at the old shared bridge7 spacer/core
boundary; its local 98-base strand now matches original Table5. NanoSeq's separate
end-repair view moved to its own directory unchanged, and `lib/duplex_error.py`
was removed. This resolves that module's ownership, not every shared endpoint or
either protocol's complete source scope. dm-tRNA similarly owns its actual fusion
RT, RNA/DNA starter, RNA hydrolysis, DNA circularization and PCR; the other
structured-RNA callers remain migration work.

TGIRT now owns its original 2016 starter and single-index route. Its later 2021
starter, optional UMI and dual-index choices remain separate source versions. The
unused shared group-II starter constructor was retired; its blocked-RNA and
single-overhang restrictions do not establish universal enzyme requirements.

Duplex Sequencing owns its corrected-SI adapter strands, six-position index and
CAGTA spacer. ChEC likewise removed the shared Henikoff and TruSeq recipes after
checking the actual 2011/2013 predecessors and separate Orsi procedure. Its
illustrated Orsi endpoint preserves the whole 58-base Universal Adapter without
the old duplicated ACAC overlap. Source-owned steps distinguish negative ChEC
selection from positive library cleanups. An illustrated adapter variant is not
an identified original sample assignment; unresolved historical sources remain
open independently of these ownership fixes.

ORGANIC and CUT&Tag now own their separate source-qualified preparation routes,
adapter states and product fates. The unused Henikoff recipe was removed after
ChEC/ORGANIC migration and removal of NICE's unused imports. The 16S declaration
now extends its actual PCR1 product, retaining both locus-primer regions; the old
unused `targeted_ngs` shell was retired. Canonical PCR2 caps remain inferred until
the historical whole-primer identity is assigned. These changes resolve specific
ownership and model errors, not the protocols' remaining linked-source audits.

The queue also identifies `nextera.TAGMENTATION_OUTCOMES` and `amplifiable`: their
blanket rejection of same-end products needs caller-specific review. PCR product
formation, flow-cell compatibility and later filtering are separate claims. This
recipe-level issue does not make canonical Nextera sequence constants invalid.

NicE-seq and CUT&RUN now keep their unassigned study-specific adapter and index
assemblies local, without filling them from a generic TruSeq endpoint. NicE-seq
records the source aliquot and unresolved remainder; CUT&RUN separates soluble
library DNA, optional pellet QC and other extraction variants. CUTAC instead
declares the cited in-house adapter pool and a conditional custom Buenrostro PCR
family. Its warning against commercial Nextera/NEBNext PCR primers is retained;
shared adapter bases alone do not identify the ordered reagent. The source
reviews remain open for identified linked material and experimental assignments.

The nascent-RNA audit also tracks `rna_cell_assays.token/custom_library`: their
placeholder length defaults can appear to assign measured lengths to unspecified
RNA and library arms. ChRO now owns the source-specific reaction route and
unknown-length roles; its six-base UMI is a designated-library variant. The other
static importers still need actual call-site and original-source review. PRO-seq's
2013 original supplement is a distinct source from its 2016 detailed protocol;
same-name primer families do not establish an unchanged historical reagent.

GRO-seq now owns the recovered 2008 source's direct adapter-ligation route, with
5′ ligation first, three anti-BrU captures and 33-base reads. PRO-cap and mNET-seq
also removed the shared eight-base index and RT/PCR finish recipes. The original
PRO-cap AP/TAP route differs from its later exonuclease variant; mNET-seq uses
on-bead phosphatase-minus PNK before RNA extraction. Those distinctions remain
local even where the named downstream kit families overlap. Source reviews remain
open for identified originals and linked material; ownership is not a full audit.

The inventory separately tracks the former `nascent_protocols.FourSUReadout` and
`NascentEndCapture`. These constructors selected molecular states while leaving
population versus RNA-tract meaning implicit; a shared U-to-C readout also assigned
a generic C product. The original TimeLapse
source describes a cytidine analogue; matching acquired bases do not establish
matching chemical products. Each caller needs its own population members,
modified sites and end states. The shared RNA-library `molecule` constructor's
eighteen-base defaults are also queued separately from its reaction recipes.
These are implementation risks to review against each source, not evidence that
every caller's chemistry is wrong.

CAGEscan, Start-seq and TimeLapse now own their reaction routes and endpoints.
CAGEscan retains its eight printed oligos, inline three/six-base TSO barcode,
two pilot qPCRs and separate production PCRs. Start-seq keeps independently
deposited 5′- and 3′-reading preparations distinct. TimeLapse represents its K562
total-RNA route with a source-qualified, conditional dated-kit continuation;
the modified C* site is distinct from both canonical C and the sequencing
readout. The unused `NascentEndCapture` constructor was retired.
These migrations leave the protocols' outstanding source obligations open.

The remaining sf4sU, SNU and TUC callers now own their declarations too. sf4sU
preserves HPDP capture/release, gel recovery and a conditional NEB guide route,
including its anti-dimer hybridization and the study-observed adapter-dimer class.
SNU distinguishes the original 2021 Ion route from later Ion/FWD implementations,
with RNA removal and random second-strand synthesis before PCR. TUC retains its
source-positive native-C conversion while its unrecovered original handling and
library details remain explicit warnings. No later DUAL method is substituted.

All users of the shared nascent-library and readout recipes have migrated. The
unused `lib/nascent_protocols.py` was removed after checking Python references;
its last neutral presentation fields moved to GRO-cap with byte-identical HTML.
These two queue entries are resolved for ownership. Their source audits and the
other shared recipe entries remain open. RNA-class labels describe pulse
incorporation: an intact labelled transcript may contain a pre-existing 5′ tract,
so capture does not date every base of the molecule.

The same boundary applies within `protocol_model`: its legacy `role_library`
selects Illumina endpoint roles, and `da_tailed` uses an implicit end
representation. These are separately queued for their actual caller paths and
source-defined topology. The generic module name does not establish neutrality;
`role`, supplied-part drawing and graph continuity remain useful shared operations.

A 2026-10-17 caller check retired shared code that no protocol still imported:
`lib/spatial_epigenome.py` (its tagged fragment omitted the Tn5 transferred strand),
`lib/antibody_tags.py`, `damid.damid_adapter`, `rare_cell_and_affinity._finish_pcr`,
`crispr.SCAFFOLD_FE` with its unused `scaffold_of`, the shared ONT cDNA primers,
`molecular_consensus.SingleStrandUID`, the two-step `pooled_crispr` intermediate and
graph, the legacy ATAC and GeoMx `workflow()` recipes, and seven uncalled telomere
graph recipes. Spatial-epigenome and the rare-cell PCR finish are resolved. SPATAC now
supplies its GEO nuclei route and PE150 acquisition to the split-pool graph. The
telomere adapter drawing now marks the adapter 3-prime end at the native C-rich
5-prime end. Removing dead code does not complete any caller's source audit.


## Initial validation migration snapshot — 2026-10-18

The latest approach combines protocol-owned typed provenance graphs with the
independent narrow molecular checker described in DESIGN.md. Passing graph
completeness is not evidence that molecular transitions are validated. Conversely,
a deliberately unresolved relation is not automatically an architecture defect.

Three TSV snapshots list protocol names, directories, code locations, queued
dependencies, graph/check coverage and next actions:

- [Migration worklist](validation_migration_pending.tsv): **182 protocols**.
  **181** reference currently queued shared recipes/defaults in non-test code,
  directly or through the inspected TruSeq PCR-free wrapper; **CRISPR-MIP**
  additionally needs a native typed graph. Priorities inherited from the ownership
  queue are P0: 122, P1: 46, P2: 13; the graph-less CRISPR-MIP is listed first.
- [Molecular-coverage follow-ups](validation_migration_coverage_followup.tsv):
  **190 graph-backed protocols** whose current checker emits no established or
  conditional claim. This overlaps the ownership worklist; it is not another 190
  disjoint migrations. Missing source evidence, unsupported molecular states and
  absent checker primitives require different follow-up decisions.
- [Whole-catalogue inventory](validation_migration.tsv): **410 protocols**,
  including six ledger-only ownership entries that needed reconciliation (migrated
  2026-10-19, round 2026-10-19-14),
  notes-only entries, non-sequencing assays and partial molecular coverage.

Method: loaded each catalogue protocol in an isolated Python process and ran the
existing chemistry checker; no contradictions or import errors occurred. Compared
open `chemistry_decoupling.tsv` rows against non-test Python imports/references,
with file/line pointers. The SMART-seq3 predecessor import was inspected separately;
its queued alias has no callers listed in the old ledger. TruSeq PCR-free reaches
`batch_ngs.truseq_pcr_free_library` through the vendor wrapper and still selects
full-length-adapter defaults locally. Static references are evidence of code
dependencies, not proof that every referenced branch executes or that its chemistry
is erroneous. The snapshot is seeded by known queue entries, not a proof that all
unlisted shared operations are neutral. It does not authorize source-independent
recipe substitutions. Commit identity and observation date are columns in each TSV.

The old queue names 186 distinct protocols. The reconciled migration worklist
retains 180 of them, adds SMART-seq3, and adds graph-less CRISPR-MIP. The six old
entries without a current non-test reference to their queued symbols are
Chromium Single Cell ATAC v1, Duplex-Repair, LamPORE, Olink Explore, PEM-seq and
TIME-seq. Their stale import entries need reconciliation; this does not declare
those protocols fully validated. Historical ownership rows are preserved.

Runtime coverage: 402 native graphs; 201 have some established/conditional local
checks but no declared `MolecularRelation` transition; 190 have only unresolved
claims; 11 declare at least one checked transition relation. Those 11 are
SMART-seq2, in situ Hi-C, Nagano single-cell Hi-C, Stevens single-cell Hi-C,
FISSEQ, ASTAR-seq, SPLiT-seq, microSPLiT, sci-ATAC-seq (2018), MALBAC and smMIP.
All 402 graphs still contain unresolved claims. Exact PCR binding or paired-column
checks are legitimate typed checks without a `MolecularRelation` flag, so flag
absence alone is not a migration verdict. No protocol is marked fully validated
by this inventory.

CRISPR-MIP is the single documented NGS page without a native graph. The other
seven graph-less entries are three notes-only protocols (CapSeq, SELECT-seq and
SIMPLE-seq) and four non-sequencing assays (superRCA, RNAscope, PrimeFlow and
PLAYR); keep their respective modelling/assay scopes separate. This snapshot
changes neither source-review acceptance nor the capped-curation census.


## Current migration progress — 2026-10-18

The user-selected cohort is the initial 182-protocol worklist.
[Pending migrations](validation_migration_pending.tsv) now contains **0**
protocols; [completed migrations](validation_migration_completed.tsv) records
**189** verified migrations (the 182-protocol cohort plus seven later entries, see
below) with their former dependencies and current scoped
coverage. Entries leave the
pending file only after their source-qualified declarations and affected checks
are verified; a temporarily removed import alone does not complete a migration.

The completed record contains all 177 protocol names and per-protocol coverage.
The batches cover CRISPR-MIP, three Chromium protocols, five multimodal protocols,
PLAC/SNuBar, three editing-outcome methods, two nuclease protocols, four extended
RNA workflows, three break-capture methods, four RNA-end methods, three ATAC
methods, DMS-MaPseq and PacBio Iso-Seq. Additional migrations cover Multiome, Flex,
Visium, BotSeqS, ACE-seq, DRIP-seq, AlkAniline-Seq and CirSeq. Further completed
entries are BLISS, CRISPR CIRCLE-seq, CITE-seq, Direct-seq, cfMeDIP-seq, CODEC,
ChAR-seq and ChIL-seq. Further migrations cover Digenome-seq, END-seq, eccDNA Circle-Seq,
DBiT-seq, DR-seq, ECCITE-seq, DART-seq, FragSeq, HRF-seq, E/L Repli-seq
and DRIPc-seq. Batch02 adds Drop-BS, GAM, Long-read ChIA-PET, easySHARE-seq,
G&T-seq, GAGE-seq, MATQ-seq, mim-tRNAseq, Nano3P-seq, LIGR-seq and MARIO.
Batch03 adds DNA Prep shotgun metagenomics, NanoSeq, NG Capture-C, INVADEseq,
OAK, Paired-seq, OTTR-seq, PANDORA-seq, PARS, oxBS-seq and Ovation RNA-Seq V2.
Batch04 adds MALBAC-DT, Paired-tag, PARIS, R-ChIP, RAGE-seq, RAISIN RNA-seq,
RamDA-seq, ReDeeM and RIC-seq.
Batch05 adds scCARE-seq, sci-CAR-seq, scRICA-seq, scM&T-seq, scMT-seq,
scNMT-seq, SHERRY, SLAM-seq and snRandom-seq.
Batch06 adds scSNV-seq, scTrio-seq, SHARE-seq, snmC-seq2, snmC2T-seq,
sn-m3C-seq, SMART-seq, SMART-seq3, Structure-seq2, scSPLAT, Slide-tags and SPLASH.

The migrations include actual molecular/source corrections, not just relocated
helpers: declared precursor strands and PCR templates, index/UMI provenance,
source-specific recovery and cleanup steps, and hydrolysis versus physical
discard. CRISPR-MIP now has a native graph built from the existing published-data
constructors. Its initial P7 copying and subsequent P5 copying are distinguished
from later duplex PCR within the same reaction. Missing graph/checker machinery
is an implementation gap, not an assertion of missing publication information.

Independent review corrected premature template-switch/barcode installation in
the VDJ precursor and false initial duplex primer availability in CRISPR-MIP.
The neutral PCR renderer now preserves a source-limited single-stranded Scene
without drawing an invented complementary strand. Further independent review
separated cellular containers and independent molecular populations that had
been falsely concatenated on one strand. BotSeqS now propagates an unresolved
Y-adapter precursor into both PCR branches, rather than a completed duplex.
Source-review targets are explicitly curated rather than generated from graph nodes.
The next reviews corrected CODEC SID7 orientation against its source bridge and
the retained nine-base bead-handle suffix in Direct-seq. ChAR-seq now preserves
RNA through contact ligation and follows the actual USER/PCR/selection states;
ChIL-seq uses its actual single-strand RT precursor. CIRCLE-seq hydrolysis
products remain molecular material rather than invented physical discards.
Further review corrected copied ECCITE tag-barcode polarity, separated the later
DBiT implementation from the defining experiment's read acquisition, removed
unsupported reuse of DR-seq's initial RT enzyme in its CEL-Seq finish, and restored
FragSeq's published PNK/gel order, post-PCR selection and 50-base read statement.
RT RNA-template fates remain explicit where physical recovery is unresolved.
Batch02 restores MARIO's source-selected ES-1 DSN branch and actual PCR inputs,
retains physical blocked nucleotides, corrects GAGE primer/index transcriptions,
and derives Drop-BS's blunt duplex from paired strands. LIGR's published SPRI
size wording and final PAGE threshold remain visible without a guessed repair.
A scoped MARIO circle-order relation preserves source-defined part identity; it
does not establish complete circularization, restriction or PCR chemistry.

Batch03 separates source-version reference implementations from historical executed runs,
removes unsupported canonical adapter shells/index widths, preserves identifier lineage
through intermediate states and supplies actual native PCR inputs. Matched ordinary BS
is represented in the explicitly selected 2013 oxBS reference; original 2012 oxRRBS
construction remains deferred. Independent review corrected dA overhang geometry,
control continuity and displayed primer transcription; physical blocked nucleotides
and mixed off-target populations are retained. Removing unsupported exact assumptions
can increase unresolved molecular coverage even as ownership migration completes.

Batch04 preserves native RNA/template populations, physical terminal modifications,
source-defined library routes and actual ss/circular/capture-bound PCR substrates.
RamDA models coupled RT/nicking/displacement and the selected standard plate finish;
MALBAC-DT separates original Nextera from automated TruePrep without inheriting unnamed
formulations. R-ChIP selects the detailed 2019 reference, with its actual partly duplex
ligation precursor and separate test/final PCR aliquots; its original 2017 execution
is not inferred from that reference. ReDeeM's known RNA reaction sequence is explicit.
Four unchanged DNA entries remain deferred in the batch04 ledger after worker tool
filter failures; they are retained in the pending file, not counted as completions.

Batch05 replaces completed-library intermediates with actual ss-cDNA, converted DNA,
native tagged fragments and copied-product populations. Separate RNA/DNA routes,
source-specific oligo modifications, native recovery stages and sequenced primer-derived
products remain explicit. Unverified kit geometry and index acquisition remain bounded;
SLAM's unread QuantSeq construction and scRICA's delegated low-throughput ATAC preparation
remain genuine source boundaries. Parent review restored scRICA's full structured RT
primer/RNA-tailed TSO and scCARE's native tagged cDNA end. The batch05 worker filter stop
left Seq-Scope, SHAPE-MaP and SHAPE-seq deferred; together with batch04, seven entries
remain pending under task-execution deferrals, not disclosure findings. No filtered task
was retried or bypassed.

Batch06 preserves separate sampled populations in scSNV, independent native RNA/DNA
and original/copy lineages, exact modified oligos and native ss/circle PCR precursors.
SHARE seals the three preassembled barcode rounds once and separates analytical qPCR
from production libraries. snmC-seq2 retains source ten-base index definitions without
canonical arm boundaries; snmC2T separates chemical methylation from copied-base signal;
sn-m3C retains DpnII cohesive ends through inactivation and physical converted U.
SMART-seq's published single-end acquisition, SMART-seq3 retained-TSO/internal classes,
Structure-seq2's retained RT hexamer and physical index versus acquired cycles are local.
scSPLAT and SPLASH include source-known successive recovery stages and explicit original
strand/circle fates. Slide-tags uses separate experimental allocations and actual mixed
SPRI supernatants; known CG000338 operations are explicit, while unread CG000317 details
remain bounded. Independent implementation reviews found and corrected these issues;
implementation-only readers are not counted as original-source absence readers.
No new tool-filter deferral occurred in batch06; the earlier seven remain pending.

The initial snapshot above is historical: runtime inventories now contain 403
graph-backed protocols and 193 with only unresolved claims. Completed ownership
entries may remain in molecular/source-coverage follow-up. The original goal
remains active until every selected pending migration is completed and verified.
Details and source scopes: [migration round](migration_rounds/2026-10-18-01.tsv).


### Batch07 native graph migration

Twelve more protocols have verified native, source-qualified declarations: stLFR,
TAB-seq, TSA-seq, TrAEL-seq, Tri-C, ZipSeq, NEBNext Enzymatic Methyl-seq v2,
NEBNext Low-bias Small RNA, NEBNext Single Cell/Low Input RNA, TruSeq DNA
PCR-Free, SurePlex and YAMAT-seq. See the [batch07 ledger](migration_rounds/2026-10-18-07.tsv)
for actual capped source sections and separate independent implementation reviews.

Critical changes retain physical forked PCR-free adapters, RNA U and dideoxy
terminal residues, carrier/barcode identity through intermediate states, original
versus copied converted bases and actual native PCR substrates. Source-known
recoveries are explicit, including the two SPB rounds, sequential RNA/cDNA
bead bindings, SurePlex ss normalization and YAMAT native PAGE. stLFR physical
circle/DNB material stays distinct from its readout reference. ZipSeq owns native
GEX/zipcode branches with an explicit unread post-GEM-recovery boundary. Native
index widths and acquisition remain unassigned where the inspected originals do
not support a choice; this does not establish publication absence.

All12 focused suites/pages/scoped chemistry/source guards/review schemas pass
after critical reviews. Pairing and P5 lint pass, as do shared1277 and catalogue83
checks. Site builds407/407 diagrams, publishes406 schematics, satisfies402/402
documented NGS graph contracts and resolves all relative links. Current runtime
coverage is403 graph-backed protocols,199 with only unresolved claims and no
contradictions. The original182 cohort partitions into123 completed and59 pending;
the seven earlier filter deferrals remain pending. Ownership completion does not
close molecular/source follow-ups. The goal remains active.


### Batch08 native graph migration

Twelve more migrations are verified: Methyl-HiC, Nagano single-cell Hi-C, Stevens
single-cell Hi-C, CEL-seq, MARS-seq, MARS-seq2.0, NEBNext Ultra II Directional RNA,
NEBNext Ultra II FS DNA PCR-free, Nextera XT, ChEC-seq, NicE-seq and ORGANIC
profiling. The [batch08 ledger](migration_rounds/2026-10-18-08.tsv) distinguishes capped
source readings from independent implementation-only critical reviews.

Native IVT graphs retain paired T7 promoters, initialGGG, physicalRNAU, actual
RT2/PCR inputs and original-versus-copy fates. Hi-C graphs retain their different
restriction enzymes/order, native tagged Y substrates and dT/A junctions; primer
phosphorothioate bonds remain distinct from template copying. Commercial graphs
retain dU/USER and source forks/ss normalization, with separately chosen oligo-kit
and source-revision limits. Source-positive adapter installation is explicit.

ChEC/ORGANIC native graphs now keep the later Orsi oligos as a separate reference,
without assigning historical native reagents or index acquisition. ChEC's cited
library preparation remains a bounded unread-source stage; its RNA/carrier fates
stay unresolved. ORGANIC follows rendered defining FigS1a/c, including S1/S2 and
sequential native recoveries. NicE uses actual bead-bound original templates and
keeps internally biotinylated originals distinct from PCR copies.

All12 focused suites/pages/chemistry/source guards/review schemas pass. Guards
check13 strands with zero failures, two historical narrowly scoped exceptions
and24 inputs without exact states. Pairing/P5 lint pass. Runtime coverage remains
403 graph-backed protocols, with201 emitting only unresolved claims and no
contradictions. The increase from199 follows removal of unsupported native Orsi
assignments, not new evidence of source absence. Ownership completion stays
separate from molecular/source follow-up. The original182 cohort is135 completed
and47 pending; the seven prior tool-filter deferrals remain pending. Whole-site
verification passes:407/407 diagrams,406 published schematics,402/402 documented
NGS contracts,813 files and no broken relative links; catalogue83 passes. See
AUDIT.md for builder-discovery corrections and compatibility checks.


### Batch09 native graph migration

Twelve migrations are verified: ChIPmentation, CoTECH, MulTI-Tag, HDST,
Slide-seq, sci-Space, FIND-seq, GoT, TARGET-seq and the three Oxford Nanopore V14
cDNA-PCR, ligation and native-barcode kits. The [batch09 ledger](migration_rounds/2026-10-18-09.tsv) records each capped
source pass separately from independent implementation reviews that read zero
additional originals.

Native graphs retain actual incoming PCR material, independently allocated
RNA/DNA/control populations, identifiers and physical versus copied modifications.
Critical corrections include Slide decoding-primer inputs and native initial
extensions; HDST dU/IVT/RNA and copied VN-anchor handling; sci-Space independent
hash/transcript lineages and selected nuclear cohort; FIND matrix-linked originals
versus free WTA copies; TARGET mixed aliquots and the explicit PCR1-to-PCR2 sample;
and the three Nanopore kits' real adapter geometry, allocations, digestion products
and associated motor boundaries. Ligation leaves unadapted molecular classes
unresolved through later recovery instead of implying complete adaptation.

CoTECH does not promote superseded preprint order: published Table1 components
are conditional native declarations with unread execution/allocation boundaries.
MulTI-Tag keeps its source final joint tagmentation and excludes future conjugate
lineages from raw input. GoT's selected linear graph remains distinct from its
source-known circular/Gibson variants, whose complete native topology is a named
follow-up. Later Nanopore sequence references do not prove historical reagent-lot
identity. No missing-information claim or source-review status was upgraded.

All 12 focused suites/pages/chemistry/source guards/review schemas pass following
critical reviews. Source guards check 53 strands with zero failures,
3 narrowly recorded transcription exceptions and 35 inputs without
exact states. Shared 1277, source-fetcher 226 and catalogue 83 checks pass. The exact
cached sci-Space grid PDF URL is now registered alongside its cited DOI; its printed
DOI/ID/date, rather than its normalized download slug, establish the inspected
revision. No oligo exemption was needed for that locator repair.

All 410 runtime entries reconcile: 403 graph-backed protocols, 193 with only
unresolved claims, no errors or contradictions; 390 still lack declared transition
relations. These are scoped model checks, not whole-protocol validation. The 57-row
ownership queue remains intact. The original 182 cohort partitions into 147 completed
and 35 pending, including the seven prior tool-filter deferrals. The overall goal
remains active; nothing staged or committed. Whole-site verification passes: 407/407 diagrams build, 406 finished schematics
publish, 402/402 documented NGS graph contracts pass, and all 813 site files have
resolving relative links. Pairing and P5 lint have zero findings (the existing
source-known PEM exemption remains). All 777 note pages are current.


### Batch10 native graph migration

Twelve more protocols now have verified native migrations: Tang2009, Quartz-Seq,
Quartz-Seq2, Fiber-seq, SAMOSA-Tag, ScISOr-seq2, HiRES, LiMCA, iMARGI,
sciPlex-ATAC-seq, spatial-ATAC-RNA-seq and spatial-CUT&Tag-RNA-seq. The
[batch10 ledger](migration_rounds/2026-10-18-10.tsv) records bounded source passes separately from independent
implementation reviews that read zero additional originals.

Critical review corrected Tang's native P2 orientation and retained distal TT
before PCR; Quartz-Seq2's opposite-strand adapter ligation and original terminal
dA; SAMOSA's post-repair molecular mixture; ScISOr's capture-reagent timing;
and the spatial assays' free molecular classes after bead removal. Native graphs
keep actual incoming PCR material, original versus copied modifications, independent
RNA/DNA/hash populations, physical allocations and source-known recovery steps.
Completed product references never supply an unknown incoming PCR substrate.

The spatial DNA reference uses source N501/N701 members conditionally; native RNA
library ends and index acquisition remain unassigned. sciPlex's published capture
LNA, gel-window and hash-read discrepancies stay public. Fiber-seq does not inherit
a later SMRTbell kit. LiMCA preserves the source 13bp duplex/6nt overhang rather
than repairing it by appending bases. The explicitly linked Wu2019 iMARGI original
now enters the source guard from its existing canonical DOI cache; no source copies,
recursive lookup or oligo exemption was added. Synthetic tests also exclude
undeclared neighboring documents.

All twelve focused suites, pages, chemistry checks, source guards and review schemas
pass after critical review. Source guards check 73 strands with zero failures or
exceptions, one absent dedicated bundle skipped (Quartz-Seq2 checks the actual
shared original in its own suite), and 16 inputs without exact states. Shared 1280
and source-fetcher 226 checks pass. All 410 runtime entries reconcile: 403 graph-backed,
193 with only unresolved claims, no errors or contradictions; 390 still lack declared
transition relations. These scoped results do not establish whole-protocol validation
or accept a missing-information claim. The 57-row ownership queue stays intact.
The original 182 cohort now partitions into 159 completed and 23 pending, including
seven prior tool-filter deferrals. The overall goal remains active; nothing staged
or committed. Whole-site integration passes: 407/407 diagrams build, 406 finished schematics
publish, 402/402 documented NGS graph contracts pass and all 813 site files have
resolving relative links. Pairing and P5 lint have zero findings (the existing
source-known PEM exemption remains). Catalogue 83 checks pass; all 777 note pages
are regenerated after the final checkpoint update.


### Batch11 native graph migration

Ten further migrations are verified: CAP-seq, icSHAPE, LAST-seq, PAIso-seq,
Prox-seq, scGESTALT, scGRO-seq, SDR-seq and lentiCRISPR v1/v2. The
[batch11 ledger](migration_rounds/2026-10-18-11.tsv) preserves each capped source pass and independent implementation
review separately. Forty source-transcription checks pass; source guards check 31
strands with zero failures, skips or exceptions and 20 inputs without exact states.

Critical corrections retain icSHAPE's physical terminalddC and multipart degenerate
UMI; LAST's original RNA/DNA/IVT lineage and unassigned index acquisition; PAIso's
native single-stranded PCR precursor and distinct size-selection portions; Prox's
actual plate-specific DNA-only preparation and constrained UMI encoding; scGESTALT's
real 5µl allocation and duplex T7 input; and scGRO's source chemical residue versus
ordinary copiedDNA, original populations and single-stranded native PCR material.
SDR preserves its source V2 read allocation and public source discrepancies.

The v1 screen no longer receives a later compatible Joung library as its own native
endpoint: its original unread SupplementaryMethods remain an explicit boundary.
The v2 source plasmid-construction and independent duplicate validation-PCR route
retains actual synthesis-pool inputs, 5OH backbone/filler lineage, the 0.5µl assembly
aliquot and independent technical QC barcodes. Guide identifiers persist through
unplaced native inventories. A 150-cycle kit does not assign acquired read/index
cycles. These graphs expose unresolved source/molecular transitions rather than
claiming complete chemistry validation or accepting publication absence.

The batch selected 16 protocols before reading. The worker handling four protocols
was stopped by the content filter: Digital Telomere Measurement, scRRBS, meta3C and scDam&T-seq.
Independent review of the partial Joung/Moffat implementations was also stopped.
Those six remain pending; no retry, reassignment or migration-completion claim is
recorded. Partial code is retained with compatible existing build entry points.
Together with the seven earlier tool-filter deferrals, 13 entries remain pending.
The original 182 cohort partitions into 169 completed and 13 pending, with all 57
ownership-queue rows retained. All 410 runtime entries reconcile: 403 graph-backed,
193 with only unresolved claims, 390 without declared transition relations and no
errors or contradictions. Completion remains unproven: thirteen original-cohort migrations/reviews are
incomplete. The same automatic content-filter blocker persists across the Batch09,
Batch10 and Batch11 checkpoints; no remaining migration can proceed without an
external change. The goal is blocked, not complete. Nothing staged or committed.
Whole-site verification passes: 407/407 diagrams build, 406 finished schematics publish,
402/402 documented NGS graph contracts pass and all 813 site files have resolving
relative links. Pairing and P5 lint report zero findings (the source-known PEM
exception remains); catalogue 83 checks pass. All 777 generated note pages are
current after this final checkpoint update.


### Batch12 migration checkpoint

Eight further migrations are verified: Joung one-PCR, Moffat two-step, LAMP-Seq,
Seq-Scope, SHAPE-MaP, SHAPE-Seq, scDam&T-seq and Digital Telomere Measurement.
The [batch12 ledger](migration_rounds/2026-10-18-12.tsv) distinguishes bounded
source reading from independent implementation reviews, which read no originals.
Native incoming substrates, independent molecular populations, modification and
identifier lineage, actual sampling/pooling, and physical cleanup fates were checked.
Review fixes removed nonphysical identifier strands and duplicate original-population
outputs, retained possible residual classes in carried mixtures, and qualified
unassigned read acquisition and source-version boundaries. Unknown chemistry remains
explicit rather than supplied by a canonical endpoint.

The original 182 cohort now partitions into **177 completed and 5 pending** (three
implementations awaiting independent review). PTA, scRepli-seq and scRRBS have native
implementations and passing author checks; the combined independent-review task hit
a content filter before any per-protocol review result. Its triggering protocol was
not identified, so none of those three is counted complete. REPLI-g Single Cell is
blocked on exact historical August 2014 RG35 and Nextera 15027987 Rev B October 2012 sources;
the available 2018 Advanced DNA kit is not a revision substitute. The meta3C worker hit
the content filter and returned no protocol changes. These stops do not establish
missing published information. Work moved to other protocols after each stop.

All 410 runtime entries reconcile: 403 graph-backed, 196 with only unresolved claims,
390 without declared transition relations, and no errors or contradictions. All 57
ownership-queue rows remain. Completed migration means verified local ownership and
scoped molecular checks; whole-protocol validation and source-absence obligations
remain separate. Current pending and completed files retain the original 182 scope.

Final batch12 verification: focused suites/pages, source guards and review schemas pass; catalogue 83 checks pass. Whole-site build passes 407/407 diagrams, 406 finished schematics, 402/402 documented NGS graph contracts and all 813 site-file relative links. Pairing and P5 lint have zero findings beyond the existing source-known PEM exemption. Changes remain unstaged and uncommitted.

## Final ownership-queue entries — 2026-10-19

After the cohort closed, six protocols outside it still sat in the ownership queue:
10x Chromium Single Cell ATAC v1 (`nextera` outcomes), Duplex-Repair, LamPORE, Olink
Explore, TIME-seq (`batch_ngs`), and PEM-seq (`batch_ngs`, `protocol_model.da_tailed`).
An indirect-import scan added UMI-4C. It executed `da_tailed` through the
`conformation_replication_specs` re-export, which the static direct-caller count missed.
All seven are migrated, each with an independent implementation review
([round 2026-10-19-14](migration_rounds/2026-10-19-14.tsv)). The `batch_ngs`,
`nextera` and `protocol_model.da_tailed` rows are now `ownership_resolved` with no
direct protocol callers; no queue row remains `queued` or `in_progress`.

A runtime trace over every protocol suite showed that no protocol executes the
`batch_ngs` library constructors or `da_tailed`. `nextera.amplifiable` is called only
by the shared `checks.run_common` invariant test. Several legacy `lib/` modules still
import these constructors, so they are not retired.


## Post-migration code cleanup (2026-10-20)

Retired 21 unused shared modules: the 18 initial no-import candidates, their orphaned
`pooled_crispr` dependency, and `rna_library_operations` / `replication_labeling` after
the unused fallback and imports were removed. Reference checks include Python
imports/re-exports, dependency chains and executable Markdown expressions. The
Markdown loader's generic discovery is not evidence of an actual consumer.

Pruned unused library constructors, adapter assemblies, graph/PCR finishers and
legacy end/library shortcuts from live helpers. `batch_ngs` retains supplied-state
drawing primitives; `protocol_specs` retains the metadata dataclass still declared
locally by SPLASH; `protocol_build` accepts the module/spec objects used by all its
callers and no longer offers an unused string-key lookup into an empty registry.
`tenx_graph` retains `partition_material`. Removed HiChIP's unused `sections()`
fallback, which the graph renderer never selected; its rebuilt schematic is
byte-for-byte identical. No source claims, molecular states or read layouts were
changed. Updated current maintainer pointers, the Tang constructor reference and
TODOFINAL's stale pending wording; historical audit records and all 57 queue rows remain.

Verification: all 407 protocol suites pass, along with shared 1281, source-tool 226
and catalogue 83 checks. All 407 diagram pages build; 402 documented NGS graph
contracts pass; all 813 published site files have resolving relative links. Pairing
and P5 lint have zero findings beyond the existing source-known PEM exemption.
The static dependency/executable-Markdown scan finds no remaining orphan shared
modules. In retained modules, 50 unused top-level definitions were removed.
All 57 resolved ownership rows and the empty migration pending list are preserved.
The cleanup is uncommitted.

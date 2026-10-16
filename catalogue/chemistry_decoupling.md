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

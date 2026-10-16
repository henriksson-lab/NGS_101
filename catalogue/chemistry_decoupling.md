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

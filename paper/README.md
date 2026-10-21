# Paper figures

Architecture figures for the protocol-definition and validation system. The final
artifacts are editable SVGs with vector shapes and live text, designed at 180 mm
width. There are no embedded raster images or external assets.

- [Figure 1 — protocol definitions](figures/figure_1_protocol_definition.svg):
  declaration fields, an illustrative provenance hypergraph, and consumers of the
  same declarations.
- [Figure 2 — validation architecture](figures/figure_2_validation_architecture.svg):
  operation labels, explicit molecular data, reusable checks, and per-claim results.
- [Figure 3 — sequence annotation and rendering](figures/figure_3_sequence_annotation.svg):
  segment metadata, identifier semantics, computed alignment and actual renderer output.
- [Figure 4 — constraint-checked folding](figures/figure_4_constraint_checked_folding.svg):
  topology constraints, geometric routing and actual synthetic padlock/hairpin views.
- [Draft captions](figure_captions.md).

## What the operation mapping actually means

The current implementation does **not** lower every named protocol operation into
a smaller executable instruction set. `StepKind` supplies operation labels and
structural contracts. A step can separately declare zero or more
`MolecularRelation` tags. The checker dispatches on those tags and the supplied
molecular object types; it does not infer a relation from an operation's name.
A step without a relation receives an unresolved transition claim even if another
check, such as primer binding, succeeds.

The snapshot used here has **41 operation labels**, **7 molecular-relation types**,
and **4 claim statuses**. The 41 include generic and explicitly unmodelled labels.
The seven relation types are not the whole validation framework: the chemistry
checker also inspects paired columns, PCR binding geometry, restriction cuts and
rendered/chemical RCA agreement directly. Constructors, graph contracts,
sequencing/read checks and source review provide additional, separate checks.
No count in these figures measures how many protocols or reactions are validated.

Suggested paper wording:

> Protocols are represented as typed provenance hypergraphs with source-qualified
> molecular states and role-labelled inputs. A reusable checker evaluates explicit
> molecular relations and typed-state geometry, reporting the status and scope of
> each claim independently of the operation label.

The arrows in Figure 2 indicate declaration and analysis flow. They are not an
operation-to-relation equivalence map. The Figure 1 graph is illustrative; it is
not a source-derived experimental protocol or a claim of complete molecular purity.
The reagent capsules group input categories for readability; the actual graph uses
separate materials and role-labelled ports.

## Rebuild

From the repository root:

```sh
python3 paper/tools/build_figures.py
```

The generator uses Python's standard library only. It reads the operation,
relation and status enums from the code. Relation explanations are deliberately
written for the paper; adding a new relation requires updating those explanations.
Paper labels use the Arial/Helvetica/sans-serif font stack; molecular text uses
a monospace stack. All text remains editable in a vector editor. Figures 3 and 4
embed output from the actual project renderers, with standalone styling, viewport
cropping and consistent paper legends. The embedding preserves molecular text
whitespace and per-segment hover titles. Rebuilding overwrites the generated figures, so edit the generator
for persistent changes, or save manual variants under a separate name.

The generator also exports the complete [operation vocabulary](figures/operation_vocabulary.tsv),
[relation vocabulary](figures/molecular_relations.tsv), and an
[implementation snapshot](figures/implementation_snapshot.json). The snapshot records
the repository commit and SHA-256 fingerprints of the relevant files, so the source
state remains identifiable even when the working tree differs from the commit.

## Synthetic molecular examples

[figure_examples.py](tools/figure_examples.py) declares the synthetic annotation,
padlock and hairpin examples. It uses the same `Construct`, `Scene`,
`PadlockCaptureTopology` and `HairpinTopology` APIs as protocol modules. Sequences,
complements, paired columns and folded nucleotide paths are generated from these
objects. The examples are explanatory software fixtures, not published oligos or
experimental recipes.

The folding figure describes a geometric layout of a declared pairing topology,
not thermodynamic structure prediction. In `topology.py`, `_choose_route` scores
candidate paths by length and bends for some complex assemblies. Padlock capture
and hairpin classes instead construct topology-specific routes directly. Both
are validated by `FoldLayout`. The scoring function is an aesthetic geometry cost,
not a free-energy model.

## Implementation references

| Figure element | Current implementation |
|---|---|
| Materials, role-labelled inputs and operation steps | `Material`, `Input`, `Step`, `ProtocolGraph.react` in [protocol_graph.py](../lib/protocol_graph.py) |
| Hypergraph construction and source-scope contracts | `add_step`, `_validate_step_shape`, `publication_problems` in [protocol_graph.py](../lib/protocol_graph.py) |
| Sequencing acquisition and explicit limitations | `SequencingPlan` in [protocol_graph.py](../lib/protocol_graph.py); `verify`, `read_contents`, `feature_spans` in [seqprimers.py](../lib/seqprimers.py) |
| Molecular states and structured identifiers | `Construct`, `Scene`, `MolecularFeature` in [chemdraw.py](../lib/chemdraw.py) |
| Typed PCR geometry | `PCRGeometry` in [pcr.py](../lib/pcr.py) |
| Seven relation tags | `MolecularRelation` in [protocol_graph.py](../lib/protocol_graph.py) |
| Relation dispatch and four direct check families | `_transition`, `_pairing`, `_pcr`, `_cut`, `check_graph` in [chemistry_validation.py](../lib/chemistry_validation.py) |
| Sequence annotations and SVG geometry | `Segment`, `MolecularFeature`, `Scene.anneal`, typed `*Visual` records and `panel_svg` in [chemdraw.py](../lib/chemdraw.py) |
| Fold constraints and route validation | `FoldStrand`, `FoldPairing`, `StraightHint`, `FoldLayout`, `_sample_polyline`, `_choose_route` in [topology.py](../lib/topology.py) |
| Padlock and hairpin layouts | `PadlockCaptureTopology`, `HairpinTopology` in [topology.py](../lib/topology.py) |
| Four statuses and scoped claim records | `Status`, `Claim` in [chemistry_validation.py](../lib/chemistry_validation.py) |
| Visible unresolved/contradictory warnings | `warnings` in [chemistry_validation.py](../lib/chemistry_validation.py) |
| Source fidelity, separate from molecular consistency | [source oligo guard](../tools/source_oligo_check.py), source-review records, and human reading; see [DESIGN.md](../DESIGN.md) |

These figures describe the architecture, not a coverage survey. A successful local
claim does not prove an entire reaction, source fidelity, experimental efficiency,
or complete chemistry validation. The checked graph contracts establish whether a
page has the required declarations; they do not establish that an original source
was correctly interpreted. Source inventories must be manually reviewed.

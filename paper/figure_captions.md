# Draft figure captions

## Figure 1. Source-qualified protocol definitions and their derived views

**A**, Each protocol owns its molecular materials, operation choices, products,
readout and evidence boundaries. Materials can carry exact or symbolic states,
structured identifier features, modifications and inference metadata. Steps declare
role-labelled inputs, outputs, geometry and optional molecular relations. Explicit
outcomes distinguish retained material, physical discard, QC and unresolved fates;
sequencing terminals carry a plan or an explicit limitation. Source inventories
record the reviewed scope rather than being inferred from the graph itself.
**B**, Illustrative representation of a typed provenance hypergraph, shown using
material nodes and operation nodes. Multiple reagent inputs feed a PCR step, and a
physical cleanup creates retained and removed fractions. Original polymers and
possible residual classes remain distinct within a mixture; amplification alone
does not imply their removal. Capsules grouping reagent categories are a visual
abstraction of separate role-labelled inputs. This example specifies no experimental
recipe. **C**, Shared readers use these declarations to generate schematics,
interpret sequencing geometry and audit graph structure and local molecular claims.

## Figure 2. Explicit molecular declarations connect a broad operation vocabulary to reusable local checks

**A**, The implementation provides 41 operation labels, including generic and
explicitly unmodelled labels; representative operations are shown. These labels
also define structural contracts for a graph step. **B**, Protocol authors supply
molecular input/output states, geometry and optional relation tags. Operation labels
do not automatically select a molecular relation. **C**, Seven explicit relation
types support narrow tests of circle closure, restriction-end fill, rolling-circle
copying, derived contact junctions, hairpin-stem pairing, strand-order preservation
and end combinations with declared downstream primer sites. Additional checks
inspect paired columns, PCR primer-binding geometry, restriction-cut coordinates
and agreement between rendered and chemical RCA states. **D**, Each question yields
an established, conditional, unresolved or contradictory claim with a target and
reason. Unknown or inferred bases can make a compatible result conditional;
missing relations or insufficient states remain unresolved. An explicitly
unmodelled step remains unresolved and produces a warning. These are local
consistency checks, not proofs of entire reactions or protocols. Source fidelity,
graph completeness and experimental performance are separate obligations.

## Interpretation notes for the Methods

The relation kernel should not be described as a complete biochemical simulator,
a compiler that lowers all operations to primitive reactions, or a single
protocol-level pass/fail system. For example, the PCR geometry check recomputes
inward-facing primer sites on the declared template; it does not by itself derive
or establish the declared PCR product. A hairpin relation tests the supplied stem's
pairing, not whether it is the preferred thermodynamic fold. A circle-order check
does not establish ligation-end chemistry or yield. End-combination eligibility
is relative to the specified downstream primers, not a universal claim that
other end classes can never amplify.

Construction-time invariants and graph publication contracts are distinct from
the independent chemistry audit. The chemistry audit reports contradictions and
unresolved scope; successful page construction does not establish that every step
has a checked molecular input-to-output relation. Human source review and source
transcription guards remain necessary to establish what the defining authority
actually reports.

## Figure 3. Sequence annotations drive semantic SVG rendering

**A**, Named `Segment` objects describe sequence content, placeholders, inference
and optional known lengths. Strand-end modifications belong to the strand declaration.
A separate `MolecularFeature` describes identifier meaning, including stable identity,
role, encoding and multipart grouping. Colour tags are presentation metadata and do
not establish biological meaning. **B**, `Construct` supplies ordered segments and
computed offsets. `Scene` receives strands in molecular 5′→3′ order, locates named
annealing spans and verifies overlapping paired columns. Typed row metadata describes
strands, binding spans, process arrows and comments; the shared renderer converts
that metadata into SVG geometry and selectable text. **C**, Actual renderer output
for a synthetic example containing a fixed handle, symbolic UMI, insert and inferred
tail. The complementary strand and the marked spans are computed by the model.
The example does not claim experimental provenance; its inference flag is included
to demonstrate rendering. **D**, Span positions and tooltip content are derived from
the same molecular objects. Melting temperature is included by the existing tooltip
helper for sufficiently long, concrete DNA segments; it is not a secondary-structure
prediction. Molecular base text and paired columns in this figure are generated by
`Scene` and `panel_svg`, not drawn or aligned independently for the paper.

## Figure 4. Constraint-checked geometric views of declared topologies

**A**, Topology objects supply sequence-resolved strands, paired intervals, optional
straightness hints and covalent bonds. Topology-specific rules place paired regions
and assign coordinates to flexible spacers. Some complex assemblies choose among
deterministic candidate routes using a geometric cost based on path length and bends.
`FoldLayout` checks declared pairing, coordinate bounds, collinearity hints and
unintended crossings or overlaps. **B**, Actual model-rendered synthetic padlock and
hairpin examples. `PadlockCaptureTopology` locates unique reverse-complementary arm
sites on a supplied target and routes the intervening probe backbone. `HairpinTopology`
uses two declared, reverse-complementary stem spans and routes the unpaired loop.
These two topology classes construct their routes directly; they do not use the
candidate-route search illustrated by more complex assemblies. Each nucleotide
retains a point in molecular 5′→3′ order, with pairing rungs and end labels derived
from the same objects. **C**, These checks establish consistency of a supplied
topology and its drawing. They do not minimize physical free energy, predict the
preferred fold, or establish binding occupancy, reaction efficiency or yield. Unknown
bases permitted by a declared constraint remain unknown; their compatibility is not
evidence of a known exact pair. The figures use synthetic sequences solely to explain
the software, with no experimental recipe or source-derived probe design.

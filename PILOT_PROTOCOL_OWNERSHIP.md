# Protocol ownership pilot

This pilot tests whether a named protocol can own its source-specific molecular claims
while sharing chemistry operations and page rendering. It covers exactly ten documented
entries in `catalogue/ours.tsv`. This file records the migration and its results before
any wider rollout.

| Protocol | Starting layout | Ownership after the first round |
|---|---|---|
| SMART-seq2 | Family graph and page in SMART-seq; named chemistry in `lib/smartseq2.py` | Local oligos, AC/ACAT evidence, drawn cDNA and final library constructors, and graph entry; shared family reactions remain in SMART-seq |
| PAIso-seq | Local graph; Smart-seq2 bases imported from `lib/` | Local PAIso-seq claims with explicit import of the Smart-seq2 subprotocol |
| scCAT-seq | Local graph; Smart-seq2 bases imported from `lib/` | Local graph explicitly establishes the linker its paper prints |
| scM&T-seq | Local graph; Smart-seq2 bases imported from `lib/`; local inference wrapper | Local graph imports the subprotocol constructor, whose inherited linker is inferred by default |
| scMT-seq | Same pattern as scM&T-seq | Same boundary; redundant inference wrapper removed |
| scNMT-seq | Same pattern as scM&T-seq | Same boundary; redundant inference wrapper removed |
| CEL-Seq2 | Model, graph and page delegated to CEL-Seq directory | Local primer design, randomhexRT sequence, final construct, read recipe and graph entry |
| sci-ATAC-seq (2018) | Model, graph and page delegated to original sci-ATAC-seq directory | Local DOI, PCR-index choice, final construct, read recipe and graph entry |
| SPLASH | Named definition and graph selected in `lib/protocol_specs.py`; oligos in `lib/splash_chem.py` | Local definition, oligos, psoralen choices and graph entry; shared builder takes those choices |
| sf4sU-Seq | Named definition and 4sU choices in `lib/nascent_protocols.py` | Local definition, 4sU readout, gel-selection choice and graph entry; shared builder takes those choices |

## Acceptance checks

- Each selected directory contains its source-specific definitions, chosen variant,
  evidence and uncertainty decisions, reaction graph, and final construct, or imports an
  explicitly named published subprotocol for a branch. `lib/` retains generic operations
  and renderers, with no selected protocol's named definition or oligo constants.
- `build_page.py` remains the executable entry point discovered by the site build.
- Affected protocol self-tests, source guards where local sources exist, documentation
  build, pairing lint and site build pass. Compare rendered pages against the pre-pilot
  revision and explain intentional changes.
- Record exceptions and lessons here before deciding whether to extend the pattern.

## Findings and limits

The Smart-seq2 group exposed an evidence problem beyond file placement. The ACAT linker
is printed by scCAT-seq, while scM&T-seq, scMT-seq and scNMT-seq inherit Smart-seq2
without printing their selected TSO. Their former copies of an inference wrapper are
replaced by a constructor that marks the inherited linker inferred by default; scCAT-seq
explicitly opts out. PAIso-seq prints its own TSO and uses the published subprotocol
sequence. The Smart-seq2 directory now also constructs its source-marked diagram cDNA
and final Nextera library. This removes the second independently assembled final library
from the SMART-seq family module.

SPLASH now owns its defining citation, source claims, oligos, selection reagents and
graph invocation. Its shared builder accepts those inputs instead of selecting SPLASH by
name. sf4sU-Seq now owns its defining citation, 4sU readout and gel-selection route;
the shared 4sU builder accepts that route. Neither selected definition remains in `lib/`.
All ten pilot pages are byte-identical to pages built from the pre-pilot revision. The
original SMART-seq, CEL-Seq and sci-ATAC-seq pages and their sibling variant pages also
render identically after the family-builder changes.

The first round exposed a remaining boundary problem: SMART-seq2, CEL-Seq2 and
sci-ATAC-seq (2018) had local graph entries whose reaction sequence still lived in an
older family member's `tools/` directory. SPLASH and sf4sU-Seq similarly selected a
named route through a graph assembler in `lib/`. The second round moved the selected
reaction sequences into those five protocol modules. All ten now assemble their graph
locally. The other five already did so in the first round.

The local graphs still reuse bounded chemistry operations and molecular scenes:
CEL-Seq2 uses `barcoded_ivt_front`, SPLASH uses `finish_rt_pcr`,
sf4sU-Seq uses `ligate_rna_adapter` and `finish_cdna`, and the two versioned methods import
published predecessor molecular drawings. The parent family graphs remain for their
own pages and other siblings. Cross-directory imports therefore remain for those
published predecessor states and renderers; moving a parent directory would require
updating the child imports. The reaction order and source decisions for the selected
ten are visible in their own modules.

A third pass moved the inherited SMART-handle scenes into the predecessor's
`smartseq_states.py` and the two-round sci-ATAC molecular states and pooling operation
into its `sciatac.py`. The selected graph modules now import those model interfaces
directly; page renderers import the same functions. The generic RNA library finishing
operations used by SPLASH and sf4sU-Seq now live in `lib/rna_library_operations.py`,
which contains no named protocol definitions. The source guard gained `--details` to
list the material IDs behind its coverage count. This makes gaps actionable without
asserting that an unrepresented oligo has a published exact sequence.

The source-text guard discovered all ten local graph entries. It checked five exact
oligo inputs and reported 53 inputs without exact states; those unrepresented inputs
remain a chemistry-completeness limitation of the existing models, not evidence that
the source texts establish them. Pairing lint found zero problematic duplexes across the
ten rendered pages.

## Verification

The first-round checks were:

- Built the pre-pilot revision in an isolated `/tmp` archive and compared the ten pilot
  schematic HTML files byte for byte: all ten are identical. The three family parents
  and their other checked sibling pages also match their baselines.
- Ran the ten source-oligo guards: ten graphs discovered, five exact inputs checked,
  zero failures, and 53 inputs without exact states.
- Ran `tools/lint_pairing.py` on the ten pages: zero problematic duplexes.
- Ran `build_docs.py --check`: all 765 generated note pages are current.
- Ran `build_site.py`: all 405 diagram pages built, all 400 documented NGS schematics
  passed the typed graph completeness contract, and every assembled site link resolved.

After moving the five delegated graphs into their protocol directories:

- The ten selected protocol suites, the three predecessor suites and 1,134 shared
  checks passed.
- The ten source-oligo guards found ten graphs, checked five exact inputs, and
  reported zero failures and 53 inputs without exact states.
- The five rebuilt schematics were byte-identical to the pre-pilot revision; pairing
  lint reported zero problematic duplexes.
- `build_site.py` built 405 diagram pages, checked the completeness contract for all
  400 documented NGS schematics, and resolved every site link.

After separating reusable molecular states and bounded RNA operations:

- The ten protocol suites, 1,134 shared checks and 80 source-tool checks passed.
- `tools/selftest_source_oligos.py --details` listed the unrepresented material IDs;
  the aggregate remained five exact inputs checked, zero failures and 53 inputs without
  exact states across ten graphs.
- All ten pilot schematic HTML files remained byte-identical to the pre-pilot build;
  pairing lint reported zero problematic duplexes and `build_docs.py --check` found all
  765 generated note pages current.
- `build_site.py` again built all 405 diagram pages, passed the 400 documented NGS
  graph contracts and resolved every site link.

## Chemistry-checking follow-up

`python3 tools/check_chemistry.py --details` now analyzes the ten local graph
declarations and five boundary protocols (smMIP, FISSEQ, Small-seq, MALBAC and in situ
Hi-C). These add hybridization, circle closure and rolling-circle copying, RNA ligation
and dephosphorylation, symbolic folding, and restriction fill/contact junctions to the
78-step pilot sample. The checker consumes the existing Python objects; FISSEQ and
MALBAC expose optional `chemical_state` objects where the diagram state cannot retain
the circle or fold topology. In situ Hi-C passes its existing restriction digest and
fill objects through the graph.

The first audit of 15 protocols reports 25 established, 29 conditional, 114 unresolved,
and zero contradictory *scoped claims*. It can recheck exact local pairing and primer
binding, compare a circle to its linear precursor, derive FISSEQ repeat sequence from
that circle, check the in situ Hi-C filled ends and junction, and classify MALBAC's
unknown handle pairing as conditional. It cannot yet derive most full products from
reactants. Small-seq's two RNA ligations and smMIP's target-gap closure remain
unresolved because their graph declarations do not contain the necessary exact
reactant/product chemistry. No status here verifies the source transcription or a
reaction's yield. A publishable graph and an established local claim are different
contracts.

## 2026-10-15 — Primary-source recovery after the first audit

The initial unresolved count was a worklist, not a count of missing publications.
Source review of the same 15 protocols recovered exact vendor PCR primers for
CEL-Seq2 (Illumina TruSeq Small RNA) and SPLASH (NEB E7335L), adapter and primer
oligos plus decapping/ligation order for sf4sU-Seq (NEB E7300), and the Small-seq
SR001 linker, RA3/RTP duplex, and smMIP probe phosphorylation and 112-nt captured
gap. The source transcription guard exposed an extra C in the modeled Small-seq
SR001 linker; the corrected 22-nt-miRNA library is 152 bp. SPLASH's detailed JoVE
author protocol adds RNA end repair, RNA-linker ligation, RT, and cDNA circle
closure to its graph. scCAT-seq, scM&T-seq and scNMT-seq now declare source-known
primer inputs. The oligo audit ledger records source scope and remaining gaps.

The current chemistry audit reports 27 established, 32 conditional, 119 unresolved,
and zero contradictory scoped claims. Its larger unresolved total includes newly
declared reactions, so it does not measure regression. Across these 15 directories,
the local source guard found 14 loadable graphs and checked 18 exact inputs with
zero failures; 52 inputs have no exact sequence state. The residual work separates
source limits (for example, the MALBAC Science supplement and smMIP backbone
table) from absent intermediate molecular states and checker transition rules.
The inventor MALBAC patent prints a candidate 27-nt handle but describes a
different primer-foot architecture, so it does not establish the Science paper's
handle bases. In situ Hi-C's original procedure still names an Illumina adapter
and primers without fixing their exact sequence or member choice.

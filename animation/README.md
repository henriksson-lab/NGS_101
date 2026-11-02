# Native chemistry animations

Presentation movies rendered by Python/Cairo directly into FFmpeg. No Manim,
browser, screenshots or intermediate frame directory is used for MP4 generation.
The [implementation plan](../plans/chemistry-animation.md) records the accepted scope.

## Install and run

Use Python3.10+, pycairo and the native Cairo library. MP4 export requires FFmpeg
with libx264; suite verification also requires ffprobe. Use your platform's native
Cairo/FFmpeg packages and a dedicated Python environment:

```sh
python3 -m pip install -r animation/requirements.txt
python3 -m animation list
python3 -m animation smartseq2 --output /tmp/smartseq2.mp4
python3 -m animation crispr-mip --shot gap --output /tmp/gap.mp4
python3 -m animation boloramis --output /tmp/rca-draft.mp4 --size 960x540
python3 -m animation smartseq2 --list-shots
```

Commands run from repository root. The default is opaque-background16:9 H.264 MP4,
1920×1080,30fps,yuv420p and fast-start metadata. Playback defaults to **1.5×** the
authored timeline; use `--speed 1` for the original pace or e.g. `--speed 2` for twice the pace. Keep generated artifacts outside the
checkout. Existing CLI outputs and MP4 manifests are refused rather than overwritten.
Rendering never downloads dependencies or source documents. Cairo/FFmpeg were already
installed in the development environment; dependency ranges are not a claim that every
platform/version combination has been verified.

The normal chemistry library and offline website build do not import animation
rendering dependencies. Animation dependencies are an explicit AGENTS.md exception.

## Examples

| CLI name | Selected mechanism |
|---|---|
| `smartseq2` | Substrate-driven RT, template switching, PCR, two Tn5 insertions, original-end repair, index PCR and complete library anatomy; 2014 TSO and source-era N501/N701 example indexes |
| `smartseq` | Original2012 SMART-seq: poly(A)/RT-primer binding, retained RNA through template switching, first complementary DNA, PCR sites, Tn5 dimers/cuts, index-primer tails, final functional library anatomy (cleanup omitted) |
| `crispr-mip` | Probe folding, gap copying, phosphorylated closure and end-dependent exonuclease selection |
| `boloramis` | Source probe/primer co-hybridization, folding, RNA-splinted nick sealing, substrate-driven RCA and conserved-product overview |
| `atac` | Loaded Tn5 homodimers, staggered strand transfer and both original-end repairs; conditional ME′ release |
| `smallseq` | Preadenylated3′ adapter ligation, free-adapter processing/degradation and5′ adapter ligation |
| `gibson` | Source method with the repository's invented DNA example: resection, complementary overlap annealing, nick sealing |
| `hic` | MboI cuts, both biotin-bearing end fills and both strands of a selected proximity junction |
| `dropseq` | Poly(A) capture, full bead-primer barcode/UMI anatomy and tethered RT |
| `scrrbs` | Prepared genomic window: C→physicalU, U→A complementary copying, then A→T copying |
| `celseq` | T7 transcription and successive initiations on the retained DNA template |
| `golden-gate` | BsaI Type IIS cuts, cohesive annealing, all four nick seals and site-free synthetic circle; all flanking pieces retained |
| `reference` | Synthetic exact primer/template-copy fixture, not a published specimen |
| `rnaseh` | H1-class cleavage of RNA in a supplied RNA–DNA hybrid; RNA fragments and intact DNA retained |

Whole movies have multiple narration-paced scenes. Export one shot for a short
presentation clip; `--list-shots` shows the narration and estimated durations. The
batch command also writes one focused clip per example.

**SMART-seq2, BOLORAMIS and the ten additional examples** use nucleotide identities
and reaction-generated bonds. The original2012 `smartseq` adapter remains a historical
prototype. The CLI registry selects the current implementations; older layouts in
`protocols.py` remain regression fixtures and are not the current ten-example movies.
The [ten-example plan and audit](../plans/chemistry-animation-ten-more.md) records
sources, explicit scopes, review iterations and exported artifacts.

The step counter and progress bar are hidden by default; `--timeline` enables them.
An enzyme enters when its reaction starts. `persistent=True` keeps it through the
result hold; `enzyme_prebound=True` separately declares an enzyme already present
from a preceding reaction, such as RT during template switching.

## Narration and presentation settings

Default scenes show short titles, sparse molecule/enzyme labels and necessary scope
qualifiers. Sources, complete sequences, narrator text and scientific limits go in the
MP4's neighboring `.manifest.json`, not a wall of on-screen prose. No voice synthesis,
audio recording or audio muxing is required by this implementation.

```sh
python3 -m animation smartseq2 --shot rt --narration-seconds 18 --output /tmp/rt.mp4
python3 -m animation smartseq2 --shot rt --duration 24 --output /tmp/rt-slow.mp4
python3 -m animation smartseq2 --output /tmp/rehearsal.mp4 --speed 1.5
python3 -m animation smartseq2 --output /tmp/original-pace.mp4 --speed 1
python3 -m animation reference --output /tmp/detail.mp4 --bases --directions
python3 -m animation reference --output /tmp/detail.svg --time 4 --coordinates --reagents
```

Each scene establishes its substrate, performs an action and holds the result.
The default125words/minute estimates how long narration will take. Explicit
`--narration-seconds` uses a known recording/speaking duration. Long narration extends
the final hold instead of slowing each base addition indefinitely. `--duration`
rejects truncation of narration or minimum establishing/action/hold intervals.
`--speed` (default1.5) scales the entire timeline, including its narration allowance.
Timing-file values and `--duration` / `--narration-seconds` use the authored timeline;
`--list-shots`, still timestamps, MP4s and exported manifests use the selected playback rate.
Use `--speed 1` with a fixed-length narration recording unless it is retimed too. Word counts are estimates, not audio alignment.

Nucleotide letters and strand arrows are enabled by default; `--bases` and
`--directions` remain accepted. Unknown bases remain literal placeholder classes,
without becoming established sequences. SMART-seq2's final molecule wraps onto two
continued rows to preserve legibility. `--coordinates`, `--reagents` and `--evidence`
add optional information. `--font` chooses a locally installed label font; nucleotide
letters use DejaVu Sans Mono. Font availability depends on the local Cairo setup.

A timing file applies scene-specific edits without altering molecular events:

```json
{
  "rt": {"narration_seconds": 18, "hold": 4},
  "switch": {"action": 7, "narration": "Your own spoken explanation."}
}
```

Pass `--timings /path/to/timing.json`. Supported keys are `duration`,
`narration_seconds`, `establish`, `action`, `hold`, `narration` and `cue`; unknown scene
ids/keys are errors. Molecular framing is controlled separately by each shot's
`camera` center x/y/zoom and its strand paths/labels, not by changing chemistry. The palette is defined separately
in `animation/render.py` as `PALETTE`; an authoring script can customize those colors
without changing protocol events or source declarations.
The native composition uses a1600×900 logical canvas with a protected title/footer.

## Stills, interchange and reproducible suite export

```sh
python3 -m animation smartseq2 --output /tmp/smartseq2.json
python3 -m animation json --input /tmp/smartseq2.json --output /tmp/replay.mp4
python3 -m animation crispr-mip --shot closure --time 5 --output /tmp/closure.svg
python3 -m animation crispr-mip --shot closure --time 5 --output /tmp/closure.png
python3 -m animation.batch --directory /tmp/chemistry-movies
python3 -m animation.batch --directory /tmp/chemistry-movies --resume
python3 -m animation.selftest --native
```

SVG and PNG use exactly the same evaluated scene as MP4 at the requested timestamp.
The versioned `ngs101.animation/4` JSON contains explicit molecular/timeline data,
including nucleotide identities and reaction deltas for generated views. Unsupported
versions are rejected. Versions1–3 are historical; regenerate those scenes from the
adapter. Backbone bonds, duplex contacts and nucleotide births are separate from
coordinates, labels and enzyme glyphs.

The batch command exports movies, sourceJSON, SVG/PNG action stills, per-shot contact
sheets, actual decoded-video stills, focused clips and `verification.json`. It checks
codec, pixel format, dimensions, rate and frame count with ffprobe. Focus clips are
frame-accurate re-encodes of their native-rendered parent, with source intervals and
narration sidecars. `--examples` selects a subset; `--size`/`--fps` control drafts.
`--resume` reuses only movies whose serialized scene/settings match the manifest.
This is artifact generation, not automatic website publication.

## Add a protocol adapter

`protocols.py` contains legacy adapters and shared loading/provenance support;
`smartseq2.py` is the separate current SMART-seq2 adapter importing its protocol owner.
Each declares the selected steps, actual sequence/end values, unresolved regions,
reaction events, narration and directing choices. `qualify()` checks that each shot
maps to actual step ids and snapshots the protocol module hash and source authority.
A step's type alone never generates detailed chemistry automatically.

- `Strand`: physical5′→3′ sequence, sugar class/intervals, end modifications, ordered
  path, initial/preserved intervals, nicks, marks, displaced pieces and retained-molecule lineage.
- `Region`: model-coordinate interval, functional label/color and callout placement.
  Regions remain colored while copying; unknown native values still stay schematic.
- `Pairing`: explicit antiparallel coordinates; exact supplied complementarity is checked.
- `Event`: copying coordinates, selected nick chemistry, cleavage, folding, resection,
  conversion, transfer, span highlights or physical motion. `displace` moves retained
  intervals without deleting cut pieces. `complex="homodimer"` draws linked identical
  subunits; `persistent=True` retains an enzyme through a multi-shot copying step.
  A fold never implicitly seals a circle.
- `Junction`: explicit3′OH and5′P/App donor geometry for a selected successful join.
- `Shot`: source scope, narration/timing, molecule inventory, geometry and camera.
  Generated shots also declare the actual starting atoms/bonds and event deltas.
- `Movie`: ordered shots, original protocol identity and provenance/limits.

An extension event declares which template interval supplies new bases and preserves
the original primer identity. `exact=True` requires matching explicit bases; unknown
schematic classes cannot establish exact copying. RCA additionally requires a closed
template and checks the wrapped copy. The displaced product3′end follows the enzyme
instead of growing disconnected from it. Retained lineage intervals must match their
source molecules; joining a nick does not delete or duplicate nucleotide intervals.

## Scientific limits

These are explanatory mechanism animations, not molecular dynamics, protein
structures, measured trajectories, thermodynamic folding predictions, reaction-rate
models or estimates of successful yield. Native kit geometry remains visibly
source-qualified. Unknown molecule lengths in coarse tracks are graphic subdivisions,
not inferred measurements. Reagents and off-products outside a selected view are not
asserted to have vanished; sidecars describe omitted populations/workups.

In particular, original SMART-seq native kit/first-cycle geometry is not replaced with
SMART-seq2. BOLORAMIS's displayed RNA and Small-seq's masking target are explicitly
conditional complementarity examples. Gibson's exact toy is zero-gap: the movie does
not invent polymerase synthesis there. scRRBS's sequence/methylation occupancy is
illustrative and distinguishes physical DNAU from a subsequent copied/reported base.
Legacy ATAC transfer remains coarse. SMART-seq2 uses explicit staggered strand
transfer and repair, with no invented ligase step. Complete
NGS read landing and all reaction species are not automatically validated by animation.

The simulator and scene validator check their stated molecular relations; source
reading remains authoritative. Do not narrow a protocol's evidence requirements to
obtain a more detailed movie. Legacy adapters without generated reaction deltas still
have weaker endpoint/geometry checks and must not be described as fully simulated.


## SMART-seq2 presentation revision (2026-10-31)

`smartseq2_simulation.py` supplies molecular inputs and explicit enzyme capabilities
to `simulator.py`; polymerase products are generated by walking actual templates.
It never supplies a desired product to an extension operation. `smartseq2.py` projects
those nucleotide identities onto twelve presentation scenes. `molecular_frame.py`
replays actual bond changes, independently of geometry. Labels follow oligos;
bound strands remain tightly spaced and incoming primers enter from outside the view.

PCR uses `Simulator.extension_sites()` and `extend_available()` to discover all
eligible3′OH ends in the supplied molecular complexes. Template and oligo ends
receive the same test: pairing, accessible next template base, sugar capabilities,
protein occupancy and strand displacement policy. The projector creates an actor
for every generated extension and grows the existing template tracks too. In the
first two index-PCR cycles this gives four polymerases; after the tails are complete,
only the new primer ends can extend. No fixed actor count determines the chemistry.
Polymerase and ligase glyphs use half the original width and height.

The movie shows oligo(dT)/poly(A) binding and the unpaired ISPCR handle, RT with RNA
retained through template switching, actual first complementary-strand synthesis,
PCR binding, two Tn5 homodimers and staggered strand transfer, departing flanks,
original-end repair, and successive indexed PCR copies. The final slide shows every
base and all nine component groups on both strands. SPRI cleanup is omitted.

The TSO is explicitly the2014 detailed-protocol form. N501/N701 are source-era
available-kit examples, not recovered study-well assignments. A synthetic short
transcript, selected cut sites and terminal-C extent make the illustration reproducible.
The repair path pauses at bound ME complements, releases them as separate molecules,
and continues synthesis; the source does not assign a unique release mechanism.
No rate, yield or whole-population prediction follows. See the
[sequence source review](smartseq-sequence-research.md) and
[Tn5 mechanism review](tn5-mechanism-research.md).

The other protocol recipes/views are at different migration stages. Their existence
or a green endpoint check is not evidence that their public adapters already use
substrate-generated molecular histories. Original2012 `smartseq` remains a separate
historical prototype; it is outside the current SMART-seq2 presentation goal.

## BOLORAMIS padlock presentation

```sh
python3 -m animation boloramis --output /tmp/boloramis.mp4
python3 -m animation boloramis --shot rca --output /tmp/rca.mp4
python3 -m animation boloramis --shot product --output /tmp/product.svg
```

Five scenes show the source ACTN1A1 probe and RCA primer co-hybridizing, folding of
the same probe, SplintR sealing an actual nick, Phi29 extending around the closed
template, and a conserved-product overview. Source RNA is a conditional matching
window. Both primer and RNA persist through ligation; no extra primer or RNase is
inserted. The nucleotide pool enters at the source's post-ligation RCA stage.

The shared layout follows actual duplex contacts. RCA retains its full bound suffix,
while displaced DNA accumulates in separated coils and unpaired RNA moves aside.
The final overview lays out the same nucleotide IDs and bonds, including the paired
growing end. Geometry never determines whether the nick is sealed or RNA is displaced.

The displayed trajectory selects successful RNA-hybrid passage. Disabling that
explicit Phi29 capability stalls copying at the occupied hybrid. Applicability of
the mechanistic study to the source NxGen37°C reaction remains inferred. New DNA
T-equivalents retain unresolved dT/aminoallyl-dU occupancy annotations; the movie
shows backbone lettering without inventing modification or crosslink sites. The
three-traversal budget illustrates one product rather than predicting yield. See
[the source review](boloramis-rna-followup.md).

BOLORAMIS reframing uses `Shot.transition` during the establish hold. Coordinates
are matched by nucleotide identity, rigid circles rotate without collapsing, and
the camera eases between views; the renderer
rejects transitions whose molecular frames differ. Chemistry begins after reframing.
The source probe's closure is marked only when its actual backbone bond exists.
SplintR is labelled as a ligase and leaves after sealing; the retained final enzyme
is explicitly labelled Phi29 polymerase. `Strand.sealed_marks` highlights the joined
bases and bond in later views without implying a new reaction.


## Auditing exported movies

```sh
python3 -m animation.batch --directory /tmp/movies --examples golden-gate rnaseh
python3 -m animation.audit --directory /tmp/movies --examples golden-gate rnaseh
python3 -m animation.selftest --native
```

The audit compares a canonicalized scene manifest with current definitions, verifies
codec/settings/frame counts and decodes actual MP4 reaction/transition frames into
review strips. It leaves the visual verdict **pending** until a reader inspects them.
SVG final frames are separate presentation assets, not proof of the MP4 content.
Released-nucleotide display is opt-in per shot (`product_pool=True`) and derives its
count and letters from actual released atoms. Intact reagent arrival uses a generated
`introduce` event, distinct from nucleotide incorporation. Unborn molecules have no
captions or physical end markers. Canonical serialization sorts unordered released
pools, allowing reproducible manifests across processes.

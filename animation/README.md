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
python3 -m animation smartseq --output /tmp/smartseq.mp4
python3 -m animation crispr-mip --shot gap --output /tmp/gap.mp4
python3 -m animation boloramis --output /tmp/rca-draft.mp4 --size 960x540
python3 -m animation smartseq --list-shots
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
| `smartseq` | Original2012 SMART-seq: poly(A)/RT-primer binding, retained RNA through template switching, first complementary DNA, PCR sites, Tn5 dimers/cuts, index-primer tails, final functional library anatomy (cleanup omitted) |
| `crispr-mip` | Continuous padlock folding, target gap copying, phosphorylated nick closure |
| `boloramis` | Conditional RNA-splinted closure and derived three-repeat RCA with displaced product |
| `atac` | Coarse loaded-Tn5 adapter transfer and extension boundary |
| `smallseq` | Conditional rRNA masking,3′ adapter ligation, free-adapter degradation and5′ adapter ligation |
| `gibson` | Source method with the repository's invented DNA example: resection, complementary overlap annealing, nick sealing |
| `hic` | MboI staggered cuts, biotin fill-in, selected proximity junction and physical enrichment |
| `dropseq` | Droplet capture, full bead-primer barcode/UMI anatomy and tethered RT |
| `scrrbs` | MspI-end processing/A-tailing; illustrative C→chemicalDNAU and later complementary copying |
| `celseq` | T7 transcription and successive initiations on the retained DNA template |
| `reference` | Synthetic exact primer/template-copy fixture, not a published specimen |

Whole movies have multiple narration-paced scenes. Export one shot for a short
presentation clip; `--list-shots` shows the narration and estimated durations. The
batch command also writes one focused clip per example.

## Narration and presentation settings

Default scenes show short titles, sparse molecule/enzyme labels and necessary scope
qualifiers. Sources, complete sequences, narrator text and scientific limits go in the
MP4's neighboring `.manifest.json`, not a wall of on-screen prose. No voice synthesis,
audio recording or audio muxing is required by this implementation.

```sh
python3 -m animation smartseq --shot rt --narration-seconds 18 --output /tmp/rt.mp4
python3 -m animation smartseq --shot rt --duration 24 --output /tmp/rt-slow.mp4
python3 -m animation smartseq --output /tmp/rehearsal.mp4 --speed 1.5
python3 -m animation smartseq --output /tmp/original-pace.mp4 --speed 1
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

`--bases`, `--directions`, `--coordinates`, `--reagents` and `--evidence` add optional
information. Letters stay hidden on source-unresolved schematic regions. Long known
molecules use a backbone view instead of unreadably tiny base glyphs. `--font` chooses
a locally installed font; the fallback depends on Cairo/system font availability.

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
python3 -m animation smartseq --output /tmp/smartseq.json
python3 -m animation json --input /tmp/smartseq.json --output /tmp/replay.mp4
python3 -m animation crispr-mip --shot closure --time 5 --output /tmp/closure.svg
python3 -m animation crispr-mip --shot closure --time 5 --output /tmp/closure.png
python3 -m animation.batch --directory /tmp/chemistry-movies
python3 -m animation.batch --directory /tmp/chemistry-movies --resume
python3 -m animation.selftest --native
```

SVG and PNG use exactly the same evaluated scene as MP4 at the requested timestamp.
The versioned`ngs101.animation/3` JSON contains explicit molecular/timeline data, not
Python object dumps or SVG screenshots. Unsupported versions are rejected. Versions1/2 are historical and are not silently treated as version3. Version3 adds
functional regions, timed highlights, preserved fragment displacements and explicit
homodimer glyphs; regenerate older scene JSON from the adapter.

The batch command exports movies, sourceJSON, SVG/PNG action stills, per-shot contact
sheets, actual decoded-video stills, focused clips and `verification.json`. It checks
codec, pixel format, dimensions, rate and frame count with ffprobe. Focus clips are
frame-accurate re-encodes of their native-rendered parent, with source intervals and
narration sidecars. `--examples` selects a subset; `--size`/`--fps` control drafts.
`--resume` reuses only movies whose serialized scene/settings match the manifest.
This is artifact generation, not automatic website publication.

## Add a protocol adapter

`protocols.py` contains explicit adapters importing existing protocol constructors.
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
Tn5 transfer is coarse and not misrepresented as ligase-mediated insertion. Complete
NGS read landing and all reaction species are not automatically validated by animation.

The renderer checks the declared local relations; source reading remains authoritative.
Do not narrow a protocol's evidence requirements to obtain a more detailed movie.


## Presentation revision (2026-10-31)

SMART-seq uses the published CDS oligo as an **oligo(dT) RT primer** on screen. Its
known 5′ handle and poly(dT) tract are separate colored regions; the displayed GA
anchor is one illustrative VN member. RNA coordinates remain fixed through RT,
terminal addition and template switching. The complementary DNA starts as a short
primer and grows before a later duplex PCR cycle is shown. Both index-PCR tails
remain distinct from their 3′ binding regions. Tn5 is a linked pair of identical
subunits carrying adapters; two insertion events bracket a selected fragment,
with flanking pieces moving apart after cleavage. SPRI is omitted from this movie.

The final slide labels **functional** library anatomy. Original SMART-seq's kit TSO,
IS-PCR footprint and Epicentre index-adapter ordering are not fully curated in the
source model. The movie therefore does not assert modern Nextera/SMART-seq2 oligos,
native index widths, or a sequence-validated final construct. The end-subregion order
and drawn counts on that slide remain schematic. This limitation is also in the
movie manifest and the slide's scope line.

Other adapters now label capture/ligation arms, RNA-adapter end chemistry, UMI and
barcode inheritance, Gibson overlaps, Hi-C recognition/contact regions and the T7
promoter. Both Hi-C restriction products and Tn5 cleavage products remain present
when separated. These are selected mechanism illustrations, not whole-protocol
population/kinetics simulations or newly verified source chemistry.

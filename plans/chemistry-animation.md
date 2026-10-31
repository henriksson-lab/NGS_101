# Chemistry animation tool

Status: SMART-seq2 and BOLORAMIS selected presentation examples complete (2026-10-31); other protocol
simulator migrations remain in progress.
Prototype completion/export records below are historical; they do not establish
that the revised simulator or all protocol adapters are complete.

## Accepted direction

Produce attractive, concise MP4 chemistry animations for live presentations. Use
Python for chemistry, timeline evaluation and initial rendering. Use a native 2D
backend (initial candidate: Cairo) and stream frames directly to FFmpeg. No Manim,
headless browser, screenshots or browser-based frame export. Optional SVG/PNG stills
should use the same evaluated scene. Rust is deferred until a measured bottleneck
justifies it; a versioned animation interchange format preserves that option.

Animation dependencies and external tools are explicitly allowed by AGENTS.md.
Keep them isolated from ordinary chemistry tests and the offline website build.
Document installation, check dependencies before export, and keep binaries, rendered
videos and frame buffers outside version control. No automatic site publication,
commit or push is implied by this goal.

## Presentation and timing requirements

- Default to clean 16:9 composition, opaque background, readable typography,
  consistent strand/enzyme colors, smooth motion and sparse labels. Make enzymes
  visibly bind, act and release; do not simply crossfade before/after diagrams.
- Avoid paragraphs, permanent legends, complete reagent lists or audit tables on
  screen. Put sources, evidence and narration in a sidecar manifest. Necessary
  uncertainty remains visible through a short qualifier or visual convention.
- Provide optional flags for base letters, coordinates, direction labels, reagent
  details and evidence annotations; these must not be required for a useful movie.
- Give each shot a narration cue and editable duration. Compose an establishing
  hold, reaction interval and concluding hold. The viewer needs time to recognize
  the substrate, hear the explanation and inspect the product.
- A configured speaking rate can suggest a minimum shot duration from narration
  word count; explicit narration duration overrides that estimate. Extend holds
  when speech is long rather than slowing every nucleotide addition indefinitely.
  Do not claim word-count estimates align recorded speech automatically.
- Support an overall speed multiplier and per-shot duration overrides. Motion is
  explanatory, not calibrated kinetics or evidence of reaction probability/yield.
- Default export: H.264 MP4, 1920x1080, 30fps, yuv420p, fast-start metadata. Provide
  quick low-resolution drafts, selected-shot export, deterministic timestamp stills
  and configurable resolution/fps. No voice synthesis or audio recording required.

## Architecture

1. Protocol adapters import existing Python definitions/constructors. Do not fork
   chemistry into the renderer or substitute SMART-seq2 for original SMART-seq.
2. An animation model describes strand identities, sequence intervals, pairing,
   covalent connectivity, modified ends and selected reaction correspondence.
   Retained bases preserve identity; newly synthesized bases identify their source
   template interval. Use interval mappings instead of heavyweight objects per base.
3. Explicit events describe annealing, separation, extension, tailing, nick sealing,
   cleavage, folding/circle closure, strand displacement, binding/release and physical
   selection. A StepKind label alone cannot authorize a detailed molecular event.
4. A versioned JSON document holds molecular data, checked event arguments, source
   pointers, selected/inferred/unresolved scope and shot/narration information.
   Serialize explicit data, not arbitrary Python objects, HTML or SVG screenshots.
5. The director controls layout, camera/framing, emphasis and timing separately
   from chemistry. Evaluate frame(t) directly so seeking/export order cannot change
   the result. No accumulated integration or random movement without a fixed seed.
6. A native drawing backend draws paths, bases, bonds, enzymes and labels. A molecule
   may bend while preserving sequence order and connectivity. Base letters use level
   of detail; long molecules must not force unreadably tiny text or giant canvases.
7. An exporter streams raw frames to FFmpeg, checks encoder failure and writes a
   provenance/timing/options sidecar. Stills and video share the evaluator.

Reuse Scene/padlock geometry and validation where applicable. New correspondence
checks should benefit chemistry validation where appropriate, while drawing/layout
logic stays in the animation package. No broad migration of existing protocols.

## Scientific boundaries

Illustrative enzyme silhouettes and paths are not structures, measured trajectories,
kinetics or thermodynamic fold predictions. Preserve 5'/3' polarity and distinguish
RNA, DNA and chemical modifications. Ligase seals an explicitly compatible selected
junction; polymerase extends a supplied extendable end 5' to 3' using the declared
substrate. Never silently repair missing chemistry. Keep catalysts and original
strands present unless an actual physical separation or declared fate removes them.
An animation shows a selected illustrative event, not all molecules succeeding.

Original SMART-seq's native kit bases/first-copy geometry are partly unresolved in
current definitions. Use source-qualified coarse regions for those stages and do not
present invented bases as published. Exact base-level reference examples must be
labeled illustrative and connected to their actual protocol scope. Keep predecessor,
commercial-kit, protocol-version and selected-branch boundaries intact.

## Agreed implementation order and coverage (user accepted items 1–6)

1. **SMART-seq**: principal example. Primer annealing, RT, untemplated addition,
   template switching, retained RNA, DNA amplification/strand separation and selected
   native tagmentation/index-PCR boundaries. Check original versus SMART-seq2 scope
   explicitly; user did not authorize silently changing the named protocol.
2. **CRISPR-MIP**: target/probe arms, folding, polymerase gap fill and ligase closure.
3. **BOLORAMIS**: RNA-splinted DNA ligation, circle, rolling-circle copying and strand
   displacement. Sequencing-by-ligation cycles are an extension once core scenes work.
4. **ATAC-seq**: loaded Tn5, target engagement, staggered adapter transfer, release and
   source-qualified gap filling. Do not depict ligase-mediated insertion.
5. **Small-seq and Gibson assembly**: sequential RNA-adapter ligation, blocked ends,
   unligated-adapter removal and masking-oligo binding; exonuclease
   recession, overlap annealing, filling and nick sealing. Gibson's repository example
   is illustrative constructed DNA, not a fabricated published specimen.
6. **Additional coverage scenes**: in situ Hi-C (restriction/fill-in/proximity ligation
   and affinity selection), Drop-seq (tethered bead oligos/droplets and sample-to-molecule
   scale), scRRBS (end preparation and chemical versus copied-base conversion), and
   CEL-Seq (T7 transcription/linear amplification). These may be short focused clips,
   not complete movies of every workup or every alternative branch.

Correction to the initial chat suggestion: the actual Small-seq module uses sequential
3′/5′ RNA-adapter ligation, not poly(A) tailing. Do not invent a tailing reaction in
Small-seq. A separate explicitly illustrative tailing primitive can be tested without
claiming it is part of that protocol.

CUT&Tag and SLAM-seq remain optional future extensions, not required for this goal.

## Milestones and acceptance

- [x] Native backend/FFmpeg prototype: actually export and inspect an MP4, an SVG
      still and a PNG. Record rendering throughput; speed is measured, not assumed.
- [x] Versioned animation data, deterministic evaluator, timeline/narration holds,
      sparse default visuals, detail flags and clear invalid-event diagnostics.
- [x] SMART-seq reference movie showing several real enzyme actions. Visually inspect
      establishing/action/product frames and verify protocol/version provenance.
- [x] Focused 10–30s clips for the remaining agreed coverage, longer only when the
      configured narration needs it. Reuse event primitives rather than copy scenes.
- [x] Meaningful checks: direction/complement/correspondence and modification/end
      compatibility where events claim them; circle continuity; nick versus gap;
      source-scope/unknown handling; deterministic seek; narration holds; FFmpeg
      failure reporting; independent native exports without a browser.
- [x] Export each accepted example, probe codec/dimensions/duration/frame count,
      inspect representative frames/contact sheets, correct clipping/illegibility.
      Keep artifacts outside checkout and record output paths in handoff.
- [x] Document CLI, dependencies, styling/timing customization, event-adapter authoring,
      limits and next steps. Verify chemistry/site import paths have no new mandatory
      dependency. Update this checklist to actual evidence and outstanding limits.

Deliverable: runnable tool plus presentation-ready example MP4s, sources for the scenes,
useful still exports, tests and a reproducible handoff. A JSON exporter alone or a mock
animation without native MP4 output does not complete the goal.

## Initial prototype checkpoint — 2026-10-31 (historical)

Active goal was created from this plan. It remains incomplete. Current optional
package: `animation/` (`model.py`, `render.py`, `export.py`, `examples.py`, CLI,
self-tests and README). Cairo native drawing and direct raw-frame FFmpeg encoding
work; no browser is involved. No dependencies were installed in this session because
pycairo/Cairo and FFmpeg were already present.

Implemented so far: versioned JSON round-trip; explicit strand/event data; checked
primer+new-base complement mapping and extendable ends; direct timestamp evaluation;
narration-paced establishing/action/final holds; timing overrides/detail switches;
native PNG/SVG/MP4; encoder errors and atomic video replacement; MP4 provenance
sidecars. Exact ligation is intentionally rejected until two-end chemistry is modeled.
The event/layout vocabulary is provisional, not a finished molecular animation model.

SMART-seq imports the original protocol module (not SMART-seq2) and uses coarse native
kit regions in four initial RT/TS/later-PCR shots. This is a partial movie: native
tagmentation/index-PCR/workup and terminal-addition geometry remain outstanding. It
does not reconstruct a source-unresolved first PCR conversion or assert track lengths
as actual native nucleotide counts.

Verification: `python3 -m animation.selftest --native` passes12 checks, including
incorrect copied bases, blocked3prime ends, template coordinates, unknown/exact scope,
narration duration, schema version, seek determinism, native stills and real encoder
output/error reporting. FFprobe confirmed H.264/yuv420p/30fps and actual frame counts.
Representative reference/SMART frames were visually inspected; label collisions in
the reference were corrected. Further visual work is required for polished movies.

Artifacts (temporary, not committed):

- `/tmp/chem-animation-reference-1080.mp4`:1920×1080,362frames,12.067s, finaldetail layout.
  Render+encode10.99s (about32.95frames/s) on this machine, not a general benchmark.
- `/tmp/chem-animation-smartseq-draft.mp4`:960×540,1964frames,65.467s; partial4shot draft.
- `/tmp/chem-animation-reference.svg`, `/tmp/chem-animation-reference-detail.png`,
  `/tmp/chem-animation-smartseq-mid.png`, `/tmp/chem-animation-smartseq-switch.png`.
- Each MP4 has a matching `.manifest.json`; benchmark metadata is
  `/tmp/chem-animation-benchmark.json`.

Next: introduce explicit checked covalent connectivity/junction and interval
correspondence (including circle closure versus a mere bent path), then finish
SMART-seq and implement CRISPR-MIP/BOLORAMIS in the accepted order. Do not satisfy
the remaining coverage merely by adding enzyme labels to generic before/after tracks.
No commit/push or staging performed for this animation work.

## Completion evidence — 2026-10-31

The native tool and all agreed coverage adapters are implemented. The initial
prototype checkpoint above is historical; its partial-movie/unsupported-ligation
statements are superseded by this section. The current format is
`ngs101.animation/2`, with checked explicit junctions, ordered routes, retained
interval lineage, paired spans, copy mappings and separate camera/timing controls.

- `animation/model.py` validates declared molecular relations and timing.
  `evaluate.py` computes a frame directly from its timestamp. `render.py` uses
  native Cairo; `export.py` streams directly to FFmpeg.
- `protocols.py` explicitly adapts the existing ten protocol constructors, recording
  original module/source identity. It does not infer a mechanism from StepKind.
- The SMART-seq principal movie has nine shots through the selected original
  Epicentre Nextera/index-PCR route,142.934s at the default narration estimate.
  Its native unknowns and later-PCR-cycle qualification remain visible.
- CRISPR-MIP checks the19nt landing of the long probe separately from its retained
  unpaired backbone, derives the copied gap and checks the source circle rotation.
  BOLORAMIS checks the actual derived repeat and keeps the displaced product3′end
  beside the polymerase. Folding does not imply covalent closure.
- Remaining coverage includes loaded-adapter transfer, sequential RNA-adapter
  ligation, masking and adapter degradation, resection/overlap/nick sealing,
  restriction/fill-in/proximity enrichment, tethered capture/extension, chemical
  DNAU versus copied sequence, and successive RNA transcription on retained DNA.
- Default titles/labels are sparse; detail flags are optional. A narration estimate,
  explicit speech duration, per-shot timing file and global playback speed control
  the establishing/action/product intervals. No audio synthesis/alignment claimed.

Run `python3 -m animation.selftest --native`:30checks pass. These include wrong
base/coordinate/end chemistry rejection, gap-before-seal rejection, noncovalent
folding versus closure, circular template/repeat checks, preserved molecule identity,
protectedC/physicalU distinction, source-version/unknown handling, deterministic
seeking, narration pacing, native stills/video and encoder failure handling. CLI
timing/JSON/SVG/detail flags and invalid-duration errors were also exercised.
Ordinary source-model loading works under `python3 -S` without Cairo. No new
mandatory dependency was added to chemistry or site code.

Artifacts: `/tmp/chemistry-animations-2026-10-31/`. Eleven full movies (ten protocols
plus the synthetic reference),eleven focused clips,sourceJSON,manifests,PNG/SVG
stills,per-shot contact sheets and decoded-video frames are present. All22MP4s
are H.264/yuv420p,1920×1080,30fps. Full movies total15499frames/516.635s; focused
clips are12.067–16.367s. `verification.json` records actual codec/dimensions/frame
counts,source-scene consistency and successful resume verification. Each native
movie was probed and an actual decoded frame visually reviewed; all35shot layouts
were sampled at start/mid/end, including optional-detail bounds checks.

Tested environment: Python3,pycairo1.29.0,Cairo1.16.0,FFmpeg4.4.2 with libx264.
The full142.934sSMART-seq export/verification pipeline took107.59s on this machine;
this includes encoding/probing/clip overhead and is not a universal throughput claim.
No Rust acceleration is justified by that initial measurement.

Reproduce with `python3 -m animation.batch --directory /tmp/new-chemistry-movies`.
README documents dependencies,CLI,adapter authoring,timing/styling and scientific
limits. The separate Spatial VDJ reader-page improvement is recorded in HANDOFF.md.
No staging,commit,push or site deployment was performed. Optional future extensions
include sequencing-by-ligation cycles,more native exact geometry where sources
support it,audio muxing and additional protocols; none replaces a missing source.


## 2026-10-31 — presentation revision after movie review

The first exported set and its completion counts above are historical. User review
identified essential teaching gaps. The revision uses schema `ngs101.animation/3`:
functional region intervals/callouts, timed highlights, retained cut-piece motion,
linked identical enzyme subunits and persistent enzymes across RT steps.

SMART-seq now labels the oligo(dT) RT primer, poly(A) tail and IS-PCR handle. RNA
remains paired through terminal addition and template switching. The first
complementary DNA grows from a short PCR primer before a duplex PCR cycle is shown.
PCR landing sites and both index-primer 5′ tails are highlighted. Loaded Tn5
homodimers dock, cut, transfer, release and expose separated fragments. The SPRI
presentation scene is omitted. A final slide labels functional library anatomy.
The native original-kit TSO/IS-PCR footprint and Epicentre end ordering remain
unassigned; no later SMART-seq2/modern dual-index architecture is substituted.

The complementary examples received the same region/role treatment: MIP arms and
gap, BOLORAMIS ligation arms/RCA primer, RNA-adapter ends/UMI, Gibson overlaps,
Hi-C cut pieces and contact junction, Drop-seq full barcode/UMI primer inheritance,
scRRBS chemical C→U label and CEL-Seq's actual T7 promoter region.

Playback defaults to1.5× throughout CLI/export/batch/stills. `--speed` remains a
positive global multiplier; `--speed 1` restores authored timing. A1.25× native
batch smoke test also verified the override, manifest, focused clip and still path.
Thirty-six checks pass, including six presentation regression checks and two native
render/encoder checks. Final full-resolution artifacts are in
`/tmp/chemistry-animations-presentation-final/`: 11 full movies and eleven
focused clips, all1080p30H.264/yuv420p, 37 shots, 11315 full-movie frames /
377.171s total. SMART-seq is 111.334s with ten shots. All eleven exports pass
source-scene/settings resume checks and ffprobe. Revised contacts and actual decoded
SMART-seq switching, post-cleavage, primer-binding and final-library frames were
visually inspected. Documentation freshness and diff whitespace checks pass.
No staging, commit or deployment performed.


## Pending sequence-rendering revision (user review, 31 October 2026)

The user rejected the current line-based presentation. The exported v3 movies above
are review checkpoints, not accepted presentation assets. Before the next export:

- Draw nucleotide sequences by default, with direction arrows matching the site.
- Attach region labels to oligo coordinates so they move with the molecule.
- Tighten duplex spacing; enzyme overlap with the opposite strand is acceptable.
- Model pairing explicitly, independently of geometric distance and transitions.
- Bring incoming oligos from outside the viewport; remove arbitrary primer-tail bends.
- Read primary Tn5/transposon literature and inspect its mechanism illustrations.
- Correct strand-transfer products, move irrelevant fragments fully off-screen, and
  animate actual gap repair instead of replacing them with a magically filled duplex.
- Show every source-defined final-library component and sequence; keep unresolved
  original-kit values explicit rather than substituting a later protocol.

This revision has not been implemented at the requested stage/commit checkpoint.
Keep the default1.5× playback and tunable `--speed` multiplier.

## Substrate-driven revision and RNase H — 31 October 2026

User clarification: **SMART-seq2 is the requested animation**. Finishing the original
2012 SMART-seq adapter is outside this goal; keep it distinct and label its current
prototype status accurately. The source-era Nextera details relevant to SMART-seq2
are integrated into its note and source-qualified example.

The revised engine preserves nucleotide identities, covalent bonds, duplex contacts,
terminal chemistry and enzyme capabilities independently of drawing coordinates.
Reaction products derive from these substrates; an expected final sequence is not
an input to polymerase extension. SMART-seq2 now exercises this path through RT,
template switching, Tn5 transfer, repair and successive index-PCR copies. It remains
a separately named example, rather than replacing original SMART-seq.

`python3 -m animation rnaseh --output /tmp/rnaseh.mp4` renders a synthetic H1-class
mechanism example. It cuts eligible RNA backbone bonds in a supplied RNA–DNA hybrid,
retains all RNA fragments and DNA, and preserves duplex contacts. Selected sites
illustrate cleavage rather than predicting specificity, rate or fragment sizes.
The four-RNA footprint and continuous antiparallel DNA partner are engine
preconditions. H2 and native CEL-Seq kit identities are not inferred. See
[the source review](../animation/rnaseh-research.md).

Outside the current SMART-seq2 goal, remaining integration includes wiring and
visually reviewing other generated protocol views, explicit RNA-primer polymerase
capability, readable displaced RCA tails, and correction of the CRISPR-MIP universal
linear-DNA exonuclease wording. Existing legacy protocol adapters must not be
described as fully substrate-driven merely because their endpoint checks pass.

## SMART-seq2 completion audit — 31 October 2026

Historical checkpoint: subsequent user review found that this revision selected
only newly added PCR primers for extension. Its endpoint and causal-copy checks
did not detect omitted extension of the original recessed3′ ends. The all-ends
correction below supersedes this movie and its reaction-coverage claim.

The user's clarified target is SMART-seq2. Original2012 SMART-seq remains a separate
prototype; its further migration is outside this completed target.

| Requirement | Current evidence |
|---|---|
| Native MP4, presentation pace | `/tmp/chemistry-animation-smartseq2-reviewed/smartseq2.mp4`: H.264/yuv420p,1920×1080,30fps,4368frames,145.6s, default1.5×. Batch ffprobe and source-scene/settings checks pass. |
| Actual sequences, attached labels, arrows and compact duplexes | Twelve base-resolved scenes; off-screen primer entry; all-shot contact sheet and decoded mechanism frames inspected. Edge-label regression checks cover both polymerases in all three index-PCR scenes. |
| Substrate-driven RT, switching and PCR | `smartseq2_simulation.py` uses the shared connectivity engine; actual RNA persists through switching and explicit thermal unpairing precedes second-strand synthesis. Input mutations change the final product; cutting a template changes where extension stops. |
| Tn5 transfer and repair | Two explicit homodimers,9nt stagger, conserved target/adapter IDs, only18inter-fragment duplex contacts removed, flanks leave the viewport, original3′ ends extend9nt then pause at occupied ME complements before release and adapter copying. Causal cut-position tests and pairing-spacing checks pass. |
| Indexed PCR and final molecule | First copies stop at physical template ends; subsequent copies obtain both tails from existing templates. The final169nt duplex has all nine feature groups and every nucleotide shown on two continued rows. Endpoint declaration is an independent assertion, never the simulator's product input. |
| Source qualifications | Source-era Nextera component orders/index widths incorporated into the SMART-seq2 note and model annotation. Movie/manifest identify the2014 TSO, synthetic transcript, N501/N701 example pair and unassigned ME-release mechanism. Actual study-well assignments are not invented. |
| Verification and handover |52animation tests including native export;4protocol checks; source-oligo guard;227offline fetcher checks; protocol page build and0problem duplexes; documentation freshness and whitespace checks. Plan/state/README/handoff distinguish this completion from unfinished protocol migrations. |

The artifact directory also contains the scene JSON, export/source manifest,
`verification.json`, a focused RT clip, SVG/PNG action stills, an all-shot contact
sheet,13decoded-video frames with timestamps, and a separate final-library SVG/PNG.
Decoded frames were inspected after the final label-layout correction. Prior
`smartseq2-simulated-review`, `smartseq2-presentation` and `smartseq2-final` exports
are superseded review checkpoints. No video, third-party original or generated HTML
is checked in. No staging, commit, push or deployment was performed.

## PCR eligibility correction — 31 October 2026

User review identified the missing counter-extension of original strands along
index-primer5′ overhangs. The recipe now binds both oligos, discovers every eligible
3′ end in the selected molecular complexes and extends until no eligible site remains.
Pre-amplification uses the same discovery operation. This is a reaction-selection
correction; adding two decorative enzyme glyphs would leave the chemistry wrong.

The projector consumes every returned operation, updates original strand sequences
and nucleotide identities, and derives template coordinates and binding spans from
actual atoms/bonds. The first two index-PCR scenes have four extensions. Later
fully paired template ends do not provide further template bases, leaving two new
primer extensions. Circular strand-displacing reactions still require budgeted RCA;
RNA primers require an explicit polymerase capability and retain their RNA sugars.
Reaction scope, blocked ends and nucleotide-pool failures remain explicit.

The polymerase and ligase drawings are half their previous width and height. Tn5
homodimers retain their existing size. Updated captions describe both directions
of extension. Generic discovery and current four-end projection have regression
checks; the final library sequence remains unchanged.

Updated artifact: `/tmp/chemistry-animation-smartseq2-all-ends/smartseq2.mp4`,
149.134s,4474frames,1920×1080,30fps,H.264/yuv420p,default1.5×. Its batch verification
checks codec/settings and includes a full contact sheet and focused RT clip. Actual
decoded frames from both four-extension scenes were inspected; all four smaller
polymerases and their labels are visible. The ligase size was checked on a CRISPR-MIP
closure still. The previous SMART-seq2 movie is superseded.

58animation tests pass, including end discovery without oligo-role selection,
blocked ends, reaction-vessel scope, multi-end rollback, RNA-primer capability,
mixed RNA/DNA projection and four-extension PCR scenes. Documentation and handover
now record the correction and preserve the wider migration limits. No files staged,
committed or deployed.

## BOLORAMIS padlock example — completion audit,31October2026

The selected source ACTN1A1 example now has five substrate-driven scenes: source
probe/primer co-hybridization, folding and short-arm binding, SplintR nick sealing,
Phi29 RCA, and a conserved-product overview. It uses the existing simulator and
protocol-owned sequence definitions. Obsolete BOLORAMIS animation recipes that
introduced the RCA primer after ligation have been retired; the CLI and public
adapter use the same corrected implementation.

| Requirement | Evidence |
|---|---|
| Source chemistry and eligibility | Published article and Supplementary TableS3 checked; probe and primer arrive together; nucleotide pool enters after ligation. The closed complex has one eligible3′ end. Missing phosphate or unsupported pairing is rejected. |
| Folding and ligation | Covalent connectivity is unchanged by folding; the bound primer follows the backbone with24logical-pixel separation. The short arm approaches above the existing hybrid. Ligation adds only the actual last-to-first bond and consumes the terminal groups. |
| RCA and displacement | Actual template traversal generates180new nucleotides after the23nt stock primer. Each circle position is copied three times. Only paired suffix bases remain on the template; older DNA spreads into legible coils. RNA contacts are displaced while its sequence, sugars and backbone remain unchanged. |
| Causal controls | Changing the probe code changes the product; shortening the primer or synthesis budget changes actual output. Disabling RNA-displacement capability stalls at the occupied hybrid and retains RNA contacts. No expected product is supplied to polymerase. |
| Conserved overview | The final re-layout retains exactly the same atoms, backbone and duplex contacts as the preceding RCA snapshot. All203product bases and the original60nt circle remain represented; RNA is displayed separately. |
| Presentation | Native1080p30H.264/yuv420p MP4,76.934s,2308frames,default adjustable1.5×. Half-size polymerase/ligase glyphs, attached region labels, nucleotide letters and strand arrows. Counter/progress bar are opt-in. |
| Verification |72 animation tests,26protocol checks and35source-oligo checks pass. Contact sheet and decoded folding/ligation/RCA/product frames inspected. Source/scene/settings consistency and ffprobe verified. Documentation freshness and whitespace checked. |

Artifacts are under `/tmp/chemistry-animation-boloramis-presentation/`: full movie,
focused RCA clip, scene/export manifests, contact sheet, decoded review frames,
verification record and a final-product SVG/PNG with explicit figure-level scope.

Scope remains selected and qualified: the RNA is a conditional25nt match, commercial
NxGen37°C RNA-hybrid passage is inferred from mechanistic evidence rather than
independently measured here, and new T-equivalents retain unresolved dT/aminoallyl-dU
occupancy. No modification or crosslink sites, native yields, mismatch-tolerance
probabilities or full-transcript geometry are invented. Other multiplex members,
imaging readout cycles and unrelated protocol migrations remain outside this example.
The source review records original locators and inspected sections. No absence claim
or publication-status promotion was introduced.

The concurrent SMART-seq2 refinement is exported at
`/tmp/chemistry-animation-smartseq2-no-timeline/smartseq2.mp4`: PCR enzymes enter at
the extension interval, while explicit prebound RT persists across switching.
`--timeline` enables the otherwise hidden step counter/progress bar. Batch resume
now compares the full rendering style so old visible-footer videos cannot be reused
under the new defaults. No staging, commit, push or deployment performed.


### BOLORAMIS transition refinement (2026-10-31)

All incoming views now ease camera and conserved nucleotide paths during the
establish hold, including the final product rearrangement. Transitions require
identical molecular frames and cannot conceal bond changes. SplintR ligase leaves
when sealing finishes; the joined bases and bond retain a pink mark through RCA
and the product view. The final retained enzyme is labelled Phi29 polymerase.
Revised artifacts: `/tmp/chemistry-animation-boloramis-refined/`.
76 animation checks pass, including endpoint continuity and rejection of a
transition across differing molecular states. Native decoded frames reviewed.

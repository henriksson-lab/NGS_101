# Ten additional chemistry animations

Goal started 1 November 2026. SMART-seq2 and BOLORAMIS are the completed visual baseline.

| Example | Distinct reaction focus | Owner |
|---|---|---|
| CRISPR-MIP | Gap fill, ligation, end-dependent exonuclease cleanup | capture worker |
| Small-seq | Activated RNA adapter ligation and unused-adapter processing | capture worker |
| In-situ Hi-C | Restriction, labelled fill-in, selected proximity ligation | capture worker |
| CEL-Seq | Promoter-dependent de novo T7 transcription | synthesis worker |
| Drop-seq | Bead-bound barcoded RT | synthesis worker |
| scRRBS | Bisulfite conversion and copied versus physical bases | synthesis worker |
| Gibson | Resection, overlap binding, synthesis and nick sealing | assembly worker |
| ATAC-seq | Tn5 dimer strand transfer and repair | assembly worker |
| RNase H | RNA-specific endonuclease action on a hybrid | assembly worker |
| Golden Gate | Type IIS cleavage, cohesive assembly and site removal | root |

Each factory must consume actual molecular substrates and operation deltas; explicit
synthetic choices stay qualified. Molecular specifications are improved when needed,
without importing chemistry from a related protocol on apparent similarity alone.
Workers inspect primary originals and record locators. Root owns shared engine,
public adapter integration, native exports and cross-example audit.

Acceptance: chemically correct generated frames, conserved identities/off-products,
clear enzyme roles, attached readable sequence annotations, concise presentation,
no unmotivated geometry cuts. Incoming transitions preserve chemistry; reagent or
molecular additions need an explicit view/event or an explained boundary. Inspect
actual decoded MP4 frames as well as contact sheets, including reaction fronts and
transition midpoints. Iterate on defects. Default1080p30 H.264, adjustable1.5×,
footer off. Source limits, tests, artifact hashes and reviews recorded per example.

Status: complete — ten movies exported, audited, refined and accepted. No publication, push or new commit requested.


## Implementation and source review

Canonical CLI factories are in `capture_ligation.py`, `synthesis_examples.py`,
`assembly_transposition.py` and `golden_gate.py`, under `animation/`. Shared reaction
primitives remain substrate-driven. Independent review workers crossed groups;
root separately inspected Golden Gate and all initial contact sheets. Each source
review records primary originals, exact inspected sections and explicit limits:

- `animation/capture-ligation-source-review.md`
- `animation/synthesis-source-review.md`
- `animation/assembly-transposition-source-review.md`
- `animation/golden-gate-source-review.md`

Golden Gate's original Figure1 was inspected and its protocol research note updated.
No protocol was promoted based on animation quality. Gibson retains a main-article
access limit and uses source-qualified generic enzyme roles. ATAC explicitly bounds
native ME′ release chemistry; scRRBS uses synthetic primers for its selected copying
example. These are declared scopes, not accepted missing-information claims.

## Review iterations

1. Worker source reading and native small-frame reviews, then root cross-example
   contact review. Removed a redundant Small-seq caption and enlarged Hi-C's contact.
2. All generated definitions validate and account for every active atom; 27 incoming
   boundaries interpolate identical molecular frames by nucleotide identity. Frames
   involving new stock oligos use an explicit atomic reagent-introduction operation.
3. Cross-process repeatability exposed unordered released-pool and end-exchange
   metadata. Canonicalized unordered released IDs and stabilized Golden Gate input
   traversal; causal chemistry was unchanged. Added a cross-process regression.
4. Independent review of actual decoded1080p MP4s exposed unborn CEL-Seq RNA captions
   and an oversized ATAC ME/ME′ spacing. Suppressed captions/termini until atoms exist;
   derived all four ME′ positions from their transferred partners at24logical pixels.
   Regression checks every paired base across docking phases. Both movies regenerated.
5. Refined remaining CRISPR-MIP, Small-seq and Hi-C caption crowding. Biotin marks
   remain on physical bases while their labels use anchored offsets/leader lines.

## Export and acceptance records

Artifacts live under `/tmp/chemistry-animations-ten-more-20261101/`, divided into
`capture/`, `synthesis/`, `assembly/` and `golden-gate/`. Every full movie has a focused
clip, scene/export manifest, native contact sheet, verification record, decoded review
frames and a final SVG. Per-example audit records check scene correspondence,
codec/dimensions/frame rate/frame count and capture code/movie hashes. Independent
review records and root final acceptance are retained with the artifacts.

The ordinary website/chemistry build remains offline and standard-library-only.
Animation-only Cairo/FFmpeg dependencies are unchanged. No generated videos, source
originals, runtime binaries or site output are checked in.


## Final verification

112 animation checks pass, including enzymatic causal controls, complete active-atom
track coverage, all27equal-state transitions, cross-process determinism, actual native
encoding, unborn-output suppression and allfour tightly paired ATAC adapter complexes.
Ten1080p30H.264/yuv420p movies contain37shots and13,164frames (438.8seconds total).
Every current scene matches its export manifest after canonicalizing unordered pools.
Actual decoded reaction and transition frames were reviewed and findings resolved.
Final artifact index: `/tmp/chemistry-animations-ten-more-20261101/README.md`.
`summary.json`, `audit.json` and per-group independent reviews retain exact hashes.
Focused clips and finalSVGfigures are included. No videos/sources/runtime binaries
were checked in; no staging, commit, push or deployment was performed for this goal.

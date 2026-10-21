# Annealing tab: shared Python prototype

Status: implemented and verified; browser-runtime exception authorized by the user.
Date: 2026-10-21.

## Decision being evaluated

Run the same Python candidate search, pairing checks and SVG renderers in a Pyodide
Web Worker. No TypeScript chemistry implementation or Rust migration is required.
The normal Python build remains standard-library-only. The approved exception is a
pinned, on-demand browser runtime for this tab; it is now used by the annealing tab.

The user authorized the browser-only exception on 2026-10-21. [AGENTS.md](../AGENTS.md)
records its narrow scope: pinned Pyodide for this tab; Python/site builds stay offline and
standard-library-only, with no additional Python packages, npm dependencies or Rust toolchain.

For production, the runtime is Pyodide 314.0.7 loaded on demand from its pinned
manufacturer-documented jsDelivr URL. That avoids checking in runtime binaries or adding
a network dependency to the Python/site build. The prototype can serve a separately downloaded core locally, or test the proposed CDN
path with `--cdn-runtime`; neither mode changes the normal site build.
Submitted sequences are data passed to a fixed Python API, not code, and are not uploaded.

## Implementation

- `lib/annealing.py`: single shared exact candidate search for two actual ssDNA strands
  in 5′→3′ order. Linear candidates use maximal complementary spans. Padlock candidates
  locate both terminal arms on the other supplied strand and preserve alternate target
  assignments, repeated sites and competing maximal arm boundaries. N is allowed in unpaired regions such as padlock backbones and interrupts exact
  matching; it is never counted as a known pair. RNA/other ambiguity symbols are rejected.
- `LinearCandidate.scene()` uses the existing `Scene.anneal()` check and strand layout.
  Full strands and their overhangs are retained; other regions are explicitly unpaired.
- `PadlockCandidate.topology()` uses the existing `PadlockCaptureTopology` and `FoldLayout`
  for pairing and planarity checks and rendering. Optional explicit target sites now let
  this shared constructor select among repeated matches; omitted sites still require a
  unique match, preserving existing caller behavior. Supplied coordinates are checked
  against the actual complementary arm sequence.
- `assets/anneal_worker.mjs` loads the four existing/shared Python modules and dispatches
  only fixed search/render operations. Input strings are passed through JSON, never
  interpolated into executable Python. The worker serializes requests.
- `assets/anneal_ui.js` provides two sequence inputs, linear/padlock/both modes, minimum
  span/arm length, candidate selection, cancellation, synthetic examples and SVG export.
  Editing clears stale results while retaining an idle runtime. Cancelling a calculation
  terminates the worker. Calculations have an explicit timeout.
- `tools/build_anneal_prototype.py` builds outside the checkout from an explicitly supplied
  runtime directory; it never fetches or installs runtime packages. Its manifest records
  core file sizes and SHA-256 hashes. The standalone prototype shares the production page/model builder and stays outside
  the checkout. `build_anneal.py` supplies the normal site page and source bundle; the
  site assembler explicitly includes the worker and JSON fetch dependencies.

Limits are 3,000 bases per input, 250 returned candidates and 500,000 counted search
operations. Truncation is explicit. Candidates describe exact binding geometry, not a
preferred thermodynamic fold, binding probability, polymerase activity or successful
ligation/circularization. Geometry can be possible without the corresponding reaction.
Non-planar selected layouts fail visibly through the shared renderer rather than being
silently repaired. Backbone and target gaps remain source-input sequence, not inferred
reaction products. Positions are zero-based half-open in the API and 1-based in the UI. Render responses
also expose original-strand paired spans and overhang lengths, or target/probe arm,
backbone and intervening target-gap intervals, without asserting a reaction product.

## Feasibility evidence

The separately downloaded pinned runtime core has five files totalling 13,531,207 bytes
before HTTP compression. These files and downloaded binaries live only in `/tmp`.

A real isolated browser run on the development machine (HeadlessChrome 104.0.5112.79):

| Measurement | Observed |
|---|---:|
| Cold Python runtime startup | 5.409 s |
| Model import/bundle loading | 0.714 s |
| Synthetic padlock candidate search | 7.3 ms |
| First padlock SVG render | 67.8 ms |
| Subsequent linear SVG render | 10.7 ms |
| Main-thread timer beats during startup/render | 259 |
| Maximum observed 25 ms timer gap | 290 ms |

A preliminary smaller smoke test took 4.9 s including startup/render and had a 31 ms
maximum timer gap. These are individual local-file-served observations, not a guaranteed
loading time, distribution benchmark or measurement of CDN transfer. Node 20 also loaded
this runtime, but its first startup was about 31 s; browser measurements are the relevant
UI evidence and startup varies materially with environment.

A separate fresh-profile run using the proposed pinned CDN path also passed all browser
checks: runtime startup 4.413 s, model import 0.990 s, search 4.9 ms, padlock render 56.9 ms
and subsequent linear render 8.0 ms. It had 228 main-thread timer beats and a 263 ms
maximum timer gap. This confirms the cross-origin worker/runtime loading path works in
the tested browser. It remains a single measurement on this machine/network, not a
guaranteed first-visit transfer or startup time.

The real-browser test verifies padlock and linear results, selecting another candidate,
invalid sequence handling, an unknown N-containing padlock backbone, reuse of the loaded
runtime across edits, proper SVG namespace
and the SVG download action, cancellation and recovery from a transient model-bundle load
failure without reloading the page. The Python suite independently checks coordinates, pairing,
repeated sites, ambiguity, truncation and explicit-site validation.

## Reproduce

For the production page, run `python3 build_site.py`, then
`python3 tools/benchmark_anneal_browser.py --site _site --url-prefix /NGS_101`.
The following standalone option remains useful for comparing local and CDN runtimes.

Use a separately downloaded, version-matched Pyodide core containing `pyodide.mjs`,
`pyodide.asm.mjs`, `pyodide.asm.wasm`, `python_stdlib.zip` and `pyodide-lock.json`.
See the official deployment documentation below; do not put binaries in the checkout.

```sh
python3 tools/build_anneal_prototype.py --runtime /tmp/chem-anneal-pyodide --output /tmp/chem-anneal-review
python3 tools/benchmark_anneal_browser.py --prototype /tmp/chem-anneal-review
# Optional: rebuild the opt-in prototype with --cdn-runtime to test CDN loading.
python3 tools/selftest_annealing.py
```

The benchmark uses an isolated temporary Chrome profile and localhost server. It writes
its report and screenshot under `/tmp` by default. A browser that can start and bind a
local debugging port is required; restrictive sandboxes may prevent this manual test.

## Production verification

The normal page uses the same pinned CDN path already tested in the prototype. The site
build includes `anneal.html`, `assets/anneal_worker.mjs`, `assets/anneal_python.json`, UI/CSS
and shared navigation. Final checks passed against the assembled site under `/NGS_101/`, including relative
worker/model URLs, navigation, startup recovery, cancellation, coordinates and SVG export.
The full site built all 407 diagrams, included 406 published schematics and passed its
link audit. Offline verification passed 47 annealing checks (including packaging), 1,288
shared checks and the existing sequence/sample-sheet JavaScript suites. The generated
source bundle contains only our four shared Python modules and is ignored build output.
No runtime binaries are included. The final assembled-page browser run observed 7.346 s
for runtime startup, 1.035 s for model import, 11.6 ms for search, 92.1 ms for padlock
rendering and 25.1 ms for a subsequent linear render. Startup varies with machine/network;
these are measured observations rather than a latency guarantee.

## Runtime references

- [Pyodide: deploying on static hosting](https://pyodide.org/en/stable/usage/downloading-and-deploying.html)
- [Pyodide: Web Worker execution](https://pyodide.org/en/stable/usage/webworker.html)
- [Pyodide: JavaScript API](https://pyodide.org/en/stable/usage/api/js-api.html)

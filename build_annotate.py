#!/usr/bin/env python3
"""Build the static, browser-local sequence annotation tool."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'lib'))
from page import head
from illumina import P5, TRUSEQ_READ1
OUT = ROOT / 'annotate.html'

def render_page():
    example = P5 + 'ACGTACGT' + TRUSEQ_READ1
    return '<!doctype html>\n<html lang="en">\n' + head('Annotate sequence · NGS 101', home='index.html') + '''
<link rel="stylesheet" href="assets/sequence_annotation.css">
<link rel="stylesheet" href="assets/sequence_alignment.css">
<main class="wrap annotation-tool">
<div class="annotation-intro"><h1>Annotate sequence</h1></div>
<form id="annotate-form" class="annotation-card">
<label for="sequence-input" class="field-title">Sequence <span class="muted">5′ → 3′</span></label>
<textarea id="sequence-input" rows="6" spellcheck="false" autocapitalize="characters" autocomplete="off" aria-describedby="input-help" placeholder="Paste DNA or a single FASTA record"></textarea>
<p id="input-help" class="muted">DNA or single FASTA · up to 100,000 nt</p>
<div class="annotation-options"><label><input id="convert-u" type="checkbox"> Convert U to T</label><label><input id="include-inferred" type="checkbox"> Include inferred parts</label></div>
<div class="annotation-actions"><button id="annotate-button" type="submit" class="primary-button">Annotate</button><button id="example-button" type="button" data-sequence="''' + example + '''">Set example sequence</button><button id="clear-button" type="button">Clear</button></div>
</form>
<p id="annotation-status" role="status" aria-live="polite"></p>
<noscript><p>This browser-local tool requires JavaScript.</p></noscript>
<details class="annotation-help"><summary>About matching</summary><p>Matches are exact known primer, adapter or component sequences, including reverse complements. Whitespace is removed; ambiguity letters interrupt matches. Coordinates are 1-based inclusive on the normalized query. Strand refers to the canonical component; local declarations can use the opposite orientation.</p><p>Highlighted regions show the longest matches. Select a primer or component below each region to inspect overlapping matches. A component may be part of a longer primer. Sequence equality does not establish protocol origin, modifications, binding or a preferred fold. Protocol suggestions rank shared informative parts; they are not confidence scores.</p></details>
<section id="annotation-results" hidden aria-labelledby="results-title">
<h2 id="results-title">Matches</h2><p id="result-summary"></p>
<div class="annotation-card"><div class="preview-header"><h3>Sequence overview</h3><div><button id="preview-prev" type="button">Previous</button> <span id="preview-range"></span> <button id="preview-next" type="button">Next</button></div></div><pre id="sequence-preview" tabindex="0" aria-label="Normalized sequence with exact matches"></pre><p id="region-detail" class="muted" role="status"></p><div id="region-alignments" class="region-alignments"></div></div>
<h3>Matched components</h3><p id="match-table-limit" class="muted"></p><div class="tw"><table class="match-table"><thead><tr><th scope="col">Position / strand</th><th scope="col">Component and aliases</th><th scope="col">Declaration</th><th scope="col">Protocols sharing this part</th></tr></thead><tbody id="match-rows"></tbody></table></div><button id="more-matches" type="button" hidden>Show more matches</button>
<h3>Related protocols</h3><p id="protocol-explanation" class="muted"></p><ol id="protocol-matches" class="protocol-matches"></ol>
</section></main>
<script src="assets/sequence_annotation.js" defer></script><script src="assets/sequence_annotation_ui.js" defer></script>
</html>\n'''

def main():
    OUT.write_text(render_page(), encoding='utf-8')
    print(OUT.name)
if __name__ == '__main__':
    main()

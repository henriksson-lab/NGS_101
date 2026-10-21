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
<main class="wrap annotation-tool">
<div class="annotation-intro"><p class="eyebrow">Sequence tools</p><h1>Annotate sequence</h1>
<p>Find exact matches to declared primers, adapters and their known parts. Your sequence stays in this browser.</p></div>
<form id="annotate-form" class="annotation-card">
<label for="sequence-input" class="field-title">Sequence <span class="muted">5′ → 3′</span></label>
<textarea id="sequence-input" rows="6" spellcheck="false" autocapitalize="characters" autocomplete="off" aria-describedby="input-help" placeholder="Paste DNA or a single FASTA record"></textarea>
<p id="input-help" class="muted">Up to 100,000 nucleotides. Whitespace is removed; IUPAC ambiguity letters interrupt exact matches. Coordinates refer to the normalized sequence.</p>
<div class="annotation-options"><label><input id="convert-u" type="checkbox"> Convert U to T</label><label><input id="include-inferred" type="checkbox"> Include inferred declarations</label></div>
<div class="annotation-actions"><button id="annotate-button" type="submit" class="primary-button">Annotate sequence</button><button id="example-button" type="button" data-sequence="''' + example + '''">Try a synthetic adapter example</button><button id="clear-button" type="button">Clear</button></div>
</form>
<p id="annotation-status" role="status" aria-live="polite">Paste a sequence to begin. The component catalogue loads when you annotate.</p>
<noscript><p>This browser-local tool requires JavaScript.</p></noscript>
<section id="annotation-results" hidden aria-labelledby="results-title">
<h2 id="results-title">Exact sequence matches</h2><p id="result-summary"></p>
<p class="muted">A match describes sequence equality. It does not establish protocol origin, modifications, binding or a preferred fold.</p>
<div class="annotation-card"><div class="preview-header"><h3>Sequence overview</h3><div><button id="preview-prev" type="button">Previous</button> <span id="preview-range"></span> <button id="preview-next" type="button">Next</button></div></div><p class="muted">Longest matches are highlighted; overlapping and contained parts remain in the table. Select a highlighted region to inspect its matches.</p><pre id="sequence-preview" tabindex="0" aria-label="Normalized sequence with exact matches"></pre><p id="region-detail" class="muted" role="status"></p></div>
<h3>Matched components</h3><p id="match-table-limit" class="muted"></p><div class="tw"><table class="match-table"><thead><tr><th scope="col">Position / strand</th><th scope="col">Component and aliases</th><th scope="col">Declaration</th><th scope="col">Protocols sharing this part</th></tr></thead><tbody id="match-rows"></tbody></table></div><button id="more-matches" type="button" hidden>Show more matches</button>
<h3>Protocols sharing matched parts</h3><p id="protocol-explanation" class="muted"></p><ol id="protocol-matches" class="protocol-matches"></ol>
<section class="annotation-card samplesheet-section" aria-labelledby="samplesheet-title"><h2 id="samplesheet-title">Illumina sample-sheet indexes</h2>
<p>Indexes are extracted only from matched, declared index templates. Select your instrument and sample-sheet workflow to obtain <code>index</code> (i7) and <code>index2</code> (i5). A free index-sized sequence alone is insufficient.</p>
<p id="index-summary" class="muted"></p><div class="sample-fields"><label>Instrument<select id="index-instrument"><option value="">Choose instrument</option></select></label><label>Sample-sheet workflow<select id="index-workflow"><option value="">Choose workflow</option></select></label><label>Sample_ID<input id="index-sample-id" type="text" value="Sample_1" autocomplete="off"></label><label>i7 candidate<select id="index-i7"><option value="">Choose i7</option></select></label><label>i5 candidate (optional)<select id="index-i5"><option value="">No i5 / single index</option></select></label></div>
<p class="muted">After choosing i7, i5 choices are limited to the same protocol declaration. This does not establish that independently listed index choices are a valid experimental pair.</p>
<div id="index-candidate-details"></div><ul id="index-unresolved" class="muted"></ul><p id="index-profile-note" role="status"></p><div id="index-fields"></div><div id="index-snippet" hidden><label for="index-output" class="field-title">Copyable sample-sheet text</label><textarea id="index-output" rows="6" readonly spellcheck="false"></textarea><p class="muted">This is a data-section snippet, not a complete run configuration. Preserve any mandatory orientation header shown above the section.</p><div class="annotation-actions"><button id="copy-index-tsv" type="button">Copy header + row (TSV)</button><button id="copy-index-csv" type="button">Copy CSV snippet</button></div></div><p id="copy-status" role="status" aria-live="polite"></p>
</section></section></main>
<script src="assets/sequence_annotation.js" defer></script><script src="assets/illumina_samplesheet.js" defer></script><script src="assets/sequence_annotation_ui.js" defer></script>
</html>\n'''

def main():
    OUT.write_text(render_page(), encoding='utf-8')
    print(OUT.name)
if __name__ == '__main__':
    main()

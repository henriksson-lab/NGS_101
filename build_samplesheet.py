#!/usr/bin/env python3
"""Build the browser-local primer-pair sample-sheet tool."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'lib'))
from page import head


def render_page():
    return '<!doctype html>\n<html lang="en">\n' + head('Sample-sheet indexes · NGS 101', home='index.html') + '''
<link rel="stylesheet" href="assets/sequence_annotation.css">
<main class="wrap annotation-tool">
<div class="annotation-intro"><p class="eyebrow">Sequence tools</p><h1>Illumina sample-sheet indexes</h1>
<p>Enter one primer pair per sample, in either order, written 5′ → 3′. Known template matches identify i7 and i5 automatically. Either primer can be left blank. Your sequences stay in this browser.</p></div>
<div class="annotation-card">
<div class="sample-fields"><label>Instrument<select id="batch-instrument"></select></label><label>Sample-sheet workflow<select id="batch-workflow"></select></label></div>
<p id="batch-profile" class="muted"></p>
<details><summary>Paste a list of primer pairs</summary><label for="batch-paste" class="field-title">Tab-separated: Sample_ID, primer 1, primer 2 (optional)</label>
<textarea id="batch-paste" rows="4" spellcheck="false" placeholder="Sample_1&#9;primer sequence&#9;primer sequence"></textarea><button id="batch-import" type="button">Replace table with pasted rows</button></details>
<p class="muted">Use a primer’s “Use for all i5” or “Use for all i7” button to reuse it across the table. Its role is checked automatically before copying. Identical index values across matching protocols are combined; conflicting values stay unresolved.</p>
<div class="tw"><table class="primer-table"><thead><tr><th scope="col">Sample_ID</th><th scope="col">Primer 1</th><th scope="col">Primer 2</th><th scope="col">index (i7)</th><th scope="col">index2 (i5)</th><th scope="col">Detection / actions</th></tr></thead><tbody id="primer-rows"></tbody></table></div>
<div class="annotation-actions"><button id="batch-add" type="button">Add sample</button><button id="batch-detect" type="button" class="primary-button">Extract indexes</button></div>
</div>
<p id="batch-status" role="status" aria-live="polite">NovaSeq X and manual bcl2fastq are selected by default. Enter primers to begin.</p>
<section id="batch-export" class="annotation-card" hidden><label for="batch-output" class="field-title">Copyable sample-sheet text</label><textarea id="batch-output" rows="9" readonly spellcheck="false"></textarea>
<p class="muted">This is a data-section snippet. Retain the remaining run settings and any orientation header shown here. An omitted primer produces an empty field; rows with unresolved supplied primers are excluded. An i5-only row has an empty index field and needs its i7 information before use for demultiplexing.</p>
<div class="annotation-actions"><button id="batch-copy-tsv" type="button">Copy header + rows (TSV)</button><button id="batch-copy-csv" type="button">Copy CSV snippet</button></div></section>
<p id="batch-copy-status" role="status" aria-live="polite"></p><textarea id="batch-copy-fallback" aria-label="Text to copy" readonly hidden></textarea>
<p class="muted">Only anchored, declared sample-index templates are used. A free index-sized sequence alone is insufficient. Sequence matches do not verify experimental primer pairing or protocol identity. Inferred declarations are excluded.</p>
<noscript><p>This browser-local tool requires JavaScript.</p></noscript>
</main><script src="assets/sequence_annotation.js" defer></script><script src="assets/illumina_samplesheet.js" defer></script><script src="assets/samplesheet_batch.js" defer></script><script src="assets/samplesheet_ui.js" defer></script>
</html>
'''


if __name__ == '__main__':
    (ROOT / 'samplesheet.html').write_text(render_page(), encoding='utf-8')
    print('samplesheet.html')

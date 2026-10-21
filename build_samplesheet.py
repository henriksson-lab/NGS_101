#!/usr/bin/env python3
"""Build the browser-local primer-pair sample-sheet tool."""
from pathlib import Path
import json
import sys
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'lib'))
from page import head
import nextera


def example_rows():
    """Synthetic index combinations on canonical Nextera primer arms; not a kit pair."""
    return [
        {"sampleId": "Example_1", "primer1": nextera.n5xx_primer("ACGTTGCA"),
         "primer2": nextera.n7xx_primer("AACCGGTA")},
        {"sampleId": "Example_2", "primer1": nextera.n5xx_primer("ACGTTGCA"),
         "primer2": nextera.n7xx_primer("GTCAGTAC")},
        {"sampleId": "Example_3", "primer1": "",
         "primer2": nextera.n7xx_primer("CTGACATG")},
    ]


def render_page():
    return '<!doctype html>\n<html lang="en">\n' + head('Sample-sheet indexes · NGS 101', home='index.html') + '''
<link rel="stylesheet" href="assets/sequence_annotation.css">
<main class="wrap annotation-tool samplesheet-tool">
<div class="annotation-intro"><h1>Illumina sample-sheet indexes</h1>
<p>Enter primers 5′ → 3′ in either order; either may be blank. i5/i7 are detected automatically.</p></div>
<div class="annotation-card">
<div class="sample-fields"><label>Instrument<select id="batch-instrument"></select></label><label>Sample-sheet workflow<select id="batch-workflow"></select></label></div>
<div id="batch-profile" class="muted"></div>
<details><summary>Paste a list of primer pairs</summary><label for="batch-paste" class="field-title">Tab-separated: Sample_ID, primer 1, primer 2 (optional)</label>
<textarea id="batch-paste" rows="4" spellcheck="false" placeholder="Sample_1&#9;primer sequence&#9;primer sequence"></textarea><button id="batch-import" type="button">Replace table with pasted rows</button></details>
<div class="tw"><table class="primer-table"><colgroup><col style="width:11%"><col style="width:23%"><col style="width:23%"><col style="width:10%"><col style="width:10%"><col style="width:23%"></colgroup><thead><tr><th scope="col">Sample_ID</th><th scope="col">Primer 1</th><th scope="col">Primer 2</th><th scope="col">index (i7)</th><th scope="col">index2 (i5)</th><th scope="col">Detection</th></tr></thead><tbody id="primer-rows"></tbody></table></div>
<div class="annotation-actions"><button id="batch-add" type="button">Add sample</button><button id="batch-examples" type="button">Set example primers</button><button id="batch-detect" type="button" class="primary-button">Extract indexes</button></div>
</div>
<p id="batch-status" role="status" aria-live="polite"></p>
<section id="batch-export" class="annotation-card" hidden><label for="batch-output" class="field-title">Copyable sample-sheet text</label><textarea id="batch-output" rows="9" readonly spellcheck="false"></textarea>
<p class="muted">Data section only: retain your run settings and any orientation header. Unresolved rows are excluded. Blank primers leave empty fields; i5-only rows still need i7 before demultiplexing.</p>
<div class="annotation-actions"><button id="batch-copy-tsv" type="button">Copy header + rows (TSV)</button><button id="batch-copy-csv" type="button">Copy CSV snippet</button></div></section>
<p id="batch-copy-status" role="status" aria-live="polite"></p><textarea id="batch-copy-fallback" aria-label="Text to copy" readonly hidden></textarea>
<noscript><p>This browser-local tool requires JavaScript.</p></noscript>
</main><script id="batch-example-data" type="application/json">''' + json.dumps(example_rows()) + '''</script><script src="assets/sequence_annotation.js" defer></script><script src="assets/illumina_samplesheet.js" defer></script><script src="assets/samplesheet_batch.js" defer></script><script src="assets/samplesheet_ui.js" defer></script>
</html>
'''


if __name__ == '__main__':
    (ROOT / 'samplesheet.html').write_text(render_page(), encoding='utf-8')
    print('samplesheet.html')

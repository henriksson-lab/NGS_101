#!/usr/bin/env python3
"""Build the offline restriction-analysis tab and licensed enzyme snapshot."""
import csv
import html
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'lib'))
from page import head

STATES={'unknown','unaffected','blocked','impaired','overlap_blocked','some_overlap_blocked','overlap_impaired','some_overlap_impaired','requires'}


def catalogue():
    data=json.loads((ROOT/'catalogue/restriction_enzymes.json').read_text())
    known={e['name']:e for e in data['enzymes']}
    for enzyme in known.values():
        enzyme['methylation']=dict.fromkeys(('cpg','dam','dcm'),'unknown')
    rows=list(csv.DictReader((ROOT/'catalogue/restriction_methylation.tsv').open(),delimiter='\t'))
    if len({r['enzyme'] for r in rows})!=len(rows):raise ValueError('Duplicate methylation annotations')
    for row in rows:
        if row['enzyme'] not in known:raise ValueError('Unknown methylation enzyme: '+row['enzyme'])
        if any(row[k] not in STATES for k in ('cpg','dam','dcm')):raise ValueError('Invalid methylation state')
        if any(not s.startswith('https://') for s in row['source'].split(';')):raise ValueError('Methylation source must be an identified HTTPS original')
        enzyme=known[row['enzyme']]
        enzyme['methylation']={k:row[k] for k in ('cpg','dam','dcm')}
        enzyme['methylation_sources']=row['source'].split(';')
        enzyme['methylation_checked']=row['checked']
    data['methylation_annotation_count']=len(rows)
    data['methylation_scope']='Locally curated factual annotations for exact named enzymes; source pages are referenced, not reproduced. Unannotated enzymes are unknown. Do not transfer annotations to isoschizomers or product variants.'
    return data


def methylation_select(name):
    return f'<label>{name}<select id="restriction-{name.lower()}"><option value="any">Any sensitivity</option><option value="unaffected">Not sensitive</option><option value="sensitive">Sensitive / conditional</option><option value="requires">Requires methylation</option><option value="unknown">Unknown</option></select></label>'


def render_page():
    return '<!doctype html><html lang="en">'+head('Restriction sites · NGS 101',home='index.html')+'''
<link rel="stylesheet" href="assets/sequence_annotation.css"><link rel="stylesheet" href="assets/restriction.css">
<main class="wrap annotation-tool restriction-tool"><h1>Restriction sites</h1>
<form id="restriction-form" class="annotation-card">
<label for="restriction-sequence" class="field-title">DNA sequence · 5′ → 3′</label>
<textarea id="restriction-sequence" rows="4" spellcheck="false" autocomplete="off" placeholder="DNA or one FASTA record" required></textarea>
<div class="restriction-controls"><label>Topology<select id="restriction-topology"><option value="linear">Linear</option><option value="circular">Circular</option></select></label>
<button id="restriction-analyze" type="submit" class="primary-button">Analyze</button><button id="restriction-cancel" type="button" disabled>Cancel</button>
<button id="restriction-example" type="button">Set example sequence</button><button id="restriction-clear" type="button">Clear</button>
<label class="restriction-file">Load FASTA<input id="restriction-file" type="file" accept=".fa,.fasta,.fna,.txt"></label></div>
<details class="restriction-custom"><summary>Custom enzymes</summary><p>Name, IUPAC motif, top and bottom cut offsets from the recognition start. Use <code>-</code> for a strand that is not cut (nicking); two <code>-</code> values leave cleavage geometry unassigned. Optional columns: CpG, Dam, Dcm sensitivity codes.</p>
<textarea id="restriction-custom" rows="2" spellcheck="false" placeholder="MyEnzyme GAATTC 1 5"></textarea>
<p>Codes: <code>unknown unaffected blocked impaired overlap_blocked some_overlap_blocked overlap_impaired some_overlap_impaired requires</code>. Built-in names cannot be replaced.</p></details></form>
<p id="restriction-status" role="status" aria-live="polite"></p>
<section class="annotation-card restriction-filters" aria-label="Restriction filters">
<div class="restriction-controls"><label>Enzymes<input id="restriction-search" type="search" placeholder="Name or comma-separated names"></label>
<label>Count<select id="restriction-count-mode"><option value="exact">Exactly</option><option value="any">Any number</option><option value="atleast">At least</option><option value="atmost">At most</option><option value="range">Range</option></select></label>
<label id="restriction-count-label">Number<input id="restriction-count" type="number" value="1" min="0" step="1"></label>
<label id="restriction-max-label" hidden>Maximum<input id="restriction-max" type="number" value="2" min="0" step="1"></label>
<button id="restriction-unique" type="button">Unique cutters</button><button id="restriction-all" type="button">All counts</button></div>
<details><summary>More filters</summary><div class="restriction-controls restriction-advanced">
<label>Count basis<select id="restriction-basis"><option value="cleavage">Cut sites</option><option value="recognition">Recognition sites</option></select></label>
<label>Enzyme class<select id="restriction-kind"><option value="double">Double-strand cutters</option><option value="nick">Nicking enzymes</option><option value="any">All, including unassigned cuts</option><option value="unknown">Unassigned cuts</option></select></label>
<label>Ends<select id="restriction-end"><option value="any">Any</option><option value="5prime">5′ overhang</option><option value="3prime">3′ overhang</option><option value="blunt">Blunt</option><option value="nick">Nick</option><option value="mixed">Mixed</option></select></label>
<label>Recognition length ≥<input id="restriction-min-length" type="number" value="1" min="1" max="100" step="1"></label>
<label>Recognition length ≤<input id="restriction-max-length" type="number" value="100" min="1" max="100" step="1"></label>
<label>Site context<select id="restriction-context"><option value="any">Any</option><option value="overlap_cpg">Overlaps CpG</option><option value="avoid_cpg">Avoids CpG</option><option value="overlap_dam">Overlaps Dam motif (GATC)</option><option value="avoid_dam">Avoids Dam motif</option><option value="overlap_dcm">Overlaps Dcm motif (CCWGG)</option><option value="avoid_dcm">Avoids Dcm motif</option></select></label>
''' + ''.join(methylation_select(k) for k in ('CpG','Dam','Dcm')) + '''
<label><input id="restriction-commercial" type="checkbox"> Supplier listed in snapshot</label>
<label><input id="restriction-custom-only" type="checkbox"> Custom enzymes only</label>
<label><input id="restriction-possible" type="checkbox"> Include possible matches in counts</label></div>
<button id="restriction-reset-filters" type="button">Reset filters</button></details></section>
<section id="restriction-results" hidden><div class="restriction-result-heading"><h2 id="restriction-summary">Results</h2><div><button id="restriction-export" type="button">Export summary TSV</button><button id="restriction-export-sites" type="button">Export sites TSV</button></div></div>
<p id="restriction-result-note" class="muted"></p><div class="tw"><table><thead><tr><th scope="col">Enzyme</th><th scope="col">Recognition</th><th scope="col">Cuts</th><th scope="col">Sites</th><th scope="col">Ends</th><th scope="col">CpG</th><th scope="col">Dam</th><th scope="col">Dcm</th></tr></thead><tbody id="restriction-rows"></tbody></table></div>
<button id="restriction-more" type="button" hidden>Show more enzymes</button></section>
<details class="restriction-about"><summary>Data and interpretation</summary>
<p>Counts are sequence-based predictions, without applying methylation effects or reaction conditions. The default requires exactly one certain complete cleavage and no additional possible cuts. A two-cut enzyme can contribute two cleavage sites for one recognition site. Linear-end cuts without a bond on either required strand are excluded; circular-origin sites are included.</p>
<p>Input accepts one FASTA record or IUPAC DNA, up to 100,000 bases. Ambiguous bases can yield possible sites. Site-context filters inspect recognition intervals plus overlapping flanks; they select enzymes without changing whole-sequence cut counts. Avoidance requires no certain or possible overlap.</p>
<p>Excessive ambiguity, repeated matches or run time stops analysis with an error rather than returning partial counts. Cancel stops the analysis worker.</p>
<p>Recognition starts are 1-based. Cut coordinates are boundaries after the numbered base; a circular-origin boundary is shown as the sequence length. Fragment lengths use top-strand reference boundaries, with no methylation correction. Display and position exports retain up to 1,000 recognition matches per enzyme; all counts remain complete and truncation is labelled.</p>
<p>Methylation sensitivity is separate from sequence motifs. Unknown sensitivity is never treated as resistance. Curated annotations cover exact named enzymes; other product variants and isoschizomers do not inherit them. Conditional sensitivity can depend on flanking sequence. Source links accompany annotated results.</p>
<div id="restriction-data-info"></div><details><summary>Catalogue license</summary><pre id="restriction-license" class="restriction-license"></pre></details></details>
<noscript><p>This tool requires JavaScript.</p></noscript></main>
<script src="assets/sequence_annotation.js" defer></script><script src="assets/restriction_analysis.js" defer></script><script src="assets/restriction_ui.js" defer></script></html>
'''


def build(output=ROOT):
    assets=output/'assets';assets.mkdir(parents=True,exist_ok=True)
    data=catalogue()
    (assets/'restriction_catalogue.json').write_text(json.dumps(data,ensure_ascii=False,separators=(',',':'))+'\n')
    (output/'restriction.html').write_text(render_page())
    return output/'restriction.html'


if __name__=='__main__':print(build().name)

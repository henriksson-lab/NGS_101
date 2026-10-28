#!/usr/bin/env python3
"""Build the offline primer-Tm tab; data are published facts, not a vendor engine."""
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'lib'))
from page import head
from chemdraw import tm_model_parameters
ROCHE='https://elabdoc-prod.roche.com/eLD/api/downloads/20bd78d2-bc2e-f011-2f91-005056a71a5d?countryIsoCode=XG'


def model_data():
    return {'version':1,'thermodynamics':tm_model_parameters(),
            'preset':{'id':'kapa_hifi_hotstart_readymix_kr0370_v16_25','name':'KAPA HiFi HotStart ReadyMix PCR Kit',
                      'mg_mM':2.5,'dntpEach_mM':.3,'primer_nM':300,'monovalent_assumption_mM':50,
                      'source':ROCHE,'revision':'KR0370 v16.25 · May 2025'}}


def concentration_fields():
    preset=model_data()['preset']
    fields=[('primer','Primer concentration (nM)','primer_nM',.001,1000000),
            ('mg','Mg²⁺ (mM)','mg_mM',0,100),
            ('dntp','Each dNTP (mM)','dntpEach_mM',0,25),
            ('monovalent','Monovalent ions (mM, assumed)','monovalent_assumption_mM',0,1000)]
    return ''.join(f'<label>{label}<input id="tm-{name}" type="number" min="{minimum}" max="{maximum}" step="any" value="{preset[key]}" required></label>' for name,label,key,minimum,maximum in fields)


def render_page():
    return '<!doctype html><html lang="en">'+head('Primer Tm · NGS 101',home='index.html')+'''
<link rel="stylesheet" href="assets/sequence_annotation.css"><link rel="stylesheet" href="assets/primer_tm.css">
<main class="wrap annotation-tool tm-tool"><h1>Primer Tm</h1>
<form id="tm-form" class="annotation-card">
<div class="tm-primer"><label for="tm-first">Primer 1 · 5′ → 3′</label><input id="tm-first" type="text" autocomplete="off" spellcheck="false" required placeholder="DNA sequence"><label>5′ tail to exclude (nt)<input id="tm-tail-first" type="number" min="0" step="1" value="0" required></label></div>
<div class="tm-primer"><label for="tm-second">Primer 2 · optional</label><input id="tm-second" type="text" autocomplete="off" spellcheck="false" placeholder="DNA sequence"><label>5′ tail to exclude (nt)<input id="tm-tail-second" type="number" min="0" step="1" value="0" required></label></div>
<div class="tm-actions"><button type="submit" class="primary-button">Calculate</button><button id="tm-example" type="button">Set example primers</button><button id="tm-clear" type="button">Clear</button></div>
<fieldset class="tm-buffer"><legend>KAPA HiFi HotStart ReadyMix · final 1× concentrations</legend>
<div class="tm-conditions">'''+concentration_fields()+'''</div>
<p class="tm-assumption">Monovalent salt is undisclosed. 50 mM is an editable model assumption, not Roche’s formulation.</p><button id="tm-reset" type="button">Reset KAPA conditions</button>
</fieldset></form>
<p id="tm-status" role="status" aria-live="polite"></p>
<section id="tm-results" class="annotation-card" hidden><h2>Estimated Tm</h2><div id="tm-values" class="tm-values"></div><p id="tm-difference"></p><p id="tm-conditions-summary" class="muted"></p><p id="tm-warnings" class="tm-assumption"></p><button id="tm-export" type="button">Export TSV</button><details><summary>Calculation details</summary><div id="tm-details"></div></details></section>
<details class="tm-method"><summary>Model and KAPA guidance</summary>
<p>SantaLucia (1998) DNA nearest-neighbour thermodynamics with terminal and self-complementarity corrections; Owczarzy (2008) monovalent/magnesium salt correction. Mg–dNTP binding uses a 1:1 equilibrium with an association constant of 30,000 M⁻¹. Each dNTP input is multiplied by four before this calculation.</p>
<p>The primer is assumed to be in excess over its template (template concentration neglected). Self-complementary strands use total strand concentration and the symmetry entropy correction. Perfectly matched, unmodified A/C/G/T DNA only, 2–200 binding bases. This model does not predict mismatches, folding, primer dimers, effects of additives, or polymerase activity.</p>
<p>Only the binding sequence after the specified 5′ tail enters the calculation. These are estimates under entered ionic conditions, not experimentally validated KAPA-specific predictions. Changing disclosed concentrations creates custom conditions.</p>
<p>Roche’s guidance for this ReadyMix is an initial annealing gradient of 60–72°C, or 65°C as a first approach when a gradient is impractical. This guidance is separate from calculated Tm; no fixed Tm-to-Ta offset is assumed. Fidelity, GC and Library Amplification formulations are not represented by this preset.</p>
<p><a href="'''+ROCHE+'''">Roche KR0370 v16.25 (May 2025), pp. 1–3</a> · <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC19045/">SantaLucia 1998</a> · <a href="https://doi.org/10.1021/bi702363u">Owczarzy 2008</a></p></details>
<noscript><p>This calculator requires JavaScript.</p></noscript></main>
<script type="application/json" id="tm-model">'''+json.dumps(model_data()).replace('<','\\u003c')+'''</script>
<script src="assets/sequence_annotation.js" defer></script><script src="assets/primer_tm.js" defer></script><script src="assets/primer_tm_ui.js" defer></script></html>
'''


def build(output=ROOT):
    (output/'tm.html').write_text(render_page(),encoding='utf-8')
    return output/'tm.html'


if __name__=='__main__':print(build().name)

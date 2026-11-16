#!/usr/bin/env python3
"""Meaningful failure modes for offline presentation-priority heuristics."""
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'lib')]
from checks import Check
from presentation_audit import measure,priority
from model_completeness import SIGNALS

c=Check()
def row(sequence='ACGTACGT'):
    return '<polygon class="chem-strand-box"/><text class="chem-strand-text"><tspan><title>primer 100 nt</title>'+sequence+'</tspan></text>'
def workflow(sequences=('ACGTACGT','ACGTACGT'),label='Extend'):
    return '<svg class="chem-svg chem-workflow">'+row(sequences[0])+''.join(
        '<text class="chem-reaction-label">'+label+str(i)+'</text>'+row(s) for i,s in enumerate(sequences[1:]))+'</svg>'
visible=workflow(('ACGTACGT',)*4)
hidden='<details><summary>Details</summary>'+visible+'</details>'
c('Closed workflow is not initially visible',measure(hidden)['visible_steps'],0)
c('Closed workflow stays available in separate hidden metrics',measure(hidden)['hidden_steps'],3)
c('Opening details makes its steps visible',measure(hidden.replace('<details>','<details open>'))['visible_steps'],3)
c('Open inner details do not override a closed outer details',measure('<details><details open>'+visible+'</details></details>')['visible_steps'],0)
c('Repeated display spines do not inflate unique coverage',measure(visible+visible)['visible_steps'],3)
c('SVG/CSS words do not masquerade as prose',measure('<style>word '*100+'</style>'+visible)['visible_words'],0)
c('Final-state diagrams alone are not reaction steps',measure('<svg>'+row()+'</svg>')['visible_steps'],0)
empty=measure(workflow(('[unassembled product]',)*4))
c('Placeholder text is not fixed sequence',empty['visible_sequence_steps'],0)
c('Concrete bases alongside a variable insert count',measure(workflow(('ACGTACGTNNNNNN','ACTGACTGNNNNNN')))['visible_sequence_steps'],1)
letters='<polygon class="chem-strand-box"/><text class="chem-strand-text">'+''.join('<tspan>'+b+'</tspan>'for b in 'ACGTACGT')+'</text>'
c('Nucleotide-resolved individual tspans still count',measure('<svg class="chem-workflow"><text class="chem-reaction-label">Copy</text>'+letters+'</svg>')['visible_sequence_steps'],1)
c('Descriptive text in tooltip cannot invent fixed sequence',measure(workflow(('[unknown]','[unknown]')).replace('primer 100 nt','ACGTACGTACGT'))['visible_sequence_steps'],0)
c('Hidden long sequences cannot improve visible sequence coverage',measure(hidden)['visible_sequence_steps'],0)
c('Inline-hidden diagrams are separate from visible ones',measure('<div style="display:none">'+visible+'</div>')['visible_steps'],0)
c('A visible molecular walkthrough ranks ahead of a hidden one',priority(measure(visible))['score']<priority(measure(hidden))['score'],True)
c('Many empty reaction arrows do not earn a low score',priority(empty)['score']>priority(measure(visible))['score'],True)
report={'version':1,'signals':{s:0 for s in SIGNALS}}
def metadata(value):return '<script type="application/json" data-model-completeness="1">'+json.dumps(value)+'</script>'
before=priority(measure(visible+metadata(report)))['score']
report['signals']['inferred_objects']=1000
c('Legitimate inference alone is not penalized',priority(measure(visible+metadata(report)))['score'],before)
c('Old completeness formats stay unassessed',measure(visible+metadata({'version':0}))['model_level'],'unassessed')
c('Scores never exceed their documented range',0<=priority(measure('<p>'+'words '*10000+'</p>'))['score']<=100,True)
# A missing deployment/build output is a separate repair queue, even when a
# good sibling page exists. Inspection must not execute a protocol to fill it.
sys.path[:0]=[str(ROOT),str(ROOT/'tools')]
import audit_protocol_presentation as audit
from tempfile import TemporaryDirectory
from unittest.mock import patch
fixture={'dir':'presentation_fixture','protocol':'Fixture','status':'documented'}
with TemporaryDirectory() as tmp:
    site=Path(tmp)
    with patch.object(audit.build_index,'diagram_out',return_value=ROOT/'presentation_fixture/page.html'):
        missing=audit.inspect(fixture,site=site)
        c('Missing site output is labelled unbuilt',missing['page_status'],'unbuilt')
        c('Missing output is explicitly a build queue',missing['score_components'],{'page_unavailable':100.0})
        target=site/'presentation_fixture/page.html';target.parent.mkdir();target.write_text(visible)
        built=audit.inspect(fixture,site=site)
        c('Explicit site root is actually inspected',built['visible_steps'],3)
        c('Input content is fingerprinted for later comparisons',len(built['html_sha256']),64)
c('Weak reaction labels are exposed for targeted follow-up',empty['visible_weak_step_labels'],['Extend0','Extend1','Extend2'])
c('An unavailable relevant original blocks the optional source filter',audit.review_flags({'status':'in_progress','source_checklist':[{'relevance':'relevant','availability':'unavailable'}]})['source_blocked'],True)
c('Ordinary unread originals remain a review queue, not automatic source blockers',audit.review_flags({'source_checklist':[{'relevance':'relevant','availability':'unread'}]})['source_blocked'],False)
c('An unavailable irrelevant original does not block available-source triage',audit.review_flags({'source_checklist':[{'relevance':'irrelevant','availability':'unavailable'}]})['source_blocked'],False)
wide='<svg class="chem-workflow" viewBox="0 0 2400 600"><polygon class="chem-strand-box"/><text class="chem-strand-text">A lengthy label 5′-<tspan><title>Tooltip must not count as label</title>ACGTACGT</tspan>-3′</text></svg>'
c('Native SVG width is exposed as an advisory metric',measure(wide)['visible_max_diagram_width'],2400.0)
c('Only the prefix before sequence contributes to the label metric',measure(wide)['visible_max_strand_prefix_columns'],len('A lengthy label 5′-'))
native_prefix=wide.replace('A lengthy label 5′-','<tspan>A lengthy label 5′-</tspan>')
c('Native Scene prefixes in a tspan remain measurable',measure(native_prefix)['visible_max_strand_prefix_columns'],len('A lengthy label 5′-'))
c('Hidden wide diagrams do not dominate the visible geometry metric',measure('<details>'+wide+'</details>')['visible_max_diagram_width'],0.0)
c('Long biological diagrams are not automatically penalized',priority(measure(wide))['score'],priority(measure(wide.replace('2400','400')))['score'])
c.report()

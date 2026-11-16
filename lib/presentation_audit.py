"""Offline triage of rendered molecular explanations, never quality certification.

Counts native diagram semantics rather than '<svg' strings or words in source
code. Closed details are separate from the initially visible page. No browser,
protocol imports, source downloads or chemistry execution are needed.
"""
from html.parser import HTMLParser
import json
import re

from model_completeness import level, SIGNALS

VOID={'area','base','br','col','embed','hr','img','input','link','meta','param','source','track','wbr'}
NONPROSE={'head','script','style','nav','header','footer','svg','pre','code','seq'}


class Presentation(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.stack=[];self.diagrams=[];self.words=[];self.reports=[];self.legacy_panels=0

    def handle_starttag(self,tag,attrs):
        attrs=dict(attrs);classes=set(attrs.get('class','').split())
        closed=(tag=='details' and 'open' not in attrs)
        hidden=any(f['hidden'] for f in self.stack) or closed or 'hidden' in attrs
        hidden=hidden or bool(re.search(r'(?:display\s*:\s*none|visibility\s*:\s*hidden)',attrs.get('style',''),re.I))
        f={'tag':tag,'hidden':hidden,'classes':classes,'text':[]}
        svg=next((f['svg'] for f in reversed(self.stack) if 'svg' in f),None)
        if tag=='svg':
            svg={'hidden':hidden,'workflow':'chem-workflow' in classes,
                 'groups':[{'label':None,'strands':0,'fixed_bases':0}],
                 'width':0.0,'label_columns':0}
            try:
                svg['width']=float(attrs.get('viewbox','').split()[2])
            except (ValueError,IndexError):
                try:svg['width']=float(attrs.get('width','0'))
                except ValueError:pass
            self.diagrams.append(svg);f['svg']=svg
        group=svg['groups'][-1] if svg else None
        if svg and tag=='text' and 'chem-reaction-label' in classes:
            group={'label':'','strands':0,'fixed_bases':0};svg['groups'].append(group)
            f['reaction']=group
        if group is not None and (('chem-strand-box' in classes and tag=='polygon')
                                   or ('fold-strand' in classes and tag=='g')):
            group['strands']+=1
        row=next((f['row'] for f in reversed(self.stack) if 'row' in f),None)
        if group is not None and 'chem-strand-text' in classes:
            row={'group':group,'svg':svg,'text':[]};f['row']=row
        if row is not None and tag=='tspan':
            f['bases']=row['group']
        if group is not None and classes & {'fold-base','chem-routed-base'}:f['fold_bases']=group
        if tag=='script' and 'data-model-completeness' in attrs:f['report']=True
        if tag=='pre' and not hidden:self.legacy_panels+=1
        if tag not in VOID:self.stack.append(f)

    def handle_startendtag(self,tag,attrs):
        self.handle_starttag(tag,attrs)
        if tag not in VOID:self.handle_endtag(tag)

    def handle_data(self,data):
        # Ignore tooltips when counting sequence characters and reaction labels.
        tooltip=any(f['tag']=='title' for f in self.stack)
        if not tooltip:
            row=next((f['row']for f in reversed(self.stack)if 'row'in f),None)
            if row is not None:row['text'].append(data)
            for frame in self.stack:
                if any(k in frame for k in ('reaction','bases','fold_bases','report')):
                    frame['text'].append(data)
        if (self.stack and not any(f['hidden'] or f['tag'] in NONPROSE for f in self.stack)):
            self.words.extend(re.findall(r"[A-Za-z][A-Za-z’'-]*",data))

    def handle_endtag(self,tag):
        index=next((i for i in range(len(self.stack)-1,-1,-1) if self.stack[i]['tag']==tag),None)
        if index is None:return
        closing=self.stack[index:];del self.stack[index:]
        for frame in reversed(closing):
            text=''.join(frame['text']).strip()
            if 'row' in frame:
                row=frame['row']
                # Scene puts its label/end marker in a first tspan, while
                # strand_row may use raw prefix text. Accept both layouts.
                marker=re.search(r"[35]['′][-−–]",''.join(row['text']))
                if marker is not None:
                    row['svg']['label_columns']=max(row['svg']['label_columns'],marker.end())
            if 'reaction' in frame:frame['reaction']['label']=' '.join(text.split())
            if 'bases' in frame:
                # A literal adapter alongside N/X insert placeholders still counts.
                # Descriptive role tokens and tooltip lengths never become bases.
                if re.fullmatch('[ACGTUNRYSWKMBDHVXacgtunryswkmbdhvx]+',text):
                    frame['bases']['fixed_bases']+=sum(c in 'ACGTUacgtu' for c in text)
            if 'fold_bases' in frame:
                frame['fold_bases']['fixed_bases']+=sum(c in 'ACGTUacgtu' for c in text)
            if 'report' in frame:
                try:self.reports.append(json.loads(text))
                except (ValueError,TypeError):self.reports.append(None)

    def metrics(self):
        result={'visible_words':len(self.words),'legacy_panels':self.legacy_panels}
        for hidden,prefix in [(False,'visible'),(True,'hidden')]:
            diagrams=[d for d in self.diagrams if d['hidden']==hidden]
            # A repeated display spine is not extra reaction coverage. For each
            # identical label retain its most informative drawn occurrence.
            steps={}
            for d in diagrams:
                if not d['workflow']:continue
                for g in d['groups']:
                    if g['label'] is None:continue
                    old=steps.get(g['label'])
                    if old is None or (g['fixed_bases'],g['strands'])>(old['fixed_bases'],old['strands']):steps[g['label']]=g
            result.update({prefix+'_diagrams':len(diagrams),prefix+'_steps':len(steps),
                prefix+'_molecular_steps':sum(g['strands']>0 for g in steps.values()),
                prefix+'_sequence_steps':sum(g['strands']>0 and g['fixed_bases']>=6 for g in steps.values()),
                prefix+'_strand_states':sum(g['strands']>0 for d in diagrams for g in d['groups']),
                prefix+'_sequence_states':sum(g['strands']>0 and g['fixed_bases']>=6 for d in diagrams for g in d['groups'])})
            # Advisory geometry counters: real long sequences legitimately need
            # scrolling, so these do not change the reaction/sequence priority.
            result[prefix+'_max_diagram_width']=max((d['width']for d in diagrams),default=0.0)
            result[prefix+'_max_strand_prefix_columns']=max((d['label_columns']for d in diagrams),default=0)
            result[prefix+'_weak_step_labels']=[label for label,g in steps.items()
                if not g['strands'] or g['fixed_bases']<6]
        result['model_level']='unassessed'
        result['model_signals']={k:0 for k in SIGNALS}
        levels=[level(r) for r in self.reports]
        for candidate in ('unassessed','unresolved','limited','inferred','declared'):
            if candidate in levels:result['model_level']=candidate;break
        for report in self.reports:
            if level(report)!='unassessed':
                for key in SIGNALS:result['model_signals'][key]+=report['signals'][key]
        return result


def measure(text):
    p=Presentation();p.feed(text);p.close();return p.metrics()


def priority(metrics):
    """A bounded 0–100 review-need heuristic with additive explanations.

    Raw model gaps can add only 12 points; invisible / non-molecular steps
    dominate. Unknown biological sequences and inference alone are not faults.
    """
    parts={};reasons=[]
    def add(key,points,reason):
        if points>0:parts[key]=round(points,1);reasons.append(reason)
    steps=metrics['visible_steps'];drawn=metrics['visible_molecular_steps'];seq=metrics['visible_sequence_steps']
    if not steps:
        add('no_visible_steps',30,'no initially visible native reaction steps')
    elif steps==1:add('single_step',12,'only one initially visible reaction step')
    if not metrics['visible_strand_states']:
        add('no_strand_states',25,'no initially visible native strand states')
    if steps:
        add('undrawn_steps',25*(1-drawn/steps),f'{steps-drawn}/{steps} visible steps have no native strand drawing')
        add('sequence_sparse_steps',15*(1-seq/steps),f'{steps-seq}/{steps} visible steps have no detected fixed sequence')
    elif not metrics['visible_sequence_states']:
        add('no_sequence_states',15,'no initially visible sequence-bearing native states')
    if metrics['hidden_steps'] and steps<3:
        add('collapsed_walkthrough',10,'reaction workflow is mainly inside closed details')
    density=metrics['visible_words']/max(seq,1)
    add('prose_density',min(10,max(0,(density-250)/75)),
        f'{metrics["visible_words"]} visible prose words for {seq} sequence-bearing reaction steps')
    signals=metrics['model_signals']
    add('unmodelled_steps',min(6,signals['unmodelled_steps']),
        str(signals['unmodelled_steps'])+' explicitly unmodelled graph steps')
    add('model_limits',min(6,signals['unrepresented_molecules']/5+signals['limited_geometry']/3+signals['limited_readouts']),
        'declared molecular / geometry / readout limits; inspect source availability')
    return {'score':round(min(100,sum(parts.values())),1),'score_components':parts,'reasons':reasons}

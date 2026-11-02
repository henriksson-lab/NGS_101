"""Ancestry coverage does not confuse missing links, blockers and inference."""
from copy import deepcopy
from pathlib import Path
import sys
sys.path[:0] = [str(Path(__file__).resolve().parents[1]/'lib')]
from protocol_ancestry import ancestry_audit, protocols_without_predecessors, validate_reviews
from protocol_lineage import build_payload, local_view


def run(check):
    rows=[dict(dir=p,protocol=p,status='documented',section='published',representation='published_instance',modality='DNA')for p in 'abcd']
    def edge(parent,child,basis):
        return dict(parent=parent,child=child,relation='adapted_component',basis=basis,scope='Specified component',evidence=[dict(url='https://example.org/paper',locator='Methods attribution')],rationale='Located source or scoped design comparison',limits='No chemistry inheritance',reviewed='2026-11-02')
    registry={'schema_version':1,'relations':[edge('a','b','documented'),edge('b','c','inferred')]}
    def review(formal='cannot_formally_assign',inference='no_defensible_candidate'):
        return dict(status='reviewed',formal_assignment=formal,reason='No formal attribution established in the bounded source scope',sources=[dict(url='https://example.org/paper',locator='Original methods',identity='Identified original',availability='read_original')],external_predecessors=[],inference_status=inference,inference_reason='Compared earlier designs within scope',blockers=[],limits='Not proof of independent origin',reviewed='2026-11-02')
    reviews={'schema_version':1,'reviews':{'b':review('internal','not_needed'),'c':review(inference='supported'),'d':review()}}
    external=dict(title='Older method outside catalogue',url='https://example.org/old',relation='adapted_component',basis='documented',scope='One component',evidence=[dict(url='https://example.org/paper',locator='Methods explicitly cite predecessor')],rationale='Named component adaptation',limits='Not a registered molecular producer',reviewed='2026-11-02')
    reviews['reviews']['d'].update(formal_assignment='external',inference_status='not_needed',external_predecessors=[external])
    audit=ancestry_audit(rows,registry,reviews,visible_ids=set('bcd'))
    by={r['id']:r for r in audit['protocols']}
    check('a root with descendants lacks a parent but is not isolated',(by['a']['without_internal_parent'],by['a']['unlinked']),(True,False))
    check('inferred parent does not become formal ancestry',(by['c']['without_internal_parent'],by['c']['without_formal_parent']),(False,True))
    check('a hidden registered parent is still counted',(by['b']['without_internal_parent'],by['b']['without_visible_parent']),(False,True))
    check('external predecessor is distinguishable from no ancestry',(by['d']['without_internal_parent'],by['d']['without_formal_parent'],by['d']['formal_assignment']),(True,False,'external'))
    check('external source-linked ancestry is not visually isolated',(by['d']['unlinked'],by['d']['without_visible_parent']),(False,False))
    check('unrecorded review remains unreviewed',by['a']['review_status'],'unreviewed')
    check('without predecessor function excludes inferred parents',[r['id']for r in protocols_without_predecessors(rows,registry,reviews)],['a','d'])
    check('documented-only query retains inference-only protocols',[r['id']for r in protocols_without_predecessors(rows,registry,reviews,include_inferred=False)],['a','c','d'])
    bad=deepcopy(reviews);bad['reviews']['a']=review();bad['reviews']['a'].update(status='source_blocked',inference_status='blocked')
    check.raises('source blocker needs actual unavailable original target',lambda:validate_reviews(bad,set('abcd'),registry['relations']))
    bad['reviews']['a']['blockers']=[dict(url='https://example.org/missing',reason='Publisher source unavailable')]
    check('explicit source blocker is a retained outcome',ancestry_audit(rows,registry,bad)['protocols'][0]['review_status'],'source_blocked')
    comparison=deepcopy(reviews);comparison['reviews']['d']['blockers']=[dict(url='https://example.org/older-comparison',reason='An additional earlier comparator original is unavailable')]
    check('a completed primary-source pass can retain comparison blockers',next(x for x in ancestry_audit(rows,registry,comparison)['protocols']if x['id']=='d')['has_source_blocker'],True)
    bad=deepcopy(reviews);bad['reviews']['c']['formal_assignment']='internal'
    check.raises('inferred link cannot satisfy formal internal assignment',lambda:validate_reviews(bad,set('abcd'),registry['relations']))
    bad=deepcopy(reviews);bad['reviews']['a']=review(inference='supported')
    check.raises('supported likely ancestor requires a real inferred edge',lambda:validate_reviews(bad,set('abcd'),registry['relations']))
    eligible=[dict(id=p,url=p+'/page.html')for p in 'abcd']
    payload=build_payload(rows,eligible,registry,{'components':[]},reviews)
    externals=[node for node in payload['nodes']if node.get('external')]
    check('external node never becomes an in-catalogue protocol',(len(externals),payload['stats']['protocols']),(1,4))
    check('external original URL is retained without a relative prefix','href="https://example.org/old"'in local_view('d',payload),True)
    check('external ancestry is not used as sequence evidence',payload['similarities'],[])
    check('inability to assign formally remains visible beside inferred link','No formal ancestor assigned'in local_view('c',payload),True)
    filtered=build_payload(rows,[e for e in eligible if e['id']!='d'],registry,{'components':[]},reviews)
    check('an external predecessor is exposed only for an eligible child',any(n.get('external')for n in filtered['nodes']),False)
    withheld=build_payload(rows,[dict(id='b',url='b/page.html')],registry,{'components':[]},reviews)
    check('hidden curated parent is explained without publishing its private node','links are withheld'in local_view('b',withheld)and'No curated ancestry links yet'not in local_view('b',withheld),True)
    guessed=deepcopy(registry);guess=edge('a','d','guessed');guess['chronology']={'parent_year':2010,'child_year':2020,'locator':'Inspected article dates'};guessed['relations'].append(guess)
    guessed_reviews=deepcopy(reviews);guessed_reviews['reviews']['d']=review(inference='guessed_candidate')
    guessed_report=ancestry_audit(rows,guessed,guessed_reviews)
    d=next(x for x in guessed_report['protocols']if x['id']=='d')
    check('a loose guessed analogue does not establish an internal or formal predecessor',(d['without_internal_parent'],d['without_formal_parent']),(True,True))
    check('guessed candidates are excluded from predecessor coverage by default',[r['id']for r in protocols_without_predecessors(rows,guessed,guessed_reviews)],['a','d'])
    check('caller may explicitly count a guessed candidate',[r['id']for r in protocols_without_predecessors(rows,guessed,guessed_reviews,include_guessed=True)],['a'])
    bad_guess=deepcopy(guessed);bad_guess['relations'][-1]['chronology']['parent_year']=2025
    check.raises('a later paper cannot be guessed as a previous method',lambda:ancestry_audit(rows,bad_guess,guessed_reviews))
    bad_guess=deepcopy(guessed);bad_guess['relations'][-1]['chronology']['parent_year']=2020
    check.raises('equal-year guesses need exact date ordering',lambda:ancestry_audit(rows,bad_guess,guessed_reviews))
    guessed['relations'][-1]['chronology'].update(parent_year=2020,parent_date='2020-02-01',child_date='2020-03-01')
    check('a verified earlier date can order same-year chemistry analogues',ancestry_audit(rows,guessed,guessed_reviews)['stats']['protocols'],4)
    bad=deepcopy(reviews);bad['reviews']['d']['external_predecessors'][0]['url']='javascript:alert(1)'
    check.raises('unsafe external source destination rejected',lambda:build_payload(rows,eligible,registry,{'components':[]},bad))

    from protocol_ancestry_similarity import rank_candidates, publication_years, cited_protocols
    candidates=[dict(dir=p,protocol=p,doi='fixture/'+p,section='published',representation='published_instance')for p in ['old','recent','same','future','unknown','example']]
    candidates[-1]['representation']='own_design'
    props={p:{'adapter_installation':('splint ligation',),'topology':('circular',),'amplification':('rolling-circle amplification',)}for p in ['old','recent','same','future','unknown','example']}
    # Another chemically distinct method makes shared rare features informative.
    candidates.append(dict(dir='other',protocol='other',doi='fixture/other',section='published',representation='published_instance'))
    props['other']={'fragmentation':('Tn5 / tagmentation',)}
    years={'old':2010,'recent':2020,'same':2020,'future':2025,'example':2005,'other':2000}
    suggestions=rank_candidates('recent',candidates,props,years)
    check('suggestions exclude future, same-year, unknown-date and invented methods',[r['protocol']for r in suggestions['candidates']],['old'])
    check('unknown sequence channel is an uncertainty interval, not known agreement',(suggestions['candidates'][0]['sequence_score'],suggestions['candidates'][0]['score'],suggestions['candidates'][0]['score_upper_bound']),(None,.75,1.0))
    check('a missing target date is reported rather than inferred',rank_candidates('unknown',candidates,props,years)['status'],'unknown_target_publication_year')
    check('discovery output contains no automatic ancestry edge','basis'in suggestions['candidates'][0],False)
    metadata=[dict(our_dir='recent',doi='fixture/recent',year='2020'),dict(our_dir='recent',doi='unrelated/old',year='2010')]
    check('grouped source bibliography does not supply another papers year',publication_years(candidates,metadata),{'recent':2020})
    metadata.append(dict(our_dir='recent',doi='fixture/recent',year='2021'))
    check('conflicting defining-year metadata stays unknown',publication_years(candidates,metadata),{})
    citation_rows=[dict(dir='target',protocol='target',doi='10.1/target'),dict(dir='old',protocol='old',doi='10.1/method'),dict(dir='other',protocol='other',doi='10.1/other')]
    matches=cited_protocols('target',citation_rows,'10.1/target\nReference DOI: 10.1/METHOD.\nNot the same DOI: 10.1/other-version')
    check('citation helper preserves DOI identity and source line pointers',[(x['protocol'],x['lines'])for x in matches],[('old',[2])])
    check('citation helper does not promote a citation into ancestry','basis'in matches[0],False)
    candidates.append(dict(dir='illustrated_method',protocol='An original published method',doi='fixture/method',section='published',representation='illustrative_example'))
    years['illustrated_method']=2005;props['illustrated_method']=props['old']
    seq='ACGTGCTACAGTCGATGCTA'
    toy_components={'components':[{'sequence':seq,'occurrences':[{'protocol':'recent','inferred':False,'name':'fixed'},{'protocol':'illustrated_method','inferred':False,'name':'invented example'}]}]}
    method=next(x for x in rank_candidates('recent',candidates,props,years,toy_components)['candidates']if x['protocol']=='illustrated_method')
    check('an illustrated page may identify its original method, never its invented sequence',(method['method_concept_only'],method['coverage']['sequence_components'],method['sequence_score']),(True,0,None))


if __name__=='__main__':
    from checks import Check
    check=Check();run(check);check.report()

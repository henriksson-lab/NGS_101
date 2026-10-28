"""Offline relationship registry, public eligibility and sequence-similarity checks."""
from copy import deepcopy
import json
from pathlib import Path
import sys
sys.path[:0]=[str(Path(__file__).resolve().parents[1]),str(Path(__file__).resolve().parents[1]/'lib')]
from protocol_lineage import (build_payload,enrich_page,local_view,sequence_signatures,
                              sequence_similarity,validate_registry)
ROOT=Path(__file__).resolve().parents[1]


def run(check):
    check.section('protocol relationships, provenance and independent similarity')
    def edge(a,b,basis='documented'):
        return {'parent':a,'child':b,'relation':'derived_from','basis':basis,'scope':'Specified component',
                'evidence':[{'url':'https://example.org/original','locator':'Methods, specified passage'}],
                'rationale':'Source relation or explicitly qualified design inference.',
                'limits':'No chemistry inheritance.','reviewed':'2026-10-28'}
    def registry(edges):return {'schema_version':1,'relations':edges}
    def rejects(data,known):
        try:validate_registry(data,known)
        except ValueError:return True
        return False
    family=registry([edge('a','c'),edge('b','c','inferred'),edge('c','d')])
    check('one protocol can have multiple separately qualified parents',
          [(r['parent'],r['basis']) for r in validate_registry(family,set('abcd')) if r['child']=='c'],
          [('a','documented'),('b','inferred')])
    check('ancestry rejects self links',rejects(registry([edge('a','a')]),{'a'}),True)
    check('ancestry rejects unknown protocol IDs',rejects(registry([edge('a','z')]),set('abc')),True)
    check('ancestry rejects duplicate parent/child entries',rejects(registry([edge('a','b'),edge('a','b','inferred')]),set('ab')),True)
    check('ancestry rejects directed cycles',rejects(registry([edge('a','b'),edge('b','c'),edge('c','a')]),set('abc')),True)
    bad=edge('a','b','inferred');bad['evidence']=[]
    check('inferred links still require an identified evidence basis',rejects(registry([bad]),set('ab')),True)
    bad=edge('a','b');bad['limits']=''
    check('scope limits cannot be silently omitted',rejects(registry([bad]),set('ab')),True)
    chain=registry([edge(str(i),str(i+1)) for i in range(1000)])
    check('future thousand-protocol chains do not depend on recursion depth',
          len(validate_registry(chain,{str(i) for i in range(1001)})),1000)
    def comp(seq,members,inferred=False):
        return {'sequence':seq,'occurrences':[{'protocol':p,'inferred':inferred,'name':'Part '+seq[:5]} for p in members]}
    common='ACGTCAGTGCTACGATCGTA'
    rare='GACTGCTACAGTCGATGCTA'
    second='TCGACATGGTACCTAGCTGA'
    unique='CTAGTCGATGACCTGAGTCA'
    reverse=rare.translate(str.maketrans('ACGT','TGCA'))[::-1]
    components={'components':[comp(common,'abc'),comp(rare,'ab'),comp(reverse,'ab'),
                              comp(second,'ab'),comp(unique,'c'),comp('A'*20,'ab'),
                              comp('ACGTNNNNACGTACGTACGT','ab'),comp('ACGTACGT','ab'),
                              comp('GCTAGCATGACTCGATCGTA','ab',True)]}
    signatures,metadata=sequence_signatures(components,set('abc'))
    check('exact reverse complements deduplicate instead of doubling similarity evidence',
          (len(signatures['a']),len(metadata)),(3,4))
    scores,metric=sequence_similarity(components,set('abc'),minimum=.01)
    check('ubiquitous platform components do not manufacture similarity neighbours',
          [(e['first'],e['second'],e['score'],e['shared_count']) for e in scores],[('a','b',1.0,2)])
    check('sequences unique to one protocol remain denominator evidence',len(signatures['c']),2)
    check('unknown sequence coverage is not asserted as a zero score',
          sequence_similarity({'components':[]},set('abc'))[0],[])
    rows=[dict(dir=p,protocol=p,status='documented',section='published',representation='published_instance',modality='RNA') for p in 'abcd']
    rows.extend([dict(dir='private',protocol='private',status='documented',section='wip',representation='own_design',modality='DNA'),
                 dict(dir='example',protocol='example',status='construct',section='published',representation='illustrative_example',modality='DNA'),
                 dict(dir='notes',protocol='notes',status='notes',section='published',representation='published_instance',modality='RNA')])
    eligible=[{'id':r['dir'],'url':r['dir']+'/diagram.html'} for r in rows]
    payload=build_payload(rows,eligible,family,components)
    check('private, notes-only and illustrative references stay out of the public graph',
          [n['id'] for n in payload['nodes']],list('abcd'))
    check('similarities never become curated ancestry',
          len(payload['relations']),3)
    check('unavailable public ancestors are not exposed by relation endpoints',
          len(build_payload(rows,[e for e in eligible if e['id']!='a'],family,components)['relations']),2)
    local=local_view('c',payload)
    check('local graph includes both parents, descendants and evidence legend',
          all(s in local for s in ['a/diagram.html','b/diagram.html','d/diagram.html','Solid: source-documented','dashed: inferred']),True)
    page='<html><h1>Protocol</h1><p>Body</p></html>'
    enriched=enrich_page(page,'c',payload)
    check('assembly enrichment adds exactly one local disclosure before body',
          enriched.count('class="protocol-lineage-local"'),1)
    try:enrich_page(enriched,'c',payload);duplicate=False
    except ValueError:duplicate=True
    check('duplicate page enrichment is rejected',duplicate,True)
    check('uncurated is explicit rather than a claim of independent origin',
          'does not establish independent origins' in local_view('private',{'nodes':[{'id':'private','title':'<unsafe>','url':'private/diagram.html'}],'relations':[]}),True)
    import csv
    known={r['dir'] for r in csv.DictReader((ROOT/'catalogue/ours.tsv').open(),delimiter='\t')}
    checked=validate_registry(json.loads((ROOT/'catalogue/protocol_lineage.json').read_text()),known)
    check('tracked registry validates against current protocol identities',len(checked)>=8,True)
    from build_lineage import render_page
    check('graph tab has independent evidence/similarity controls and shared navigation',
          all(s in render_page() for s in ['lineage-documented','lineage-inferred','lineage-similar','href="lineage.html"']),True)


if __name__=='__main__':
    from checks import Check
    c=Check();run(c);c.report()

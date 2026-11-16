#!/usr/bin/env python3
"""Rank protocol pages needing a step-by-step molecular presentation review.

Read generated HTML, without executing protocol modules or fetching sources.
Run build_site.py first for a complete current snapshot. Scores prioritize
inspection, never certify scientific accuracy, completeness or source absence.
"""
import argparse
import csv
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'lib')]
import build_index
from presentation_audit import measure,priority

FINISHED={'documented','construct','assay'}
COLUMNS=('score','protocol','dir','catalogue_status','page_status','page','model_level',
    'visible_steps','visible_molecular_steps','visible_sequence_steps','hidden_steps',
    'visible_strand_states','visible_sequence_states','visible_words','legacy_panels',
    'visible_max_diagram_width','visible_max_strand_prefix_columns',
    'source_review_status','unread_relevant_sources','unavailable_relevant_sources','historical_disposition','reasons')


def review_flags(data):
    relevant=[s for s in data.get('source_checklist',[]) if s.get('relevance')=='relevant']
    unavailable=sum(s.get('availability')=='unavailable' for s in relevant)
    status=data.get('status','unreviewed')
    return {'source_review_status':status,
            'unread_relevant_sources':sum(s.get('availability')!='read' for s in relevant),
            'unavailable_relevant_sources':unavailable,
            'source_blocked':status=='blocked' or unavailable>0}


def inspect(record,*,site=None):
    name=record['dir']
    row={'dir':name,'protocol':record['protocol'],'catalogue_status':record['status']}
    review=ROOT/'catalogue/source_reviews'/ (name+'.json')
    row.update(review_flags(json.loads(review.read_text()) if review.exists() else {}))
    try:page=build_index.diagram_out(name)
    except (ValueError,SyntaxError) as error:
        page=None;row['discovery_error']=str(error)
    path=(site/page.relative_to(ROOT)) if site and page else page
    row['page']=str(page.relative_to(ROOT)) if page else ''
    row['page_status']='built' if path and path.is_file() else 'unbuilt' if page else 'no_builder'
    text=path.read_text() if row['page_status']=='built' else ''
    row['html_sha256']=hashlib.sha256(text.encode()).hexdigest() if text else ''
    row.update(measure(text))
    row.update(priority(row))
    if row['page_status']!='built':
        row.update(score=100.0,score_components={'page_unavailable':100.0},
                   reasons=['generated schematic unavailable; build or inspect the entry before grading presentation'])
    return row


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('protocols',nargs='*',help='Exact catalogue directories; default: all published finished schematics')
    p.add_argument('--include-unfinished',action='store_true',help='Also inventory notes/draft entries; missing pages are a separate build queue')
    p.add_argument('--exclude-source-blocked',action='store_true',help='Exclude blocked reviews or explicitly unavailable relevant originals; ordinary unread sources stay reported')
    p.add_argument('--built-only',action='store_true',help='Exclude missing pages from the presentation ranking')
    p.add_argument('--site',type=Path,help='Read pages from an assembled site directory instead of sibling build outputs')
    p.add_argument('--format',choices=['table','tsv','json'],default='table')
    p.add_argument('--limit',type=int,default=None,help='Default: 30 rows for table, all rows for TSV/JSON; 0 means all')
    p.add_argument('--output',type=Path,help='Write report to this path; no files changed by default')
    a=p.parse_args(argv)
    if a.limit is not None and a.limit<0:p.error('--limit must be nonnegative')
    if a.site is not None and not a.site.is_dir():p.error('--site must be an existing assembled-site directory')
    records=list(build_index.cat.ours());known={r['dir']:r for r in records}
    unknown=set(a.protocols)-known.keys()
    if unknown:p.error('Unknown protocol directories: '+', '.join(sorted(unknown)))
    selected=[known[n] for n in a.protocols] if a.protocols else [r for r in records if r['section']=='published' and (a.include_unfinished or r['status'] in FINISHED)]
    historical={}
    with (ROOT/'catalogue/molecular_state_audit.tsv').open() as handle:
        historical={r['slug']:r for r in csv.DictReader(handle,delimiter='\t')}
    rows=[]
    for record in selected:
        row=inspect(record,site=a.site)
        row['historical_disposition']=historical.get(row['dir'],{}).get('disposition','')
        row['historical_finding']=historical.get(row['dir'],{}).get('finding','')
        if a.built_only and row['page_status']!='built':continue
        if a.exclude_source_blocked and row['source_blocked']:continue
        rows.append(row)
    rows.sort(key=lambda r:(-r['score'],r['protocol'].casefold(),r['dir']))
    total=len(rows);limit=(30 if a.format=='table' else 0) if a.limit is None else a.limit
    if limit:rows=rows[:limit]
    import io
    out=io.StringIO()
    if a.format=='json':
        json.dump({'schema_version':1,'generated_at':datetime.now(timezone.utc).isoformat(),
                   'scope':'Heuristic rendered-page review queue, not scientific validation or a human visual audit.',
                   'input':'site:'+str(a.site) if a.site else 'sibling generated protocol HTML',
                   'matching_protocols':total,'rows':rows},out,indent=2);out.write('\n')
    elif a.format=='tsv':
        writer=csv.DictWriter(out,fieldnames=COLUMNS,delimiter='\t',extrasaction='ignore',lineterminator='\n');writer.writeheader()
        for row in rows:writer.writerow({**row,'reasons':'; '.join(row['reasons'])})
    else:
        out.write('Priority = review need (0–100), not a quality verdict. Rebuild pages before comparing changes.\n')
        out.write(f'{len(rows)} of {total} entries; larger score = earlier inspection.\n\n')
        for r in rows:
            out.write(f'{r["score"]:5.1f}  {r["protocol"]}  [{r["page_status"]}; model {r["model_level"]}]\n')
            out.write(f'       steps {r["visible_steps"]} visible / {r["hidden_steps"]} hidden; '
                      f'{r["visible_sequence_steps"]} sequence-bearing; {r["visible_words"]} prose words\n')
            out.write('       '+('; '.join(r['reasons']) or 'No weighted flags; still requires visual/source review.')+'\n')
    if a.output:
        a.output.write_text(out.getvalue());print(f'Wrote {len(rows)} of {total} entries to {a.output}')
    else:print(out.getvalue(),end='')


if __name__=='__main__':main()

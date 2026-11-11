#!/usr/bin/env python3
"""Audit kit-based index extraction and declared source primer sequences in Node.

The reference trials use kit binding arms without protocol or library-feature
annotations. Source primer trials instantiate unknown bases only as synthetic
queries. Coverage is recognition of binding geometry, not evidence of experimental
barcode use. Optional generated JSON reports must remain outside the checkout.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]

BROWSER = r"""
const fs=require('fs'),api=require('./assets/samplesheet_batch.js');
const matcher=require('./assets/sequence_annotation.js');
const catalogue=JSON.parse(fs.readFileSync(0,'utf8'));
const geometries=catalogue.index_geometries;
if (!Array.isArray(geometries) || !geometries.length)
  throw new Error('Rebuild the catalogue: index_geometries is missing or empty.');
if (!Array.isArray(catalogue.primer_queries))
  throw new Error('Rebuild the catalogue: primer_queries is missing.');
const failures=[],inventory=[],counts={kit:0,source:0};
function check(query,expected,blocked,metadata,data=catalogue) {
  counts[metadata.kind]++;
  const result=api.analyze([{sampleId:'Audit_sample',primer1:query,primer2:metadata.companion||''}],data);
  const row=result.results[0],actual=result.valid.length?row.values:null;
  const ok=blocked?!actual:actual&&['i7','i5'].every(role=>(actual[role]||'')===(expected[role]||''));
  if(!ok)failures.push({...metadata,query,expected,blocked,actual,error:row.error});
}
function valuesFor(result) {
  if(result.truncated || !result.candidates.length)return null;
  const values={};
  for(const candidate of result.candidates) {
    if(candidate.status!=='ready' || !candidate.forwardSequence)return null;
    if(values[candidate.role]&&values[candidate.role]!==candidate.forwardSequence)return null;
    values[candidate.role]=candidate.forwardSequence;
  }
  return values;
}
const referenceSeed='ACGTTGCAAAGCCATG';
for(const geometry of geometries) for(const length of [6,8,10,12,16]) {
  for(const mode of ['source','substituted']) {
    const raw=referenceSeed.slice(0,length).replace(/./g,base=>mode==='source'?base:{A:'C',C:'G',G:'T',T:'A'}[base]);
    const expected={[geometry.role]:geometry.to_forward==='forward'?raw:matcher.revcomp(raw)};
    const query=geometry.outer+raw+geometry.inner;
    const ambiguous=geometry.outer+'N'+raw.slice(1)+geometry.inner;
    for(const reverse of [false,true]) for(const scope of ['family','all_kits']) {
      const data=scope==='family'?{index_geometries:[geometry],nonindexed_primers:catalogue.nonindexed_primers||[]}:catalogue;
      const metadata={kind:'kit',family:geometry.id,length,mode,reverse,scope};
      check(reverse?matcher.revcomp(query):query,expected,false,metadata,data);
      check(reverse?matcher.revcomp(ambiguous):ambiguous,{},true,{...metadata,ambiguous:true},data);
    }
  }
}
const supportedSequences=new Set(),supportedProtocols=new Set(),queryCache=new Map(),universalQueries=new Set();
for(const source of catalogue.primer_queries||[]) {
  const item={protocol:source.protocol,label:source.label,locator:source.locator,sequence:source.sequence};
  if(source.inferred) {inventory.push({...item,status:'inferred'});continue;}
  const query=source.sequence.replace(/N/g,'A');
  const result=matcher.indexGeometryCandidates(query,geometries),expected=valuesFor(result);
  if(!expected) {
    let status;
    if(result.candidates.length || result.truncated)status='ambiguous_geometry';
    else if(api.detector(catalogue)(query).nonindexed)status='known_universal';
    else if(geometries.some(g=>query.includes(g.outer)||matcher.revcomp(query).includes(g.outer)))status='canonical_outer_without_supported_inner';
    else status='no_kit_binding';
    inventory.push({...item,status});
    if(status==='known_universal' && !universalQueries.has(query)) {
      universalQueries.add(query);
      const g=geometries.find(g=>g.role==='i7'),raw='ACGTTGCA';
      const companion=g.outer+raw+g.inner,values={i7:matcher.revcomp(raw)};
      for(const reverse of [false,true]) {
        const input=reverse?matcher.revcomp(query):query;
        const metadata={kind:'source',protocol:source.protocol,label:source.label,mode:'universal',reverse};
        check(input,{},true,metadata);
        check(input,values,false,{...metadata,companion});
      }
    }
    continue;
  }
  supportedSequences.add(source.sequence);supportedProtocols.add(source.protocol);
  inventory.push({...item,status:'supported'});
  // The source inventory retains duplicate occurrences, but identical DNA needs
  // only one set of engine trials. No protocol annotation enters extraction.
  if(queryCache.has(query))continue;
  queryCache.set(query,true);
  const positions=new Set();
  for(const candidate of result.candidates)for(let i=candidate.start;i<candidate.end;i++)positions.add(i);
  let variant=Array.from(query);
  for(const i of positions)variant[i]={A:'C',C:'G',G:'T',T:'A'}[variant[i]];
  variant=variant.join('');
  for(const [mode,input] of [['source',query],['substituted',variant]]) {
    const direct=matcher.indexGeometryCandidates(input,geometries),values=valuesFor(direct);
    for(const reverse of [false,true])check(reverse?matcher.revcomp(input):input,values||{},!values,
      {kind:'source',protocol:source.protocol,label:source.label,mode,reverse});
    const ambiguous=Array.from(input);ambiguous[result.candidates[0].start]='N';
    for(const reverse of [false,true])check(reverse?matcher.revcomp(ambiguous.join('')):ambiguous.join(''),{},true,
      {kind:'source',protocol:source.protocol,label:source.label,mode,reverse,ambiguous:true});
  }
}
process.stdout.write(JSON.stringify({counts,failures,inventory,
  supported_source_sequences:supportedSequences.size,supported_protocols:supportedProtocols.size,
  distinct_instantiated_source_sequences:queryCache.size}));
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalogue', type=Path, default=ROOT / 'assets/sequence_components.json')
    parser.add_argument('--report', type=Path, help='Optional JSON audit report outside the checkout')
    args = parser.parse_args()
    catalogue = json.loads(args.catalogue.read_text())
    result = subprocess.run(['node', '-e', BROWSER], cwd=ROOT, text=True,
                            input=json.dumps(catalogue), capture_output=True)
    if result.returncode:
        print(result.stderr)
        return result.returncode
    report = json.loads(result.stdout)
    report['catalogue'] = str(args.catalogue)
    counts = Counter(item['status'] for item in report['inventory'])
    print('Source occurrence inventory: ' + ', '.join(f'{key}={value}' for key, value in sorted(counts.items())))
    print(f"Supported source sequences: {report['supported_source_sequences']} across {report['supported_protocols']} protocols")
    print(f"Kit trials: {report['counts']['kit']}; source trials: {report['counts']['source']}; failures: {len(report['failures'])}")
    for failure in report['failures'][:30]:
        print('FAIL ' + failure['kind'] + ' / ' + failure.get('family', failure.get('protocol', '')) +
              ' / ' + failure['mode'] + (' / reverse' if failure['reverse'] else '') +
              ': ' + (failure.get('error') or f"got {failure['actual']}, expected {failure.get('expected')}"))
    if len(report['failures']) > 30:
        print('Additional failures omitted; use --report for all records.')
    if args.report:
        target = args.report.resolve()
        if target.is_relative_to(ROOT):
            parser.error('--report must be outside the checkout (generated analysis is not tracked)')
        target.write_text(json.dumps(report, indent=2) + '\n')
    return bool(report['failures'])


if __name__ == '__main__':
    raise SystemExit(main())

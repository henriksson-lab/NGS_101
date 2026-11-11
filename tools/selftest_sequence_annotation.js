#!/usr/bin/env node
'use strict';
const assert = require('node:assert/strict');
const SA = require('../assets/sequence_annotation.js');
let count = 0;
function check(label, fn) { fn(); count++; console.log('PASS ' + label); }
const seq = 'ACGTCAGTTCGAGATC';
function component(id,sequence,protocols,inferred=false) {
  return {id,sequence,aliases:[id],protocols,occurrences:protocols.map(protocol => ({protocol,name:id,source:'fixture',inferred,orientation:'+'}))};
}
const protocols = Array.from({length:30},(_,i) => ({id:'p'+i,title:'Protocol '+i,url:'p'+i+'.html'}));
const data = {schema_version:1,protocols,components:[component('distinct',seq,['p1']),component('nested',seq.slice(2,14),['p1']),component('common','AAGTCGTACGAT',protocols.map(p=>p.id)),component('inferred','GGTACAGTCGAT',['p2'],true)],templates:[]};
check('single FASTA and explicit RNA conversion', () => {
  assert.equal(SA.normalizeInput(' >x\n acgu \nRYn'), 'ACGURYN');
  assert.equal(SA.normalizeInput('acgu',{convertUToT:true}), 'ACGT');
});
check('reject multi FASTA, punctuation, empty and excessive length', () => {
  for (const value of ['>x\nACTG\n>y\nA','AC-TG','  ','A'.repeat(100001)]) assert.throws(() => SA.normalizeInput(value));
});
check('both orientations and ambiguity as a barrier', () => {
  const result = SA.annotate(seq+'N'+SA.revcomp(seq),data);
  const hits = result.hits.filter(h=>h.component==='distinct');
  assert.equal(hits.length,2);
  assert.deepEqual(hits.map(h=>[h.start,h.end,h.strand]),[[0,16,'+'],[17,33,'-']]);
  assert.equal(SA.annotate(seq.slice(0,8)+'N'+seq.slice(9),data).hits.length,0);
});
check('overlapping hits retained, long component in overview', () => {
  const result = SA.annotate(seq,data);
  assert.equal(result.hits.length,2);
  assert.deepEqual(result.overview.map(h=>h.component),['distinct']);
  assert.deepEqual(result.protocols[0].support,['distinct']);
});
check('palindrome appears once with both orientations', () => {
  const palindrome='ACGTACGTACGT';
  const result=SA.annotate(palindrome,{protocols,components:[component('pal',palindrome,['p1'])]});
  assert.equal(result.hits.length,1); assert.equal(result.hits[0].strand,'both');
});
check('inferred occurrences opt in; repeated evidence counts once', () => {
  assert.equal(SA.annotate('GGTACAGTCGAT',data).hits.length,0);
  assert.equal(SA.annotate('GGTACAGTCGAT',data,{includeInferred:true}).hits.length,1);
  assert.equal(SA.annotate(seq+seq,data).protocols[0].score,SA.annotate(seq,data).protocols[0].score);
});
check('ubiquitous parts do not provide informative ranking; caps deterministic', () => {
  const common=SA.annotate('AAGTCGTACGAT',data);
  assert.equal(common.protocols.length,20);
  assert.ok(common.protocols.every(p=>!p.informative && p.score===0));
  assert.equal(SA.annotate('AAGTCGTACGAT'+seq,data).protocols[0].id,'p1');
  const result=SA.annotate(seq.repeat(5),data,{maxHits:3});
  assert.equal(result.hits.length,3); assert.equal(result.truncated,true);
});
const left='AGCTACGTCAGTTCGA', right='CGATGCTACAGTGCAT';
const template={id:'t1',protocol:'p1',label:'adapter example',source:'fixture',sequence:left+'N'.repeat(8)+right,index_spans:[{start:left.length,end:left.length+8,feature_id:'index7',role:'i7',to_forward:'reverse_complement',basis:'explicit fixture orientation'}]};
const indexed={protocols,components:[],templates:[template]};
const query=left+'AACCGTTA'+right;
check('index windows require literal flanks and declared orientation', () => {
  const hit=SA.annotate(query,indexed).indexCandidates[0];
  assert.equal(hit.sequence,'AACCGTTA'); assert.equal(hit.forwardSequence,'TAACGGTT');
  assert.equal(hit.role,'i7'); assert.equal(hit.status,'ready');
  assert.deepEqual([hit.start,hit.end],[16,24]);
  assert.equal(SA.annotate('T'+query.slice(1),indexed).indexCandidates.length,0);
});
check('reverse query preserves index orientation and query coordinates', () => {
  const hit=SA.annotate('GG'+SA.revcomp(query),indexed).indexCandidates[0];
  assert.equal(hit.sequence,'AACCGTTA'); assert.equal(hit.forwardSequence,'TAACGGTT');
  assert.equal(hit.strand,'-'); assert.deepEqual([hit.start,hit.end],[18,26]);
});
check('ambiguous indexes and unproven orientation are unresolved', () => {
  assert.equal(SA.annotate(left+'AANCGTTA'+right,indexed).indexCandidates[0].status,'unresolved');
  const noBasis=JSON.parse(JSON.stringify(indexed)); delete noBasis.templates[0].index_spans[0].basis;
  assert.equal(SA.annotate(query,noBasis).indexCandidates[0].forwardSequence,null);
});
check('concrete source sample index supplies a variable window without changing literal annotations', () => {
  const source=left+'AACCGTTA'+right, changed=left+'GTCAGTAC'+right;
  const concrete={protocols,components:[component('source-primer',source,['p1'])],
    templates:[Object.assign({},template,{sequence:source})]};
  for (const query of [changed,SA.revcomp(changed)]) {
    const result=SA.annotate(query,concrete),hit=result.indexCandidates[0];
    assert.equal(hit.sequence,'GTCAGTAC'); assert.equal(hit.forwardSequence,'GTACTGAC');
    assert.equal(hit.status,'ready'); assert.equal(result.hits.length,0);
  }
  assert.equal(concrete.templates[0].sequence,source);
  assert.equal(SA.annotate(source,concrete).hits[0].component,'source-primer');
  assert.equal(SA.annotate(left.slice(0,-1)+'C'+'GTCAGTAC'+right,concrete).indexCandidates.length,0);
  const ambiguous=SA.annotate(left+'GTNAGTAC'+right,concrete).indexCandidates[0];
  assert.equal(ambiguous.status,'unresolved'); assert.equal(ambiguous.forwardSequence,null);
});
check('all structured sample-index windows vary independently in the same template', () => {
  const middle='GACCTAGTCGATGCAT', secondStart=left.length+8+middle.length;
  const concrete={protocols,components:[],templates:[Object.assign({},template,{
    sequence:left+'AACCGTTA'+middle+'TTCAGGAC'+right,
    index_spans:[template.index_spans[0],{start:secondStart,end:secondStart+8,
      feature_id:'index5',role:'i5',to_forward:'forward',basis:'explicit fixture orientation'}]
  })]};
  const query=left+'GTCAGTAC'+middle+'AGCTTCAG'+right;
  for (const input of [query,SA.revcomp(query)]) {
    const candidates=SA.annotate(input,concrete).indexCandidates;
    assert.equal(candidates.length,2);
    assert.deepEqual(Object.fromEntries(candidates.map(c=>[c.role,c.forwardSequence])),{i7:'GTACTGAC',i5:'AGCTTCAG'});
  }
});
check('short flanks never yield a guessed index', () => {
  const short=JSON.parse(JSON.stringify(indexed)); short.templates[0].sequence='AAA'+'N'.repeat(8)+right; short.templates[0].index_spans[0].start=3;short.templates[0].index_spans[0].end=11;
  const candidate=SA.annotate('AAA'+'AACCGTTA'+right,short).indexCandidates[0];
  assert.equal(candidate.status,'unresolved'); assert.equal(candidate.sequence,null);
});
check('U remains distinct unless conversion explicitly requested', () => {
  assert.equal(SA.annotate(seq.replace('T','U'),data).hits.filter(h=>h.component==='distinct').length,0);
  assert.equal(SA.annotate(seq.replace('T','U'),data,{convertUToT:true}).hits.filter(h=>h.component==='distinct').length,1);
});
check('successive start positions retain genuine overlapping hits', () => {
  const repeated={protocols,components:[component('repeat','A'.repeat(12),['p1'])]};
  assert.deepEqual(SA.annotate('A'.repeat(14),repeated).hits.map(h=>h.start),[0,1,2]);
});
check('inferred aliases do not leak protocol associations by default', () => {
  const mixed=component('mixed',seq,['p1']);
  mixed.occurrences.push({protocol:'p2',inferred:true}); mixed.protocols.push('p2');
  assert.deepEqual(SA.annotate(seq,{protocols,components:[mixed]}).hits[0].protocols,['p1']);
});
check('unrelated templates never create unresolved index records', () => {
  const short=JSON.parse(JSON.stringify(indexed)); short.templates[0].sequence='AAA'+'N'.repeat(8)+right; short.templates[0].index_spans[0].start=3;short.templates[0].index_spans[0].end=11;
  assert.equal(SA.annotate('C'.repeat(50),short).indexCandidates.length,0);
});
check('reverse-complement extraction preserves U as a mismatch', () => {
  const reverse=SA.revcomp(query);
  const bad=reverse.replace('T','U');
  assert.equal(SA.annotate(bad,indexed).indexCandidates.length,0);
});
check('index candidates are capped with an explicit truncation flag', () => {
  const many=SA.annotate(query.repeat(1100),indexed);
  assert.equal(many.indexCandidates.length,1000); assert.equal(many.indexCandidatesTruncated,true);
});
check('aliases follow selected provenance, including inference opt in', () => {
  const mixed=component('mixed',seq,['p1']);
  mixed.aliases=['inferred alias','real alias'];
  mixed.occurrences[0].name='real alias';
  mixed.occurrences.push({protocol:'p2',name:'inferred alias',inferred:true});
  const catalogue={protocols,components:[mixed]};
  assert.deepEqual(SA.annotate(seq,catalogue).hits[0].aliases,['real alias']);
  assert.deepEqual(SA.annotate(seq,catalogue,{includeInferred:true}).hits[0].aliases,['inferred alias','real alias']);
});
check('core result limits cannot be disabled or exceed protocol maximum', () => {
  assert.equal(SA.annotate('AAGTCGTACGAT',data,{maxProtocols:1000}).protocols.length,20);
  assert.equal(SA.annotate('AAGTCGTACGAT',data,{maxProtocols:3}).protocols.length,3);
  assert.equal(SA.annotate(seq.repeat(20),data,{maxHits:'invalid'}).hits.length,40);
});
check('symmetric template anchors cannot assign a nonpalindromic index orientation', () => {
  const anchor='AGCTACGTCAGTTCGA';
  // A concrete source index must not make symmetric flanks appear directional.
  const symmetric={protocols,components:[],templates:[Object.assign({},template,{sequence:anchor+'TTCAGGAC'+SA.revcomp(anchor)})]};
  const candidates=SA.annotate(anchor+'AACCGTTA'+SA.revcomp(anchor),symmetric).indexCandidates;
  assert.equal(candidates.length,1); assert.equal(candidates[0].strand,'both');
  assert.equal(candidates[0].status,'unresolved'); assert.equal(candidates[0].forwardSequence,null);
  const capped=SA.annotate((anchor+'AACCGTTA'+SA.revcomp(anchor)).repeat(1100),symmetric);
  assert.equal(capped.indexCandidatesTruncated,true);
  assert.ok(capped.indexCandidates.every(candidate=>candidate.status==='unresolved' && !candidate.forwardSequence));
});
const geometries=[{id:'canonical.i7',role:'i7',outer:left,inner:right,
  to_forward:'reverse_complement',basis:'Synthetic canonical geometry fixture',min_anchor:12}];
check('canonical arms directly resolve variable-width indexes without a protocol catalogue',()=>{
  for(const width of [6,8,10,12,16]) {
    const index='GTCAGTACAGTTCGAC'.slice(0,width),primer=left+index+right+'GACT';
    for(const input of [primer,SA.revcomp(primer)]) {
      const result=SA.indexGeometryCandidates(input,geometries);
      assert.equal(result.truncated,false);assert.equal(result.candidates.length,1);
      const hit=result.candidates[0];assert.equal(hit.forwardSequence,SA.revcomp(index));
      assert.equal(hit.status,'ready');assert.equal(hit.protocol,undefined);
      assert.equal(hit.end-hit.start,width);
    }
  }
});
check('canonical arms accept exact terminal inner prefixes and preserve strand coordinates',()=>{
  const index='AGCTTCAG';
  for(const width of [12,13,15,16]) for(const strand of ['+','-']) {
    const primer=left+index+right.slice(0,width),input=strand==='+'?primer:SA.revcomp(primer);
    const hit=SA.indexGeometryCandidates(input,geometries).candidates[0];
    assert.equal(hit.forwardSequence,SA.revcomp(index));assert.equal(hit.strand,strand);
    assert.deepEqual([hit.start,hit.end],strand==='+'?[left.length,left.length+8]:[width,width+8]);
  }
});
check('canonical geometry retains literal arms and ambiguous observed indexes cannot export',()=>{
  for(const primer of [left.slice(0,-1)+'C'+'AGCTTCAG'+right,left+'AGCTTCAG'+right.slice(0,11),left+'AGCTTCAG'+right.slice(0,12)+'A',left+right])
    assert.equal(SA.indexGeometryCandidates(primer,geometries).candidates.length,0);
  const hit=SA.indexGeometryCandidates(left+'AGNTTCAG'+right,geometries).candidates[0];
  assert.equal(hit.status,'unresolved');assert.equal(hit.forwardSequence,null);
});
check('duplicate canonical families agree on one index span without increasing candidate counts',()=>{
  const duplicate={...geometries[0],id:'same-arms'},shorter={...geometries[0],id:'shorter-arm',inner:right.slice(0,12)};
  const result=SA.indexGeometryCandidates(left+'AGCTTCAG'+right,[...geometries,duplicate,shorter]);
  assert.equal(result.candidates.length,1);assert.equal(result.candidates[0].forwardSequence,'CTGAAGCT');
});
check('nested PE-style full and suffix arms preserve competing boundaries',()=>{
  const rules=[{...geometries[0],id:'full-PE-fixture',inner:'CGGTC'+right},
    {...geometries[0],id:'PE-suffix-fixture',inner:right}];
  const result=SA.indexGeometryCandidates(left+'AGCTTCAG'+'CGGTC'+right,rules);
  assert.deepEqual(result.candidates.map(hit=>hit.end-hit.start),[8,13]);
  assert.equal(new Set(result.candidates.map(hit=>hit.forwardSequence)).size,2);
});
check('canonical geometry candidates are bounded and explicitly report truncation',()=>{
  const result=SA.indexGeometryCandidates((left+'AGCTTCAG'+right).repeat(1100),geometries);
  assert.equal(result.candidates.length,1000);assert.equal(result.truncated,true);
});
console.log(count+' checks passed');

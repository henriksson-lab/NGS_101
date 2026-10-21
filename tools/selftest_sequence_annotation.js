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
  const symmetric={protocols,components:[],templates:[Object.assign({},template,{sequence:anchor+'N'.repeat(8)+SA.revcomp(anchor)})]};
  const candidates=SA.annotate(anchor+'AACCGTTA'+SA.revcomp(anchor),symmetric).indexCandidates;
  assert.equal(candidates.length,1); assert.equal(candidates[0].strand,'both');
  assert.equal(candidates[0].status,'unresolved'); assert.equal(candidates[0].forwardSequence,null);
  const capped=SA.annotate((anchor+'AACCGTTA'+SA.revcomp(anchor)).repeat(1100),symmetric);
  assert.equal(capped.indexCandidatesTruncated,true);
  assert.ok(capped.indexCandidates.every(candidate=>candidate.status==='unresolved' && !candidate.forwardSequence));
});
console.log(count+' checks passed');

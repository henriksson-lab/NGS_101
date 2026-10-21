#!/usr/bin/env node
/* Browser-local annotation interaction tests with the real exact matcher. */
'use strict';
const fs = require('fs'), vm = require('vm'), assert = require('assert'), path = require('path');
const root = path.resolve(__dirname,'..');
class Element {
  constructor(tag) { this.tag=tag; this.children=[]; this.listeners={}; this.dataset={}; this.value=''; this.checked=false; this.hidden=false; this._text=''; this.focused=false; }
  set textContent(value) { this._text=value; this.children=[]; }
  get textContent() { return this._text+this.children.map(child=>child.textContent || '').join(''); }
  append(...children) { this.children.push(...children); }
  replaceChildren(...children) { this.children=children; this._text=''; }
  addEventListener(name,fn) { this.listeners[name]=fn; }
  setAttribute(name,value) { this[name]=value; }
  get options() { return this.children; }
  focus() { this.focused=true; }
}
const ids = ['annotate-form','sequence-input','convert-u','include-inferred','annotate-button','example-button','clear-button','annotation-status','annotation-results','result-summary','sequence-preview','preview-range','preview-prev','preview-next','region-detail','region-alignments','match-rows','match-table-limit','more-matches','protocol-matches','protocol-explanation'];
const pageSource=fs.readFileSync(path.join(root,'build_annotate.py'),'utf8');
for (const id of ids) assert(pageSource.includes('id="'+id+'"'),'Page declares '+id);
assert(!pageSource.includes('id="index-i7"'),'Sample-sheet indexes have moved out of the detect view');
assert(pageSource.includes('sequence_alignment.css'),'Page includes alignment styles');
const elements = Object.fromEntries(ids.map(id=>[id,new Element('div')]));
const get = id => elements[id];
const left='ACGTCAGATCGT', middle='GCTACGATCGTA', right='TGCATACGTGCA';
const full = left+middle, prefix='TTGGA', suffix='CCATT';
const catalogue={
  protocols:[{id:'p',title:'Synthetic protocol',url:'synthetic/diagram.html'}],
  components:[
    {id:'part',sequence:left,aliases:['INFERRED_ONLY_ALIAS'],occurrences:[{name:'<img src=x onerror=evil()>',protocol:'p',role:'PCR primer',source:'<script>evil</script>',orientation:'-'}]},
    {id:'whole',sequence:full,occurrences:[{name:'Full primer',protocol:'p',role:'PCR primer',orientation:'+'}]},
    {id:'nested',sequence:middle,occurrences:[{name:'Nested component',protocol:'p',role:'oligo',orientation:'+'}]},
    {id:'overlap',sequence:middle+right,occurrences:[{name:'Overlapping primer',protocol:'p',role:'PCR primer',orientation:'+'}]},
    {id:'inferred',sequence:'GGTTACGTCGAT',occurrences:[{name:'Inferred primer',protocol:'p',inferred:true,role:'PCR primer'}]}
  ].concat(Array.from({length:105},(_,i)=>({id:'capped_'+i,sequence:'ACTGGTCACTGAGT',occurrences:[{name:'Capped primer '+i,protocol:'p',role:'PCR primer'}]}))),templates:[]
};
const api=require(path.join(root,'assets/sequence_annotation.js'));
const context={document:{getElementById:get,createElement:tag=>new Element(tag),createTextNode:text=>({textContent:text})},window:{SequenceAnnotation:api},fetch:async()=>({ok:true,json:async()=>catalogue}),setTimeout,Map,Promise,console};
vm.runInNewContext(fs.readFileSync(path.join(root,'assets/sequence_annotation_ui.js'),'utf8'),context);
const descendants=el=>[el,...(el.children || []).flatMap(descendants)];
const byTag=(el,tag)=>descendants(el).filter(child=>child.tag===tag);
async function submit(value) { get('sequence-input').value=value; await get('annotate-form').listeners.submit({preventDefault(){}}); }
function pick(select,text) { const option=select.options.find(o=>o.textContent.includes(text)); assert(option,'Selector contains '+text); select.value=option.value; select.listeners.change(); }
(async()=>{
  await submit(prefix+full+right+suffix);
  assert.equal(get('annotation-results').hidden,false);
  assert.equal(get('match-rows').children.length,4);
  assert(get('match-rows').textContent.includes('<img src=x onerror=evil()>'));
  assert(get('match-rows').textContent.includes('PCR primer'));
  assert(!get('protocol-matches').textContent.includes('INFERRED_ONLY_ALIAS'),'Protocol support labels use filtered match metadata');
  assert(descendants(get('match-rows')).some(el=>el.href==='synthetic/diagram.html'),'Same-site protocol href is usable');
  const regions=get('region-alignments').children;
  assert.equal(regions.length,1,'Longest overview region has one selector');
  const select=byTag(regions[0],'select')[0]; assert.equal(select.options.length,4,'Overlap and contained components remain selectable');
  pick(select,'<img');
  const aligned=byTag(regions[0],'pre')[0];
  assert.equal(aligned.children[0].textContent,'Query    '+prefix+full+right+suffix,'Query keeps context on both sides');
  assert.equal(aligned.children[1].textContent,'Selected '+' '.repeat(prefix.length)+left,'Exact component occupies its query coordinates');
  assert(regions[0].textContent.includes(api.revcomp(left)),'Declaration preserves its local reverse orientation');
  assert.equal(byTag(regions[0],'img').length,0,'Catalogue names remain literal text');
  pick(select,'Overlapping primer');
  assert.equal(byTag(regions[0],'pre')[0].children[1].textContent,'Selected '+' '.repeat(prefix.length+left.length)+middle+right,'Overlapping component can extend beyond selected overview region');
  const mark=byTag(get('sequence-preview'),'mark')[0]; mark.listeners.click();assert(select.focused,'Highlighted region focuses its selector');
  let prevented=false;mark.listeners.keydown({key:'Enter',preventDefault(){prevented=true;}});assert(prevented,'Highlight is keyboard accessible');
  await submit(prefix+api.revcomp(full)+suffix);
  const reverseRegion=get('region-alignments').children[0]; const reverseSelect=byTag(reverseRegion,'select')[0];pick(reverseSelect,'Full primer');
  assert(reverseRegion.textContent.includes('reverse-complemented'),'Reverse strand orientation is explicit');
  assert.equal(byTag(reverseRegion,'pre')[0].children[1].textContent,'Selected '+' '.repeat(prefix.length)+api.revcomp(full),'Reverse hit aligns exact query-oriented sequence');
  await submit('GGTTACGTCGAT');assert.equal(get('region-alignments').children.length,0,'Inferred declarations are excluded by default');
  get('include-inferred').checked=true; await submit('GGTTACGTCGAT');assert(get('region-alignments').textContent.includes('Includes inferred declaration'));
  await submit('ACTGGTCACTGAGT');
  assert.equal(byTag(get('region-alignments'),'select')[0].options.length,100,'Region selector remains bounded for common sequences');
  assert(get('region-alignments').textContent.includes('capped at 100'),'Additional region candidates are disclosed');
  await submit('INVALID!');assert.equal(get('annotation-status').dataset.error,'true');assert.equal(get('annotation-results').hidden,true);
  await submit('A'.repeat(995)+full+'T'.repeat(1000));
  assert.equal(get('sequence-preview').textContent.length,1000);assert.equal(get('preview-next').disabled,false);
  assert(get('region-alignments').textContent.includes('996–1019'),'Page-crossing region retains full coordinates');
  get('preview-next').listeners.click();assert(get('preview-range').textContent.startsWith('1001'));
  assert(get('region-alignments').textContent.includes('996–1019'),'Region selection survives across preview boundaries');
  await submit('A'.repeat(100000));assert.equal(get('sequence-preview').textContent.length,1000);assert.equal(get('region-alignments').children.length,0,'New query clears stale alignments');
  get('clear-button').listeners.click();assert.equal(get('sequence-input').value,'');assert.equal(get('annotation-results').hidden,true);
  console.log('PASS annotation UI: literal text, roles, protocol links, per-region selectors, nested and overlapping exact alignments, declaration orientation, reverse complements, keyboard access, inference opt-in, preview boundaries, invalid input and 100kb input');
})().catch(error=>{console.error(error);process.exitCode=1;});

#!/usr/bin/env node
/* UI boundary tests using the real matcher and index-orientation policy, without a browser dependency. */
'use strict';
const fs = require('fs'), vm = require('vm'), assert = require('assert'), path = require('path');
const root = path.resolve(__dirname,'..');
class Element {
  constructor(tag) { this.tag=tag; this.children=[]; this.listeners={}; this.dataset={}; this.value=''; this.checked=false; this.hidden=false; this._text=''; }
  set textContent(value) { this._text=value; this.children=[]; }
  get textContent() { return this._text+this.children.map(child=>child.textContent || '').join(''); }
  append(...children) { this.children.push(...children); }
  replaceChildren(...children) { this.children=children; this._text=''; }
  addEventListener(name,fn) { this.listeners[name]=fn; }
  setAttribute(name,value) { this[name]=value; }
  get options() { return this.children; }
  focus() {} select() {}
}
const ids = ['annotate-form','sequence-input','convert-u','include-inferred','annotate-button','example-button','clear-button','annotation-status','annotation-results','result-summary','sequence-preview','preview-range','preview-prev','preview-next','region-detail','match-rows','match-table-limit','more-matches','protocol-matches','protocol-explanation','index-i7','index-i5','index-instrument','index-workflow','index-sample-id','index-summary','index-unresolved','index-fields','index-candidate-details','index-profile-note','index-snippet','index-output','copy-status','copy-index-tsv','copy-index-csv'];
const pageSource=fs.readFileSync(path.join(root,'build_annotate.py'),'utf8');
for (const id of ids) assert(pageSource.includes('id="'+id+'"'),'Page declares '+id);
const elements = Object.fromEntries(ids.map(id=>[id,new Element('div')]));
const get = id => elements[id];
for (const id of ['index-i5','index-i7','index-instrument','index-workflow']) { const option = new Element('option'); option.value=''; get(id).append(option); }
const left='ACGTCAGATCGT', middle='GCTACGATCGTA', right='TGCATACGTGCA';
const sequence=left+'AACCGGTT'+middle+'ACGTTGCA'+right;
const catalogue={
  protocols:[{id:'p',title:'Synthetic protocol',url:'synthetic/diagram.html'}],
  components:[{id:'part',sequence:left,aliases:['INFERRED_ONLY_ALIAS'],occurrences:[{name:'<img src=x onerror=evil()>',protocol:'p',role:'pcr_primer',source:'<script>evil</script>'}]}],
  templates:[{id:'template',protocol:'p',label:'Synthetic dual-index template',sequence:left+'NNNNNNNN'+middle+'NNNNNNNN'+right,source:'Synthetic fixture',index_spans:[
    {feature_id:'i7',role:'i7',start:12,end:20,to_forward:'forward',basis:'Explicit synthetic forward basis'},
    {feature_id:'i5',role:'i5',start:32,end:40,to_forward:'forward',basis:'Explicit synthetic forward basis'}]}]
};
const context={document:{getElementById:get,createElement:tag=>new Element(tag),createTextNode:text=>({textContent:text})},window:{SequenceAnnotation:require(path.join(root,'assets/sequence_annotation.js')),IlluminaSampleSheet:require(path.join(root,'assets/illumina_samplesheet.js'))},navigator:{},fetch:async()=>({ok:true,json:async()=>catalogue}),setTimeout,Map,Promise,console};
vm.runInNewContext(fs.readFileSync(path.join(root,'assets/sequence_annotation_ui.js'),'utf8'),context);
const descendants=el=>[el,...(el.children || []).flatMap(descendants)];
const change=(id,value)=>{ get(id).value=value; get(id).listeners[id==='index-sample-id'?'input':'change'](); };
async function submit(value) { get('sequence-input').value=value; await get('annotate-form').listeners.submit({preventDefault(){}}); }
(async()=>{
  await submit(sequence);
  assert.equal(get('annotation-results').hidden,false);
  assert.equal(get('match-rows').children.length,1);
  assert(get('match-rows').textContent.includes('<img src=x onerror=evil()>'));
  assert(get('match-rows').textContent.includes('pcr_primer'));
  assert(!get('protocol-matches').textContent.includes('INFERRED_ONLY_ALIAS'),'Protocol support labels use filtered match metadata');
  assert(descendants(get('match-rows')).some(el=>el.href==='synthetic/diagram.html'),'Same-site protocol href is usable');
  get('index-sample-id').value='Synthetic_1';
  change('index-instrument','nextseq500'); change('index-workflow','manual-bcl2fastq');
  change('index-i7','0'); change('index-i5','1');
  assert(get('index-fields').textContent.includes('TGCAACGT'),'Manual NextSeq i5 is reverse-complemented');
  assert(get('index-output').value.startsWith('[Data]\nSample_ID,index,index2\nSynthetic_1,AACCGGTT,TGCAACGT'));
  assert(descendants(get('index-profile-note')).some(el=>el.href && el.href.startsWith('https://knowledge.illumina.com/')),'Orientation guidance href is usable');
  change('index-workflow','bssh-bclconvert-v2');
  assert(get('index-output').value.startsWith('[Header]\nIndexOrientation,Forward\n\n[BCLConvert_Data]\n'),'Mandatory orientation flag is inside Header');
  assert(get('index-output').value.endsWith('Synthetic_1,AACCGGTT,ACGTTGCA'));
  change('index-workflow','run-planning');
  assert.equal(get('index-snippet').hidden,true); assert(get('index-profile-note').textContent.includes('Run Planning inputs only'));
  assert.equal(get('index-fields').children.length,2);
  change('index-instrument','hiseq4000');change('index-workflow','standalone-bclconvert');
  assert.equal(get('index-snippet').hidden,true);assert.equal(get('index-fields').children.length,0);
  change('index-instrument','nextseq500');change('index-workflow','manual-bcl2fastq');
  change('index-sample-id','invalid sample');assert.equal(get('index-fields').children.length,0);assert.equal(get('index-output').value,'');assert.equal(get('index-snippet').hidden,true);
  change('index-sample-id','Synthetic_1');
  await get('copy-index-csv').listeners.click(); await new Promise(resolve=>setTimeout(resolve,0));
  assert(get('copy-status').textContent.includes('Select and copy'));
  change('index-i7',''); change('index-i5','1');
  assert.equal(get('index-fields').children.length,1);assert.equal(get('index-snippet').hidden,true);
  catalogue.templates[0].inferred=true;
  get('include-inferred').checked=true;await submit(sequence);
  assert(get('index-i7').textContent.includes('INFERRED'),'Inferred candidate is labelled before selection');
  change('index-i7','0');assert(get('index-candidate-details').textContent.includes('INFERRED'));
  await submit('INVALID!'); assert.equal(get('annotation-status').dataset.error,'true');assert.equal(get('annotation-results').hidden,true);
  await submit('A'.repeat(100000)); assert.equal(get('sequence-preview').textContent.length,1000);assert.equal(get('preview-next').disabled,false);
  get('preview-next').listeners.click();assert(get('preview-range').textContent.startsWith('1001'));
  get('clear-button').listeners.click();assert.equal(get('sequence-input').value,'');assert.equal(get('annotation-results').hidden,true);
  console.log('PASS annotation UI: matching, literal text, role labels, links, real orientation policy, required Header, planning-only inputs, unsupported profiles, invalid IDs, clipboard fallback, i5-only, error and bounded 100kb preview');
})().catch(error=>{console.error(error);process.exitCode=1;});

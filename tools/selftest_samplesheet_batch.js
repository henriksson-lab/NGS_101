#!/usr/bin/env node
'use strict';
const assert=require('node:assert/strict'),fs=require('fs'),vm=require('vm');
const api=require('../assets/samplesheet_batch.js'),matcher=require('../assets/sequence_annotation.js'),policy=require('../assets/illumina_samplesheet.js');
const left='ACGTCAGATCGT',right='GCTACGATCGTA',l2='TTGCGATCCGTA',r2='CACTGATCGTGC';
const p7=left+'AACCGGTA'+right,p5=l2+'ACGTTGCA'+r2;
function template(id,role,l,r,protocol='p'){return {id,protocol,label:id,sequence:l+'NNNNNNNN'+r,index_spans:[{feature_id:role,role,start:12,end:20,to_forward:'forward',basis:'Synthetic fixture'}]};}
const catalogue={protocols:[{id:'p',title:'<b>Synthetic</b>',url:'p/page.html'}],components:[],templates:[template('seven','i7',left,right),template('five','i5',l2,r2)]};
let count=0;function check(label,fn){fn();count++;console.log('PASS '+label);}
const row={sampleId:'Sample_1',primer1:p7,primer2:p5};
check('NovaSeq X bcl2fastq defaults and automatic roles in either primer order',()=>{
 for(const r of [row,{...row,primer1:p5,primer2:p7}]){const result=api.analyze([r],catalogue);assert.equal(result.profile.i5Orientation,'reverse_complement');assert.deepEqual(result.values,[['Sample_1','AACCGGTA','TGCAACGT']]);}
});
check('reverse-complement primer input preserves exported orientation',()=>assert.deepEqual(api.analyze([{...row,primer1:matcher.revcomp(p7)}],catalogue).values,[['Sample_1','AACCGGTA','TGCAACGT']]));
check('either primer can be omitted and blank input is an empty field',()=>{
 assert.deepEqual(api.analyze([{...row,primer2:''}],catalogue).values,[['Sample_1','AACCGGTA','']]);
 assert.deepEqual(api.analyze([{...row,primer1:''}],catalogue).values,[['Sample_1','','TGCAACGT']]);
 assert.equal(api.analyze([{sampleId:'',primer1:'',primer2:''}],catalogue).valid.length,0);
});
check('duplicate template matches do not require protocol/candidate selection',()=>{
 const c={...catalogue,templates:[...catalogue.templates,template('seven-copy','i7',left,right,'other')]};assert.equal(api.analyze([row],c).valid.length,1);
});
check('unresolved alternate role is preserved as a warning without overriding resolved geometry',()=>{
 const other=template('unknown','i7',left,right);delete other.index_spans[0].role;const c={...catalogue,templates:[...catalogue.templates,other]};
 const result=api.analyze([row],c);assert.equal(result.valid.length,1);assert(result.results[0].results[0].warnings.length);
});
check('conflicting index values and unknown supplied primers never export',()=>{
 const c={...catalogue,templates:[...catalogue.templates,{...template('conflict','i7',left,right),index_spans:[{...catalogue.templates[0].index_spans[0],to_forward:'reverse_complement'}]}]};
 assert.equal(api.analyze([row],c).valid.length,0);assert.match(api.analyze([row],c).results[0].error,/Conflicting/);
 assert.equal(api.analyze([{...row,primer2:'ACGTACGTACGT'}],catalogue).valid.length,0);
});
check('disagreeing paired primers, duplicate samples and invalid IDs are excluded',()=>{
 assert.equal(api.analyze([{...row,primer2:left+'TTTTGGTA'+right}],catalogue).valid.length,0);
 assert.equal(api.analyze([row,row],catalogue).valid.length,0);
 assert.equal(api.analyze([{...row,sampleId:'=bad'}],catalogue).valid.length,0);
});
check('reuse targets detected roles, preserves partner and accepts blank slots',()=>{
 const input=[row,{sampleId:'Sample_2',primer1:p5,primer2:left+'TTTTGGTA'+right},{sampleId:'Sample_3',primer1:p5,primer2:''}];
 const output=api.reuse(input,0,1,'i7',catalogue);assert.equal(output[1].primer1,p5);assert.equal(output[1].primer2,p7);assert.equal(output[2].primer2,p7);
 assert.equal(api.reuse(input,0,2,'i5',catalogue)[1].primer1,p5);
 assert.throws(()=>api.reuse(input,0,1,'i5',catalogue));
});
check('paste list accepts header and omitted first or second primer; bounds enforced',()=>{
 assert.deepEqual(api.parseTSV('Sample_ID\tPrimer 1\tPrimer 2\nS1\t'+p7+'\nS2\t\t'+p5),[{sampleId:'S1',primer1:p7,primer2:''},{sampleId:'S2',primer1:'',primer2:p5}]);
 assert.throws(()=>api.parseTSV('S1 '+p7));assert.throws(()=>api.analyze(Array(201).fill(row),catalogue));
});
check('workflow changes apply RC and required Header exactly once',()=>{
 assert.match(api.analyze([row],catalogue).csv,/\[Data\]\nSample_ID,index,index2\nSample_1,AACCGGTA,TGCAACGT/);
 assert.equal(api.analyze([row],catalogue,'novaseqx','standalone-bclconvert').values[0][2],'ACGTTGCA');
 assert.match(api.analyze([row],catalogue,'nextseq500','bssh-bclconvert-v2').csv,/^\[Header\]\nIndexOrientation,Forward\n\n\[BCLConvert_Data\]/);
 assert(api.analyze([row],catalogue,'novaseqx','run-planning').inputOnly);
});
check('inferred template is excluded',()=>assert.equal(api.analyze([row],{...catalogue,templates:catalogue.templates.map(t=>({...t,inferred:true}))}).valid.length,0));

// Small DOM boundary fixture uses the real batch detector and orientation formatter.
class Element{
 constructor(tag){this.tag=tag;this.children=[];this.listeners={};this.dataset={};this.value='';this.hidden=false;this._text='';}
 set textContent(v){this._text=v;this.children=[];}get textContent(){return this._text+this.children.map(c=>c.textContent||'').join('');}
 append(...items){for(const item of items){item.parent=this;this.children.push(item);}}replaceChildren(...items){this.children=[];this._text='';this.append(...items);}
 addEventListener(name,fn){this.listeners[name]=fn;}setAttribute(k,v){this[k]=v;}focus(){}select(){}remove(){this.parent.children=this.parent.children.filter(c=>c!==this);}
}
const ids=['batch-instrument','batch-workflow','batch-profile','batch-paste','batch-import','primer-rows','batch-add','batch-examples','batch-example-data','batch-detect','batch-status','batch-export','batch-output','batch-copy-tsv','batch-copy-csv','batch-copy-status','batch-copy-fallback'];
const elements=Object.fromEntries(ids.map(id=>[id,new Element('div')])),get=id=>elements[id];
get('batch-example-data').textContent=JSON.stringify([{sampleId:'Example_1',primer1:p5,primer2:p7},{sampleId:'Example_2',primer1:'',primer2:p7}]);
const source=fs.readFileSync(require('path').join(__dirname,'../build_samplesheet.py'),'utf8');for(const id of ids)assert(source.includes('id="'+id+'"'));
const context={document:{getElementById:get,createElement:tag=>new Element(tag),createTextNode:text=>({textContent:text})},window:{SampleSheetBatch:api,IlluminaSampleSheet:policy},navigator:{},fetch:async()=>({ok:true,json:async()=>catalogue}),setTimeout,Map,Promise};
vm.runInNewContext(fs.readFileSync(require('path').join(__dirname,'../assets/samplesheet_ui.js'),'utf8'),context);
(async()=>{
 const controls=row=>row.children.slice(1,3).flatMap(cell=>cell.children[1].children);
 const compactRows=()=>{
   get('primer-rows').children.forEach((tr,index)=>{
     assert.equal(controls(tr).length,index===0?4:0,'Reuse buttons occur only on the first current row');
     for(const cell of tr.children.slice(1,3)){assert.equal(cell.children[0].tag,'input');assert.equal(cell.children[0].type,'text');}
   });
 };
 compactRows();
 get('primer-rows').children[0].children[5].children[1].listeners.click();compactRows();
 assert.equal(get('primer-rows').children.length,2,'Removal promotes reuse controls to the next row');
 get('batch-add').listeners.click();compactRows();
 assert.equal(get('batch-instrument').value,'novaseqx');assert.equal(get('batch-workflow').value,'manual-bcl2fastq');assert.equal(get('primer-rows').children.length,3);
 get('batch-paste').value='Sample_ID\tPrimer1\tPrimer2\nOne\t'+p5+'\t'+p7+'\nTwo\t'+p7;
 get('batch-import').listeners.click();compactRows();await get('batch-detect').listeners.click();
 assert.equal(get('batch-export').hidden,false);assert.match(get('batch-output').value,/One,AACCGGTA,TGCAACGT/);assert.match(get('batch-output').value,/Two,AACCGGTA,/);
 const tr=get('primer-rows').children[0];assert.match(tr.children[5].textContent,/<b>Synthetic<\/b>/,'Catalogue names remain text');
 get('batch-copy-csv').listeners.click();await new Promise(r=>setTimeout(r,0));assert.equal(get('batch-copy-fallback').hidden,false);assert.match(get('batch-copy-fallback').value,/\[Data\]/);
 // Use the first row's i5 primer for every sample; second row's i7 remains intact.
 await tr.children[1].children[1].children[0].listeners.click();
 await get('batch-detect').listeners.click();assert.match(get('batch-output').value,/Two,AACCGGTA,TGCAACGT/);
 const active=get('primer-rows').children[0];active.children[1].children[0].value='BAD!';active.children[1].children[0].listeners.input();
 assert.equal(get('batch-export').hidden,true);assert.equal(get('batch-output').value,'');await get('batch-detect').listeners.click();assert.match(get('batch-status').textContent,/1 unresolved/);
 get('batch-workflow').value='run-planning';get('batch-workflow').listeners.change();await get('batch-detect').listeners.click();assert.equal(get('batch-copy-csv').hidden,true);assert.match(get('batch-output').value,/Sample_ID\tindex\tindex2/);
 get('batch-examples').listeners.click();compactRows();assert.equal(get('batch-export').hidden,true,'Loading examples clears stale output');
 assert.equal(get('primer-rows').children.length,2);assert.equal(get('primer-rows').children[0].children[0].children[0].value,'Example_1');
 assert.equal(get('primer-rows').children[1].children[1].children[0].value,'','Example may omit a primer');
 assert.match(get('batch-status').textContent,/Synthetic example/);await get('batch-detect').listeners.click();assert.match(get('batch-output').value,/Example_1/);
 console.log('PASS compact table DOM: single-line inputs, first-row reuse after add/import/removal, example loading and stale output clearing');
 console.log('PASS primer table DOM: defaults, list import, automatic roles, missing primer, literal text, per-row reuse, stale export clearing, clipboard fallback and planning-only inputs');
 console.log(count+' batch checks passed');
})().catch(e=>{console.error(e);process.exitCode=1;});

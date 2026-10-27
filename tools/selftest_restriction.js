'use strict';
const assert=require('assert'),fs=require('fs'),path=require('path');
const api=require('../assets/restriction_analysis.js');
const catalogue=JSON.parse(fs.readFileSync(path.join(__dirname,'../catalogue/restriction_enzymes.json'),'utf8'));
const annotations=fs.readFileSync(path.join(__dirname,'../catalogue/restriction_methylation.tsv'),'utf8').trim().split('\n').slice(1);
for(const line of annotations){const [name,dam,dcm,cpg]=line.split('\t');catalogue.enzymes.find(e=>e.name===name).methylation={dam,dcm,cpg};}
let checks=0;
function check(name,fn){fn();checks++;console.log('PASS '+name);}
const enzyme=name=>catalogue.enzymes.find(e=>e.name===name);
const analyze=(sequence,name,options)=>api.analyze(sequence,[typeof name==='string'?enzyme(name):name],options).results[0];
const eq=assert.deepStrictEqual;
check('single FASTA and lowercase normalize without losing query positions',()=>eq(api.normalize('>query\naCg t\nCC'),'ACGTCC'));
for(const bad of ['','>a\nACGT\n>b\nACGT','ACGU','AC-GT','ACGT123','A'.repeat(100001)])check('invalid input rejected: '+bad.slice(0,12),()=>assert.throws(()=>api.normalize(bad)));
check('palindromic EcoRI counts once, with separate top/bottom coordinates',()=>{const r=analyze('TGAATTCA','EcoRI');eq([r.cutCount,r.siteCount,r.sites[0].cuts[0].top,r.sites[0].cuts[0].bottom],[1,1,2,6]);eq(r.fragments.lengths,[2,6]);});
check('PstI retains a three-prime overhang rather than reversing its interval',()=>{const r=analyze('TCTGCAGA','PstI');eq([r.end,r.sites[0].cuts[0].top,r.sites[0].cuts[0].bottom],['3prime',6,2]);});
check('type IIS cutting is outside recognition and transforms both strands on reversal',()=>{const seq='AGGTCTCAAAAAAAA',forward=analyze(seq,'BsaI'),reverse=analyze(require('../assets/sequence_annotation.js').revcomp(seq),'BsaI');eq([forward.sites[0].cuts[0].top,forward.sites[0].cuts[0].bottom],[8,12]);eq([reverse.sites[0].cuts[0].top,reverse.sites[0].cuts[0].bottom],[seq.length-12,seq.length-8]);});
check('linear sites with out-of-bounds cuts remain recognition matches but not complete cuts',()=>{const r=analyze('GGTCTCA','BsaI');eq([r.cutCount,r.siteCount,r.outside],[0,1,1]);});
check('a cut at an existing linear end has no bond and is excluded',()=>{const r=analyze('AGAATTC','EcoRI');eq([r.cutCount,r.siteCount],[1,1]);const e={name:'Terminal',motifs:['GAATTC'],cuts:[[0,4]]};eq(analyze('GAATTCA',e).cutCount,0);});
check('circular geometry normalizes outside cuts to inter-base coordinates',()=>{const r=analyze('GGTCTCA','BsaI',{circular:true});eq([r.cutCount,r.sites[0].cuts[0].top,r.sites[0].cuts[0].bottom],[1,0,4]);eq(r.fragments.lengths,[7]);});
check('recognition across circular origin is found only in circular mode',()=>{eq(analyze('TTCAGAA','EcoRI').siteCount,0);const r=analyze('TTCAGAA','EcoRI',{circular:true});eq([r.siteCount,r.sites[0].start,r.sites[0].cuts[0].top,r.sites[0].cuts[0].bottom],[1,4,5,2]);});
check('motifs longer than the circular molecule do not match by repeatedly traversing it',()=>eq(analyze('GA','EcoRI',{circular:true}).siteCount,0));
check('two-cut enzymes expose two cleavage sites for one recognition site',()=>{const e=enzyme('BcgI'),r=analyze('T'.repeat(30)+e.motifs[0].replace(/N/g,'A')+'T'.repeat(35),e);eq([r.siteCount,r.cutCount],[1,2]);assert(!api.filter([r],{minimum:1,maximum:1}).length);assert(api.filter([r],{basis:'recognition',minimum:1,maximum:1}).length);});
check('overlapping recognition matches are not skipped',()=>eq(analyze('AAAAA',{name:'Overlap',motifs:['AAA'],cuts:[[1,2]]}).cutCount,3));
check('degenerate catalogue motif recognizes all allowed concrete letters',()=>eq(analyze('AGACTCA','HinfI').siteCount,1));
check('query ambiguity inside a motif wildcard can still be a certain match',()=>eq(analyze('AGANTCA','HinfI').siteCount,1));
check('ambiguous query bases outside motif wildcards are possible, not certain',()=>{const r=analyze('NNNNNN','EcoRI');eq([r.cutCount,r.possibleCuts],[0,1]);eq(api.filter([r],{minimum:0,maximum:0}).length,0);eq(api.filter([r],{minimum:1,maximum:1}).length,0);eq(api.filter([r],{minimum:1,maximum:1,includePossible:true}).length,1);});
check('nicking geometry preserves its single strand without inventing a duplex fragment',()=>{const r=analyze('TAAAT',{name:'Nick',motifs:['AAA'],cuts:[[1,null]]});eq([r.kind,r.cutCount,r.fragments],['nick',1,null]);});
check('unknown cut geometry is unassigned, never a zero-cut declaration',()=>{const r=analyze('TGAATTCA',{name:'Unknown',motifs:['GAATTC'],cuts:[]});eq([r.cutCount,r.siteCount],[null,1]);eq(api.filter([r],{minimum:0,maximum:0}).length,0);eq(api.filter([r],{basis:'recognition',minimum:1,maximum:1}).length,1);});
check('custom recognition alternatives with unknown cuts retain both patterns',()=>{const r=analyze('AGACTATGAT',{name:'Alternatives',motifs:['GACT','TGAT'],cuts:[]});eq(r.siteCount,2);});
check('cut range and enzyme search intersect rather than replace one another',()=>{const r=api.analyze('TGAATTCAGGAATTCATGGATCCAA',[enzyme('EcoRI'),enzyme('BamHI')]).results;eq(api.filter(r,{minimum:1,maximum:1,search:'EcoRI,BamHI'}).map(r=>r.enzyme.name),['BamHI']);eq(api.filter(r,{minimum:2}).map(r=>r.enzyme.name),['EcoRI']);});
check('CpG content is separate from sensitivity of isoschizomers',()=>{const rows=api.analyze('AACCGGTT',[enzyme('HpaII'),enzyme('MspI')]).results;eq(api.filter(rows,{context:'overlap_cpg'}).length,2);eq(api.filter(rows,{cpg:'unaffected'}).map(r=>r.enzyme.name),['MspI']);eq(api.filter(rows,{cpg:'sensitive'}).map(r=>r.enzyme.name),['HpaII']);});
check('Dam context includes an overlapping recognition flank',()=>{const r=analyze('GATCTAGATT','XbaI');eq(r.contexts.dam.certain,1);eq(api.filter([r],{context:'avoid_dam'}).length,0);});
check('Dcm motifs include either CCAGG or CCTGG',()=>{eq(api.contextAt('CCAGG',0,5,false).dcm,2);eq(api.contextAt('CCTGG',0,5,false).dcm,2);});
check('context filters can find motifs wrapping the circular origin',()=>eq(api.contextAt('GAAAAC',0,1,true).cpg,2));
check('uncertain CpG overlap cannot qualify as CpG-free',()=>{const r=analyze('NTGAATTCA','EcoRI');r.contexts.cpg={certain:0,possible:1};eq(api.filter([r],{context:'avoid_cpg'}).length,0);});
check('methylation requirement is distinct from insensitive and blocked',()=>{const r=analyze('AAGATCAA','DpnI');eq(api.filter([r],{dam:'requires'}).length,1);eq(api.filter([r],{dam:'unaffected'}).length,0);});
check('unannotated methylation sensitivity is never inferred from recognition',()=>{const r=analyze('ACCGGTT',{name:'NoAnnotation',motifs:['CCGG'],cuts:[[1,3]]});eq(api.filter([r],{cpg:'unaffected'}).length,0);eq(api.filter([r],{cpg:'unknown'}).length,1);});
check('position retention caps do not truncate counts or zero/unique filtering',()=>{const r=analyze('A'.repeat(30),{name:'Dense',motifs:['AAA'],cuts:[[1,2]]},{positionLimit:2});eq([r.cutCount,r.sites.length,r.positionsTruncated],[28,2,true]);eq(api.filter([r],{minimum:28,maximum:28}).length,1);assert(api.exportTSV([r],true).includes('\t2\t28'));});
check('custom enzyme offsets and optional sensitivity states parse explicitly',()=>{const [r]=api.parseCustom('Mine GAATTC 1 5 blocked unaffected unknown');eq([r.cuts,r.methylation.cpg,r.custom],[[[1,5]],'blocked',true]);});
for(const bad of ['Bad AXX 1 2','EcoRI GAATTC 1 5','Bad GAATTC 1.2 5','Bad GAATTC - - blocked unaffected nope','<script> ACGT 1 2'])check('invalid/duplicate custom declaration rejected: '+bad,()=>assert.throws(()=>api.parseCustom(bad,catalogue.enzymes)));
check('empty custom set is valid',()=>eq(api.parseCustom(''),[]));
check('invalid count ranges fail explicitly',()=>assert.throws(()=>api.filter([],{minimum:3,maximum:2})));
check('TSV export does not turn a custom name into a spreadsheet formula',()=>{const r=analyze('TAAAT',{name:'+Custom',motifs:['AAA'],cuts:[[1,2]]});assert(api.exportTSV([r],false).includes("'+Custom\t"));});
check('catalogue validates completely without a runtime Biopython import',()=>catalogue.enzymes.forEach(api.validateEnzyme));
check('a match limit fails explicitly instead of returning misleading partial counts',()=>assert.throws(()=>api.analyze('A'.repeat(30),[{name:'Dense',motifs:['AAA'],cuts:[[1,2]]}],{matchLimit:5}),/no partial counts/));
check('possible extra cuts prevent a certain unique-cutter claim',()=>{const r=analyze('AGAATTCAANNNNNNA','EcoRI');assert(r.cutCount===1 && r.possibleCuts>0);eq(api.filter([r],{minimum:1,maximum:1}).length,0);eq(api.filter([r],{minimum:1}).length,1);});
check('exhaustive short concrete sequences match an independent two-orientation cut oracle',()=>{
 const e={name:'Oracle',motifs:['RAN'],cuts:[[0,2]]};
 for(let code=0;code<1024;code++){
  let value=code,sequence='';for(let j=0;j<5;j++){sequence+='ACGT'[value%4];value=Math.floor(value/4);}
  for(const circular of [false,true]){
   const cuts=new Set(),hay=circular?sequence+sequence.slice(0,2):sequence;
   for(let p=0;p<(circular?5:3);p++){
    const s=hay.slice(p,p+3),pairs=[];
    if(/^[AG]A[ACGT]$/.test(s))pairs.push([p,p+2]);
    if(/^[ACGT]T[CT]$/.test(s))pairs.push([p+1,p+3]);
    for(let [top,bottom] of pairs)if(circular || top>0 && bottom<5){if(circular){top%=5;bottom%=5;}cuts.add(top+':'+bottom);}
   }
   eq(analyze(sequence,e,{circular}).cutCount,cuts.size);
  }
 }
});
console.log(`All ${checks} restriction checks passed.`);

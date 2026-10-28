'use strict';
// Independent numerical examples: Owczarzy 2008 worked example and Biopython 1.88 docs.
// No Biopython package is installed or imported; the shared repo exports its own constants.
const assert=require('assert'),cp=require('child_process'),path=require('path');
const api=require('../assets/primer_tm.js');
const py=cp.spawnSync('python3',['-c',"import json,sys;sys.path.insert(0,'lib');from chemdraw import tm_model_parameters;print(json.dumps(tm_model_parameters()))"],{cwd:path.join(__dirname,'..'),encoding:'utf8'});
assert.strictEqual(py.status,0,py.stderr);const model=JSON.parse(py.stdout);
let count=0;function check(label,fn){fn();count++;console.log('PASS '+label);}function near(got,want,epsilon=.005){assert(Math.abs(got-want)<epsilon,`${got} vs ${want}`);}
const conditions={target_nM:0,primer_nM:300,monovalent_mM:50,mg_mM:2.5,dntpTotal_mM:1.2},s='ACGTCAGTGCATGCTAGTCA';
const calc=(seq=s,c=conditions,tail=0)=>api.calculate(seq,c,model,tail);
check('published Owczarzy magnesium example predicts 67.9 C from 76.3 C reference',()=>near(api.correctTemperature(76.3+273.15,'AAGGCGAGTCAGGCTCAGTG',{monovalent_mM:5,mg_mM:1.5,dntpTotal_mM:0}).tm_C,67.9,.06));
const ref='CGTTCCAAAGATGTGGGCATGAGCTTAC',base={primer_nM:25,target_nM:25,monovalent_mM:50,mg_mM:0,dntpTotal_mM:0};
check('independent documented monovalent NN example gives 59.78 C',()=>near(calc(ref,base).tm_C,59.78));
check('independent documented mixed-ion NN example gives 66.81 C',()=>near(calc(ref,{...base,monovalent_mM:55,mg_mM:1.5}).tm_C,66.81));
check('independent documented dNTP-chelated mixed-ion example gives 66.04 C',()=>near(calc(ref,{...base,monovalent_mM:55,mg_mM:1.5,dntpTotal_mM:.6}).tm_C,66.04));
check('magnesium/dNTP binding conserves mass and satisfies association equilibrium',()=>{const free=api.freeMagnesium(2.5,1.2)/1000,bound=.0025-free,freeDntp=.0012-bound;assert(free>0&&bound>0&&freeDntp>0);near(bound/(free*freeDntp),30000,.0001);});
check('zero magnesium stays zero and absent dNTP leaves magnesium free',()=>{near(api.freeMagnesium(0,1.2),0);near(api.freeMagnesium(2.5,0),2.5);});
check('dNTP excess produces a small positive free Mg, not negative salt',()=>{const free=api.freeMagnesium(.5,5);assert(free>0&&free<.01);assert(Number.isFinite(calc(s,{...conditions,mg_mM:.5,dntpTotal_mM:5}).tm_C));});
check('more primer raises predicted Tm',()=>assert(calc(s,{...conditions,primer_nM:1000}).tm_C>calc().tm_C));
check('higher disclosed magnesium raises Tm at fixed other conditions',()=>assert(calc(s,{...conditions,mg_mM:2.5}).tm_C>calc(s,{...conditions,mg_mM:1.5}).tm_C));
check('monovalent assumption has a measurable effect and is not ignored',()=>assert(Math.abs(calc(s,{...conditions,monovalent_mM:150}).tm_C-calc().tm_C)>.2));
check('5prime tail does not enter thermodynamics or binding GC',()=>{const withTail=calc('GGGGGG'+s,conditions,6),plain=calc();near(withTail.tm_C,plain.tm_C);near(withTail.gc_percent,plain.gc_percent);assert.strictEqual(withTail.sequence,s);assert.strictEqual(withTail.fullSequence,'GGGGGG'+s);});
check('primer-excess concentration and equimolar concentration remain distinct',()=>assert(calc(ref,{...base,target_nM:0}).tm_C>calc(ref,base).tm_C));
check('self-complementary strands use symmetry and total concentration',()=>{const r=calc('GAATTC');assert(r.selfComplementary);assert.strictEqual(r.effective_nM,300);assert.strictEqual(calc('GAATTC',{...conditions,target_nM:300}).effective_nM,600);});
check('pure magnesium, mixed salts and monovalent regimes are explicit',()=>{assert.strictEqual(calc().regime,'mixed');assert.strictEqual(calc(s,{...conditions,monovalent_mM:0}).regime,'magnesium');assert.strictEqual(calc(s,{...conditions,mg_mM:0}).regime,'monovalent');});
check('short binding lengths produce a scope warning',()=>assert(calc('ACGT').warnings.some(w=>w.includes('length'))));
for(const bad of ['ACGNACGT','ACGUACGT','ACG-T','/5Phos/ACGT','>a\nACGT\n>b\nACGT','A','A'.repeat(201)])check('unsupported sequence rejected: '+bad.slice(0,18),()=>assert.throws(()=>calc(bad)));
for(const c of [{...conditions,primer_nM:0},{...conditions,primer_nM:NaN},{...conditions,monovalent_mM:-1},{...conditions,mg_mM:-1},{...conditions,dntpTotal_mM:-1},{...conditions,target_nM:301},{...conditions,monovalent_mM:0,mg_mM:0}])check('invalid concentration rejected '+JSON.stringify(c),()=>assert.throws(()=>calc(s,c)));
check('invalid excluded tail cannot silently create a partial result',()=>{for(const tail of [-1,2.5,s.length-1])assert.throws(()=>calc(s,conditions,tail));});
console.log(`all ${count} primer Tm checks passed`);

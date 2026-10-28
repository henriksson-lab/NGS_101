/* Primer Tm form. Shared model and published preset metadata are supplied by the build. */
(function(){
  'use strict';
  const $=id=>document.getElementById('tm-'+id),api=window.PrimerTm,data=JSON.parse($('model').textContent),preset=data.preset;
  let submitted=false,last;
  const node=(tag,text,cls)=>{const e=document.createElement(tag);if(text!==undefined)e.textContent=text;if(cls)e.className=cls;return e;};
  function number(id){const raw=$(id).value,value=Number(raw);if(raw.trim()===''||!Number.isFinite(value))throw new Error('Enter a number for '+id.replace(/-/g,' ')+'.');return value;}
  function reset(){for(const [id,key]of [['primer','primer_nM'],['mg','mg_mM'],['dntp','dntpEach_mM'],['monovalent','monovalent_assumption_mM']])$(id).value=preset[key];}
  function calculate(){
    last=null;$('results').hidden=true;$('status').textContent='';
    try{
      const conditions={target_nM:0,primer_nM:number('primer'),mg_mM:number('mg'),dntpTotal_mM:number('dntp')*4,monovalent_mM:number('monovalent')};
      const results=[{name:'Primer 1',result:api.calculate($('first').value,conditions,data.thermodynamics,number('tail-first'))}];
      if($('second').value.trim())results.push({name:'Primer 2',result:api.calculate($('second').value,conditions,data.thermodynamics,number('tail-second'))});
      else if(number('tail-second')!==0)throw new Error('Enter primer 2 or set its excluded tail to zero.');
      const cards=results.map(({name,result:r})=>{const card=node('div',undefined,'tm-value');card.append(node('h3',name),node('strong',r.tm_C.toFixed(1)+' °C'),node('p',`${r.length} binding bases · ${r.gc_percent.toFixed(1)}% GC`),node('code',r.sequence));if(r.tail)card.append(node('small',`${r.tail} nt 5′ tail excluded`));return card;});
      $('values').replaceChildren(...cards);
      $('difference').textContent=results.length===2?'ΔTm: '+Math.abs(results[0].result.tm_C-results[1].result.tm_C).toFixed(1)+' °C':'';
      const r=results[0].result,custom=conditions.mg_mM!==preset.mg_mM||conditions.dntpTotal_mM!==preset.dntpEach_mM*4;
      $('conditions-summary').textContent=`${custom?'Custom conditions':'KAPA ReadyMix disclosed Mg²⁺/dNTP concentrations'} · assumed monovalent ions ${conditions.monovalent_mM} mM · Mg²⁺ ${conditions.mg_mM} mM (${r.freeMg_mM.toFixed(3)} mM estimated free) · dNTPs ${(conditions.dntpTotal_mM/4).toFixed(3)} mM each, ${conditions.dntpTotal_mM.toFixed(3)} mM total · primer ${conditions.primer_nM} nM.`;
      $('warnings').textContent=[...new Set(results.flatMap(({name,result:r})=>r.warnings.map(w=>name+': '+w)))].join(' ');
      const table=node('table'),head=node('tr');for(const label of ['Primer','ΔH (kcal/mol)','ΔS (cal/mol/K)','Salt regime','Symmetry'])head.append(node('th',label));table.append(head);
      for(const {name,result:r}of results){const row=node('tr');for(const value of [name,r.dh_kcal.toFixed(2),r.ds_cal.toFixed(2),r.regime,r.selfComplementary?'Self-complementary':'Non-self-complementary'])row.append(node('td',value));table.append(row);}const scroll=node('div',undefined,'tw');scroll.append(table);$('details').replaceChildren(scroll);
      last={conditions,results};$('results').hidden=false;
    }catch(error){$('status').textContent=error.message;}
  }
  $('form').addEventListener('submit',event=>{event.preventDefault();submitted=true;calculate();});
  $('form').addEventListener('input',()=>{if(submitted)calculate();});
  $('example').addEventListener('click',()=>{$('first').value='ACGTCAGTGCATGCTAGTCA';$('second').value='GCTACGATCGTACCTGACGA';$('tail-first').value=$('tail-second').value='0';submitted=true;calculate();});
  $('clear').addEventListener('click',()=>{submitted=false;last=null;$('first').value=$('second').value='';$('tail-first').value=$('tail-second').value='0';$('results').hidden=true;$('status').textContent='';});
  $('reset').addEventListener('click',()=>{reset();if(submitted)calculate();});
  $('export').addEventListener('click',()=>{
    if(!last)return;const c=last.conditions;
    const rows=[['primer','full_sequence_5to3','binding_sequence_5to3','excluded_5prime_nt','binding_length','gc_percent','estimated_tm_C','primer_nM','target_nM','assumed_monovalent_mM','mg_total_mM','mg_free_estimate_mM','dntp_each_mM','dntp_total_mM','model','conditions']];
    for(const {name,result:r}of last.results)rows.push([name,r.fullSequence,r.sequence,r.tail,r.length,r.gc_percent.toFixed(2),r.tm_C.toFixed(3),c.primer_nM,c.target_nM,c.monovalent_mM,c.mg_mM,r.freeMg_mM.toFixed(6),c.dntpTotal_mM/4,c.dntpTotal_mM,'SantaLucia1998+Owczarzy2008',c.mg_mM===preset.mg_mM&&c.dntpTotal_mM===preset.dntpEach_mM*4?preset.id:'custom']);
    const blob=new Blob([rows.map(row=>row.join('\t')).join('\n')+'\n'],{type:'text/tab-separated-values'}),url=URL.createObjectURL(blob),a=node('a');a.href=url;a.download='primer-tm.tsv';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
  });
}());

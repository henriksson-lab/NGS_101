/* Primer-pair table; sequences and declaration strings remain literal text. */
(function(){
  'use strict';
  const $=id=>document.getElementById(id),api=window.SampleSheetBatch,policy=window.IlluminaSampleSheet;
  const node=(tag,text)=>{const el=document.createElement(tag);if(text!==undefined)el.textContent=text;return el;};
  let rows=[],catalogue,loading,output,revision=0;
  function status(text,error=false){$('batch-status').textContent=text;$('batch-status').dataset.error=String(error);}
  function invalidate(){revision++;output=null;$('batch-export').hidden=true;$('batch-output').value='';$('batch-copy-status').textContent='';$('batch-copy-fallback').hidden=true;for(const row of rows){row.i7.replaceChildren();row.i5.replaceChildren();row.note.textContent='';}}
  function inputs(){return rows.map(r=>({sampleId:r.id.value,primer1:r.primers[0].value,primer2:r.primers[1].value}));}
  async function load(){
    if(catalogue)return catalogue;
    if(!loading)loading=fetch('assets/sequence_components.json').then(r=>{if(!r.ok)throw new Error('The component catalogue could not be loaded. Please try again.');return r.json();}).then(data=>catalogue=data).catch(error=>{loading=null;throw error;});
    return loading;
  }
  function addRow(value={}){
    if(rows.length>=api.MAX_ROWS){status('The table is limited to '+api.MAX_ROWS+' samples.',true);return;}
    invalidate();const tr=node('tr'),id=node('input'),primers=[node('textarea'),node('textarea')],i7=node('td'),i5=node('td'),note=node('div'),actions=node('td');
    let next=1;while(rows.some(r=>r.id.value==='Sample_'+next))next++;
    id.type='text';id.value=value.sampleId===undefined?'Sample_'+next:value.sampleId;id.setAttribute('aria-label','Sample_ID');id.addEventListener('input',invalidate);
    const first=node('td');first.append(id);tr.append(first);
    const row={tr,id,primers,i7,i5,note};
    primers.forEach((input,index)=>{
      input.value=value['primer'+(index+1)]||'';input.rows=3;input.spellcheck=false;input.setAttribute('aria-label','Primer '+(index+1));input.addEventListener('input',invalidate);
      const cell=node('td'),buttons=node('div');buttons.className='primer-reuse';cell.append(input,buttons);
      for(const role of ['i5','i7']){const button=node('button','Use for all '+role);button.type='button';button.addEventListener('click',()=>reuse(row,index+1,role));buttons.append(button);}
      tr.append(cell);
    });
    const remove=node('button','Remove');remove.type='button';remove.addEventListener('click',()=>{invalidate();rows=rows.filter(r=>r!==row);tr.remove();});
    actions.append(note,remove);tr.append(i7,i5,actions);rows.push(row);$('primer-rows').append(tr);
  }
  function replaceRows(values){invalidate();rows=[];$('primer-rows').replaceChildren();values.forEach(addRow);}
  function profile(){
    const p=policy.describeProfile($('batch-instrument').value,$('batch-workflow').value),box=$('batch-profile');
    box.textContent=p.note+' i7 → index; i5 → index2 ('+(p.i5Orientation==='reverse_complement'?'reverse complement':'forward')+').';
    for(const url of p.sources||[])if(/^https:\/\/[a-zA-Z0-9.-]*illumina\.com\//.test(url)){const a=node('a',' Illumina guidance');a.href=url;box.append(a);}
    return p;
  }
  async function reuse(row,slot,role){
    invalidate();const token=revision,values=inputs(),index=rows.indexOf(row);status('Checking primer role…');
    try{const data=await load();if(token!==revision)return;const replaced=api.reuse(values,index,slot,role,data);replaceRows(replaced);status('Applied the detected '+role+' primer to all '+rows.length+' rows.');}
    catch(error){if(token===revision)status(error.message,true);}
  }
  function field(cell,value,column){
    if(!value){cell.textContent='—';return;}
    const code=node('code',value),copy=node('button','Copy '+column);copy.type='button';copy.addEventListener('click',()=>copyText(value));cell.append(code,copy);
  }
  async function detect(){
    invalidate();const token=revision,values=inputs(),instrument=$('batch-instrument').value,workflow=$('batch-workflow').value;
    $('batch-detect').disabled=true;status('Loading templates and detecting indexes…');
    try{
      const data=await load();if(token!==revision)return;
      await new Promise(resolve=>setTimeout(resolve,0));if(token!==revision)return;
      const analysis=api.analyze(values,data,instrument,workflow);output=analysis;
      const protocols=new Map(data.protocols.map(p=>[p.id,p]));
      analysis.results.forEach((result,index)=>{
        const row=rows[index];
        if(result.empty){row.note.textContent='No primers supplied.';return;}
        if(result.error){row.note.textContent=result.error;return;}
        for(const role of ['i7','i5'])field(row[role],result.values[role]?policy.formatIndex(result.values[role],role,instrument,workflow).value:'',role==='i7'?'index':'index2');
        row.note.textContent=Object.keys(result.values).length===1?'One primer/index omitted. ':'Indexes detected. ';
        for(const detail of result.results)for(const warning of detail.warnings||[])row.note.append(document.createTextNode(warning+' '));
        row.note.append(document.createTextNode('Matching templates: '));
        for(const id of result.protocols.slice(0,5)){
          const p=protocols.get(id);if(!p)continue;
          const a=node('a',p.title+' ');if(typeof p.url==='string'&&!p.url.includes('..')&&/^[a-zA-Z0-9_+./%~-]+\.html$/.test(p.url)&&!p.url.startsWith('/'))a.href=p.url;row.note.append(a);
        }
        if(result.protocols.length>5)row.note.append(document.createTextNode('('+result.protocols.length+' total)'));
      });
      $('batch-export').hidden=!analysis.valid.length;
      $('batch-copy-csv').hidden=analysis.inputOnly;
      $('batch-output').value=analysis.inputOnly?analysis.tsv:analysis.csv;
      const skipped=analysis.results.filter(r=>!r.empty&&r.error).length;
      status(analysis.valid.length+' sample rows resolved; '+skipped+' unresolved rows excluded.'+(analysis.inputOnly?' These are Run Planning inputs. Use the sheet exported by that tool.':''),!!skipped);
    }catch(error){if(token===revision)status(error.message,true);}
    finally{$('batch-detect').disabled=false;}
  }
  async function copyText(text){
    const token=revision;
    try{if(!navigator.clipboard)throw new Error('clipboard unavailable');await navigator.clipboard.writeText(text);if(token===revision)$('batch-copy-status').textContent='Copied.';}
    catch(_){if(token!==revision)return;const fallback=$('batch-copy-fallback');fallback.hidden=false;fallback.value=text;fallback.focus();fallback.select();$('batch-copy-status').textContent='Select and copy the highlighted text using your browser’s Copy command.';}
  }
  for(const [id,options,selected]of [['batch-instrument',policy.instruments,'novaseqx'],['batch-workflow',policy.workflows,'manual-bcl2fastq']]){
    for(const item of options){const option=node('option',item.label);option.value=item.id;$(id).append(option);}$(id).value=selected;
    $(id).addEventListener('change',()=>{invalidate();profile();status('Settings changed. Extract indexes to update the output.');});
  }
  $('batch-add').addEventListener('click',()=>addRow());
  $('batch-detect').addEventListener('click',detect);
  $('batch-import').addEventListener('click',()=>{try{replaceRows(api.parseTSV($('batch-paste').value));status('Imported '+rows.length+' sample rows.');}catch(error){status(error.message,true);}});
  $('batch-copy-tsv').addEventListener('click',()=>{if(output&&output.valid.length)copyText(output.tsv);});
  $('batch-copy-csv').addEventListener('click',()=>{if(output&&output.valid.length&&!output.inputOnly)copyText(output.csv);});
  profile();for(let i=0;i<3;i++)addRow();
})();

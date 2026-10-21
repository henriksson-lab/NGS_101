/* Thin browser interface; candidate discovery and SVG come from Python. */
(function(){
  'use strict';
  const $=id=>document.getElementById(id);
  let worker,serial=0,revision=0,pending=new Map(),result,currentSVG='';
  function status(text){$('anneal-status').textContent=text;}
  function stop(){
    revision++;if(worker)worker.terminate();worker=null;
    for(const task of pending.values()){clearTimeout(task.timer);task.reject(new Error('Cancelled.'));}pending.clear();
    $('anneal-search').disabled=false;$('anneal-cancel').disabled=true;
  }
  function request(operation,payload){
    if(!worker){
      worker=new Worker('assets/anneal_worker.mjs',{type:'module'});
      worker.onmessage=({data})=>{const task=pending.get(data.id);if(!task)return;if(data.startupError){stop();status(data.error);return;}clearTimeout(task.timer);pending.delete(data.id);data.error?task.reject(new Error(data.error)):task.resolve(data);};
      worker.onerror=event=>{const error=event.message || 'The Python worker could not start.';stop();status(error);};
    }
    const id=++serial;
    return new Promise((resolve,reject)=>{
      const timer=setTimeout(()=>{stop();status('Calculation timed out. Shorten the sequences and try again.');},45000);
      pending.set(id,{resolve,reject,timer});worker.postMessage({id,operation,payload});
    });
  }
  function invalidate(){if(pending.size)stop();else revision++;clearOutput();$('anneal-results').hidden=true;}
  function clearOutput(){currentSVG='';$('anneal-diagram').replaceChildren();$('anneal-download').hidden=true;$('anneal-summary').textContent='';$('anneal-coordinate-text').textContent='';}
  function label(c){
    if(c.kind==='linear')return `Linear · ${c.length} bp · sequence 1: ${c.first_start+1}–${c.first_start+c.length}; sequence 2: ${c.second_start+1}–${c.second_start+c.length}`;
    return `Padlock sequence ${3-c.target} on ${c.target} · ${c.ligation_length} + ${c.extension_length} arm bases · target starts ${c.ligation_start+1}, ${c.extension_start+1}`;
  }
  function fitDiagram(svg){
    // Measure the shared renderer's actual geometry/text; keep its scale and trim
    // unused canvas space in both the displayed diagram and the exported SVG.
    const box=svg.getBBox(),padding=8;
    if(!box.width||!box.height)return;
    svg.setAttribute('viewBox',`${box.x-padding} ${box.y-padding} ${box.width+2*padding} ${box.height+2*padding}`);
    svg.setAttribute('width',String(Math.ceil(box.width+2*padding)));
    svg.setAttribute('height',String(Math.ceil(box.height+2*padding)));
    currentSVG=new XMLSerializer().serializeToString(svg);
  }
  async function show(){
    clearOutput();const token=++revision,candidate=result.candidates[Number($('anneal-candidate').value)];if(!candidate)return;
    $('anneal-cancel').disabled=false;status('Rendering…');
    try{
      const response=await request('render',{first:result.first,second:result.second,candidate});if(token!==revision)return;
      const documentSVG=new DOMParser().parseFromString(response.result.svg,'image/svg+xml');
      if(documentSVG.querySelector('parsererror')||documentSVG.documentElement.localName!=='svg')throw new Error('The model returned an invalid SVG.');
      const svg=document.importNode(documentSVG.documentElement,true);
      $('anneal-diagram').replaceChildren(svg);currentSVG=response.result.svg;fitDiagram(svg);
      if(document.fonts)document.fonts.ready.then(()=>{if(token===revision&&svg.isConnected)fitDiagram(svg);});
      $('anneal-summary').textContent=response.result.summary;
      const coordinates=response.result.coordinates || {},names={first:'Sequence 1 paired span',second:'Sequence 2 paired span',first_overhangs:'Sequence 1 overhangs (5′ / 3′)',second_overhangs:'Sequence 2 overhangs (5′ / 3′)',target_sequence:'Target sequence',probe_sequence:'Padlock sequence',target_extension_arm:'Target extension arm',target_ligation_arm:'Target ligation arm',probe_extension_arm:'Padlock extension arm',probe_ligation_arm:'Padlock ligation arm',probe_backbone:'Padlock backbone',target_gap:'Target gap'};
      $('anneal-coordinate-text').textContent=Object.entries(coordinates).map(([key,value])=>{
        if(!Array.isArray(value))return `${names[key] || key}: ${value}`;
        if(key.endsWith('_overhangs'))return `${names[key]}: ${value[0]} / ${value[1]} bases`;
        return `${names[key] || key}: ${value[0]===value[1]?'boundary after base '+value[0]:value[0]+1+'–'+value[1]} (${value[1]-value[0]} bases)`;
      }).join('\n');$('anneal-download').hidden=false;status('');
    }catch(error){if(token===revision)status(error.message);}
    finally{if(token===revision)$('anneal-cancel').disabled=true;}
  }
  $('anneal-form').addEventListener('submit',async event=>{
    event.preventDefault();const token=++revision;clearOutput();$('anneal-results').hidden=true;$('anneal-search').disabled=true;$('anneal-cancel').disabled=false;status('Starting Python and finding complementary spans…');
    try{
      const response=await request('search',{first:$('anneal-first').value,second:$('anneal-second').value,mode:$('anneal-mode').value,minimum:Number($('anneal-minimum').value)});if(token!==revision)return;
      result=response.result;const select=$('anneal-candidate');select.replaceChildren();
      result.candidates.forEach((candidate,i)=>{const option=document.createElement('option');option.value=String(i);option.textContent=label(candidate);select.append(option);});
      $('anneal-results').hidden=!result.candidates.length;$('anneal-count').textContent=result.candidates.length+' candidates'+(result.truncated?' · Search limit reached; alternatives are incomplete.':'');
      if(result.candidates.length){select.value='0';await show();}else status(result.truncated?'Search limit reached; shorten the sequences.':'No exact complementary spans met the minimum length.');
    }catch(error){if(token===revision)status(error.message);}
    finally{$('anneal-search').disabled=false;$('anneal-cancel').disabled=true;}
  });
  $('anneal-candidate').addEventListener('change',show);
  $('anneal-cancel').addEventListener('click',()=>{stop();status('Cancelled.');});
  for(const id of ['anneal-first','anneal-second','anneal-mode','anneal-minimum'])$(id).addEventListener('input',()=>{invalidate();status('');});
  $('anneal-example').addEventListener('click',()=>{invalidate();const examples=JSON.parse($('anneal-example').dataset.examples),value=examples[$('anneal-mode').value==='linear'?'linear':'padlock'];$('anneal-first').value=value.first;$('anneal-second').value=value.second;status('Synthetic example loaded.');});
  $('anneal-download').addEventListener('click',()=>{
    if(!currentSVG)return;const url=URL.createObjectURL(new Blob([currentSVG],{type:'image/svg+xml'})),link=document.createElement('a');link.href=url;link.download='annealing.svg';link.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
  });
})();

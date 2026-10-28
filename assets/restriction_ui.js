/* Local restriction analysis: worker search, inexpensive result filtering, safe DOM. */
(function () {
  'use strict';
  const api=window.RestrictionAnalysis,$=id=>document.getElementById('restriction-'+id);
  const node=(tag,text,cls)=>{const e=document.createElement(tag);if(text!==undefined)e.textContent=text;if(cls)e.className=cls;return e;};
  let data,loading,result,visible=[],shown=50,worker,generation=0;
  const ends={blunt:'Blunt','5prime':'5′ overhang','3prime':'3′ overhang',mixed:'Mixed',nick:'Nick',unknown:'Unassigned'};
  function status(text,error){$('status').textContent=text;$('status').className=error?'restriction-error':'';}
  function busy(value){$('analyze').disabled=value;$('cancel').disabled=!value;}
  function cancel(){generation++;if(worker)worker.terminate();worker=null;busy(false);}
  function invalidate(){viewer.clear();cancel();result=null;visible=[];$('results').hidden=true;status('');}
  async function catalogue(){
    if(data)return data;
    if(!loading)loading=fetch('assets/restriction_catalogue.json').then(r=>{if(!r.ok)throw new Error('Could not load the enzyme catalogue. Try again.');return r.json();}).then(value=>{
      if(value.version!==1 || !Array.isArray(value.enzymes))throw new Error('Unsupported enzyme catalogue.');
      value.enzymes.forEach(api.validateEnzyme);data=value;
      $('license').textContent=value.copyright+'\n\n'+value.license_text;
      const p=node('p',`${value.enzymes.length.toLocaleString()} named enzymes · ${value.enzymes.filter(e=>e.cuts.length).length} cut definitions · ${value.methylation_annotation_count} exact-name methylation annotations. ${value.release}, REBASE ${value.rebase_release}. `);
      const a=node('a','Pinned data source');a.href=value.source;p.append(a,document.createTextNode(' · '));
      const license=node('a',value.license_name);license.href=value.license_url;p.append(license);
      $('data-info').replaceChildren(p);return value;
    }).catch(e=>{loading=null;throw e;});
    return loading;
  }
  function number(id,min,max){const value=Number($(id).value);if($(id).value.trim()==='' || !Number.isInteger(value)||value<min||value>max)throw new Error('Enter a valid integer for '+id.replace(/-/g,' ')+'.');return value;}
  function options(){
    const mode=$('count-mode').value,n=mode==='any'?0:number('count',0,1000000);
    const minimum=['any','atmost'].includes(mode)?0:n;
    const maximum=['any','atleast'].includes(mode)?null:mode==='range'?number('max',n,1000000):n;
    const minLength=number('min-length',1,100),maxLength=number('max-length',minLength,100);
    return {minimum,maximum,minLength,maxLength,search:$('search').value,basis:$('basis').value,kind:$('kind').value,end:$('end').value,
            context:$('context').value,cpg:$('cpg').value,dam:$('dam').value,dcm:$('dcm').value,
            includePossible:$('possible').checked,commercial:$('commercial').checked,customOnly:$('custom-only').checked};
  }
  function countControls(){const mode=$('count-mode').value;$('count-label').hidden=mode==='any';$('max-label').hidden=mode!=='range';}
  function counts(n,possible){return n===null?'—':String(n)+(possible?` + ${possible} possible`:'');}
  function cutCoordinate(p){return p===null?'—':String(result.circular && p===0?result.length:p);}
  function download(text,name,type){const url=URL.createObjectURL(new Blob([text],{type}));const a=node('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
  const viewer=window.RestrictionView.create({$,node,download});
  function svgMap(row){
    const ns='http://www.w3.org/2000/svg',svg=document.createElementNS(ns,'svg');
    svg.setAttribute('viewBox','0 0 840 86');svg.setAttribute('class','restriction-map');svg.setAttribute('role','img');svg.setAttribute('aria-label',row.enzyme.name+' cut positions');
    const make=(tag,attrs,text)=>{const e=document.createElementNS(ns,tag);for(const [k,v]of Object.entries(attrs))e.setAttribute(k,v);if(text!==undefined)e.textContent=text;svg.append(e);return e;};
    make('rect',{x:0,y:0,width:840,height:86,fill:'#ffffff'});
    make('line',{x1:25,x2:815,y1:42,y2:42,stroke:'#243d52','stroke-width':2});
    make('text',{x:25,y:18,fill:'#16202b','font-family':'sans-serif','font-size':12},row.enzyme.name+' · '+result.length+' bp · '+(result.circular?'circular':'linear'));
    const coordinates=[...new Set(row.sites.flatMap(s=>s.cuts.filter(c=>c.valid).map(c=>c.top===null?c.bottom:c.top)))].sort((a,b)=>a-b);
    for(const p of coordinates.slice(0,250)){
      const x=25+790*p/result.length;const line=make('line',{x1:x,x2:x,y1:31,y2:54,stroke:'#a31b32','stroke-width':1.4});
      const title=document.createElementNS(ns,'title');title.textContent='Cut boundary '+cutCoordinate(p);line.append(title);
      if(coordinates.length<=16)make('text',{x,y:73,'text-anchor':'middle',fill:'#16202b','font-family':'sans-serif','font-size':11},cutCoordinate(p));
    }
    make('text',{x:25,y:84,fill:'#566879','font-family':'sans-serif','font-size':10},row.positionsTruncated?'Map uses retained positions; counts are complete.':'Top-strand boundaries; nicking enzymes use their specified strand.');
    return svg;
  }
  function contextLine(site){
    const n=result.length,flank=12,start=result.circular?site.start-flank:Math.max(0,site.start-flank),end=result.circular?site.start+site.length+flank:Math.min(n,site.start+site.length+flank);
    const seq=(a,b)=>{let s='';for(let i=a;i<b;i++)s+=result.sequence[((i%n)+n)%n];return s;};
    const pre=node('span',undefined,'restriction-context');pre.append(document.createTextNode(seq(start,site.start)),node('mark',seq(site.start,site.start+site.length)),document.createTextNode(seq(site.start+site.length,end)));return pre;
  }
  function details(row){
    const box=node('div',undefined,'restriction-sites'),map=svgMap(row);box.append(map);
    const save=node('button','Export map SVG');save.type='button';save.addEventListener('click',()=>download(new XMLSerializer().serializeToString(map),row.enzyme.name+'-restriction-map.svg','image/svg+xml'));box.append(save);
    if(row.fragments){const f=row.fragments;box.append(node('p','Reference fragment lengths: '+f.lengths.join(', ')+' bp'+(f.truncated?` (${f.lengths.length} of ${f.total} retained)`:'')));}
    if(row.positionsTruncated)box.append(node('p',`First ${row.sites.length} of ${row.sitesTotal} recognition matches retained for display/export; all counts are complete.`,'restriction-error'));
    if(row.outside)box.append(node('p',`${row.outside} candidate cut pair(s) at or beyond a linear end are excluded from cut counts.`,'muted'));
    const list=node('ol');
    for(const site of row.sites){
      const li=node('li'),wrap=site.start+site.length>result.length;
      li.append(document.createTextNode(`Site ${site.start+1}, ${site.strand==='both'?'palindromic':site.strand+' strand'}${wrap?', wraps origin':''}${site.certain?'':', possible'} · `));
      li.append(document.createTextNode(site.cuts.length?site.cuts.map(c=>`top/bottom ${cutCoordinate(c.top)} / ${cutCoordinate(c.bottom)}${c.valid?'':' (at/outside end)'}`).join('; '):'cut geometry unassigned'));
      const contexts=Object.entries(site.context).filter(([,v])=>v).map(([k,v])=>k.toUpperCase()+(v===1?' possible':''));
      if(contexts.length)li.append(document.createTextNode(' · '+contexts.join(', ')));
      li.append(contextLine(site));list.append(li);
    }
    box.append(list);
    if(!row.sites.length)box.append(node('p','No recognition matches.'));
    const sources=row.enzyme.methylation_sources||[];
    if(sources.length){const p=node('p','Methylation annotations: ');sources.forEach((url,i)=>{if(i)p.append(document.createTextNode(' · '));const a=node('a','Source '+(i+1));if(/^https:\/\//.test(url))a.href=url;p.append(a);});box.append(p);}
    if(Number.isInteger(row.enzyme.rebase_id)){const a=node('a','REBASE entry');a.href='https://identifiers.org/rebase:'+row.enzyme.rebase_id;box.append(a);}
    return box;
  }
  function appendRows(start){
    for(const row of visible.slice(start,shown)){
      const tr=node('tr'),name=node('td'),button=node('button',row.enzyme.name,'restriction-name');button.type='button';button.setAttribute('aria-expanded','false');name.append(button);if(row.enzyme.custom)name.append(node('small',' Custom'));tr.append(name,node('td',row.enzyme.motifs.join(' | '),'restriction-motif'),node('td',counts(row.cutCount,row.possibleCuts)),node('td',counts(row.siteCount,row.possibleSites)),node('td',ends[row.end]));
      for(const key of ['cpg','dam','dcm'])tr.append(node('td',api.LABELS[(row.enzyme.methylation||{})[key]||'unknown'],'restriction-sensitivity'));
      let extra;
      button.addEventListener('click',()=>{if(extra){extra.remove();extra=null;button.setAttribute('aria-expanded','false');return;}extra=node('tr',undefined,'restriction-details');const cell=node('td');cell.colSpan=8;cell.append(details(row));extra.append(cell);tr.after(extra);button.setAttribute('aria-expanded','true');});
      $('rows').append(tr);
    }
    $('more').hidden=shown>=visible.length;
  }
  function render(){
    countControls();if(!result)return;
    try {
      visible=api.filter(result.results,options());shown=50;$('rows').replaceChildren();$('results').hidden=false;
      $('summary').textContent=`${visible.length} of ${result.results.length} enzymes`;
      $('result-note').textContent='';
      const capped=visible.filter(r=>r.positionsTruncated).length;
      if(capped)$('result-note').textContent+=` Position display/export is capped for ${capped} enzyme(s); all counts are complete.`;
      $('export').disabled=!visible.length;$('export-sites').disabled=!visible.length;
      viewer.render(result,visible);appendRows(0);status('');
    }catch(error){$('rows').replaceChildren();$('results').hidden=true;visible=[];status(error.message,true);}
  }
  $('form').addEventListener('submit',async event=>{
    event.preventDefault();viewer.clear();cancel();const id=generation;result=null;$('results').hidden=true;busy(true);status('Loading catalogue…');
    try {
      const sequence=api.normalize($('sequence').value),catalog=await catalogue();if(id!==generation)return;
      const custom=api.parseCustom($('custom').value,catalog.enzymes),enzymes=catalog.enzymes.concat(custom);
      worker=new Worker('assets/restriction_worker.js');
      worker.onmessage=({data:message})=>{
        if(id!==generation)return;
        if(message.type==='progress')status(`Analyzing ${message.done} of ${message.total} enzymes…`);
        else if(message.type==='result'){result=message.result;worker.terminate();worker=null;busy(false);render();}
        else if(message.type==='error'){cancel();status(message.error,true);}
      };
      worker.onerror=event=>{event.preventDefault();if(id===generation){cancel();status('Analysis worker could not start. Try again.',true);}};
      worker.postMessage({sequence,enzymes,circular:$('topology').value==='circular'});
    }catch(error){if(id===generation){cancel();status(error.message,true);}}
  });
  $('sequence').addEventListener('input',invalidate);$('custom').addEventListener('input',invalidate);$('topology').addEventListener('change',invalidate);
  $('cancel').addEventListener('click',()=>{cancel();status('Cancelled.');});
  $('clear').addEventListener('click',()=>{invalidate();$('sequence').value='';$('custom').value='';$('file').value='';});
  $('example').addEventListener('click',()=>{invalidate();$('sequence').value='ACGTACGAGGAATTCTACGATCGGATCCAGTCCCGGGTACGACTGCAGATGCCGGCATAGCGGCCGCTAGCTTCTAGAGATCAGCTAGGGTCTCAACGTACTGACTAGT';});
  $('file').addEventListener('change',async()=>{const file=$('file').files[0];if(!file)return;invalidate();const id=generation;try{if(file.size>500000)throw new Error('The input file is too large.');const text=await file.text();if(id===generation){$('sequence').value=text;api.normalize(text);}}catch(error){if(id===generation)status(error.message,true);}});
  const filters=['search','count-mode','count','max','basis','kind','end','min-length','max-length','context','cpg','dam','dcm','possible','commercial','custom-only'];
  for(const name of filters)$(name).addEventListener(['search','count','max','min-length','max-length'].includes(name)?'input':'change',render);
  $('unique').addEventListener('click',()=>{$('basis').value='cleavage';$('kind').value='double';$('count-mode').value='exact';$('count').value='1';$('possible').checked=false;render();});
  $('all').addEventListener('click',()=>{$('count-mode').value='any';render();});
  $('reset-filters').addEventListener('click',()=>{$('search').value='';$('count-mode').value='exact';$('count').value='1';$('max').value='2';$('basis').value='cleavage';$('kind').value='double';$('end').value='any';$('min-length').value='1';$('max-length').value='100';for(const name of ['context','cpg','dam','dcm'])$(name).value='any';for(const name of ['possible','commercial','custom-only'])$(name).checked=false;render();});
  $('more').addEventListener('click',()=>{const start=shown;shown+=50;appendRows(start);});
  $('export').addEventListener('click',()=>download(api.exportTSV(visible,false),'restriction-enzymes.tsv','text/tab-separated-values'));
  $('export-sites').addEventListener('click',()=>download(api.exportTSV(visible,true),'restriction-sites.tsv','text/tab-separated-values'));
  countControls();
}());

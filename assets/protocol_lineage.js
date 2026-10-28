/* Curated ancestry and computed similarity are separate visual layers. No chemistry input. */
(function(){
  'use strict';
  const $=id=>document.getElementById('lineage-'+id),ns='http://www.w3.org/2000/svg';
  const node=(tag,text,cls)=>{const e=document.createElement(tag);if(text!==undefined)e.textContent=text;if(cls)e.className=cls;return e;};
  let data,nodes=new Map(),selected=null,focused=false,svg,world,positions=new Map(),layout={width:1,height:1},scale=1,pan={x:0,y:0},drag=null,shown=[];
  function element(tag,attrs,text,parent=world){const e=document.createElementNS(ns,tag);for(const [key,value]of Object.entries(attrs))e.setAttribute(key,value);if(text!==undefined)e.textContent=text;parent.append(e);return e;}
  function sourceLink(text,url){const a=node('a',text);a.href=url;return a;}
  function protocolButton(id){const b=node('button',nodes.get(id).title,'lineage-detail-button');b.type='button';b.addEventListener('click',()=>choose(id,true));return b;}
  function relations(){return data.relations.filter(e=>$(e.basis).checked);}
  function adjacency(edges,similarities){
    const map=new Map();for(const n of nodes.keys())map.set(n,new Set());
    for(const e of edges){map.get(e.parent).add(e.child);map.get(e.child).add(e.parent);}
    for(const e of similarities){map.get(e.first).add(e.second);map.get(e.second).add(e.first);}return map;
  }
  function choose(id,focus){
    if(!nodes.has(id))return;selected=id;if(focus!==undefined)focused=focus;
    $('search').value=nodes.get(id).title;$('focus').disabled=false;
    history.replaceState(null,'','?protocol='+encodeURIComponent(id));draw();details(id);
  }
  function details(id){
    const item=nodes.get(id),box=$('detail');box.replaceChildren(node('h2',item.title),sourceLink('Open protocol schematic',item.url));
    const edges=data.relations.filter(e=>id===e.parent||id===e.child);
    if(!edges.length)box.append(node('p','No curated ancestry links yet. This does not establish independent origins.','lineage-detail-muted'));
    for(const [heading,field]of [['Parents','child'],['Descendants','parent']]){
      const list=edges.filter(e=>e[field]===id);if(!list.length)continue;box.append(node('h3',heading));const ul=node('ul');
      for(const e of list){
        const li=node('li'),other=e.parent===id?e.child:e.parent;
        li.append(protocolButton(other),node('div',e.basis==='documented'?'Source-documented':'Inferred','lineage-basis'),node('div',e.scope));
        const more=node('button','Evidence','lineage-detail-button');more.type='button';more.addEventListener('click',()=>edgeDetails(e));li.append(more);ul.append(li);
      }box.append(ul);
    }
    const peers=data.similarities.filter(e=>e.first===id||e.second===id).sort((a,b)=>b.score-a.score);
    box.append(node('h3','Sequence similarity · not ancestry'));
    if(!peers.length)box.append(node('p','No retained similarity neighbours; sequence coverage may be incomplete.','lineage-detail-muted'));
    else{
      const list=node('ul');for(const e of peers){
        const li=node('li');li.append(protocolButton(e.first===id?e.second:e.first),node('div',Math.round(e.score*100)+'% weighted overlap · '+e.shared_count+' shared components'));
        const button=node('button','Shared sequences','lineage-detail-button');button.type='button';button.addEventListener('click',()=>similarityDetails(e));li.append(button);list.append(li);
      }box.append(list);
    }
  }
  function edgeDetails(edge){
    const box=$('detail');
    box.replaceChildren(node('h2','Method relationship'),protocolButton(edge.parent),node('p','↓'),protocolButton(edge.child),node('p',edge.basis==='documented'?'Source-documented':'Inferred','lineage-basis'),node('h3','Scope'),node('p',edge.scope),node('p',edge.rationale),node('h3','Limits'),node('p',edge.limits));
    for(const source of edge.evidence)box.append(sourceLink(source.url,source.url),node('p',source.locator));
    box.append(node('small','Reviewed '+edge.reviewed));
  }
  function similarityDetails(edge){
    const box=$('detail');
    box.replaceChildren(node('h2','Sequence similarity'),protocolButton(edge.first),node('p','↔'),protocolButton(edge.second),node('p',Math.round(edge.score*100)+'% weighted overlap · '+edge.shared_count+' shared components'),node('p','Computed similarity, not a curated ancestry link.','lineage-basis'));
    const list=node('ul');for(const part of edge.shared){
      const li=node('li');li.append(node('div',part.names.join(' / ')),node('code',part.sequence),node('div','Observed in '+part.protocol_count+' eligible protocols'));list.append(li);
    }box.append(list,node('p',data.similarity_method.metric),node('p',data.similarity_method.scope));
  }
  function transform(){if(world)world.setAttribute('transform','translate('+pan.x+' '+pan.y+') scale('+scale+')');}
  function fit(){
    const box=$('canvas').getBoundingClientRect();
    scale=Math.min(1.1,(box.width-36)/layout.width,(box.height-36)/layout.height);scale=Math.max(.04,scale);
    pan={x:(box.width-layout.width*scale)/2,y:(box.height-layout.height*scale)/2};transform();
  }
  function zoom(factor,x,y){
    const box=$('canvas').getBoundingClientRect();x=x===undefined?box.width/2:x;y=y===undefined?box.height/2:y;
    const next=Math.max(.04,Math.min(3,scale*factor)),ratio=next/scale;
    pan={x:x-(x-pan.x)*ratio,y:y-(y-pan.y)*ratio};scale=next;transform();
  }
  function draw(){
    if(!data)return;const edges=relations(),similarities=$('similar').checked?data.similarities:[];
    const adjacent=adjacency(edges,similarities),included=new Set();
    if(focused&&selected){included.add(selected);for(const n of adjacent.get(selected))included.add(n);}
    else{
      for(const e of edges){included.add(e.parent);included.add(e.child);}
      for(const e of similarities){included.add(e.first);included.add(e.second);}
      if($('isolated').checked)for(const id of nodes.keys())included.add(id);if(selected)included.add(selected);
    }
    const visibleEdges=edges.filter(e=>included.has(e.parent)&&included.has(e.child));
    const visibleSimilarities=similarities.filter(e=>included.has(e.first)&&included.has(e.second));
    shown=[...included].sort((a,b)=>nodes.get(a).title.localeCompare(nodes.get(b).title));
    const rank=new Map(shown.map(id=>[id,0])),indegree=new Map(shown.map(id=>[id,0])),children=new Map(shown.map(id=>[id,[]]));
    for(const e of visibleEdges){indegree.set(e.child,indegree.get(e.child)+1);children.get(e.parent).push(e.child);}
    const queue=shown.filter(id=>indegree.get(id)===0);
    for(let i=0;i<queue.length;i++){
      const id=queue[i];for(const child of children.get(id)){rank.set(child,Math.max(rank.get(child),rank.get(id)+1));indegree.set(child,indegree.get(child)-1);if(!indegree.get(child))queue.push(child);}
    }
    // Keep disconnected method families together instead of alphabetically
    // interleaving their generations. Similarity does not determine ancestry rank.
    const familyAdj=adjacency(visibleEdges,[]),visited=new Set(),families=[],unlinked=[];
    for(const id of shown){
      if(visited.has(id))continue;
      if(!familyAdj.get(id).size){unlinked.push(id);visited.add(id);continue;}
      const family=[id];visited.add(id);
      for(let i=0;i<family.length;i++)for(const peer of familyAdj.get(family[i])){
        if(included.has(peer)&&!visited.has(peer)){visited.add(peer);family.push(peer);}
      }
      families.push(family);
    }
    positions=new Map();let yOffset=24,width=300;
    for(const family of families){
      const levels=new Map();
      for(const id of family){const r=rank.get(id);if(!levels.has(r))levels.set(r,[]);levels.get(r).push(id);}
      const familyHeight=Math.max(...[...levels.values()].map(ids=>ids.length))*72;
      for(const [level,ids]of levels){
        ids.sort((a,b)=>nodes.get(a).title.localeCompare(nodes.get(b).title));
        for(let i=0;i<ids.length;i++)positions.set(ids[i],{x:24+level*270,y:yOffset+(familyHeight-ids.length*72)/2+i*72});
        width=Math.max(width,level*270+244);
      }
      yOffset+=familyHeight+20;
    }
    // Uncurated protocols use a compact grid, including in the similarity layer.
    if(unlinked.length){
      const columns=Math.min(4,unlinked.length);width=Math.max(width,columns*270-26);
      for(let i=0;i<unlinked.length;i++)positions.set(unlinked[i],{x:24+(i%columns)*270,y:yOffset+Math.floor(i/columns)*72});
      yOffset+=Math.ceil(unlinked.length/columns)*72;
    }
    layout={width,height:Math.max(100,yOffset+8)};
    $('canvas').replaceChildren();
    svg=element('svg',{xmlns:ns,role:'img','aria-label':'Protocol ancestry and optional sequence-similarity graph'},undefined,$('canvas'));world=element('g',{},undefined,svg);
    element('rect',{x:-10000,y:-10000,width:20000,height:20000,fill:'#fff'});
    for(const e of visibleSimilarities){
      const a=positions.get(e.first),b=positions.get(e.second);
      const path=element('path',{d:'M'+(a.x+100)+','+(a.y+22)+' Q'+((a.x+b.x)/2+140)+','+((a.y+b.y)/2-36)+' '+(b.x+100)+','+(b.y+22),fill:'none',stroke:'#9072ad','stroke-width':1.3,'stroke-dasharray':'2 5',class:'lineage-edge','data-similarity':'true',tabindex:0,role:'button','aria-label':'Sequence similarity, not ancestry'});
      path.addEventListener('click',()=>similarityDetails(e));path.addEventListener('keydown',event=>{if(event.key==='Enter')similarityDetails(e);});
    }
    for(const e of visibleEdges){
      const a=positions.get(e.parent),b=positions.get(e.child),x1=a.x+200,y1=a.y+22,x2=b.x,y2=b.y+22,color=e.basis==='documented'?'#287d89':'#a67415';
      const group=element('g',{class:'lineage-edge','data-relation':e.id,tabindex:0,role:'button','aria-label':nodes.get(e.parent).title+' to '+nodes.get(e.child).title+', '+e.basis});
      const d='M'+x1+','+y1+' C'+(x1+35)+','+y1+' '+(x2-35)+','+y2+' '+x2+','+y2;
      element('path',{d,fill:'none',stroke:'transparent','stroke-width':14},undefined,group);
      element('path',{d,fill:'none',stroke:color,'stroke-width':1.8,'stroke-dasharray':e.basis==='inferred'?'6 4':''},undefined,group);
      element('path',{d:'M'+(x2-8)+','+(y2-5)+' L'+x2+','+y2+' L'+(x2-8)+','+(y2+5),fill:'none',stroke:color,'stroke-width':1.8},undefined,group);
      group.addEventListener('click',()=>edgeDetails(e));group.addEventListener('keydown',event=>{if(event.key==='Enter')edgeDetails(e);});
    }
    for(const id of shown){
      const p=positions.get(id),item=nodes.get(id),group=element('g',{class:'lineage-node','data-protocol':id,tabindex:0,role:'button','aria-label':item.title});
      element('rect',{x:p.x,y:p.y,width:200,height:44,rx:7,fill:id===selected?'#e5f2f5':'#f6f8fa',stroke:id===selected?'#287d89':'#c8d5df','stroke-width':id===selected?2:1},undefined,group);
      element('title',{},item.title,group);
      const words=item.title.length>28?item.title.slice(0,26)+'…':item.title;
      element('text',{x:p.x+100,y:p.y+27,'text-anchor':'middle','font-family':'sans-serif','font-size':13,fill:'#20364b'},words,group);
      group.addEventListener('click',()=>choose(id));group.addEventListener('keydown',event=>{if(event.key==='Enter'||event.key===' '){event.preventDefault();choose(id);}});
    }
    if(!shown.length)element('text',{x:24,y:50,'font-family':'sans-serif','font-size':14,fill:'#526579'},'No relationships match these filters. Find a protocol or include unlinked protocols.');
    $('status').textContent=shown.length+' of '+nodes.size+' public protocols · '+visibleEdges.length+' curated links'+(visibleSimilarities.length?' · '+visibleSimilarities.length+' similarity links':'')+'. Lineage curation is partial.';
    $('export').disabled=!shown.length;fit();
  }
  $('find').addEventListener('click',()=>{
    if(!data)return;const q=$('search').value.trim().toLowerCase();
    const item=[...nodes.values()].find(n=>n.title.toLowerCase()===q||n.id===q)||[...nodes.values()].find(n=>n.title.toLowerCase().includes(q));
    if(item&&q)choose(item.id,true);else $('status').textContent='No matching public protocol.';
  });
  $('search').addEventListener('keydown',event=>{if(event.key==='Enter'){$('find').click();event.preventDefault();}});
  for(const id of ['documented','inferred','similar','isolated'])$(id).addEventListener('change',draw);
  $('focus').addEventListener('click',()=>{focused=true;draw();});$('all').addEventListener('click',()=>{focused=false;draw();});
  $('fit').addEventListener('click',fit);$('plus').addEventListener('click',()=>zoom(1.25));$('minus').addEventListener('click',()=>zoom(.8));
  $('canvas').addEventListener('wheel',event=>{if(!data)return;event.preventDefault();const box=$('canvas').getBoundingClientRect();zoom(event.deltaY<0?1.12:1/1.12,event.clientX-box.left,event.clientY-box.top);},{passive:false});
  $('canvas').addEventListener('pointerdown',event=>{if(event.target.closest('.lineage-node,.lineage-edge'))return;drag={x:event.clientX,y:event.clientY,pan:{...pan}};$('canvas').setPointerCapture(event.pointerId);});
  $('canvas').addEventListener('pointermove',event=>{if(drag){pan={x:drag.pan.x+event.clientX-drag.x,y:drag.pan.y+event.clientY-drag.y};transform();}});
  $('canvas').addEventListener('pointerup',()=>{drag=null;});$('canvas').addEventListener('pointercancel',()=>{drag=null;});
  $('canvas').addEventListener('keydown',event=>{
    if(event.target!==$('canvas'))return;
    if(event.key==='+'||event.key==='=')zoom(1.25);else if(event.key==='-')zoom(.8);
    else if(['ArrowLeft','ArrowRight','ArrowUp','ArrowDown'].includes(event.key)){event.preventDefault();pan.x+=event.key==='ArrowLeft'?30:event.key==='ArrowRight'?-30:0;pan.y+=event.key==='ArrowUp'?30:event.key==='ArrowDown'?-30:0;transform();}
  });
  new ResizeObserver(()=>{if(svg)fit();}).observe($('canvas'));
  $('export').addEventListener('click',()=>{
    if(!svg)return;const copy=svg.cloneNode(true);copy.setAttribute('viewBox','0 0 '+layout.width+' '+layout.height);copy.setAttribute('width',layout.width);copy.setAttribute('height',layout.height);copy.firstChild.removeAttribute('transform');
    const url=URL.createObjectURL(new Blob([new XMLSerializer().serializeToString(copy)],{type:'image/svg+xml'})),a=node('a');a.href=url;a.download='protocol-relationships.svg';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
  });
  fetch('assets/protocol_lineage.json').then(response=>{if(!response.ok)throw new Error('Could not load relationships. Reload to retry.');return response.json();}).then(value=>{
    if(value.schema_version!==1||!Array.isArray(value.nodes)||!Array.isArray(value.relations)||!Array.isArray(value.similarities))throw new Error('Unsupported relationship data.');
    data=value;nodes=new Map(data.nodes.map(n=>[n.id,n]));
    for(const n of [...nodes.values()].sort((a,b)=>a.title.localeCompare(b.title))){const o=node('option');o.value=n.title;$('protocols').append(o);}
    const id=new URLSearchParams(location.search).get('protocol');if(id&&nodes.has(id))choose(id,true);else draw();
  }).catch(error=>{$('status').textContent=error.message;});
}());

/* Full, wrapped sequence viewer. Scene supplies aligned strand columns and polarity. */
(function () {
  'use strict';
  window.RestrictionView={create({$,node,download}){
    const ns='http://www.w3.org/2000/svg',model=JSON.parse($('duplex-model').textContent),port=$('duplex'),popup=$('hover');
    const cell=10,gap=7,left=48,rowHeight=138,topY=72,bottomY=98,midY=81;
    let result,groups=[],byRow=[],columns=80,content,active=null,pinned=false,frame,drawn=new Map(),lastPointer={x:0,y:0},dismissTimer;
    const mod=(p,n)=>((p%n)+n)%n;
    function svg(w,h,label){const s=document.createElementNS(ns,'svg');for(const [k,v]of Object.entries({xmlns:ns,viewBox:`0 0 ${w} ${h}`,width:w,height:h,role:'img','aria-label':label}))s.setAttribute(k,v);return s;}
    function add(s,tag,attrs,text){const e=document.createElementNS(ns,tag);for(const [k,v]of Object.entries(attrs))e.setAttribute(k,v);if(text!==undefined)e.textContent=text;s.append(e);return e;}
    function text(s,x,y,value,attrs={}){return add(s,'text',{x,y,fill:'#25344b','font-family':'sans-serif','font-size':12,...attrs},value);}
    const baseX=i=>left+i*cell+Math.floor(i/10)*gap;
    const boundaryX=i=>baseX(i)-(i>0&&i%10===0?gap/2:0);
    function coord(p){return p===null?'—':String(result.circular&&p===0?result.length:p);}
    function title(group){return group.rows.map(r=>r.enzyme.name).join(', ');}
    function cutDescription(group){return group.site.cuts.length?group.site.cuts.map(c=>`Top / bottom: ${coord(c.top)} / ${coord(c.bottom)}${c.valid?'':' · incomplete at linear end'}`).join('; '):'Cut geometry unassigned';}
    function ranges(group){const {start,length}=group.site,n=result.length;return start+length>n?[[start,n],[0,start+length-n]]:[[start,start+length]];}
    function clear(){active=null;pinned=false;popup.hidden=true;$('unpin').hidden=true;redrawActive();}
    function deferClear(){clearTimeout(dismissTimer);if(!pinned)dismissTimer=setTimeout(clear,120);}
    function show(group,event,candidates){
      clearTimeout(dismissTimer);if(pinned&&active!==group)return;
      active=group;lastPointer={x:event.clientX||lastPointer.x,y:event.clientY||lastPointer.y};
      popup.replaceChildren(node('strong',title(group)),node('div',`Recognition: ${group.site.start+1}–${mod(group.site.start+group.site.length-1,result.length)+1} · ${group.site.strand==='both'?'palindromic':group.site.strand+' strand'}${group.site.certain?'':' · possible'}`),node('div',cutDescription(group)));
      if(candidates&&candidates.length>1){const list=node('div',undefined,'restriction-hover-choices');list.append(node('small','Overlapping motifs'));for(const g of candidates){const button=node('button',title(g));button.type='button';button.addEventListener('pointerenter',e=>{if(!pinned)activateOnly(g,e);});button.addEventListener('focus',e=>activateOnly(g,e));button.addEventListener('click',e=>{pinned=false;show(g,e);pin();});list.append(button);}popup.append(list);}
      popup.hidden=false;placePopup();redrawActive();
    }
    function activateOnly(group){active=group;popup.querySelector('strong').textContent=title(group);popup.children[2].textContent=cutDescription(group);popup.children[1].textContent=`Recognition: ${group.site.start+1} · ${group.site.certain?'certain':'possible'}`;redrawActive();}
    function placePopup(){const box=popup.getBoundingClientRect();popup.style.left=Math.max(8,Math.min(window.innerWidth-box.width-12,lastPointer.x+12))+'px';const above=lastPointer.y-box.height-18;popup.style.top=Math.max(8,Math.min(window.innerHeight-box.height-12,above>=8?above:lastPointer.y+32))+'px';}
    function pin(){pinned=true;$('unpin').hidden=false;}
    function intersects(group,a,b){return ranges(group).some(([s,e])=>s<b&&e>a);}
    function paintActive(s,row){
      s.querySelector('[data-active-layer]')?.remove();if(!active)return;
      const a=row*columns,b=Math.min(a+columns,result.length),layer=add(s,'g',{'data-active-layer':'true','pointer-events':'none','data-active-enzyme':title(active)});
      for(const [start,end]of ranges(active)){const lo=Math.max(a,start),hi=Math.min(b,end);if(hi>lo)add(layer,'rect',{x:boundaryX(lo-a),y:topY-15,width:boundaryX(hi-a)-boundaryX(lo-a),height:bottomY-topY+20,fill:'#f8cc49',opacity:.26});}
      for(const c of active.site.cuts){
        const raw={top:c.rawTop,bottom:c.rawBottom};
        // Preserve the unwrapped cut path at the circular origin; modulo coordinates
        // alone lose which end of the displayed molecule the path continues through.
        const points=Object.values(raw).filter(p=>p!==null),lo=Math.min(...points),hi=Math.max(...points);
        const shifts=result.circular?[Math.floor((a-hi)/result.length),Math.ceil((b-lo)/result.length)]:[0,0];
        for(let turn=shifts[0];turn<=shifts[1];turn++){
          const shift=result.circular?turn*result.length:0,pos={};
          for(const [strand,y]of [['top',topY],['bottom',bottomY]]){
            const p=raw[strand];if(p===null||p+shift<a||p+shift>b)continue;
            const x=boundaryX(p+shift-a);pos[strand]=x;
            add(layer,'path',{d:`M ${x} ${strand==='top'?y-15:midY} V ${strand==='top'?midY:y+5}`,fill:'none',stroke:c.valid?'#c83445':'#818894','stroke-width':2,'stroke-dasharray':active.site.certain&&c.valid?'':'4 3','data-cut-strand':strand,'data-cut-boundary':c[strand]});
          }
          // Vertical strand marks join through a horizontal middle segment.
          if(raw.top!==null&&raw.bottom!==null){
            const start=Math.max(a,Math.min(raw.top,raw.bottom)+shift),end=Math.min(b,Math.max(raw.top,raw.bottom)+shift);
            if(end>=start)add(layer,'path',{d:`M ${boundaryX(start-a)} ${midY} H ${boundaryX(end-a)}`,fill:'none',stroke:c.valid?'#c83445':'#818894','stroke-width':2,'stroke-dasharray':active.site.certain&&c.valid?'':'4 3','data-cut-connector':'horizontal'});
          }
        }
      }
    }

    function redrawActive(){for(const [row,s]of drawn)paintActive(s,row);}
    function rowSVG(row,interactive=true){
      const a=row*columns,b=Math.min(a+columns,result.length),len=b-a,w=baseX(columns)+40;
      const s=svg(w,rowHeight,`DNA bases ${a+1}–${b}`);s.setAttribute('data-sequence-row',row);
      add(s,'rect',{width:w,height:rowHeight,fill:'#fff'});
      const matching=byRow[row]||[],laneEnds=[-Infinity,-Infinity],omitted=[];
      for(const group of matching){
        for(const [start,end]of ranges(group)){
          const lo=Math.max(a,start),hi=Math.min(b,end);if(hi<=lo)continue;
          add(s,'rect',{x:boundaryX(lo-a),y:topY-15,width:boundaryX(hi-a)-boundaryX(lo-a),height:bottomY-topY+20,fill:'#257c86',opacity:.065,'pointer-events':'none'});
          const name=group.rows[0].enzyme.name+(group.rows.length>1?` +${group.rows.length-1}`:''),labelWidth=name.length*6.6+12;
          const x=Math.min(boundaryX(lo-a),w-labelWidth-10),lane=laneEnds.findIndex(e=>e+8<x);
          if(lane<0){omitted.push(group);continue;}laneEnds[lane]=x+labelWidth;
          const g=add(s,'g',{'data-enzyme':group.rows[0].enzyme.name,'data-site':group.site.start+1,'data-motif':group.id,tabindex:interactive?0:-1,role:'button','aria-label':title(group)+', recognition starts at '+(group.site.start+1)});
          add(g,'rect',{x,y:7+lane*20,width:labelWidth,height:18,rx:3,fill:'#eef5f5'});
          text(g,x+5,20+lane*20,name,{fill:'#226c78','font-size':11});
          add(g,'path',{d:`M ${boundaryX(lo-a)} ${49} H ${boundaryX(hi-a)}`,stroke:'#8cb7bd','stroke-width':2,'stroke-dasharray':group.site.certain?'':'3 2'});
          if(interactive){g.addEventListener('pointerenter',e=>show(group,e));g.addEventListener('pointerleave',deferClear);g.addEventListener('focus',()=>{const box=g.getBoundingClientRect();show(group,{clientX:box.left,clientY:box.bottom});});g.addEventListener('blur',deferClear);g.addEventListener('click',e=>{if(pinned&&active===group){clear();return;}pinned=false;show(group,e);pin();});g.addEventListener('keydown',e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();g.dispatchEvent(new MouseEvent('click',{bubbles:true}));}});}
        }
      }
      const sequence=result.sequence.slice(a,b),strands=[sequence,window.SequenceAnnotation.revcomp(sequence).split('').reverse().join('')];
      model.strands.forEach((strand,j)=>{
        const y=j?bottomY:topY,xOffset=(strand.start-model.strands[0].start)*cell;
        text(s,12,y,j?'3′':'5′',{'font-size':11,fill:'#7c8798'});text(s,baseX(len)+12,y,j?'5′':'3′',{'font-size':11,fill:'#7c8798'});
        for(let i=0;i<len;i++)text(s,baseX(i)+cell/2+xOffset,y,strands[j][i],{'font-family':'monospace','font-size':14,'text-anchor':'middle','data-strand':j?'bottom':'top','data-base':a+i+1});
      });
      for(let i=0;i<len;i++)if(i===0||(a+i+1)%10===0)text(s,baseX(i)+cell/2,121,String(a+i+1),{'font-size':10,'text-anchor':'middle',fill:'#8491a3'});
      if(result.circular&&(a===0||b===result.length))text(s,a===0?left:baseX(len)-5,135,a===0?'↻ origin':'origin ↻',{'font-size':10,fill:'#8491a3','text-anchor':a===0?'start':'end'});
      if(omitted.length)text(s,w-12,45,`+${new Set(omitted).size}`,{'font-size':10,'text-anchor':'end',fill:'#8491a3'});
      if(interactive){
        s.addEventListener('pointermove',event=>{const base=event.target.closest('[data-base]');if(!base||pinned)return;const p=Number(base.getAttribute('data-base'))-1,candidates=matching.filter(g=>ranges(g).some(([a,b])=>p>=a&&p<b));if(!candidates.length){deferClear();return;}if(active!==candidates[0]||popup.hidden)show(candidates[0],event,candidates);else clearTimeout(dismissTimer);});
        s.addEventListener('pointerleave',deferClear);
      }
      paintActive(s,row);return s;
    }
    function draw(){
      if(!result||!content)return;
      const count=Math.ceil(result.length/columns),first=Math.max(0,Math.floor(port.scrollTop/rowHeight)-2),last=Math.min(count,Math.ceil((port.scrollTop+port.clientHeight)/rowHeight)+2);
      const keep=new Set();for(let i=first;i<last;i++){keep.add(i);if(!drawn.has(i)){const s=rowSVG(i);s.style.position='absolute';s.style.top=i*rowHeight+'px';content.append(s);drawn.set(i,s);}}
      for(const [i,s]of drawn)if(!keep.has(i)){s.remove();drawn.delete(i);}
    }
    function layout(){
      if(!result)return;const next=Math.max(20,Math.min(120,Math.floor((port.clientWidth-left-45)/(cell*10+gap))*10));
      const base=content?Math.floor(port.scrollTop/rowHeight)*columns:0;columns=next;
      byRow=Array.from({length:Math.ceil(result.length/columns)},()=>[]);
      for(const g of groups){const assigned=new Set();for(const [a,b]of ranges(g))for(let i=Math.floor(a/columns);i<=Math.floor((b-1)/columns);i++)if(!assigned.has(i)){byRow[i].push(g);assigned.add(i);}}
      for(const row of byRow)row.sort((a,b)=>a.site.start-b.site.start||title(a).localeCompare(title(b)));
      content=node('div',undefined,'restriction-sequence-content');content.style.height=byRow.length*rowHeight+'px';port.style.height=`min(${Math.min(620,Math.max(150,byRow.length*rowHeight))}px,70vh)`;content.style.minWidth=(baseX(columns)+40)+'px';port.replaceChildren(content);drawn.clear();port.scrollTop=Math.floor(base/columns)*rowHeight;draw();
    }
    port.addEventListener('scroll',()=>{if(!pinned)clear();cancelAnimationFrame(frame);frame=requestAnimationFrame(draw);},{passive:true});
    popup.addEventListener('pointerenter',()=>clearTimeout(dismissTimer));popup.addEventListener('pointerleave',deferClear);
    $('unpin').addEventListener('click',clear);port.addEventListener('keydown',e=>{if(e.key==='Escape')clear();});
    document.addEventListener('pointerdown',e=>{if(pinned&&!port.contains(e.target)&&!popup.contains(e.target))clear();});
    new ResizeObserver(()=>{if(result)layout();}).observe(port);
    $('export-view').addEventListener('click',()=>{
      if(!result)return;const first=Math.floor(port.scrollTop/rowHeight),last=Math.min(byRow.length,Math.ceil((port.scrollTop+port.clientHeight)/rowHeight)),w=baseX(columns)+40;
      const out=svg(w,(last-first)*rowHeight,'Visible double-stranded sequence and selected restriction cuts');
      for(let row=first;row<last;row++){const s=rowSVG(row,false);s.setAttribute('y',(row-first)*rowHeight);out.append(s);}
      download(new XMLSerializer().serializeToString(out),'restriction-dsdna.svg','image/svg+xml');
    });
    return {clear,render(value,rows){
      const fresh=result!==value;result=value;clear();const merged=new Map();
      for(const row of rows)for(const site of row.sites){const key=JSON.stringify([site.start,site.length,site.strand,site.certain,site.cuts]);let g=merged.get(key);if(!g){g={site,rows:[],id:merged.size};merged.set(key,g);}g.rows.push(row);}
      groups=[...merged.values()];if(fresh){content=null;port.scrollTop=0;}
      $('sequence-info').textContent=`${result.length.toLocaleString()} bp · ${result.circular?'circular':'linear'}`;
      $('view-note').textContent=rows.some(r=>r.positionsTruncated)?'Retained sites shown; counts include every match.':'';
      layout();
    }};
  }};
}());

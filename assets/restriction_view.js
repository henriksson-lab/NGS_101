/* Interactive coordinate viewer. Duplex columns/polarity come from a build-time Scene. */
(function () {
  'use strict';
  window.RestrictionView={create({$,node,download}){
    const ns='http://www.w3.org/2000/svg',model=JSON.parse($('duplex-model').textContent);
    let result,rows=[],entries=[],selected=null,start=0,viewSVG;
    const mod=(p,n)=>((p%n)+n)%n;
    function svg(width,height,label){const s=document.createElementNS(ns,'svg');for(const [k,v]of Object.entries({xmlns:ns,viewBox:`0 0 ${width} ${height}`,width,height,role:'img','aria-label':label}))s.setAttribute(k,v);return s;}
    function add(s,tag,attrs,text){const e=document.createElementNS(ns,tag);for(const [k,v]of Object.entries(attrs))e.setAttribute(k,v);if(text!==undefined)e.textContent=text;s.append(e);return e;}
    function text(s,x,y,value,attrs={}){return add(s,'text',{x,y,fill:'#243d52','font-family':'sans-serif','font-size':12,...attrs},value);}
    function line(s,x1,y1,x2,y2,color,width=1.5,dash=''){return add(s,'line',{x1,y1,x2,y2,stroke:color,'stroke-width':width,'stroke-dasharray':dash});}
    function interactive(e,entry){
      e.style.cursor='pointer';e.setAttribute('tabindex','0');e.setAttribute('role','button');e.setAttribute('aria-label',label(entry));
      const title=document.createElementNS(ns,'title');title.textContent=label(entry);e.append(title);
      const activate=()=>choose(entry);e.addEventListener('click',event=>{event.stopPropagation();activate();});e.addEventListener('keydown',event=>{if(event.key==='Enter'||event.key===' '){event.preventDefault();activate();}});
    }
    function label(e){return `${e.row.enzyme.name} · site ${e.site.start+1} · ${e.site.strand==='both'?'palindromic':e.site.strand+' strand'}${e.site.certain?'':' · possible'}`;}
    function coordinate(p){return p===null?'—':String(result.circular&&p===0?result.length:p);}
    function current(){const name=$('view-enzyme').value;return name?entries.filter(e=>e.row.enzyme.name===name):entries;}
    function siteOptions(){
      const matches=current(),select=$('view-site');select.replaceChildren(node('option','Select a site'));
      select.firstChild.value='';if(matches.length>1000)select.firstChild.textContent=`Select a site (first 1,000 of ${matches.length})`;
      const retained=matches.slice(0,1000);if(selected&&matches.includes(selected)&&!retained.includes(selected))retained.push(selected);
      for(const e of retained){const option=node('option',label(e));option.value=String(e.id);select.append(option);}
      select.value=selected?String(selected.id):'';
      $('site-prev').disabled=$('site-next').disabled=!matches.length;
    }
    function choose(entry){selected=entry;$('view-enzyme').value=entry.row.enzyme.name;siteOptions();start=Math.max(0,entry.site.start-Math.floor(Number($('view-width').value)/4));draw();}
    function color(e){return selected&&e.row.enzyme.name===selected.row.enzyme.name?'#ac2843':'#217d85';}
    function overlap(e,width){const n=result.length,p=result.circular?mod(e.site.start-start,n):e.site.start-start;return p<width||result.circular&&p+e.site.length>n||!result.circular&&p<0&&p+e.site.length>0;}
    function intervals(e,width){
      let p=e.site.start-start;if(result.circular)p=mod(p,result.length);
      const out=[];for(const q of result.circular?[p,p-result.length]:[p]){const a=Math.max(0,q),b=Math.min(width,q+e.site.length);if(b>a)out.push([a,b]);}return out;
    }
    function overview(){
      const width=1000,s=svg(width,92,'Recognition and strand-specific cut overview'),n=result.length,scale=900/n,x=p=>50+p*scale;
      add(s,'rect',{width,height:92,fill:'#fff'});line(s,50,48,950,48,'#243d52');line(s,50,60,950,60,'#243d52');
      text(s,12,51,'5′');text(s,12,64,'3′');text(s,962,51,'3′');text(s,962,64,'5′');
      const matches=current(),cap=2000;
      for(const e of matches.slice(0,cap)){
        const g=add(s,'g',{});interactive(g,e);const p=e.site.start,b=Math.min(n,p+e.site.length);
        add(g,'rect',{x:x(p),y:37,width:Math.max(3,(b-p)*scale),height:34,fill:color(e),opacity:.18});
        if(p+e.site.length>n)add(g,'rect',{x:50,y:37,width:Math.max(3,(p+e.site.length-n)*scale),height:34,fill:color(e),opacity:.18});
        for(const c of e.site.cuts.filter(c=>c.valid))for(const [strand,y]of [['top',48],['bottom',60]])if(c[strand]!==null)line(g,x(c[strand]),y-5,x(c[strand]),y+5,'#ac2843',1.4,e.site.certain?'':'3 2');
      }
      const size=Math.min(Number($('view-width').value),n),a=Math.min(size,n-start);
      add(s,'rect',{x:x(start),y:29,width:Math.max(2,a*scale),height:48,fill:'none',stroke:'#243d52','stroke-width':2,'pointer-events':'none'});
      if(a<size)add(s,'rect',{x:50,y:29,width:(size-a)*scale,height:48,fill:'none',stroke:'#243d52','stroke-width':2,'pointer-events':'none'});
      text(s,50,17,`${n.toLocaleString()} bp · ${result.circular?'circular, unrolled at origin':'linear'} · click a site or the backbone`);
      text(s,50,89,`1${matches.length>cap?` · first ${cap} of ${matches.length} retained matches on overview`:''}`);text(s,950,89,String(n),{'text-anchor':'end'});
      s.addEventListener('click',event=>{const box=s.getBoundingClientRect();start=Math.max(0,Math.min(n-1,Math.round(((event.clientX-box.left)*width/box.width-50)/scale)));draw();});
      $('overview').replaceChildren(s);
    }
    function draw(){
      if(!result)return;const n=result.length,width=Math.min(Number($('view-width').value),n);
      start=Math.max(0,Math.min(result.circular?n-1:Math.max(0,n-width),start));
      for(const name of ['view-start','view-slider']){$(name).max=result.circular?n:Math.max(1,n-width+1);$(name).value=start+1;}
      const cell=10,left=55,w=left+width*cell+55,ty=58,by=84;
      const matching=current().filter(e=>overlap(e,width)||e.site.cuts.some(c=>[c.top,c.bottom].some(p=>p!==null&&(result.circular?mod(p-start,n):p-start)>=0&&(result.circular?mod(p-start,n):p-start)<=width)));
      matching.sort((a,b)=>(b===selected)-(a===selected)||a.site.start-b.site.start||a.row.enzyme.name.localeCompare(b.row.enzyme.name));
      const tracks=matching.slice(0,18),s=svg(w,145+tracks.length*23,'Double-stranded DNA with recognition spans and strand-specific cuts');viewSVG=s;
      add(s,'rect',{width:w,height:145+tracks.length*23,fill:'#fff'});
      text(s,12,ty,'5′');text(s,w-30,ty,'3′');text(s,12,by,'3′');text(s,w-30,by,'5′');
      const sequence=Array.from({length:width},(_,i)=>result.sequence[mod(start+i,n)]).join('');
      // Same Scene column for both strands; the existing shared IUPAC helper supplies pairing.
      const strands=[sequence,window.SequenceAnnotation.revcomp(sequence).split('').reverse().join('')];
      model.strands.forEach((strand,j)=>{const y=j?by:ty;const x=left+(strand.start-model.strands[0].start)*cell;line(s,x,y+4,x+width*cell,y+4,'#667887',.6);
        for(let i=0;i<width;i++)text(s,x+(i+.5)*cell,y,strands[j][i],{'font-family':'monospace','font-size':14,'text-anchor':'middle','data-strand':j?'bottom':'top','data-base':mod(start+i,n)+1});});
      for(let i=0;i<width;i++){const p=mod(start+i,n);if(i===0||p%10===9||p===0){text(s,left+(i+.5)*cell,29,String(p+1),{'text-anchor':'middle','font-size':10});if(p===0&&i)line(s,left+i*cell,33,left+i*cell,95,'#576575',1,'3 3');}}
      for(let index=0;index<tracks.length;index++){
        const e=tracks[index],g=add(s,'g',{'data-enzyme':e.row.enzyme.name,'data-site':e.site.start+1}),y=128+index*23,c=color(e);interactive(g,e);
        for(const [a,b]of intervals(e,width)){
          add(g,'rect',{x:left+a*cell,y:ty-17,width:(b-a)*cell,height:by-ty+24,fill:c,opacity:e===selected?.23:.07,'pointer-events':'none'});
          add(g,'rect',{x:left+a*cell,y:y-9,width:Math.max(2,(b-a)*cell),height:10,fill:c,opacity:.65,'stroke-dasharray':e.site.certain?'':'3 2',stroke:c});
        }
        text(g,left,y+13,label(e),{'font-size':10,fill:c});
        for(const cut of e.site.cuts){
          const positions={};for(const [strand,ybase]of [['top',ty],['bottom',by]]){
            const p=cut[strand];if(p===null)continue;const offset=result.circular?mod(p-start,n):p-start;
            if(offset<0||offset>width)continue;const x=left+offset*cell;positions[strand]=x;
            const tick=line(g,x,ybase-15,x,ybase+7,cut.valid?'#ac2843':'#8a8a8a',2,e.site.certain&&cut.valid?'':'3 2');tick.setAttribute('data-cut-strand',strand);tick.setAttribute('data-cut-boundary',p);
          }
          if(positions.top!==undefined&&positions.bottom!==undefined)line(g,positions.top,ty+7,positions.bottom,by-15,cut.valid?'#ac2843':'#8a8a8a',1.2,e.site.certain&&cut.valid?'':'3 2');
        }
      }
      $('duplex').replaceChildren(s);overview();
      const capped=rows.some(r=>r.positionsTruncated);
      $('view-note').textContent=`Recognition spans: teal · cuts: red · dashed: possible or incomplete at a linear end.${matching.length>tracks.length?` Showing ${tracks.length} of ${matching.length} overlapping retained sites; select an enzyme to isolate it.`:''}${capped?' Retained positions only; full counts remain in the table.':''}`;
      const box=$('selection');box.replaceChildren();
      if(selected){const e=selected,p=node('p',label(e));p.append(document.createTextNode(' · '+(e.site.cuts.length?e.site.cuts.map(c=>`top / bottom cuts ${coordinate(c.top)} / ${coordinate(c.bottom)}${c.valid?'':' (incomplete at linear end)'}`).join('; '):'cut geometry unassigned')));box.append(p);}
      $('export-view').disabled=false;
    }
    $('view-enzyme').addEventListener('change',()=>{selected=null;siteOptions();const first=current()[0];if(first){selected=first;start=Math.max(0,first.site.start-10);$('view-site').value=String(first.id);}draw();});
    $('view-site').addEventListener('change',()=>{const e=entries.find(e=>String(e.id)===$('view-site').value);if(e)choose(e);});
    for(const [name,delta]of [['site-prev',-1],['site-next',1]])$(name).addEventListener('click',()=>{const list=current();if(!list.length)return;const index=list.indexOf(selected);choose(list[mod(index+delta,list.length)]);});
    for(const name of ['view-start','view-slider'])$(name).addEventListener(name==='view-slider'?'input':'change',()=>{const v=Number($(name).value);if(Number.isFinite(v))start=Math.round(v)-1;draw();});
    $('view-width').addEventListener('change',draw);
    for(const [name,delta]of [['pan-prev',-1],['pan-next',1]])$(name).addEventListener('click',()=>{start+=delta*Math.max(1,Math.floor(Number($('view-width').value)/2));if(result.circular)start=mod(start,result.length);draw();});
    $('export-view').addEventListener('click',()=>{if(viewSVG)download(new XMLSerializer().serializeToString(viewSVG),'restriction-dsdna.svg','image/svg+xml');});
    return {render(value,filtered){const fresh=result!==value;result=value;rows=filtered;entries=[];for(const row of rows)for(const site of row.sites)entries.push({row,site,id:entries.length});
      const previous=$('view-enzyme').value,select=$('view-enzyme');select.replaceChildren(node('option','All filtered enzymes'));select.firstChild.value='';
      for(const row of rows){const o=node('option',row.enzyme.name);o.value=row.enzyme.name;select.append(o);}select.value=rows.some(r=>r.enzyme.name===previous)?previous:'';
      selected=null;if(fresh)start=0;siteOptions();draw();}};
  }};
}());

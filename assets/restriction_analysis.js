/* Sequence-based restriction maps. Coordinates are zero-based inter-base boundaries. */
(function (root, factory) {
  'use strict';
  const api = factory(typeof module === 'object' && module.exports ? require('./sequence_annotation.js') : root.SequenceAnnotation);
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.RestrictionAnalysis = api;
}(typeof globalThis !== 'undefined' ? globalThis : this, function (sequenceTools) {
  'use strict';
  const MASK = {A:1,C:2,G:4,T:8,R:5,Y:10,S:6,W:9,K:12,M:3,B:14,D:13,H:11,V:7,N:15};
  const STATES = ['unknown','unaffected','blocked','impaired','overlap_blocked','some_overlap_blocked','overlap_impaired','some_overlap_impaired','requires'];
  const LABELS = {unknown:'Unknown',unaffected:'Not sensitive',blocked:'Blocked',impaired:'Impaired',overlap_blocked:'Blocked if overlapping',some_overlap_blocked:'Blocked for some overlaps',overlap_impaired:'Impaired if overlapping',some_overlap_impaired:'Impaired for some overlaps',requires:'Requires methylation'};
  const CONTEXTS = {cpg:'CG',dam:'GATC',dcm:'CCWGG'};
  const mod = (a,n) => ((a % n)+n)%n;
  function normalize(text) {
    const sequence = sequenceTools.normalizeInput(text);
    if (sequence.includes('U')) throw new Error('Use DNA bases: replace U with T.');
    return sequence;
  }
  function validateEnzyme(enzyme) {
    if (!enzyme || typeof enzyme.name !== 'string' || !/^[A-Za-z0-9_.+\-]{1,80}$/.test(enzyme.name)) throw new Error('Invalid enzyme name.');
    if (!Array.isArray(enzyme.motifs) || !enzyme.motifs.length || enzyme.motifs.some(s => typeof s !== 'string' || !/^[ACGTRYSWKMBDHVN]{1,100}$/.test(s))) throw new Error('Invalid recognition motif for '+enzyme.name+'.');
    if (!Array.isArray(enzyme.cuts) || enzyme.cuts.length>2 || enzyme.cuts.some(pair => !Array.isArray(pair) || pair.length!==2 || pair.every(x=>x===null) || pair.some(x=>x!==null && (!Number.isInteger(x) || Math.abs(x)>1000)))) throw new Error('Invalid cut offsets for '+enzyme.name+'.');
    if (enzyme.motifs.length>1 && enzyme.cuts.length) throw new Error('Alternative motifs with shared cut offsets are not supported.');
    for (const value of Object.values(enzyme.methylation || {})) if (!STATES.includes(value)) throw new Error('Invalid methylation annotation.');
    return enzyme;
  }
  function parseCustom(text, existing) {
    const used = new Set((existing || []).map(e=>e.name.toLowerCase())), enzymes=[];
    for (const raw of String(text || '').split(/\r?\n/)) {
      const line=raw.trim(); if (!line || line.startsWith('#')) continue;
      const fields=line.split(/[\s,;]+/);
      if (![4,7].includes(fields.length)) throw new Error('Custom format: Name Motif Top Bottom [CpG Dam Dcm].');
      const [name,site,top,bottom]=fields;
      if (used.has(name.toLowerCase())) throw new Error('Duplicate enzyme name: '+name);
      const offset=s=>s==='-' || s==='null' ? null : /^-?\d+$/.test(s) ? Number(s) : NaN;
      const pair=[offset(top),offset(bottom)];
      const e={name,motifs:[site.toUpperCase()],cuts:pair.every(x=>x===null)?[]:[pair],suppliers:[],custom:true,
               methylation:{cpg:fields[4]||'unknown',dam:fields[5]||'unknown',dcm:fields[6]||'unknown'}};
      validateEnzyme(e);used.add(name.toLowerCase());enzymes.push(e);
      if (enzymes.length>200) throw new Error('Use at most 200 custom enzymes.');
    }
    return enzymes;
  }
  function relation(query,motif) {
    let certain=true;
    for (let i=0;i<motif.length;i++) {
      const q=MASK[query[i]],m=MASK[motif[i]];
      if (!(q&m)) return 0;
      if ((q&m)!==q) certain=false;
    }
    return certain?2:1;
  }
  const patterns=new Map();
  function regex(motif) {
    if (!patterns.has(motif)) patterns.set(motif,motif.split('').map(b=>'['+Object.keys(MASK).filter(q=>MASK[q]&MASK[b]).join('')+']').join(''));
    return new RegExp('(?=('+patterns.get(motif)+'))','g');
  }
  function contextAt(sequence,start,length,circular) {
    const answer={},n=sequence.length;
    for (const [name,motif] of Object.entries(CONTEXTS)) {
      let best=0;
      if (motif.length<=n) for (let p=start-motif.length+1;p<start+length;p++) {
        if (!circular && (p<0 || p+motif.length>n)) continue;
        let sub='';for(let j=0;j<motif.length;j++)sub+=sequence[mod(p+j,n)];
        best=Math.max(best,relation(sub,motif));
        if(best===2)break;
      }
      answer[name]=best; // 0 absent,1 possible,2 certain motif overlap.
    }
    return answer;
  }
  function enzymeKind(enzyme) {
    if (!enzyme.cuts.length) return 'unknown';
    return enzyme.cuts.some(p=>p.includes(null))?'nick':'double';
  }
  function endType(enzyme) {
    const kind=enzymeKind(enzyme);if(kind!=='double')return kind;
    const kinds=new Set(enzyme.cuts.map(([t,b])=>t===b?'blunt':b>t?'5prime':'3prime'));
    return kinds.size===1?[...kinds][0]:'mixed';
  }
  function cutsAt(enzyme,start,length,strand,n,circular) {
    return enzyme.cuts.map(([top,bottom])=>{
      let t=strand==='-' ? bottom===null?null:start+length-bottom : top===null?null:start+top;
      let b=strand==='-' ? top===null?null:start+length-top : bottom===null?null:start+bottom;
      const valid=circular || [t,b].every(x=>x===null || (x>0 && x<n));
      const result={top:circular && t!==null?mod(t,n):t,bottom:circular && b!==null?mod(b,n):b,valid,rawTop:t,rawBottom:b,
                    end:t===null || b===null?'nick':t===b?'blunt':b>t?'5prime':'3prime',overhang:t===null || b===null?null:Math.abs(t-b)};
      return result;
    });
  }
  function fragments(length,boundaries,circular,limit) {
    const points=[...boundaries].sort((a,b)=>a-b);let values=[];
    if (!points.length) values=[length];
    else if (circular) values=points.map((p,i)=>i+1<points.length?points[i+1]-p:length-p+points[0]);
    else {let p=0;for(const x of points){values.push(x-p);p=x;}values.push(length-p);}
    return {lengths:values.slice(0,limit),total:values.length,truncated:values.length>limit};
  }
  function analyze(text,enzymes,options) {
    options=options||{};const sequence=normalize(text),n=sequence.length,circular=Boolean(options.circular);
    const limit=options.positionLimit===undefined?1000:options.positionLimit;
    if (!Number.isInteger(limit) || limit<1 || limit>10000) throw new Error('Invalid position display limit.');
    if (!Array.isArray(enzymes) || enzymes.length>2000) throw new Error('Invalid enzyme catalogue.');
    const matchLimit=options.matchLimit===undefined?2000000:options.matchLimit;
    if(!Number.isInteger(matchLimit)||matchLimit<1||matchLimit>2000000)throw new Error('Invalid match limit.');
    const deadline=Date.now()+15000;let work=0,ambiguousWork=0;
    const results=[],names=new Set();
    for (let index=0;index<enzymes.length;index++) {
      const enzyme=validateEnzyme(enzymes[index]);
      if(names.has(enzyme.name.toLowerCase()))throw new Error('Duplicate enzyme: '+enzyme.name);names.add(enzyme.name.toLowerCase());
      const matches=new Map();
      for (const motif of enzyme.motifs) {
        const reverse=sequenceTools.revcomp(motif),variants=reverse===motif?[[motif,'both']]:[[motif,'+'],[reverse,'-']];
        if(motif.length>n)continue;
        const hay=circular?sequence+sequence.slice(0,motif.length-1):sequence;
        for (const [pattern,strand] of variants) {
          const rx=regex(pattern);let found;
          while((found=rx.exec(hay))!==null) {
            const p=found.index;rx.lastIndex=p+1;if(p>=n)break;
            const certain=relation(found[1],pattern)===2,key=p+':'+pattern.length+':'+strand;
            if(++work>matchLimit || (!certain && ++ambiguousWork>100000) || (work%256===0 && Date.now()>deadline))throw new Error('Too many ambiguous or repetitive matches. Analyze a smaller or less ambiguous sequence; no partial counts were returned.');
            if(!matches.has(key) || certain)matches.set(key,{start:p,length:pattern.length,strand,certain});
          }
        }
      }
      const certain=new Set(),possible=new Set(),topBoundaries=new Set(),recognitionCertain=new Set(),recognitionPossible=new Set();
      const contexts={cpg:{certain:0,possible:0},dam:{certain:0,possible:0},dcm:{certain:0,possible:0}};
      const sites=[];let outside=0;
      for (const match of [...matches.values()].sort((a,b)=>a.start-b.start || a.strand.localeCompare(b.strand))) {
        const rkey=match.start+':'+match.length;
        (match.certain?recognitionCertain:recognitionPossible).add(rkey);
        const cuts=cutsAt(enzyme,match.start,match.length,match.strand,n,circular);
        const context=contextAt(sequence,match.start,match.length,circular);
        for(const name of Object.keys(contexts))if(context[name])contexts[name][match.certain && context[name]===2?'certain':'possible']++;
        for(const cut of cuts) {
          if(!cut.valid){outside++;continue;}
          const key=String(cut.top)+':'+String(cut.bottom);
          (match.certain?certain:possible).add(key);
          if(match.certain && cut.top!==null && cut.bottom!==null)topBoundaries.add(cut.top);
        }
        if(sites.length<limit)sites.push(Object.assign({},match,{cuts,context}));
      }
      for(const key of certain)possible.delete(key);
      for(const key of recognitionCertain)recognitionPossible.delete(key);
      const kind=enzymeKind(enzyme);
      results.push({enzyme,kind,end:endType(enzyme),cutCount:kind==='unknown'?null:certain.size,possibleCuts:possible.size,
                    siteCount:recognitionCertain.size,possibleSites:recognitionPossible.size,sites,contexts,outside,
                    sitesTotal:matches.size,positionsTruncated:matches.size>sites.length,
                    fragments:kind==='double' && !possible.size?fragments(n,topBoundaries,circular,limit):null});
      if (typeof options.progress==='function' && (index%50===0 || index===enzymes.length-1))options.progress(index+1,enzymes.length);
    }
    return {sequence,circular,length:n,results};
  }
  function filter(results,options) {
    const o=options||{},terms=String(o.search||'').toLowerCase().split(/[\s,;]+/).filter(Boolean);
    const minimum=o.minimum===undefined?0:Number(o.minimum),maximum=o.maximum===undefined || o.maximum===null?Infinity:Number(o.maximum);
    if(!Number.isInteger(minimum)||minimum<0||!(maximum===Infinity || Number.isInteger(maximum))||maximum<minimum)throw new Error('Use a valid non-negative cut-count range.');
    return results.filter(row=>{
      if(terms.length && !terms.some(t=>row.enzyme.name.toLowerCase().includes(t)))return false;
      if(o.kind && o.kind!=='any' && row.kind!==o.kind)return false;
      if(o.end && o.end!=='any' && row.end!==o.end)return false;
      if(o.commercial && !(row.enzyme.suppliers||[]).length)return false;
      if(o.customOnly && !row.enzyme.custom)return false;
      const size=Math.min(...row.enzyme.motifs.map(s=>s.length));
      if(size<(o.minLength||0) || size>(o.maxLength||Infinity))return false;
      const count=o.basis==='recognition'?row.siteCount:row.cutCount,uncertain=o.basis==='recognition'?row.possibleSites:row.possibleCuts;
      if(count===null){if(minimum>0 || maximum!==Infinity)return false;}
      else if(o.includePossible ? count+uncertain<minimum || count+uncertain>maximum : count<minimum || count+uncertain>maximum)return false;
      for(const name of ['cpg','dam','dcm']) {
        const chosen=o[name]||'any',state=(row.enzyme.methylation||{})[name]||'unknown';
        if(chosen==='sensitive' && !state.includes('blocked') && !state.includes('impaired'))return false;
        if(!['any','sensitive'].includes(chosen) && chosen!==state)return false;
      }
      if(o.context && o.context!=='any') {
        const [mode,name]=o.context.split('_'),c=row.contexts[name];
        if(!c || (mode==='overlap' ? c.certain===0 : mode==='avoid' ? c.certain+c.possible>0 : true))return false;
      }
      return true;
    }).sort((a,b)=>a.enzyme.name.localeCompare(b.enzyme.name,undefined,{numeric:true}));
  }
  function cell(value) {
    let s=String(value===null || value===undefined?'':value).replace(/[\t\r\n]/g,' ');
    return /^[=+\-@]/.test(s)?"'"+s:s;
  }
  function exportTSV(rows,sites) {
    const output=[sites?['enzyme','site_start_1based','site_length','strand','certainty','top_cut_boundary_0based','bottom_cut_boundary_0based','complete_cut','end','positions_retained','positions_total']:['enzyme','recognition','certain_cuts','possible_cuts','certain_recognition_sites','possible_recognition_sites','ends','CpG','Dam','Dcm','positions_retained','positions_total']];
    for(const row of rows) {
      if(sites)for(const site of row.sites)for(const cut of site.cuts.length?site.cuts:[{top:null,bottom:null,valid:false,end:'unknown'}])
        output.push([row.enzyme.name,site.start+1,site.length,site.strand,site.certain?'certain':'possible',cut.top,cut.bottom,cut.valid,cut.end,row.sites.length,row.sitesTotal]);
      else output.push([row.enzyme.name,row.enzyme.motifs.join('|'),row.cutCount,row.possibleCuts,row.siteCount,row.possibleSites,row.end,...['cpg','dam','dcm'].map(k=>(row.enzyme.methylation||{})[k]||'unknown'),row.sites.length,row.sitesTotal]);
    }
    return output.map(r=>r.map(cell).join('\t')).join('\n')+'\n';
  }
  return {normalize,analyze,filter,parseCustom,validateEnzyme,relation,contextAt,enzymeKind,endType,exportTSV,LABELS,STATES};
}));

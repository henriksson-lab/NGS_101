/* Automatic index resolution from primer pairs; no protocol/candidate selection. */
(function(root,factory){
  if(typeof module==='object'&&module.exports) module.exports=factory(require('./sequence_annotation.js'),require('./illumina_samplesheet.js'));
  else root.SampleSheetBatch=factory(root.SequenceAnnotation,root.IlluminaSampleSheet);
})(typeof globalThis!=='undefined'?globalThis:this,function(matcher,policy){
  'use strict';
  const MAX_ROWS=200;
  function detector(catalogue) {
    const cache=new Map();
    return function(text) {
      if(!text.trim()) return {empty:true,values:{},protocols:[],candidates:[]};
      const sequence=matcher.normalizeInput(text);
      if(cache.has(sequence)) return cache.get(sequence);
      const result=matcher.annotate(sequence,catalogue), candidates=result.indexCandidates || [],values={};
      const errors=[],warnings=[];
      if(result.indexCandidatesTruncated) errors.push('Index candidate limit reached. Shorten the primer input.');
      for(const role of ['i7','i5']) {
        const records=candidates.filter(c=>c.role===role), ready=records.filter(c=>c.status==='ready'&&c.forwardSequence);
        const distinct=[...new Set(ready.map(c=>c.forwardSequence))];
        if(distinct.length>1) errors.push('Conflicting '+role+' values were detected.');
        else if(distinct.length===1) values[role]=distinct[0];
        // A matching unresolved placement can contradict the resolved interpretation.
        if(records.some(c=>c.status!=='ready'||!c.forwardSequence)) errors.push('A matching '+role+' template has unresolved bases or orientation.');
      }
      if(!Object.keys(values).length) errors.push('No orientation-resolved index template matched this primer.');
      if(candidates.some(c=>!c.role)) warnings.push('Other matching templates have unresolved index roles.');
      const output={empty:false,values,protocols:[...new Set(candidates.filter(c=>c.status==='ready').map(c=>c.protocol))],candidates,warnings,error:[...new Set(errors)].join(' ')};
      cache.set(sequence,output);return output;
    };
  }
  function resolvePair(row,detect) {
    const primers=[row.primer1 || '',row.primer2 || ''];
    const results=primers.map(text=>{try{return detect(text);}catch(error){return {empty:false,values:{},protocols:[],error:error.message};}});
    const errors=results.map((r,i)=>r.error?'Primer '+(i+1)+': '+r.error:'').filter(Boolean),values={};
    for(const role of ['i7','i5']) {
      const distinct=[...new Set(results.map(r=>r.values[role]).filter(Boolean))];
      if(distinct.length>1) errors.push('The primers disagree on '+role+'.');
      else if(distinct.length===1) values[role]=distinct[0];
    }
    if(!/^[A-Za-z0-9][A-Za-z0-9_-]*$/.test(row.sampleId)||row.sampleId.toLowerCase()==='undetermined') errors.push('Use a Sample_ID with letters, numbers, underscores or hyphens (Undetermined is reserved).');
    return {...row,values,results,protocols:[...new Set(results.flatMap(r=>r.protocols))],empty:results.every(r=>r.empty),error:errors.join(' ')};
  }
  function analyze(rows,catalogue,instrument='novaseqx',workflow='manual-bcl2fastq') {
    if(rows.length>MAX_ROWS) throw new Error('At most '+MAX_ROWS+' samples can be processed at once.');
    const profile=policy.describeProfile(instrument,workflow);
    if(!profile.supported) throw new Error(profile.note);
    const detect=detector(catalogue),results=rows.map(row=>resolvePair(row,detect)),ids=new Map();
    for(const row of results) if(!row.empty) ids.set(row.sampleId,(ids.get(row.sampleId)||0)+1);
    for(const row of results) if(!row.empty&&ids.get(row.sampleId)>1) row.error+=' Duplicate Sample_ID.';
    const valid=results.filter(r=>!r.empty&&!r.error&&Object.keys(r.values).length);
    const headers=['Sample_ID','index','index2'];
    const values=valid.map(row=>[row.sampleId,...['i7','i5'].map(role=>row.values[role]?policy.formatIndex(row.values[role],role,instrument,workflow).value:'')]);
    const tsv=[headers,...values].map(row=>row.join('\t')).join('\n');
    const header=profile.headerRows.length?'[Header]\n'+profile.headerRows.map(row=>row.join(',')).join('\n')+'\n\n':'';
    const csv=header+'['+profile.section+']\n'+[headers,...values].map(row=>row.join(',')).join('\n');
    return {results,valid,headers,values,tsv,csv,profile,inputOnly:workflow==='run-planning'};
  }
  function reuse(rows,sourceIndex,slot,role,catalogue) {
    const detect=detector(catalogue),text=rows[sourceIndex]['primer'+slot]||'',source=detect(text);
    if(source.empty||source.error||!source.values[role]||Object.keys(source.values).length!==1) throw new Error('Use a primer with one unambiguous '+role+' index for this action.');
    return rows.map(row=>{
      const current=[row.primer1||'',row.primer2||''];
      const detected=current.map(value=>{try{return detect(value);}catch(_){return {values:{}};}});
      let target=detected.findIndex(r=>r.values[role]);
      if(target<0) target=current.findIndex(value=>!value.trim());
      if(target<0) target=detected.findIndex(r=>!r.values[role==='i7'?'i5':'i7']);
      if(target<0) throw new Error('A row has two primers for the opposite index. Clear a cell before applying this primer.');
      return {...row,['primer'+(target+1)]:text};
    });
  }
  function parseTSV(text) {
    let lines=text.trim().split(/\r?\n/).filter(line=>line.trim());
    if(/^sample_id\t/i.test(lines[0]||'')) lines.shift();
    if(!lines.length) throw new Error('Paste at least one sample row.');
    if(lines.length>MAX_ROWS) throw new Error('At most '+MAX_ROWS+' samples can be processed at once.');
    return lines.map(line=>{const parts=line.split('\t');if(parts.length<2||parts.length>3)throw new Error('Use tab-separated Sample_ID, primer 1, and optional primer 2.');return {sampleId:parts[0].trim(),primer1:parts[1].trim(),primer2:(parts[2]||'').trim()};});
  }
  return {analyze,detector,resolvePair,reuse,parseTSV,MAX_ROWS};
});

/* Browser-local sequence annotation. All query and catalogue text is inserted as text. */
(function () {
  'use strict';
  const $ = id => document.getElementById(id);
  const node = (tag, text, cls) => { const el = document.createElement(tag); if (text !== undefined) el.textContent = text; if (cls) el.className = cls; return el; };
  const clear = el => el.replaceChildren();
  const MAX_ROWS = 500, PAGE_ROWS = 50, PREVIEW = 1000;
  let catalogue, loading, result, previewStart = 0, rowsShown = PAGE_ROWS, candidates = [], rowOutput = null, csvOutput = '', generation = 0, sheetGeneration = 0;
  let protocolMap = new Map();
  function safeLink(url, title) {
    // Only public same-site protocol pages belong in catalogue navigation.
    const a = node('a', title);
    if (typeof url === 'string' && !url.startsWith('/') && !url.includes('..') && !url.includes(':') && /^[a-zA-Z0-9_+./%~-]+\.html(?:#[\w-]+)?$/.test(url)) a.href = url;
    return a;
  }
  function protocol(id) { return protocolMap.get(id) || {id, title:id}; }
  function title(hit) { return (hit.aliases || [])[0] || hit.component; }
  async function getCatalogue() {
    if (catalogue) return catalogue;
    if (!loading) loading = fetch('assets/sequence_components.json').then(response => { if (!response.ok) throw new Error('The component catalogue could not be loaded. Please try again.'); return response.json(); }).then(data => {
      catalogue = data;
      protocolMap = new Map((Array.isArray(data.protocols) ? data.protocols : Object.values(data.protocols || {})).map(p => [p.id,p]));
      return data;
    }).catch(error => { loading = null; throw error; });
    return loading;
  }
  function status(text, error) { $('annotation-status').textContent = text; $('annotation-status').dataset.error = String(Boolean(error)); }
  function protocolLinks(ids) {
    const box = node('div'), unique = [...new Set(ids || [])], initial = node('ul',undefined,'protocol-links');
    for (const id of unique.slice(0,8)) { const p = protocol(id), li = node('li'); li.append(safeLink(p.url,p.title || id)); initial.append(li); }
    box.append(node('span', `${Math.min(8,unique.length)} of ${unique.length} protocols`),initial);
    if (unique.length > 8) { const details = node('details'), summary = node('summary',`Show up to ${Math.min(20,unique.length)} protocols`), more = node('ul',undefined,'protocol-links'); details.append(summary,more); for (const id of unique.slice(8,20)) { const p = protocol(id), li = node('li'); li.append(safeLink(p.url,p.title || id)); more.append(li); } if (unique.length > 20) details.append(node('p','Display capped at 20 protocol links.','muted')); box.append(details); }
    return box;
  }
  function renderPreview() {
    const end = Math.min(result.sequence.length,previewStart+PREVIEW), pre = $('sequence-preview'); clear(pre);
    let cursor = previewStart;
    for (const hit of result.overview || []) {
      if (hit.end <= previewStart || hit.start >= end) continue;
      const start = Math.max(hit.start,previewStart), stop = Math.min(hit.end,end);
      if (start > cursor) pre.append(document.createTextNode(result.sequence.slice(cursor,start)));
      const mark = node('mark',result.sequence.slice(start,stop));
      mark.tabIndex = 0; mark.setAttribute('role','button'); mark.title = `${title(hit)} · ${hit.start+1}–${hit.end} (${hit.strand})`;
      const select = () => {
        const overlapping = result.hits.filter(h => h.start < stop && h.end > start);
        $('region-detail').textContent = overlapping.slice(0,20).map(h => `${title(h)} ${h.start+1}–${h.end} (${h.strand})`).join('; ') + (overlapping.length > 20 ? '; further overlapping matches are listed below.' : '');
      };
      mark.addEventListener('click',select); mark.addEventListener('keydown',event => { if (event.key === 'Enter' || event.key === ' ') { event.preventDefault(); select(); } });
      pre.append(mark); cursor = stop;
    }
    if (cursor < end) pre.append(document.createTextNode(result.sequence.slice(cursor,end)));
    $('preview-range').textContent = `${previewStart+1}–${end} of ${result.sequence.length.toLocaleString()} nt`;
    $('preview-prev').disabled = previewStart === 0; $('preview-next').disabled = end >= result.sequence.length;
  }
  function renderTable() {
    const tbody = $('match-rows'); clear(tbody);
    for (const hit of result.hits.slice(0,Math.min(rowsShown,MAX_ROWS))) {
      const tr = node('tr'), position = node('td',`${hit.start+1}–${hit.end} (${hit.strand})`), description = node('td'), evidence = node('td'), associations = node('td');
      position.append(node('div',`${hit.end-hit.start} nt`,'muted'));
      description.append(node('span',title(hit),'component-label'));
      if ((hit.aliases || []).length > 1) description.append(node('span',hit.aliases.slice(1,20).join('; ') + (hit.aliases.length > 20 ? `; ${hit.aliases.length-20} additional aliases` : ''),'component-aliases'));
      description.append(node('span',result.sequence.slice(hit.start,hit.end),'match-sequence'));
      const occurrences = hit.occurrences || [], kinds = [...new Set(occurrences.map(o => o.role || o.kind || o.usage || 'sequence component'))];
      evidence.append(node('div',kinds.slice(0,5).join('; ') || 'Sequence component'));
      if (occurrences.some(o => o.inferred)) evidence.append(node('div','Includes inferred declaration','muted'));
      if (occurrences.some(o => o.reference || o.reference_only)) evidence.append(node('div','Compatible reference; loaded reagent may be unverified','muted'));
      if (occurrences.length) {
        const d = node('details'), s = node('summary','Declaration details'); d.append(s);
        for (const o of occurrences.slice(0,8)) d.append(node('p',[o.name,o.locator,o.orientation ? 'Local component orientation: '+o.orientation : '',o.source,o.reason].filter(Boolean).map(v => typeof v === 'object' ? JSON.stringify(v) : v).join(' · '),'muted'));
        evidence.append(d);
      }
      associations.append(protocolLinks(hit.protocols)); tr.append(position,description,evidence,associations); tbody.append(tr);
    }
    const shown = Math.min(rowsShown,MAX_ROWS,result.hits.length);
    $('match-table-limit').textContent = `${shown.toLocaleString()} of ${result.hits.length.toLocaleString()} stored matches shown. Coordinates are 1-based inclusive. Strand is the matched canonical component relative to your submitted strand; local declarations can use the opposite orientation.` + (result.hits.length > MAX_ROWS ? ` Table display is capped at ${MAX_ROWS} rows.` : '');
    $('more-matches').hidden = shown >= Math.min(MAX_ROWS,result.hits.length);
  }
  function renderProtocols() {
    clear($('protocol-matches'));
    const matchedNames = new Map(result.hits.map(hit => [hit.component,title(hit)]));
    const rows = (result.protocols || []).slice(0,20), informative = rows.some(p => p.informative);
    $('protocol-explanation').textContent = informative ? 'Up to 20 protocols ranked by informative shared components and covered bases. Nested matches are not independent evidence. This is a navigation aid, not a confidence score.' : 'Only common components, or no informative components, matched. Protocol links are examples of shared parts and do not establish an origin.';
    for (const p of rows) { const li = node('li'); li.append(safeLink(p.url,p.title || p.id)); li.append(node('small', (p.informative ? 'Shared parts: ' : 'Common components: ') + (p.support || []).map(id => { return matchedNames.get(id) || id; }).join('; '))); $('protocol-matches').append(li); }
  }
  function candidateName(c) { const p = protocol(c.protocol); return `${c.inferred ? 'INFERRED · ' : ''}${c.label || c.feature_id || c.role} · ${p.title || c.protocol} · ${c.forwardSequence || 'unresolved'}`; }
  function chosen(role) { const value = $('index-'+role).value; return value === '' ? null : candidates[Number(value)]; }
  function fillCandidate(role, preserve) {
    const select = $('index-'+role), old = preserve ? select.value : '', i7 = chosen('i7');
    clear(select); const placeholder = node('option',role === 'i7' ? 'Choose i7' : 'No i5 / single index'); placeholder.value = ''; select.append(placeholder);
    candidates.forEach((c,index) => { if (c.role !== role || c.status !== 'ready' || !c.forwardSequence || (role === 'i5' && i7 && c.protocol !== i7.protocol)) return; const option = node('option',candidateName(c)); option.value = String(index); select.append(option); });
    if ([...select.options].some(o => o.value === old)) select.value = old;
  }
  function renderCandidates() {
    candidates = result.indexCandidates || []; fillCandidate('i7',false); fillCandidate('i5',false);
    const ready = candidates.filter(c => c.status === 'ready' && c.forwardSequence);
    $('index-summary').textContent = candidates.length ? `${ready.length} orientation-resolved candidates; ${candidates.length-ready.length} unresolved candidates. Select a candidate explicitly before copying.` : 'No declared index template matched. No index values can be exported from this query.';
    if (result.indexCandidatesTruncated) $('index-summary').append(document.createTextNode(' Candidate limit reached; index results are incomplete.'));
    clear($('index-unresolved'));
    for (const c of candidates.filter(c => c.status !== 'ready' || !c.forwardSequence).slice(0,20)) $('index-unresolved').append(node('li',`${candidateName(c)}: ${c.reason || 'Forward index orientation is not resolved.'}`));
    renderSampleSheet();
  }
  async function copyText(text) {
    const revision = sheetGeneration;
    try { if (!navigator.clipboard || !navigator.clipboard.writeText) throw new Error('clipboard unavailable'); await navigator.clipboard.writeText(text); if (revision !== sheetGeneration) return; $('copy-status').textContent = 'Copied.'; }
    catch (_) { if (revision !== sheetGeneration) return; const output = $('index-output'); $('index-snippet').hidden = false; output.value = text; output.focus(); output.select(); $('copy-status').textContent = 'Select and copy the highlighted text using your browser’s Copy command.'; }
  }
  function renderSampleSheet() {
    sheetGeneration++; rowOutput = null; csvOutput = ''; $('index-output').value = ''; clear($('index-fields')); clear($('index-candidate-details')); $('index-snippet').hidden = true; $('copy-status').textContent = '';
    const api = window.IlluminaSampleSheet, instrument = $('index-instrument').value, workflow = $('index-workflow').value, selected = ['i7','i5'].map(chosen).filter(Boolean), note = $('index-profile-note');
    note.dataset.error = 'false';
    for (const c of selected) {
      const detail = node('div',undefined,'index-detail');
      detail.append(node('strong',`${c.inferred ? 'INFERRED · ' : ''}${c.role}: ${c.label || c.feature_id || 'declared index'}`),node('div',`Illumina-forward index: ${c.forwardSequence} · query ${c.start+1}–${c.end} (${c.strand})`));
      detail.append(node('div',[c.source,c.basis].filter(Boolean).map(v => typeof v === 'object' ? JSON.stringify(v) : v).join(' · ')));
      $('index-candidate-details').append(detail);
    }
    if (!api) { note.textContent = 'Sample-sheet orientation support could not be loaded.'; note.dataset.error = 'true'; return; }
    if (!instrument || !workflow) { note.textContent = 'Choose both instrument and workflow. The i5 orientation depends on both.'; return; }
    try {
      const profile = api.describeProfile(instrument,workflow);
      if (!profile.supported) { note.textContent = profile.note || 'This instrument/workflow combination is not supported.'; note.dataset.error = 'true'; return; }
      note.textContent = `${profile.label || ''} · [${profile.section}] · i5: ${profile.i5Orientation === 'reverse_complement' ? 'reverse complement of Illumina-forward index' : 'Illumina-forward index'}. ${profile.note || ''}`;
      for (const url of profile.sources || []) { if (/^https:\/\/[a-zA-Z0-9.-]*illumina\.com\//.test(url)) { const a = node('a',' Illumina orientation guidance'); a.href = url; note.append(a); } }
      for (const c of selected) {
        const formatted = api.formatIndex(c.forwardSequence,c.role,instrument,workflow), box = node('div',undefined,'index-field'), copy = node('button','Copy '+formatted.column);
        copy.type = 'button'; copy.addEventListener('click',() => copyText(formatted.value));
        box.append(node('strong',formatted.column),node('code',formatted.value),copy,node('small',`${c.role} · ${formatted.orientation}. ${formatted.note || ''}`)); $('index-fields').append(box);
      }
      if (workflow === 'run-planning') { note.append(document.createTextNode(' These are Run Planning inputs only. Use the sample sheet exported by that setup tool; no manually prepared CSV snippet is offered here.')); return; }
      if (!chosen('i7')) return;
      rowOutput = api.makeRow({sampleId:$('index-sample-id').value,i7:chosen('i7').forwardSequence,i5:chosen('i5') ? chosen('i5').forwardSequence : undefined,instrument,workflow});
      csvOutput = (rowOutput.headerRows || []).map(r => Array.isArray(r) ? r.join(',') : r).join('\n');
      csvOutput = (csvOutput ? '[Header]\n' + csvOutput + '\n\n' : '') + `[${rowOutput.section}]\n` + rowOutput.csv;
      $('index-output').value = csvOutput; $('index-snippet').hidden = false;
    } catch (error) { rowOutput = null; csvOutput = ''; $('index-output').value = ''; $('index-snippet').hidden = true; clear($('index-fields')); note.textContent = error.message; note.dataset.error = 'true'; }
  }
  $('annotate-form').addEventListener('submit',async event => {
    event.preventDefault(); const token = ++generation; sheetGeneration++;
    $('annotate-button').disabled = true; $('annotation-results').hidden = true;
    try {
      if (!window.SequenceAnnotation) throw new Error('The annotation script could not be loaded. Please reload the page.');
      const sequence = window.SequenceAnnotation.normalizeInput($('sequence-input').value,{convertUToT:$('convert-u').checked});
      status('Loading the component catalogue…'); const data = await getCatalogue(); if (token !== generation) return;
      status('Finding exact matches…'); await new Promise(resolve => setTimeout(resolve,0)); if (token !== generation) return;
      result = window.SequenceAnnotation.annotate(sequence,data,{includeInferred:$('include-inferred').checked,maxHits:10000,maxProtocols:20});
      previewStart = 0; rowsShown = PAGE_ROWS;
      $('result-summary').textContent = `${sequence.length.toLocaleString()} nt · ${result.hits.length.toLocaleString()} exact matches` + ($('convert-u').checked ? ' · U converted to T' : ' · U retained') + (result.truncated ? ' · Hit limit reached: results are incomplete.' : '');
      renderPreview(); renderTable(); renderProtocols(); renderCandidates(); $('annotation-results').hidden = false;
      status(result.hits.length ? 'Annotation complete.' : 'No exact component matches found.');
    } catch (error) { if (token === generation) status(error.message || 'Annotation failed. Please try again.',true); }
    finally { if (token === generation) $('annotate-button').disabled = false; }
  });
  $('example-button').addEventListener('click',() => { $('sequence-input').value = $('example-button').dataset.sequence; $('sequence-input').focus(); });
  $('clear-button').addEventListener('click',() => { generation++; sheetGeneration++; rowOutput = null; $('sequence-input').value = ''; $('annotation-results').hidden = true; $('annotate-button').disabled = false; result = null; status('Paste a sequence to begin.'); $('sequence-input').focus(); });
  $('preview-prev').addEventListener('click',() => { previewStart = Math.max(0,previewStart-PREVIEW); renderPreview(); });
  $('preview-next').addEventListener('click',() => { previewStart += PREVIEW; renderPreview(); });
  $('more-matches').addEventListener('click',() => { rowsShown += PAGE_ROWS; renderTable(); });
  $('index-i7').addEventListener('change',() => { fillCandidate('i5',true); renderSampleSheet(); });
  for (const id of ['index-i5','index-instrument','index-workflow','index-sample-id']) $(id).addEventListener(id === 'index-sample-id' ? 'input' : 'change',renderSampleSheet);
  $('copy-index-tsv').addEventListener('click',() => { if (rowOutput) copyText(rowOutput.tsv); });
  $('copy-index-csv').addEventListener('click',() => { if (rowOutput) copyText(csvOutput); });
  const sampleApi = window.IlluminaSampleSheet;
  if (sampleApi) for (const [id,options] of [['index-instrument',sampleApi.instruments],['index-workflow',sampleApi.workflows]]) for (const option of options || []) { const el = node('option',option.label); el.value = option.id; $(id).append(el); }
}());

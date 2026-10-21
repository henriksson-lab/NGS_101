/* Browser-local sequence annotation. All query and catalogue text is inserted as text. */
(function () {
  'use strict';
  const $ = id => document.getElementById(id);
  const node = (tag, text, cls) => { const el = document.createElement(tag); if (text !== undefined) el.textContent = text; if (cls) el.className = cls; return el; };
  const clear = el => el.replaceChildren();
  const MAX_ROWS = 500, PAGE_ROWS = 50, PREVIEW = 1000;
  let catalogue, loading, result, previewStart = 0, rowsShown = PAGE_ROWS, generation = 0;
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
  function renderAlignment(container, hit) {
    clear(container);
    const flank = 35, start = Math.max(0,hit.start-flank), end = Math.min(result.sequence.length,hit.end+flank);
    const sequence = hit.strand === '-' ? window.SequenceAnnotation.revcomp(hit.sequence) : hit.sequence;
    // Align stored exact strings on query coordinates. This is sequence comparison,
    // not a model of primer annealing or a pairwise/gapped alignment.
    const heading = node('p',`${title(hit)} · query ${hit.start+1}–${hit.end} · ${hit.end-hit.start} exact bases`,'alignment-caption');
    const pre = node('pre',undefined,'region-alignment');
    const query = node('div'), component = node('div');
    query.append(node('span','Query    ','alignment-row-label'),document.createTextNode(result.sequence.slice(start,hit.start)),node('strong',result.sequence.slice(hit.start,hit.end),'alignment-match'),document.createTextNode(result.sequence.slice(hit.end,end)));
    component.append(node('span','Selected ','alignment-row-label'),document.createTextNode(' '.repeat(hit.start-start)),node('strong',sequence,'alignment-match'));
    pre.append(query,component);
    const orientation = hit.strand === '-' ? 'The selected canonical component is reverse-complemented to align with the submitted strand.' : hit.strand === 'both' ? 'This component is palindromic; either orientation has the same sequence.' : 'The canonical component is shown in the submitted strand’s orientation.';
    const details = node('p',`${orientation} Query context: ${start+1}–${end}. Rows run left to right in query order. Equality does not establish primer binding.`,'muted');
    container.append(heading,pre,details);
    const occurrences = hit.occurrences || [];
    if (occurrences.some(o => o.inferred)) container.append(node('p','Includes inferred declaration.','muted'));
    if (occurrences.some(o => o.reference_only || o.reference)) container.append(node('p','Compatible sequencing reference; loaded reagent may be unverified.','muted'));
    // A deduplicated component can be a fixed part rather than the full source oligo.
    const declarations = node('details'); declarations.append(node('summary','Primer/component declarations'));
    for (const occurrence of occurrences.slice(0,8)) {
      const row = node('p',undefined,'muted'), declared = occurrence.orientation === '-' ? window.SequenceAnnotation.revcomp(hit.sequence) : hit.sequence;
      row.append(document.createTextNode([occurrence.name,occurrence.role,protocol(occurrence.protocol).title,occurrence.locator].filter(Boolean).join(' · ')),node('code',declared,'alignment-declared'));
      declarations.append(row);
    }
    if (occurrences.length > 8) declarations.append(node('p',`${occurrences.length-8} additional declarations are available in the component table.`,'muted'));
    if (occurrences.length) container.append(declarations);
  }
  function renderPreview() {
    const end = Math.min(result.sequence.length,previewStart+PREVIEW), pre = $('sequence-preview'), regions = $('region-alignments'); clear(pre); clear(regions);
    $('region-detail').textContent = '';
    let cursor = previewStart, regionNumber = 0;
    for (const hit of result.overview || []) {
      if (hit.end <= previewStart || hit.start >= end) continue;
      const start = Math.max(hit.start,previewStart), stop = Math.min(hit.end,end);
      if (start > cursor) pre.append(document.createTextNode(result.sequence.slice(cursor,start)));
      const mark = node('mark',result.sequence.slice(start,stop));
      mark.tabIndex = 0; mark.setAttribute('role','button'); mark.title = `${title(hit)} · ${hit.start+1}–${hit.end} (${hit.strand})`;
      const region = node('section',undefined,'region-selector'), label = node('label',`Region ${hit.start+1}–${hit.end}: matching primer/component`), select = node('select'), alignment = node('div');
      const selectId = `region-select-${regionNumber++}`; select.id = selectId; label.htmlFor = selectId;
      // Use full overview coordinates, including when a region crosses a preview page.
      // Retain partial overlaps and contained components for explicit inspection.
      const overlapping = result.hits.filter(h => h.start < hit.end && h.end > hit.start);
      overlapping.sort((a,b) => (b.end-b.start)-(a.end-a.start) || a.start-b.start || a.component.localeCompare(b.component));
      const candidates = overlapping.slice(0,100), selected = Math.max(0,candidates.indexOf(hit));
      for (const [index,candidate] of candidates.entries()) {
        const option = node('option',`${title(candidate)} · ${candidate.start+1}–${candidate.end} (${candidate.strand})`); option.value = String(index); select.append(option);
      }
      select.value = String(selected);
      const update = () => { const candidate = candidates[Number(select.value)]; if (candidate) renderAlignment(alignment,candidate); };
      select.addEventListener('change',update); region.append(label,select,alignment);
      if (overlapping.length > candidates.length) region.append(node('p','Selector capped at 100 matches; further matches remain in the component table.','muted'));
      regions.append(region); update();
      const focus = () => { $('region-detail').textContent = `Inspecting region ${hit.start+1}–${hit.end}. Choose a primer/component in its selector.`; select.focus(); };
      mark.addEventListener('click',focus); mark.addEventListener('keydown',event => { if (event.key === 'Enter' || event.key === ' ') { event.preventDefault(); focus(); } });
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
  $('annotate-form').addEventListener('submit',async event => {
    event.preventDefault(); const token = ++generation;
    $('annotate-button').disabled = true; $('annotation-results').hidden = true;
    try {
      if (!window.SequenceAnnotation) throw new Error('The annotation script could not be loaded. Please reload the page.');
      const sequence = window.SequenceAnnotation.normalizeInput($('sequence-input').value,{convertUToT:$('convert-u').checked});
      status('Loading the component catalogue…'); const data = await getCatalogue(); if (token !== generation) return;
      status('Finding exact matches…'); await new Promise(resolve => setTimeout(resolve,0)); if (token !== generation) return;
      result = window.SequenceAnnotation.annotate(sequence,data,{includeInferred:$('include-inferred').checked,maxHits:10000,maxProtocols:20});
      previewStart = 0; rowsShown = PAGE_ROWS;
      $('result-summary').textContent = `${sequence.length.toLocaleString()} nt · ${result.hits.length.toLocaleString()} exact matches` + ($('convert-u').checked ? ' · U converted to T' : ' · U retained') + (result.truncated ? ' · Hit limit reached: results are incomplete.' : '');
      renderPreview(); renderTable(); renderProtocols(); $('annotation-results').hidden = false;
      status(result.hits.length ? 'Annotation complete.' : 'No exact component matches found.');
    } catch (error) { if (token === generation) status(error.message || 'Annotation failed. Please try again.',true); }
    finally { if (token === generation) $('annotate-button').disabled = false; }
  });
  $('example-button').addEventListener('click',() => { $('sequence-input').value = $('example-button').dataset.sequence; $('sequence-input').focus(); });
  $('clear-button').addEventListener('click',() => { generation++; $('sequence-input').value = ''; $('annotation-results').hidden = true; $('annotate-button').disabled = false; result = null; status('Paste a sequence to begin.'); $('sequence-input').focus(); });
  $('preview-prev').addEventListener('click',() => { previewStart = Math.max(0,previewStart-PREVIEW); renderPreview(); });
  $('preview-next').addEventListener('click',() => { previewStart += PREVIEW; renderPreview(); });
  $('more-matches').addEventListener('click',() => { rowsShown += PAGE_ROWS; renderTable(); });
}());

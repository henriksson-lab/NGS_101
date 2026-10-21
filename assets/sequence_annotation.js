/* Exact sequence comparison. Coordinates are zero-based, half-open, on the query. */
(function (root, factory) {
  'use strict';
  const api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.SequenceAnnotation = api;
}(typeof globalThis !== 'undefined' ? globalThis : this, function () {
  'use strict';
  const LIMIT = 100000;
  function revcomp(sequence) {
    const pairs = {A:'T',T:'A',C:'G',G:'C',N:'N',R:'Y',Y:'R',S:'S',W:'W',K:'M',M:'K',B:'V',V:'B',D:'H',H:'D',U:'U'};
    return Array.from(sequence).reverse().map(base => pairs[base] || base).join('');
  }
  function normalizeInput(text, options) {
    options = options || {};
    const lines = String(text).trim().split(/\r?\n/);
    if (lines[0] && lines[0].trim().startsWith('>')) lines.shift();
    if (lines.some(line => line.includes('>'))) throw new Error('Paste one sequence or a single FASTA record.');
    let sequence = lines.join('').replace(/\s/g, '').toUpperCase();
    if (!sequence) throw new Error('Enter a sequence.');
    if (sequence.length > LIMIT) throw new Error('The sequence limit is 100,000 bases.');
    if (!/^[ACGTURYSWKMBDHVN]+$/.test(sequence)) throw new Error('Use nucleotide letters only (IUPAC ambiguity letters are allowed).');
    if (options.convertUToT || options.convertU) sequence = sequence.replace(/U/g, 'T');
    return sequence;
  }
  const cache = new WeakMap();
  function compile(catalogue, includeInferred) {
    let variants = cache.get(catalogue);
    if (!variants) { variants = {}; cache.set(catalogue, variants); }
    const key = includeInferred ? 'all' : 'declared';
    if (variants[key]) return variants[key];
    const nodes = [{next:Object.create(null),fail:0,out:[]}];
    const components = new Map();
    for (const component of catalogue.components || []) {
      if (!/^[ACGT]+$/.test(component.sequence) || component.sequence.length < 12) continue;
      if (component.inferred && !includeInferred) continue;
      const occurrences = (component.occurrences || []).filter(item => includeInferred || !item.inferred);
      if ((component.occurrences || []).length && !occurrences.length) continue;
      const protocols = Array.from(new Set(occurrences.length ? occurrences.map(item => item.protocol) : component.protocols || []));
      const aliases = occurrences.length ? Array.from(new Set(occurrences.map(item => item.name).filter(Boolean))).sort() : component.aliases || [];
      const clean = Object.assign({}, component, {occurrences, protocols, aliases});
      components.set(component.id, clean);
      const reverse = revcomp(component.sequence);
      const patterns = reverse === component.sequence ? [[component.sequence,'both']] : [[component.sequence,'+'],[reverse,'-']];
      for (const [pattern,strand] of patterns) {
        let state = 0;
        for (const base of pattern) {
          if (nodes[state].next[base] === undefined) { nodes[state].next[base] = nodes.length; nodes.push({next:Object.create(null),fail:0,out:[]}); }
          state = nodes[state].next[base];
        }
        nodes[state].out.push({id:component.id,length:pattern.length,strand});
      }
    }
    const queue = Object.values(nodes[0].next);
    for (let i = 0; i < queue.length; i++) {
      const state = queue[i];
      for (const [base,next] of Object.entries(nodes[state].next)) {
        queue.push(next);
        let fail = nodes[state].fail;
        while (fail && nodes[fail].next[base] === undefined) fail = nodes[fail].fail;
        nodes[next].fail = nodes[fail].next[base] === undefined ? 0 : nodes[fail].next[base];
        nodes[next].out = nodes[next].out.concat(nodes[nodes[next].fail].out);
      }
    }
    return (variants[key] = {nodes, components});
  }
  function nonoverlapping(hits, priority) {
    const selected = [];
    // A long containing component supplies one interval of evidence, regardless of nested parts.
    const ordered = hits.slice().sort((a,b) => (priority ? priority(b)-priority(a) : 0) || (b.end-b.start)-(a.end-a.start) || a.start-b.start || a.component.localeCompare(b.component));
    const occupied = new Uint8Array(ordered.reduce((end,hit) => Math.max(end,hit.end),0));
    for (const hit of ordered) {
      let overlaps = false;
      for (let i=hit.start; i<hit.end; i++) if (occupied[i]) { overlaps=true; break; }
      if (!overlaps) { selected.push(hit); occupied.fill(1,hit.start,hit.end); }
    }
    return selected.sort((a,b) => a.start-b.start || b.end-a.end);
  }
  function rank(hits, catalogue, maxProtocols) {
    const total = Math.max(1,(catalogue.protocols || []).length);
    const byProtocol = new Map();
    for (const hit of hits) for (const id of hit.protocols) {
      if (!byProtocol.has(id)) byProtocol.set(id, []);
      byProtocol.get(id).push(hit);
    }
    const metadata = new Map((catalogue.protocols || []).map(protocol => [protocol.id,protocol]));
    const ranked = [];
    for (const [id, matches] of byProtocol) {
      const selected = nonoverlapping(matches, hit => (hit.end-hit.start)*Math.log2((total+1)/(hit.protocols.length+1)));
      const counted = new Set();
      let score = 0, informative = false;
      for (const hit of selected) {
        if (counted.has(hit.component)) continue;
        counted.add(hit.component);
        const prevalence = hit.protocols.length / total;
        const weight = Math.log2((total + 1)/(hit.protocols.length + 1));
        score += (hit.end-hit.start)*weight;
        if (prevalence <= 0.2) informative = true;
      }
      ranked.push(Object.assign({id,title:id,url:null}, metadata.get(id), {score,support:Array.from(counted),informative}));
    }
    ranked.sort((a,b) => b.score-a.score || a.title.localeCompare(b.title) || a.id.localeCompare(b.id));
    return ranked.slice(0,maxProtocols);
  }
  function indexCandidates(sequence, catalogue, includeInferred) {
    const results = [], seen = new Set();
    const candidateLimit = 1000;
    const reverseSequence = revcomp(sequence);
    function finish(truncated) {
      const grouped = new Map();
      for (const candidate of results) {
        const key = [candidate.template,candidate.feature_id,candidate.start,candidate.end].join(':');
        const existing = grouped.get(key);
        if (existing && existing.strand !== candidate.strand) {
          existing.strand = 'both';
          if (existing.sequence !== candidate.sequence || existing.forwardSequence !== candidate.forwardSequence) {
            existing.status = 'unresolved'; existing.forwardSequence = null;
            existing.reason = 'Both orientations of this template match the same index span; its orientation is ambiguous.';
          }
        } else if (!existing) grouped.set(key,candidate);
      }
      const candidates = Array.from(grouped.values()).sort((a,b) => (a.start === null ? Infinity : a.start)-(b.start === null ? Infinity : b.start) || String(a.template).localeCompare(String(b.template)));
      return {candidates,truncated};
    }
    for (const template of catalogue.templates || []) {
      if (template.inferred && !includeInferred) continue;
      const pattern = template.sequence || '';
      if (!/^[ACGTN]+$/.test(pattern)) continue;
      const runs = Array.from(pattern.matchAll(/[ACGT]+/g)).sort((a,b) => b[0].length-a[0].length || a.index-b.index);
      for (const span of template.index_spans || []) {
        const base = {template:template.id,protocol:template.protocol,label:template.label,feature_id:span.feature_id,role:span.role || null,source:template.source,locator:template.locator,inferred:!!template.inferred,basis:span.basis || null};
        if (!(Number.isInteger(span.start) && Number.isInteger(span.end) && span.start >= 0 && span.end > span.start && span.end <= pattern.length)) continue;
        if (span.inferred && !includeInferred) continue;
        const left = (pattern.slice(0,span.start).match(/[ACGT]+$/) || [''])[0].length;
        const right = (pattern.slice(span.end).match(/^[ACGT]+/) || [''])[0].length;
        const anchor = runs[0];
        if (!anchor || anchor[0].length < 12) continue;
        const insufficientAnchors = left < 12 || right < 12;
        for (const strand of ['+','-']) {
          const query = strand === '+' ? sequence : reverseSequence;
          let from = 0, position;
          while ((position = query.indexOf(anchor[0],from)) !== -1) {
            from = position + 1;
            const offset = position-anchor.index;
            if (offset < 0 || offset+pattern.length > query.length) continue;
            if (!runs.every(run => query.startsWith(run[0],offset+run.index))) continue;
            const start = strand === '+' ? offset+span.start : sequence.length-offset-span.end;
            const end = strand === '+' ? offset+span.end : sequence.length-offset-span.start;
            if (insufficientAnchors) {
              const unresolvedKey = [template.id,span.feature_id,'unanchored'].join(':');
              if (!seen.has(unresolvedKey)) {
                seen.add(unresolvedKey);
                if (results.length >= candidateLimit) return finish(true);
                results.push(Object.assign({},base,{start:null,end:null,strand:null,sequence:null,forwardSequence:null,status:'unresolved',reason:'The compatible declared template lacks 12 concrete anchoring bases on each side of this index.'}));
              }
              continue;
            }
            const unique = [template.id,span.feature_id,start,end,strand].join(':');
            if (seen.has(unique)) continue;
            seen.add(unique);
            const observed = query.slice(offset+span.start,offset+span.end);
            let reason = '', forwardSequence = null;
            if (!/^[ACGT]+$/.test(observed)) reason = 'The index contains ambiguous or RNA bases.';
            else if (!span.role || !['i5','i7'].includes(span.role)) reason = 'The sequencer index role is not declared.';
            else if (!['forward','reverse_complement'].includes(span.to_forward) || !span.basis) reason = 'The orientation relative to the Illumina forward index is not established.';
            else forwardSequence = span.to_forward === 'forward' ? observed : revcomp(observed);
            // Check the opposite placement immediately, including when the output cap
            // prevents scanning that orientation later.
            const oppositeQuery = strand === '+' ? reverseSequence : sequence;
            const oppositeOffset = strand === '+' ? sequence.length-start-span.end : start-span.start;
            const oppositeMatches = oppositeOffset >= 0 && oppositeOffset+pattern.length <= sequence.length && runs.every(run => oppositeQuery.startsWith(run[0],oppositeOffset+run.index));
            const ambiguousOrientation = oppositeMatches && oppositeQuery.slice(oppositeOffset+span.start,oppositeOffset+span.end) !== observed;
            if (ambiguousOrientation) { forwardSequence=null; reason='Both orientations of this template match the same index span; its orientation is ambiguous.'; }
            if (results.length >= candidateLimit) return finish(true);
            results.push(Object.assign({},base,{start,end,strand:ambiguousOrientation ? 'both' : strand,sequence:observed,forwardSequence,status:forwardSequence ? 'ready' : 'unresolved',reason}));
          }
        }
      }
    }
    return finish(false);
  }
  function annotate(input, catalogue, options) {
    options = options || {};
    const sequence = normalizeInput(input,options);
    const {nodes,components} = compile(catalogue,!!options.includeInferred);
    const boundedLimit = (value, fallback, maximum) => Number.isFinite(Number(value)) && value !== null && value !== undefined ? Math.max(1,Math.min(maximum,Math.floor(Number(value)))) : fallback;
    const maxHits = boundedLimit(options.maxHits,10000,100000);
    const maxProtocols = boundedLimit(options.maxProtocols,20,20);
    const hits = [];
    let state = 0, truncated = false;
    outer: for (let i=0; i<sequence.length; i++) {
      const base = sequence[i];
      while (state && nodes[state].next[base] === undefined) state = nodes[state].fail;
      state = nodes[state].next[base] === undefined ? 0 : nodes[state].next[base];
      for (const match of nodes[state].out) {
        if (hits.length >= maxHits) { truncated=true; break outer; }
        const component = components.get(match.id);
        hits.push({component:match.id,start:i+1-match.length,end:i+1,strand:match.strand,sequence:component.sequence,aliases:component.aliases || [],occurrences:component.occurrences,protocols:component.protocols});
      }
    }
    hits.sort((a,b) => a.start-b.start || b.end-a.end || a.component.localeCompare(b.component));
    const indexes = indexCandidates(sequence,catalogue,!!options.includeInferred);
    return {sequence,hits,overview:nonoverlapping(hits),protocols:rank(hits,catalogue,maxProtocols),indexCandidates:indexes.candidates,indexCandidatesTruncated:indexes.truncated,truncated};
  }
  function createMatcher(catalogue) { return {annotate:(sequence,options) => annotate(sequence,catalogue,options)}; }
  return {normalizeInput,parse:normalizeInput,annotate,createMatcher,revcomp};
}));

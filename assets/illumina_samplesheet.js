/* Sample-sheet field formatting. Inputs must already be in Illumina-forward orientation.
 * Rules checked against Illumina Knowledge #1800 and Connected Software orientation table,
 * 2026-10-21. These are data-section snippets, not complete instrument run configurations.
 */
(function (root, factory) {
  'use strict';
  const api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.IlluminaSampleSheet = api;
})(typeof globalThis !== 'undefined' ? globalThis : this, function () {
  'use strict';
  const ORIENTATION = 'https://knowledge.illumina.com/software/general/software-general-reference_material-list/000001800';
  const TABLE = 'https://help.connected.illumina.com/run-set-up/overview/index-orientation-guide/i5-index-orientation-table';
  const FORMAT = 'https://support-docs.illumina.com/SW/BCL_Convert/Content/SW/BCLConvert/SampleSheets_swBCL.htm';
  const instruments = [
    {id:'miseq',label:'MiSeq (original series)'},
    {id:'miniseq',label:'MiniSeq (standard kits)'},
    {id:'miniseq-rapid',label:'MiniSeq Rapid kits'},
    {id:'iseq100',label:'iSeq 100'},
    {id:'nextseq500',label:'NextSeq 500 / 550'},
    {id:'nextseq2000',label:'NextSeq 1000 / 2000'},
    {id:'novaseq6000',label:'NovaSeq 6000 (v1.5 reagents)'},
    {id:'novaseqx',label:'NovaSeq X / X Plus'},
    {id:'miseqi100-indexfirst',label:'MiSeq i100 — index-first sequencing'},
    {id:'miseqi100-readfirst',label:'MiSeq i100 — read-first sequencing'},
    {id:'hiseq2500',label:'HiSeq 2000 / 2500 (legacy)'},
    {id:'hiseq4000',label:'HiSeq 3000 / 4000 (legacy)'}
  ];
  const workflows = [
    {id:'standalone-bclconvert',label:'Manual sample sheet → standalone / DRAGEN BCL Convert'},
    {id:'manual-bcl2fastq',label:'Manual v1 sample sheet → bcl2fastq / FASTQ Generation'},
    {id:'run-planning',label:'BaseSpace Run Planning — enter indexes into the setup tool'},
    {id:'bssh-bclconvert-v2',label:'BaseSpace BCL Convert app — v2 sample sheet'}
  ];
  const ids = new Set(instruments.map(x => x.id)), modes = new Set(workflows.map(x => x.id));
  const reverseRead = new Set(['miniseq','iseq100','nextseq500','nextseq2000','novaseq6000','novaseqx','miseqi100-readfirst','hiseq4000']);
  const autoConvert = new Set(['nextseq2000','novaseqx','miseqi100-indexfirst','miseqi100-readfirst']);
  const legacy = new Set(['hiseq2500','hiseq4000']);
  function reverseComplement(value) { return value.split('').reverse().map(c => ({A:'T',C:'G',G:'C',T:'A'})[c]).join(''); }
  function describeProfile(instrument, workflow) {
    if (!ids.has(instrument) || !modes.has(workflow)) throw new Error('Choose a supported instrument and workflow.');
    const label = instruments.find(x=>x.id===instrument).label + ' / ' + workflows.find(x=>x.id===workflow).label;
    const result = {supported:true,label,i7Orientation:'forward',i5Orientation:'forward',
      section:workflow==='manual-bcl2fastq'?'Data':'BCLConvert_Data',headerRows:[],sources:[ORIENTATION,TABLE,FORMAT],note:''};
    if (legacy.has(instrument) && workflow!=='manual-bcl2fastq') {
      result.supported=false; result.note='This legacy instrument profile is provided for manual bcl2fastq sheets only. Choose a verified workflow rather than assuming another tool uses the same orientation.'; return result;
    }
    if (workflow==='manual-bcl2fastq') {
      result.i5Orientation=reverseRead.has(instrument)?'reverse_complement':'forward';
      result.note='Paste into a manually prepared v1 [Data] section. This is the expected index orientation, not a statement of software compatibility with every instrument.';
    } else if (workflow==='standalone-bclconvert') {
      result.i5Orientation=reverseRead.has(instrument)&&!autoConvert.has(instrument)?'reverse_complement':'forward';
      result.note=autoConvert.has(instrument)?'BCL Convert uses the instrument RunInfo.xml orientation flag; supply the forward i5 value. Preserve the instrument run files.':'These values assume manual preparation without an IndexOrientation,Forward override.';
    } else if (workflow==='run-planning') {
      // This profile is input to BaseSpace planning, not LRM/IEM or a manual sheet.
      result.note='Enter these forward values into BaseSpace Run Planning. Use the exported sheet unchanged; do not reverse-complement it a second time. This is not a Local Run Manager profile.';
      if (reverseRead.has(instrument)&&!autoConvert.has(instrument)) result.headerRows=[['IndexOrientation','Forward']];
    } else {
      if (reverseRead.has(instrument)&&!autoConvert.has(instrument)) {
        result.headerRows=[['IndexOrientation','Forward']];
        result.note='Required: keep IndexOrientation,Forward in the v2 [Header] section. Without that setting these forward i5 values are not valid for this workflow.';
      } else result.note='The v2 BCL Convert app expects forward i5 values for this instrument.';
    }
    if (instrument==='miniseq-rapid') result.note+=' Rapid-kit sequencing uses forward i5. Local Run Manager has a separate setup exception and is not covered here.';
    return result;
  }
  function bases(value) {
    if (typeof value!=='string' || !/^[ACGT]+$/i.test(value)) throw new Error('An index must contain known A/C/G/T bases only; no placeholder or unresolved index can be exported.');
    return value.toUpperCase();
  }
  function formatIndex(sequence, role, instrument, workflow) {
    const profile=describeProfile(instrument,workflow);
    if (!profile.supported) throw new Error(profile.note);
    if (role!=='i7'&&role!=='i5') throw new Error('Only resolved i7 and i5 sequencer indexes have sample-sheet columns.');
    const input=bases(sequence), orientation=role==='i7'?profile.i7Orientation:profile.i5Orientation;
    return {column:role==='i7'?'index':'index2',value:orientation==='reverse_complement'?reverseComplement(input):input,
      orientation,note:orientation==='reverse_complement'?'Reverse complement of the Illumina-forward index.':'Illumina-forward index; unchanged.'};
  }
  function makeRow({sampleId,i7,i5,instrument,workflow}) {
    if (typeof sampleId!=='string'||!/^[A-Za-z0-9][A-Za-z0-9_-]*$/.test(sampleId)||sampleId.toLowerCase()==='undetermined')
      throw new Error('Sample_ID must start with a letter or number and contain only letters, numbers, underscores or hyphens; Undetermined is reserved.');
    if (!i7) throw new Error('Choose a resolved i7 value for a sample-sheet data row. An i5 value alone can be copied separately.');
    const profile=describeProfile(instrument,workflow), first=formatIndex(i7,'i7',instrument,workflow);
    const headers=['Sample_ID','index'],values=[sampleId,first.value],fields={index:first};
    if (i5) { const second=formatIndex(i5,'i5',instrument,workflow);headers.push('index2');values.push(second.value);fields.index2=second; }
    return {headers,values,tsv:headers.join('\t')+'\n'+values.join('\t'),
      csv:headers.join(',')+'\n'+values.join(','),section:profile.section,headerRows:profile.headerRows,profile,fields};
  }
  return {instruments,workflows,describeProfile,formatIndex,makeRow};
});

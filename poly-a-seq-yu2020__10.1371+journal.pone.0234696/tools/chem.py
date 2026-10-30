"""Protocol-owned source-qualified partial chemistry. No public schematic endpoint."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
from chemdraw import Scene, Segment
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan, Terminal, Disposition
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry

def scene(name, text):
    sc = Scene(); sc.strand(name, [Segment(name, text, placeholder=True)]); return sc

def op(g, key, x, label, kind=S.MAGIC, *, reagent='Source mixture; internal chemistry unassembled', product=None, structural=()):
    catalyst = kind in {S.DIGESTION, S.END_REPAIR, S.EXTENSION, S.TAILING, S.REVERSE_TRANSCRIPTION, S.PCR, S.LIGATION, S.DEPHOSPHORYLATION, S.PHOSPHORYLATION}
    stock = g.input(key + '_reagent', reagent, K.CATALYST if catalyst else K.CONDITION)
    inputs = [(x, R.TEMPLATE), (stock, R.CATALYST if catalyst else R.CONDITION)]
    inputs.extend((a, R.STRUCTURAL) for a in structural)
    kwargs = {}
    if kind == S.MAGIC:
        kwargs['reason'] = 'Source operation established; exact molecular transition/kit internals unassembled. No accepted absence claim.'
    if kind == S.PCR:
        kwargs['geometry'] = PCRGeometry.unavailable(Scene(), reason='Exact template/primer endpoint geometry unassembled; selected oligos do not establish whole amplicon.')
        if len(structural) < 2:
            for arm in ('f', 'r'):
                inputs.append((g.input(key + '_' + arm, 'Source PCR primer; exact identity/landing unassembled', K.OLIGO, state=scene(key + '_' + arm, '[source primer]')), R.STRUCTURAL))
    if kind == S.REVERSE_TRANSCRIPTION and not structural:
        inputs.append((g.input(key + '_primer', 'Source RT primer; geometry unassembled', K.OLIGO, state=scene(key + '_primer', '[source RT primer]')), R.STRUCTURAL))
    if kind == S.LIGATION:
        kwargs['geometry'] = AdapterLigationGeometry('unavailable', 'unavailable', (scene(key + '_ends', '[source ligation; complete substrate/end geometry unassembled]'),), 'This incomplete model does not claim publication absence.')
        if not structural:
            inputs.append((g.input(key + '_adapter', 'Source ligating adapter; full geometry unassembled', K.ADAPTER, state=scene(key + '_adapter', '[source adapter]')), R.STRUCTURAL))
    y = g.product(key + '_product', product or label + '; physical companions retained until separation', K.MOLECULAR)
    g.react(key, label, kind, inputs, (y,), **kwargs); return y

def clean(g, key, x, label):
    condition = g.input(key + '_condition', label, K.CONDITION)
    kept = g.product(key + '_kept', 'Source physically recovered fraction; carryover/yield unquantified', K.MOLECULAR)
    waste = g.product(key + '_waste', 'Actual separated cleanup fractions; no perfect classification assumed', K.WASTE)
    g.react(key, label, S.CLEANUP, ((x, R.TEMPLATE), (condition, R.CONDITION)), (kept, waste))
    g.discard(waste, 'Actual physical separation rather than computational removal.'); return kept

def qc(g, key, x, label):
    kept = g.product(key + '_kept', 'Source retained preparation', K.MOLECULAR)
    assay = g.product(key + '_aliquot', 'Separate analytical aliquot; no return to library assumed', K.QC)
    g.react(key, label, S.SPLIT, ((x, R.TEMPLATE),), (kept, assay))
    g.mark_terminal(Terminal(assay.id, Disposition.QC, 'Source analytical aliquot.')); return kept

SOURCE='PublishedPMC7297374 Methods77–99/Results126–130,299;S1/S2workbookheadersonly. NextSeqScriptSeq andHiSeq/NovaSMARTer actualpreparationskeptseparate.'
def protocol_graph_for(mode='nextseq'):
    if mode not in {'nextseq','hiseq','novaseq'}:raise ValueError('Explicit sourcepreparation/readout route required')
    g=ProtocolGraph('Poly(A)-seq Yu2020 '+mode+' — privatepartialpath',authority=SOURCE,representation='published_instance')
    x=g.input('total_rna','SourcecelltotalRNAafterTrizol/RQ1DNase;5.1ugfragmentedinputNextSeqor5ugtotalinputlaterHeLaroute',K.MOLECULAR,primary=True)
    x=qc(g,'purity',x,'SeparateA260/A280andquantityspectrophotometryassay')
    x=qc(g,'integrity',x,'Separate1.5percentagaroseRNAintegrityassay')
    x=op(g,'rnase',x,'SourceRNaseT1fragmentation37C3minNextSeq;laterrouteconditionsnotassigned',S.DIGESTION,reagent='SourceRNaseT1',product='NativepolyA3prime-end-containingfragments andotherRNAfragmentswithnewendclasses; no biologicalRNAdepletionyet')
    beads=g.input('polya_beads','Sourceoligo-dT-conjugatedmagneticbeads',K.SELECTION);bound=g.product('poly_kept','SourcepolyA-associatedRNAwithunquantifiedcarryover/shorttailbias',K.MOLECULAR);w=g.product('poly_removed','Actualphysicalnonretainedfraction',K.WASTE)
    g.react('capture','PolyAcaptureAFTERpartialRNaseT1fragmentation',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(beads,R.SELECTION)),(bound,w));g.discard(w,'PhysicalcapturenotperfectpolyAclassification');x=bound
    if mode!='nextseq':
        control=g.input('spike_rna','Independentlysynthesized40A/120AspikeRNAs;source30pgor9pgsamplechoiceunassigned;upstreamPCR/BsrDI/T7/PAGEunassembled',K.MOLECULAR,primary=True)
        y=g.product('spiked','PhysicalpolyAselectedRNAplusspikeRNAs',K.MOLECULAR);g.react('spike','SourcespikeadditionAFTERoligo-dTcapture',S.POOL,((x,R.TEMPLATE),(control,R.TEMPLATE)),(y,));x=y
    x=op(g,'three_ligation',x,'SourceGnomeGensRNA-seq3primeadapterligation',S.LIGATION,reagent='HistoricalGnomeGenkit;laterrouteK02420-L;exactadapter/ligationgeometryunassembled',product='3prime-adaptedRNAclassesplusunligatedRNA/freeadapters/reagents;declaredpathdoesnotconvertcleavage3primephosphatesto3primeOH;fullreactivityunassembled')
    x=op(g,'rt',x,'SourceRTprimercomplementaryto3primeadapter',S.REVERSE_TRANSCRIPTION,reagent='HistoricalsourceRTcomponent;enzyme/conditionsunassigned',product='cDNA/RNA/adaptedtemplatesplusRTprimers/reagents;physicalRNAfateunassembled')
    x=op(g,'second',x,'ScriptSeqv2Terminal-TaggingDNA synthesis'if mode=='nextseq'else'SMARTerStrandedRNA-seq634837second-strandcomponent',reagent='HistoricalScriptSeqv2Terminal-Taggingoligo/kit'if mode=='nextseq'else'HistoricalSMARTerStrandedRNA-seqkit634837',product='SourceDNAwithsecondlibraryhandleclasses;exactterminaltaggingversusstrand-synthesismechanismunassembled;notinheritedfromSHERRYorPAT')
    if mode=='nextseq':
        x=clean(g,'cdna_clean',x,'SourcecDNApurificationbeforePCR;reagentdetailunassigned')
        x=op(g,'pcr',x,'Source18cyclePCR',S.PCR,reagent='HistoricalsourcePCRprimer/enzymestockunassigned')
        x=clean(g,'size',x,'Source250–350bpor300–500bpPCRproductpurification;per-samplewindowassignmentunresolved')
        x=qc(g,'quantify',x,'Separatepurifiedlibraryquantitationaliquot;assayvolumeunassigned')
    else:
        x=op(g,'finish',x,'SourceSMARTerlibraryfinishing/PCR/cleanupboundarynotreconstructedfromNextSeqroute',reagent='HistoricalsourceSMARTercomponentoriginalunreviewed',product='SourceactualHiSeq/NovaSeqDNA libraryclasseswithfinishgeometry/reactionorderunassembled; no18cycleassignment')
    if mode=='nextseq':x=op(g,'storage',x,'Sourcepurifiedlibraryminus80Cstorageuntiluse',S.THERMAL_TREATMENT,reagent='Sourceminus80Cstoragecondition',product='StoredphysicalNextSeqlibrarypreparation;no molecularchangeorlossassigned')
    lib=g.product('library','SourcephysicalPoly(A)-seqDNApreparation;kitends/readprimersunassembled',K.LIBRARY)
    g.react('identity','Physicalsequencinginterface',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    readout={'nextseq':'NextSeq500single-end300nt','hiseq':'HiSeqXTenpaired-end150nt','novaseq':'NovaSeq6000paired-end250nt'}[mode]
    g.sequence(lib,SequencingPlan(limitation='Published'+readout+';exactprimer/cyclemappingunassembled. ThefirstDNAstrandwithpolyTservesastemplateforpolyAsequencecalls;calledAisnotchemicalT-to-Aconversion. pA-finder/polyG/qualityfiltersoccurafteracquisition.'))
    g.add_source_item(SourceItem('pending','HistoricalGnomeGen/ScriptSeqv2/SMARTeroriginals;fulljunctions/spikepreparation/actualsamplewindow/readgeometry','Selectedpublishedsourceonly;notabsence.'));return g

def protocol_graph():return protocol_graph_for()
def selftest():
    from checks import Check
    c=Check();a=protocol_graph();b=protocol_graph_for('hiseq');d=protocol_graph_for('novaseq')
    c('PolyAcapturefollowsRNaseT1',a.steps['capture'].inputs[0].material,'rnase_product')
    c('NextSeqhasnoinventedspike','spike'in a.steps,False)
    c('LaterRNAspikeaddedaftercapture',b.steps['spike'].inputs[0].material,'poly_kept')
    c('Laterrouteavoidseighteencycleinheritance','pcr'in b.steps,False)
    c('NextSeqcDNApurifiedbeforePCR',a.steps['pcr'].inputs[0].material,'cdna_clean_kept')
    c('PhysicalnonretainedRNAremoved',a.terminals['poly_removed'].disposition,Disposition.DISCARD)
    c('Incompletecommercialchemistrypreventspublication',bool(d.publication_problems()),True)
    c.report()
if __name__=='__main__':selftest()

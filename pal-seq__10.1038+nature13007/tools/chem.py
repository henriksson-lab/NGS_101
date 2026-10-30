"""Protocol-owned source-qualified partial chemistry. No public schematic endpoint."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
from chemdraw import Scene, Segment, revcomp
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan, Terminal, Disposition
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry

def scene(name, text):
    sc = Scene(); sc.strand(name, [Segment(name, text, placeholder=True)]); return sc

def oligo(g, key, name, bases, **mods):
    parts = []
    for i, block in enumerate(bases.split('[BtndT]')):
        if i: parts.append(Segment('Internal biotin-dT', 'T', note='Source [BtndT]; physical modification not inherited by copies'))
        if block: parts.append(Segment(key + str(i), block, placeholder=bool(set(block)-set('ACGTU')), note='IUPAC synthesis pattern retained; full landing unassembled'))
    sc = Scene(); sc.strand(key, parts, **mods)
    return g.input(key, name, K.OLIGO, state=sc)

def op(g, key, x, label, kind=S.MAGIC, *, reagent='Source mixture; internal chemistry unassembled', product=None, structural=()):
    catalyst = kind in {S.DIGESTION, S.END_REPAIR, S.EXTENSION, S.TAILING, S.REVERSE_TRANSCRIPTION, S.PCR, S.LIGATION}
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

def finish(g, x, label, limitation):
    lib = g.product('library', label, K.LIBRARY)
    g.react('library_identity', 'Physical library interface; complete molecular state unassembled', S.GENERIC, ((x, R.TEMPLATE),), (lib,), changes_molecule=False)
    g.sequence(lib, SequencingPlan(limitation=limitation))
    g.add_source_item(SourceItem('pending', 'Historical exact kit, substrate and sequencing geometry', 'Private partial graph; no full source audit or accepted absence.'))
    return g

SOURCE = 'AssociatedAMPMC4086860 Methods425–482 and originalPDFp8/9 rendered. Publishedfinal first-partyPDF retrievalcurl60; final applicability unreviewed. No preprint.'
from illumina import INDEX2_PRIMER_RC,P5,TRUSEQ_READ1
CAPTURE = INDEX2_PRIMER_RC + 'AGACACATAC'
SPLINT = 'TTCCGATCTTTTTTTTT'
RT = P5 + TRUSEQ_READ1[4:18]
LEFT_DNA = 'CAAGCAGAAGACGGCATACGAGTTCAGAGTTCTA'
LEFT_RNA = 'CAGUCCGACGAUC'

@dataclass(frozen=True)
class RetainedTailCopy:
    body: str
    tail_length: int
    def __post_init__(self):
        if not self.body or set(self.body)-set('ACGU') or type(self.tail_length)is not int or self.tail_length<0:
            raise ValueError('ExplicitselectedRNAbodyandtailrequired')
    def extension(self):
        return 'T'*self.tail_length + revcomp(self.body.replace('U','T'))
    def tail_removal(self):
        raise ValueError('PAL selected preparation retains tail for copying; no3P-seqRNaseHstep')
@dataclass(frozen=True)
class MixedLeftAdapter:
    dna: str=LEFT_DNA
    rna: str=LEFT_RNA
    def scene(self):
        sc=Scene();sc.strand('left',[Segment('DNAarm',self.dna),Segment('RNAarm',self.rna,note='Physicalsource lowercaseRNA')],mod5='C3.spacer');return sc

def protocol_graph():
    g=ProtocolGraph('PAL-seq — selected associated2014manuscriptpreparation',authority=SOURCE,representation='published_instance')
    x=g.input('rna','Source totalRNAorcytoplasmicallyenrichedRNA plus tail-lengthstandards/traceisotopeRNA',K.MOLECULAR,primary=True)
    adapter=g.input('capture_adapter','Source5Pall-DNAadapter3biotin',K.ADAPTER,state=scene('capture','[exactcaptureadapterbelow]'))
    adapter.state.strands.clear();adapter.state.strand('capture',[Segment('CaptureDNA',CAPTURE)],mod5='phosphate',mod3='biotin')
    splint=oligo(g,'splint','SourceDNApolyAsplint',SPLINT)
    x=op(g,'capture_ligation',x,'SourceT4Rnl2splintedRNA3endligation18Covernight',S.LIGATION,reagent='T4Rnl2',structural=(adapter,splint),product='SelectedRNA-DNAcaptureadapterclasses retainingnativepolyAtailplusunligatedRNA/freeadapter/splint/reagents')
    x=op(g,'t1',x,'SourcepartialRNaseT1digestion',S.DIGESTION,reagent='RNaseT1; detaileddirectlycited3Pconditionsunassembled',product='SelectedpolyAtail/captureadapterRNAfragments plus otherRNAfragments/freeDNAadapter/splint/remnants; noRNAfragmentdiscard')
    x=clean(g,'rna_gel',x,'Sourcephenolchloroform/ethanolprecipitation/denaturingPAGE104–750nt; physicallyremovesresidualunreactedadapter')
    beads=g.input('beads','StreptavidinM280beads',K.SELECTION);bound=g.product('bound','BiotinboundRNA/DNA/beadcomplexesandcarryover',K.MOLECULAR);other=g.product('unbound','Separatedcaptureunboundfraction',K.WASTE)
    g.react('capture','SourceM280capture',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(beads,R.SELECTION)),(bound,other));g.discard(other,'Physicalcaptureselection; notperfectpolyAclassification')
    x=op(g,'pnk',bound,'Sourceon-beadRNA5primephosphorylation',S.END_REPAIR,reagent='3prime-phosphatase-deficientT4PNK',product='SourceRNA5Pclasseswith3endcaptureadapter/beads,companionDNA/splint/reagents retained; no3primephosphateremovalasserted')
    left=g.input('left_adapter','SourceC3.spacer mixedDNA/RNA5primeadapter',K.ADAPTER,state=MixedLeftAdapter().scene(),chemical_state=MixedLeftAdapter())
    x=op(g,'left_ligation',x,'SourceT4Rnl1RNA5primeadapterligation22Covernight',S.LIGATION,reagent='T4Rnl1',structural=(left,),product='Mixedleftadapter-RNA-polyA-DNAcaptureadapterclassesplusunligatedclasses/freeoligos/beads/reagents')
    primer=oligo(g,'rt_primer','SourceSSII RTprimer',RT)
    x=op(g,'rt',x,'SourceSuperScriptIIreverse transcription',S.REVERSE_TRANSCRIPTION,reagent='SuperScriptII; exactconditions/wholelandingunassembled',structural=(primer,),product='Tail-copyingfirststrandDNAclassesplusoriginalRNA/DNAcaptureadapter,unextendedprimers,freeleftadapter,beadsandremnants; fullnativepolyAkeptthroughcopying')
    x=op(g,'hydrolysis',x,'SourcebasehydrolysisliberatescDNAfrombeads',reagent='Sourcebasehydrolysis; exactconditionsunreviewed',product='LiberatedcDNAplusRNAhydrolysisproducts,originalDNAcaptureadapter/freeprimers/beads/remnantsbeforePAGE')
    x=clean(g,'cdna_gel',x,'Source denaturingPAGE166–790nt DNArecovery')
    x=op(g,'denature',x,'Source roomtemperature5–100mMNaOH thenHT1neutralization',reagent='NaOH/IlluminaHT1')
    x=op(g,'cluster',x,'SourcecBotclustergeneration/linearization/3endblocking/primerhybridization',reagent='HistoricalcBotflowcell/reagents/primer; fullmolecularclustergeometryunassembled',product='PhysicalflowcellDNAclusterspluswash/reagentcompanionswithfatesunassembled')
    x=op(g,'tail_extension',x,'SourceKlenowdTTP/biotin16dUTP extension37C30min',S.EXTENSION,reagent='Klenow200nMdTTP/10nMyeastor4nMotherbiotin16dUTP',product='Sourcebiotinylatedtail-copyprimerextensionclassesonclusterDNAplusunextended/partialclasses/reagents; stochasticincorporation/yieldnotmodelled')
    return finish(g,x,'PhysicalPALflowcellDNAtemplate/extendedprimerinterface; exactgeometryunassembled','SourceGenomeAnalyzerII36cycles SBS follows tail-extension. Fluorescentstreptavidin tail-length imaging occurs after SBS; fluorescence is separatequantitation,notbasecalledtailsequence. Source35ntArabidopsisbranch excludeslowqualityfirstcycle. Exactcluster/primerlandingandfullfatesunassembled.')

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();a=RetainedTailCopy('ACGU',7);c('Selected firstcopy retains fulltail length',a.extension(),'T'*7+'ACGT')
    try:a.tail_removal();rejected=False
    except ValueError:rejected=True
    c('PAL rejects inserting3Ptailremoval',rejected,True)
    g=protocol_graph();c('Onbeadphosphorylationactsaftercapture',g.steps['pnk'].inputs[0].material,'bound')
    c('NativeRNAhydrolysisonlyaftercopying',g.steps['hydrolysis'].inputs[0].material,'rt_product')
    c('Polymerase tail-extensionfollowsclustergeneration',g.steps['tail_extension'].inputs[0].material,'cluster_product')
    c('FragmentationdoesnotdiscardRNAbyitself',g.steps['t1'].outputs,('t1_product',))
    c('Partialgraphcannotpublish',bool(g.publication_problems()),True)
    src=Source(sources_root()/'pal-seq__10.1038+nature13007'/'nature13007_PMC4086860.xml.txt','get_sources --slugs pal-seq__10.1038+nature13007')
    if have(c,src,label='PALassociatedAMMethods'):
        t=src.path.read_text();c('CaptureDNAtranscription',CAPTURE in t,True);c('ActualRTprimertranscription',RT in t,True)
        c('Internaltraceasteriskis32PlabelnotPSbond','32P-label (*)' in t,True)
    c.report()
if __name__=='__main__':selftest()

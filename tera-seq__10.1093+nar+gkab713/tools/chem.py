"""Protocol-owned source-qualified partial chemistry. No public schematic endpoint."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
from chemdraw import Scene, Segment, revcomp
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

from dataclasses import dataclass
from illumina import P5, TRUSEQ_READ1

SOURCE='Published PMC8599856 Methods59–73/Results147/163; publishedTableS1 XLSXwholeadapterrows/underlinedRTAarm; supplementaryPDFidentity/Note1–2 only. HistoricalRNA001/002internalsunreviewed.'
# Source TableS1 has this canonical base order, but it is an RNA adapter for ONT,
# not an Illumina primer or sample/cell identifier.
ADAPTER5=(P5+TRUSEQ_READ1[4:]).replace('T','U')
RNA3='GUGUCAGUCACUUCCAGC'
DNA3='GG'
RTA_TOP='GGCTTCTTCTTGCTCTTAGGTAGTAGGTTC'
RTA_CORE='GGCTTCTTCTTGCTCTTAGG'
RTA_BOTTOM='GAGGCGAGCGGTCAATTTTCCTAAGAGCAAGAAGAAGCCCCGCTGGAAGTGACTGACAC'

@dataclass(frozen=True)
class MixedAdapter:
    rna:str=RNA3
    dna:str=DNA3
    def __post_init__(self):
        if not self.rna or set(self.rna)-set('ACGU') or not self.dna or set(self.dna)-set('ACGT'):raise ValueError('Explicit ordered RNA and DNA blocks required')
    def scene(self):
        sc=Scene();sc.strand('adapter3',[Segment('Source RNA block',self.rna,note='Ribonucleosides'),Segment('Source last two DNA bases',self.dna,note='Deoxyribonucleosides')],mod5='phosphate');return sc
    def flatten_as_dna(self):raise ValueError('Source adapter is mixed RNA/DNA, not ordinary DNA')

def adapter5_scene(hydroxyl=False):
    sc=Scene();sc.strand('adapter5',[Segment('Source58ntRNA',ADAPTER5,note='Source TableS1 allRNA; identical Illumina-like base order does not confer Illumina use')],mod5='photocleavable biotin',mod3='phosphate' if hydroxyl else '');return sc

def rta_scene():
    # Source-underlined20nt at the bottom3primeend pairs to the source3primeRNA/DNAadapter.
    tail=revcomp((RNA3+DNA3).replace('U','T'))
    core_partner=revcomp(RTA_CORE)
    if not RTA_TOP.startswith(RTA_CORE) or not RTA_BOTTOM.endswith(core_partner+tail):raise ValueError('Source duplex/adapter-binding arms disagree')
    lead=RTA_BOTTOM[:-len(core_partner+tail)]
    sc=Scene();sc.strand('top',[Segment('Known duplex core',RTA_CORE),Segment('Source top3primeextension',RTA_TOP[len(RTA_CORE):])],mod5='phosphate')
    sc.anneal('bottom',[Segment('Source bottom5primeextension',lead),Segment('Core complement',core_partner),Segment('Underlined adapter-binding3primearm',tail)],to='top',pair=('Core complement','Known duplex core'),unpaired=('Source bottom5primeextension','Underlined adapter-binding3primearm'))
    return sc

def affinity(g,key,x,label,*,retained):
    beads=g.input(key+'_beads',label,K.SELECTION);k=g.product(key+'_kept',retained,K.MOLECULAR);w=g.product(key+'_removed','Actual physically separated capture fraction; nonideal partition unquantified',K.WASTE)
    inputs=[(x,R.TEMPLATE),(beads,R.SELECTION)]
    if key=='deplete':inputs.append((g.input('depletion_probes','SourceantisensebiotinDNAprobepool;exactsequenceinventoryunreviewed',K.OLIGO,state=scene('probes','[source depletion probes]')),R.STRUCTURAL))
    g.react(key,label,S.AFFINITY_CAPTURE,inputs,(k,w));g.discard(w,'Physical capture separation; no perfect biological specificity assumed');return k

def protocol_graph_for(mode='both'):
    if mode not in {'5p','5oh','cap','cap_5p','tera3','both'}:raise ValueError('Explicit source TERA branch required')
    if mode in {'tera3','both'}:rta_scene()
    g=ProtocolGraph('TERA-seq '+mode+' — selected published path',authority=SOURCE,representation='published_instance')
    x=g.input('total_rna','SourceHeLatotalRNA afterTrizol/DNase; extraction/predecessorinternalsunassembled',K.MOLECULAR,primary=True)
    x=qc(g,'rna_integrity',x,'SeparateBioanalyzerRNAintegrityassay beforelibrarypreparation')
    if mode in {'tera3','both'}:
        x=affinity(g,'deplete',x,'SourceantisensebiotinDNAprobe/streptavidinrRNA/smallRNAdepletion; referencedcomponentunreviewed',retained='Source physicalunboundRNAfraction retained;boundrRNA/probesremoved without perfectdepletionclaim')
    else:
        x=affinity(g,'polya',x,'Sourceoligo-dTDynabeadspolyAselection;kitcomponentunreviewed',retained='PolyA-associatedRNA/beadfractionwithnativecap/5P/5OHclassesandcarryover')
    if mode=='cap':
        x=op(g,'cip',x,'OnbeadQuickCIPdephosphorylation',S.DEPHOSPHORYLATION,reagent='QuickCIP/historicalconditionsunreviewed',product='Noncappedaccessible5Pclassesdephosphorylated; cappedRNAandotherclassesretained')
        x=clean(g,'cip_wash',x,'ThreeWashingBufferB/oneThermopolwash')
    if mode in {'cap','cap_5p','both'}:
        x=op(g,'decap',x,'RppH37C70min200ul' if mode!='both' else 'RppH37C75min350ul',S.END_REPAIR,reagent='RppH/Thermopolbuffer;20ulstockfor5TERAor35ulforboth',product='DecappedRNAclasses pluspreexistingendclasses/cleavageproducts/reagents; perfectselectionnotpredicted')
        x=op(g,'edta',x,'Source1ul500mMEDTAquench',S.QUENCHING,reagent='SourceEDTA')
        x=clean(g,'decap_clean',x,'SourceRNACleanXP1x/20ulelution' if mode=='both' else 'ThreeWashingBufferB/oneRNL1bufferwash')
    if mode!='tera3':
        a=g.input('five_adapter','Published5primePCBio58ntRNA'+('/3PhosforRtcB' if mode=='5oh' else ';free3primeRNAterminus'),K.ADAPTER,state=adapter5_scene(mode=='5oh'))
        x=op(g,'five_ligation',x,'SourceRtcB5OHligation;historicalconditionsunreviewed' if mode=='5oh' else ('T4RNL1 37C3h7.5percentPEG' if mode=='both' else 'OnbeadT4RNL1 37C3h12.5percentPEG'),S.LIGATION,reagent='RtcB/sourcebuffer' if mode=='5oh' else 'T4RNAligase1/sourcePEGbuffer',structural=(a,),product='5prime-adaptedRNAclasses plusunligatedRNA/adapters/reagents; adaptersareRNA,notPCRprimers')
        if mode=='both':
            x=affinity(g,'five_capture',x,'MyOneC1streptavidincaptureof5prime-adaptedRNA',retained='PCBio-associatedRNA/beadswithcarryover; actualunboundmaterialseparated')
            x=clean(g,'five_wash',x,'ThreeBW/oneRNAreactionbufferwash')
        else:
            x=clean(g,'five_wash',x,'ThreeWashingBufferBwashes; retainedpolyAbeadfraction')
            x=clean(g,'five_release',x,'75C2min50ulwaterelution thenRNACleanXPrecovery;1.8xfor5Psource')
    if mode in {'tera3','both'}:
        a=g.input('three_adapter','Published5primephosphate18RNA+2DNAadapter',K.ADAPTER,state=MixedAdapter().scene(),chemical_state=MixedAdapter())
        x=op(g,'three_ligation',x,'T4RNL1 37C3h'+('12.5percentPEGonbeads' if mode=='both' else '7.5percentPEGsolution'),S.LIGATION,reagent='T4RNAligase1/sourcePEGbuffer',structural=(a,),product='3prime-adaptedRNAclasses plusnonligatedRNA/adapters/reagents; mixedRNA/DNAadapterretained')
        if mode=='both':
            x=clean(g,'three_wash',x,'ThreeBWwashes')
            x=clean(g,'release',x,'50ul95percentformamide/5mMEDTA65C5minrelease then1.8xRNACleanXP/9.5ulelution')
        else:x=clean(g,'three_clean',x,'1xRNACleanXP/9.5ulelution')
        rta=g.input('rta','PublishedcustomRTADNAduplex; fullsourceoligosretained',K.ADAPTER,state=rta_scene())
    else:rta=g.input('rta','HistoricalRNA001/002 RTAstock;wholegeometryunreviewed',K.ADAPTER,state=scene('rta','[historical RTA]'))
    x=op(g,'rta_ligation',x,'SourceT4DNAligase30C15minRTAligation',S.LIGATION,reagent='T4DNAligase/Quickligationbuffer',structural=(rta,),product='RTA-associatedadaptedRNAclasses plusotherRNA/adapters/reagents; fullRNA/RTAjunctiongeometryunassembled')
    x=op(g,'rt',x,'SourceSSIIIreverse transcription; noPCR',S.REVERSE_TRANSCRIPTION,reagent='SSIII/historicalRNA001/002reactionconditionsunreviewed',product='NativeRNAwithcDNAcompanionclasses;RNAremainssequencedmolecule,notreplacementbyDNA')
    x=clean(g,'rt_clean',x,'RNACleanXPrecovery/20ulelution')
    x=op(g,'motor_adapter',x,'T4DNAligase30C15minsequencingadapterligation',S.LIGATION,reagent='T4DNAligase/Quickligationbuffer',product='Motor-adapter-associatedRNA/cDNAclasses plusresidualsubstrates/adapters/reagents')
    x=clean(g,'final_clean',x,'RNACleanXPrecovery/21ulElutionBuffer')
    x=qc(g,'qubit',x,'Separate1ulQubitDNAHSassay;retainedlibraryfornativeRNAacquisition')
    lib=g.product('library','PhysicalnativeRNA-containingTERAlibrary;completeONTmotorgeometryunassembled',K.LIBRARY)
    g.react('identity','Physicalacquisitioninterface',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='PublishedMinION/R9.4FLO-MIN10648or72h;directRNA3prime-to5primereadout,notIllumina/PCR. AdapterdetectionandRNAselectioninsoftware occurafteracquisition.'))
    g.add_source_item(SourceItem('pending','HistoricalRNA001/RNA002/Dynabeads/predecessorprobeoriginals andfullmixedend/motorgeometry','Selectedpublishedmain/TableS1;fullsource/linkedreviewincomplete,noabsenceclaim.'));return g

def protocol_graph():return protocol_graph_for()
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();a=protocol_graph_for('5p');b=protocol_graph_for('cap');h=protocol_graph_for('5oh')
    c('Bothendcaptureonlyafter5adapterligation',g.steps['five_capture'].inputs[0].material,'five_ligation_product')
    c('5TERAdoesnotinventstreptavidin5adapterenrichment','five_capture'in a.steps,False)
    c('CaprouteCIPprecedesdecapping',b.steps['decap'].inputs[0].material,'cip_wash_kept')
    c('BothroutehasnoinventedCIP','cip'in g.steps,False)
    c('5OHroutehasnoinventeddecapping','decap'in h.steps,False)
    c('RNAretainedafterRT',g.steps['rt_clean'].inputs[0].material,'rt_product')
    try:MixedAdapter().flatten_as_dna();bad=False
    except ValueError:bad=True
    c('MixedadapterrejectsordinaryDNAcollapse',bad,True)
    c('Partialgraphcannotpublish',bool(g.publication_problems()),True)
    src=Source(sources_root()/'tera-seq__10.1093+nar+gkab713'/'nar_gkab713_supp_gkab713_supplemental_files_Supplementary_Table_1-Updated3.xlsx.txt','get_sources --slugs tera-seq__10.1093+nar+gkab713')
    if have(c,src,label='publishedTableS1'):
        t=src.path.read_text();c('WholeRTAtoptranscription',RTA_TOP in t,True);c('WholeRTAbottomtranscription',RTA_BOTTOM in t,True);c('Source3PhosRtcBvariant','/3Phos/'in t,True);c('Sourcemixedadapterliteral','5Phos/rGrUrGrUrCrArGrUrCrArCrUrUrCrCrArGrCGG'in t,True)
    c.report()
if __name__=='__main__':selftest()

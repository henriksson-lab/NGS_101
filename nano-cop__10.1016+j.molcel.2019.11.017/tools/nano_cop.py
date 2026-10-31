"""nano-COP native-RNA tail and supportive-cDNA states, qualified to inspected manuscript."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph, MaterialKind, InputRole, StepKind, Terminal, Disposition
from adapter_ligation import AdapterLigationGeometry
SOURCE='Drexler et al.2020 doi:10.1016/j.molcel.2019.11.017, accepted manuscript JATS54–56,160,176; VOR applicability unverified'

@dataclass(frozen=True)
class NativeRNA:
    artificial_tail: str='none'
    supporting_cdna: bool=False
    splint: str='none'
    def __post_init__(self):
        if self.artificial_tail not in {'none','A','I'}: raise ValueError('native tail identity must remain A or inosine, never generic DNA')
        if self.splint not in {'none','T10','C10'}: raise ValueError('unknown source splint')
        if self.splint=='C10' and self.artificial_tail!='I': raise ValueError('C10 branch requires the explicitly selected poly(I) RNA')
    @property
    def sequenced_polymer(self): return 'RNA'
    def tailed(self,kind):
        if self.artificial_tail!='none' or kind not in {'A','I'}: raise ValueError('select one explicit source tail route')
        return replace(self,artificial_tail=kind)
    def ligated(self):
        if self.artificial_tail=='none': raise ValueError('unassigned endogenous tail in no-tailing control is not a universally ligatable substrate')
        return replace(self,splint='C10' if self.artificial_tail=='I' else 'T10')
    def reverse_transcribed(self):
        if self.splint=='none': raise ValueError('source RT follows adapter ligation')
        return replace(self,supporting_cdna=True)


def graph_for_tail(tail):
    if tail not in {'A','I','none'}: raise ValueError('source branches are A,I,none')
    g=ProtocolGraph('nano-COP — '+tail+' tail branch',authority=SOURCE,representation='published_instance')
    rna=g.input('labelled_chr_rna','chromatin-associated RNA after8min4sU pulse and fractionation',MaterialKind.MOLECULAR,primary=True)
    reagent=g.input('biotin_hpdp','EZ-Link Biotin-HPDP reaction24°C1.5h',MaterialKind.CONDITION)
    tagged=g.product('biotin_mixture','biotinylated RNA plus soluble reaction products',MaterialKind.MOLECULAR)
    g.react('biotinylate','4sU-dependent biotinylation; retain reaction mixture before separation',StepKind.CHEMICAL_CONVERSION,
            ((rna,InputRole.TEMPLATE),(reagent,InputRole.CONDITION)),(tagged,))
    clean=g.product('clean_tagged_rna','RNA after chloroform/isopropanol recovery',MaterialKind.MOLECULAR)
    waste=g.product('tagging_waste','extraction phases and precipitation supernatant',MaterialKind.WASTE)
    g.react('recover_tagged','Actual extraction and isopropanol precipitation',StepKind.CLEANUP,((tagged,InputRole.TEMPLATE),),(clean,waste));g.discard(waste,'Physical separation, not reaction-generated automatic waste')
    beads=g.input('streptavidin','μMACS streptavidin beads',MaterialKind.SELECTION)
    bound=g.product('bound_rna','biotinylated4sU RNA on streptavidin',MaterialKind.MOLECULAR)
    unbound=g.product('unbound','affinity flow-through and washes',MaterialKind.WASTE)
    g.react('capture','Streptavidin affinity capture and source washes',StepKind.AFFINITY_CAPTURE,
            ((clean,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(bound,unbound));g.discard(unbound,'Physical affinity reject')
    dtt=g.input('dtt','0.1M DTT elution',MaterialKind.CONDITION)
    eluate=g.product('released_rna','DTT-eluted RNA-containing mixture; no RNA-base rewrite',MaterialKind.MOLECULAR)
    g.react('release','DTT release of selected RNA',StepKind.CHEMICAL_CONVERSION,
            ((bound,InputRole.TEMPLATE),(dtt,InputRole.CONDITION)),(eluate,))
    purified=g.product('purified_rna','miRNeasy-purified RNA with on-column DNase',MaterialKind.MOLECULAR,chemical_state=NativeRNA())
    reject=g.product('purification_waste','column and cleanup reject',MaterialKind.WASTE)
    g.react('purify','miRNeasy cleanup after elution',StepKind.CLEANUP,((eluate,InputRole.TEMPLATE),),(purified,reject));g.discard(reject,'Physical column reject')
    kit=g.input('rrna_kit','source Ribo-Zero Gold or RiboMinus v2; per-sample assignment in gated TableS1',MaterialKind.SELECTION)
    depleted=g.product('depleted','rRNA-depleted RNA population; kit molecular mechanics unreviewed',MaterialKind.MOLECULAR,chemical_state=NativeRNA())
    removed=g.product('rrna_material','depletion-rejected molecular material; physical final fate unreviewed',MaterialKind.MOLECULAR)
    g.react('rrna_deplete','Source kit rRNA depletion; molecular implementation unresolved',StepKind.MAGIC,
            ((purified,InputRole.TEMPLATE),(kit,InputRole.SELECTION)),(depleted,removed),reason='The original names two rRNA depletion kits; molecular steps and rejected-material fate remain unreviewed.')
    g.mark_terminal(Terminal(removed.id,Disposition.UNRESOLVED,'Source rRNA kit internals and rejected-material fate unreviewed.'))
    endpoint=depleted
    if tail!='none':
        enzyme=g.input('pap','yeast poly(A) polymerase37°C30min' if tail=='I' else 'E.coli PAP/NEB37°C1h (unperturbed branch)',MaterialKind.CATALYST)
        ntp=g.input('ntp','ITP' if tail=='I' else 'ATP required for poly(A) synthesis; exact buffer not curated',MaterialKind.MOLECULAR,inferred=tail=='A')
        tailed_state=NativeRNA().tailed(tail)
        tailed=g.product('tailed','native RNA with artificial poly('+tail+') tail; tail length unassigned',MaterialKind.MOLECULAR,chemical_state=tailed_state)
        g.react('tail','Source enzymatic RNA tailing',StepKind.TAILING,
                ((depleted,InputRole.TEMPLATE),(enzyme,InputRole.CATALYST),(ntp,InputRole.STRUCTURAL)),(tailed,))
        adapter=g.input('adapter','custom C10 splinted adapter (TableS4)' if tail=='I' else 'SQK-RNA001/002 T10 splinted adapter',MaterialKind.ADAPTER)
        ligase=g.input('ligase','T4 DNA ligase',MaterialKind.CATALYST)
        joined_state=tailed_state.ligated()
        joined=g.product('adapter_rna','RNA with source splinted3′adapter; full adapter sequence unreviewed',MaterialKind.MOLECULAR,chemical_state=joined_state)
        g.react('ligate','Source3′ adapter ligation15min at room temperature',StepKind.LIGATION,
                ((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(joined,),
                geometry=AdapterLigationGeometry.unavailable('Reported C10/poly(I) or T10/poly(A) recognition is known; full source oligos/end geometry not reviewed (custom TableS4 gated).'))
        rt=g.input('rt','SuperScript III',MaterialKind.CATALYST)
        rt_primer=g.input('rt_primer','source RT primer role; exact kit oligo identity unreviewed',MaterialKind.OLIGO,inferred=True)
        endpoint=g.product('supported_rna','native RNA with supportive cDNA; RNA is the sequenced polymer',MaterialKind.MOLECULAR,chemical_state=joined_state.reverse_transcribed())
        g.react('supportive_rt','Source supportive RT; resulting cDNA is not sequenced',StepKind.REVERSE_TRANSCRIPTION,
                ((joined,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(endpoint,))
    g.mark_terminal(Terminal(endpoint.id,Disposition.UNRESOLVED,'Source SQK-RNA001/002 motor-adapter attachment, complete chemical state and acquisition geometry remain unreviewed; no current ONT kit substituted.'))
    g.declare_reviewed_scope(('labelled_chr_rna','biotin_hpdp','biotinylate','recover_tagged','streptavidin','capture','dtt','release','purify','rrna_kit','rrna_deplete'),SOURCE)
    return g

def protocol_graph(): return selected_graph()

def selftest():
    from checks import Check
    c=Check();s=NativeRNA().tailed('I').ligated().reverse_transcribed()
    c('inosine tail retains C10 recognition and native RNA readout after supportive RT',
      (s.artificial_tail,s.splint,s.supporting_cdna,s.sequenced_polymer),('I','C10',True,'RNA'))
    try: NativeRNA().ligated()
    except ValueError: rejected=True
    else: rejected=False
    c('no-tail control does not invent universal adapter-compatible native ends',rejected,True)
    c('poly(A) branch retains separate T10 state',graph_for_tail('A').materials['supported_rna'].chemical_state.splint,'T10')
    c.report()



PAGE_FILE='nano_cop.html'
SELECTED_SCOPE=('labelled_chr_rna','hpdp','tag_conditions','biotinylate','tag_recovery','recover_tagged','streptavidin','capture_conditions','capture','dtt','release_conditions','release','column_reagents','column_dnase','purify','rrna_kit','rrna_deplete','pap','itp','tail_conditions','tail','ligation_allocate','custom_adapter_stock','adapter_inventory','adapter_buffer','adapter_prepare','ligase','ligation_supplies','ligation_conditions','ligate','ssiii','rt_supplies','supportive_rt','first_xp','support_clean','loaded_motor_adapter','motor_inventory','motor_ligation_supplies','motor_ligate','second_xp','motor_clean','elution_buffer','elute','running_buffer','mix_running','loading_allocate','flowcell','loaded_rna')

def selected_graph():
    from protocol_graph import MaterialKind as K,InputRole as R,StepKind as S,SequencingPlan
    from chemdraw import Scene,Segment
    g=ProtocolGraph('nano-COP — source poly(I) native-RNA sequencing',authority=SOURCE,representation='published_instance')
    def inp(k,label,kind=K.MOLECULAR,**kw):return g.input(k,label,kind,**kw)
    def clean(k,x,supply,label):
        keep=g.product(k+'_kept',label+' retained RNA-containing material; recovery/carryover unassigned',K.MOLECULAR)
        waste=g.product(k+'_waste','Actual '+label+' physical rejects',K.WASTE)
        g.react(k,label,S.CLEANUP,((x,R.TEMPLATE),(supply,R.SELECTION)),(keep,waste));g.discard(waste,'Actual source physical workup; no perfect molecular depletion asserted');return keep
    rna=inp('labelled_chr_rna','Source prepared4sU-labelled chromatin-associated RNA;8min pulse/fractionation upstream of this boundary',primary=True)
    hpdp=inp('hpdp','Actual EZ-Link Biotin-HPDP molecular contribution:20percent stock1mg/ml inDMF; occupancy unassigned')
    cond=inp('tag_conditions','Source1ugRNA/10ul,10percent biotinylation buffer100mMTris7.5/10mMEDTA;24C1.5h800rpm dark',K.CONDITION)
    tagged=g.product('tagged','Biotin-treated4sU RNA with freeHPDP/reaction products; detailed adduct occupancy unassembled',K.MOLECULAR)
    g.react('biotinylate','Source4sU-dependent biotin treatment',S.CHEMICAL_CONVERSION,((rna,R.TEMPLATE),(hpdp,R.STRUCTURAL),(cond,R.CONDITION)),(tagged,))
    x=clean('recover_tagged',tagged,inp('tag_recovery','Source1:1chloroform/isoamylalcohol24:1 phase-lock separation16000g5min thenisopropanol precipitation',K.SELECTION),'Source post-biotin recovery')
    bead=inp('streptavidin','Source microMACS streptavidin bead contribution at2:1volume relation',K.SELECTION)
    capture_conditions=inp('capture_conditions','24C15min800rpm;columnwashes65C androomtemperature3each',K.CONDITION)
    kept=g.product('bound_rna','Source streptavidin-associated biotin-treatedRNA',K.MOLECULAR);waste=g.product('capture_waste','Actual affinity flow-through and washes',K.WASTE)
    g.react('capture','Source streptavidin capture/column washing',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(bead,R.SELECTION),(capture_conditions,R.CONDITION)),(kept,waste));g.discard(waste,'Actual capture/wash fractions; selectedRNA purity unassigned')
    dtt=inp('dtt','Source0.1MDTT molecular reductant forRNA release; detailed reducedadduct unassembled')
    release_conditions=inp('release_conditions','Source0.1MDTT elution; exact time/temperature unassigned',K.CONDITION)
    released=g.product('released','DTT-releasedRNA andreaction components; no RNA-base rewrite to unmodifiedU',K.MOLECULAR)
    residual=g.product('capture_residual','Residual bead-associated material after sourceelution; amount/finalfate unassigned',K.MOLECULAR)
    g.react('release','Source DTT release from captured preparation',S.CHEMICAL_CONVERSION,((kept,R.TEMPLATE),(dtt,R.STRUCTURAL),(release_conditions,R.CONDITION)),(released,residual));g.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,'Residual affinity material outside represented RNA continuation'))
    column=inp('column_reagents','Source miRNeasyMicro cleanup supplies',K.SELECTION);dnase=inp('column_dnase','Actual source on-column DNaseI treatment',K.CATALYST)
    purified=g.product('purified','Source column-recoveredRNA; residual species unquantified',K.MOLECULAR);reject=g.product('column_reject','Source column/extraction physical rejects',K.WASTE)
    g.react('purify','Source miRNeasyMicro cleanup withon-column DNaseI',S.MAGIC,((released,R.TEMPLATE),(column,R.SELECTION),(dnase,R.CATALYST)),(purified,reject),reason='Source names actual combined column/DNase operation; detailed within-column molecular products not reconstructed.');g.discard(reject,'Source column rejects')
    rrna=inp('rrna_kit','Source Ribo-ZeroGold orRiboMinusv2 depletion supplies; per-sample selection in unreadTableS1',K.SELECTION)
    depleted=g.product('depleted','Source rRNA-depletedRNA; chosen kit/sample mapping unassigned',K.MOLECULAR);other=g.product('rrna_other','Other depletion-associated material; exact composition/finalfate unassigned',K.MOLECULAR)
    g.react('rrna_deplete','Source rRNA-depletion handoff',S.MAGIC,((purified,R.TEMPLATE),(rrna,R.SELECTION)),(depleted,other),reason='Two source named kit alternatives; no modern molecular recipe or unverified sample assignment.');g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Depletion-side material not carried into selected trace'))
    pap=inp('pap','Actual source yeast poly(A)polymerase forpoly(I)route',K.CATALYST);itp=inp('itp','Actual source inosine triphosphate substrate; dose unassigned')
    tc=inp('tail_conditions','Source37C30min; historical manufacturer-internal recipe not reconstructed',K.CONDITION)
    tailed=g.product('tailed','NativeRNA withartificialpoly(I)tail plusreactioncomponents; length/occupancy unassigned',K.MOLECULAR,chemical_state=NativeRNA().tailed('I'))
    g.react('tail','Source poly(I)tailing',S.TAILING,((depleted,R.TEMPLATE),(pap,R.CATALYST),(itp,R.STRUCTURAL),(tc,R.CONDITION)),(tailed,))
    chosen=g.product('ligation_input','Actual source500ngRNA contribution foradapterligation; amount relative totailingstock unassigned',K.MOLECULAR)
    rest=g.product('rna_stock','Untransferred RNA stock class; amount/presence/finalfate unassigned',K.MOLECULAR)
    g.react('ligation_allocate','Select source500ngRNA forligation without assumingwhole-stock transfer',S.SPLIT,((tailed,R.TEMPLATE),),(chosen,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'RNA stock outside selectedligationcontribution'))
    oligos=inp('custom_adapter_stock','Actual source ONT_oligoA/ONT_oligoB_C10 supplies; wholeorderings in unreadTableS4; stockcount/packaging unassigned',K.ADAPTER)
    ab=inp('adapter_buffer','Source adapter1.4uM,10mMTris7.5/50mMNaCl;95Cthen slowcool toroomtemperature',K.CONDITION)
    oa=g.product('oligo_a_component','Source ONT_oligoA component ofaggregate actualsupplies; wholeordering unassembled',K.OLIGO)
    ob=g.product('oligo_b_component','Source ONT_oligoB_C10 component ofaggregate actualsupplies; wholeordering unassembled',K.OLIGO)
    g.react('adapter_inventory','Named oligo components; stockcount and packagingunassigned, notextra freshaddition',S.GENERIC,((oligos,R.TEMPLATE),),(oa,ob),changes_molecule=False)
    adapter=g.product('custom_adapter','Source annealed customC10splinted adapter preparation; full duplex/end chemistry unassembled',K.ADAPTER)
    g.react('adapter_prepare','Source custom adapter annealing',S.HYBRIDIZATION,((oa,R.TEMPLATE),(ob,R.STRUCTURAL),(ab,R.CONDITION)),(adapter,))
    ligase=inp('ligase','Actual source T4DNA ligaseM0202M; exactligationstockvolume unassigned',K.CATALYST)
    ligsup=inp('ligation_supplies','Actual source ligation molecular supplies; internal recipe and stockpackaging unassigned')
    lc=inp('ligation_conditions','Source15min roomtemperature',K.CONDITION)
    sc=Scene();sc.strand('rna',[Segment('nativeRNA','[unassemblednativeRNA]',placeholder=True),Segment('polyI','[sourcepolyI;lengthunassigned]',placeholder=True)],label='Source substrate class; exact terminal/junction state unassembled')
    sc.strand('adapter',[Segment('adapter5context','[5primecontextunassembled]',placeholder=True),Segment('C10recognition','C'*10),Segment('adapter3context','[3primecontextunassembled]',placeholder=True)],label='SourceC10 recognition, not a reconstructed ordered adapter or assigned polarity')
    geometry=AdapterLigationGeometry('unavailable','unavailable',(sc,),'PolyI/C10 recognition is source-positive; actual TableS4 fulloligos, ends and occupied junction remain unassembled, not asserted absent.')
    joined=g.product('adapter_rna','Source adapter-ligated nativeRNA mixture; exact junction/unligatedfractions unassigned',K.MOLECULAR,chemical_state=NativeRNA().tailed('I').ligated())
    g.react('ligate','Source C10adapter ligation topolyI RNA',S.LIGATION,((chosen,R.TEMPLATE),(adapter,R.STRUCTURAL),(ligase,R.CATALYST),(ligsup,R.STRUCTURAL),(lc,R.CONDITION)),(joined,),geometry=geometry)
    ssiii=inp('ssiii','Actual source SuperScriptIII supportive RT',K.CATALYST)
    rtsup=inp('rt_supplies','Actual source-recommended RT molecular supplies; priming relationship within adapter preparation, individualidentities/doses and stockpackaging unassembled')
    supported=g.product('supported_rna','NativeRNA withsupportivecDNA andotherretainedreaction species; RNAis sequenced, cDNAis not',K.MOLECULAR,chemical_state=NativeRNA().tailed('I').ligated().reverse_transcribed())
    g.react('supportive_rt','Source supportive reverse transcription',S.MAGIC,((joined,R.TEMPLATE),(ssiii,R.CATALYST),(rtsup,R.STRUCTURAL)),(supported,),reason='Source explicitly reports supportive cDNA toreduceRNAsecondarystructure. Exact retainedprimer/nucleotide relation is unassembled; no new free primer addition or wholeDNAproduct asserted.')
    x=clean('support_clean',supported,inp('first_xp','Actual first AgencourtRNACleanXP purification',K.SELECTION),'Source cleanup aftersupportiveRT')
    loaded=inp('loaded_motor_adapter','Actual source sequencing-adapter preparation withpreloadedmotorprotein; RNA001/002 sampleassignment unassigned',K.ADAPTER)
    ad=g.product('motor_adapter_component','Adapter component ofsameactualpreloadedstock; wholeordering unassembled',K.ADAPTER)
    motor=g.product('motor_component','Copresent preloadedmotorprotein component ofsameactualstock; identity/stoichiometry unassembled',K.MOLECULAR)
    g.react('motor_inventory','Inventory components ofONEsourcepreloadedmotor-adapterstock, notfreshseparateproteinaddition',S.GENERIC,((loaded,R.TEMPLATE),),(ad,motor),changes_molecule=False)
    ms=inp('motor_ligation_supplies','Actual source motor-adapter ligation supplies; freshenzymeidentity/dose/packaging not assigned')
    prepared=g.product('motor_rna','Source motor-adapter-ligated nativeRNA preparation withsupportivecDNA; wholeoccupiedconstruct unassembled',K.MOLECULAR,chemical_state={'sequenced_polymer':'RNA','supporting_cDNA':'copresent,notsequenced','motor':'sourcepreloadedproteinretained;wholejunctionunassembled'})
    g.react('motor_ligate','Source attachment ofsequencingadapter withpreloadedmotor',S.MAGIC,((x,R.TEMPLATE),(ad,R.STRUCTURAL),(motor,R.STRUCTURAL),(ms,R.STRUCTURAL)),(prepared,),reason='Main176 explicitly ligates sequencingadapter withpreloadedmotor. Fulladapter/reagent/occupiedRNA geometry remains unassembled; no unreported freshligaseidentity.')
    x=clean('motor_clean',prepared,inp('second_xp','Actual second AgencourtRNACleanXP purification',K.SELECTION),'Source secondRNAcleanup')
    eb=inp('elution_buffer','Actual source ElutionBuffer contribution; exactvolume/composition unassigned')
    eluted=g.product('eluted','ElutionBuffer-resuspended sourceRNAsequencing preparation',K.MOLECULAR)
    g.react('elute','Source resuspension inElutionBuffer',S.GENERIC,((x,R.TEMPLATE),(eb,R.STRUCTURAL)),(eluted,),changes_molecule=False)
    rb=inp('running_buffer','Actual source RNA RunningBuffer; exactvolume/composition unassigned')
    library=g.product('rna_library','Source prepared motor-adapter nativeRNA library withsupportivecDNA; occupiedpolymers unassembled',K.LIBRARY,chemical_state={'polymer':'RNA','supporting_cDNA':'copresent,notsequenced','source':'AM176motoradapter/cleanup/EB+runningbuffer'})
    g.react('mix_running','Source RNA RunningBuffer mixing',S.GENERIC,((eluted,R.TEMPLATE),(rb,R.STRUCTURAL)),(library,),changes_molecule=False)
    load=g.product('loaded_rna','Actual source loading contribution ofmotor-preparedRNA library',K.LIBRARY,chemical_state=library.chemical_state)
    rem=g.product('library_stock','Untransferred sequencing-preparationstock; amount/presence/finalfate unassigned',K.LIBRARY,chemical_state=library.chemical_state)
    g.react('loading_allocate','Select actual flowcellloadingcontribution',S.SPLIT,((library,R.TEMPLATE),),(load,rem));g.mark_terminal(Terminal(rem.id,Disposition.UNRESOLVED,'Stock outside selectedload'))
    flow=inp('flowcell','Actual primedFLO-MIN106 flowcell; primingvolume/kitrecipe unassigned',K.CONDITION)
    # Flowcell is an acquisition condition, not an additional biochemical molecule or new copy.
    g.react('load_flowcell','Apply sourceRNAloadingcontribution toprimedFLO-MIN106',S.GENERIC,((load,R.TEMPLATE),(flow,R.CONDITION)),(g.product('sequencing_rna','Source nativeRNA contribution onflowcell',K.LIBRARY,chemical_state=load.chemical_state),),changes_molecule=False)
    g.sequence(g.materials['sequencing_rna'],SequencingPlan(acquisition='unspecified',limitation='Source48hMINKNOW directRNA defaults/FLO-MIN106; physicalRNAsequenced, supportivecDNA notsequenced. Exactbasecycles/primerpositions and RNA001/002 sampleassignment unassembled.'))
    g.declare_reviewed_scope(SELECTED_SCOPE+('load_flowcell','sequencing_rna'),SOURCE)
    return g.require_publishable()

def render_page():
    from page import head,info
    from protocol_views import assemble_page,WorkflowSection
    g=protocol_graph()
    return assemble_page(g,(head('nano-COP'),'<main class="wrap" data-protocol-complete="true"><h1>nano-COP</h1><p>4sU-selected chromatin RNA receives a source poly(I) tail and custom C10 adapter, supportive cDNA, and a preloaded-motor sequencing adapter. The prepared native RNA is sequenced directly.</p><p><a href="https://doi.org/10.1016/j.molcel.2019.11.017">Drexler et al. (2020)</a> · <a href="01_nano_cop.html">Research notes</a></p>',info('Accepted-author-manuscript preparation, Methods160/176; final published version equivalence is unverified. This selects the poly(I) branch without assigning a particular sample or choosing RNA001 versus RNA002 from unread TableS1. Full custom TableS4 adapter orderings, exact ligation-ready ends, occupied RNA/cDNA/motor geometry and historical kit details remain unassembled. C10 is a reported recognition segment, not a complete adapter or invented end.'),'<p>Supportive cDNA remains copresent and is not the sequenced substrate. No native RNA base is rewritten as unmodified U. The poly(A) and no-tailing alternatives remain separately qualified unresolved prefixes. No current ONT kit or Illumina comparison-library recipe is substituted.</p>','<h2>Read setup</h2><p>FLO-MIN106; 48 hours using MINKNOW direct-RNA defaults. Exact base-cycle, index and physical primer declarations are unassigned.</p>','<details><summary>Source preparation, retained stocks and physical fractions</summary>',WorkflowSection('Selected poly(I) native-RNA preparation'),'</details></main>'))

_prefix_selftest=selftest
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    _prefix_selftest()
    c=Check();g=protocol_graph()
    c('Source preparation produces an RNA library, not sequenced supportingDNA',tuple(g.materials[k].chemical_state['polymer'] for k in g.sequencing),('RNA',))
    c('Cleanup follows supportiveRT',g.steps['support_clean'].inputs[0].material,'supported_rna')
    c('First source ligation consumes selected500ngcontribution',g.steps['ligate'].inputs[0].material,'ligation_input')
    c('Motor adapter andprotein derivefromthe samepreloadedstock',g.steps['motor_inventory'].outputs,('motor_adapter_component','motor_component'))
    c('Second cleanup follows motoradapterattachment',g.steps['motor_clean'].inputs[0].material,'motor_rna')
    c('Library stock retained outside acquisition',g.terminals['library_stock'].disposition,Disposition.UNRESOLVED)
    c('Normal publicationgate',g.publication_problems(),[])
    source=Source(sources_root()/'nano-cop__10.1016+j.molcel.2019.11.017'/'j.molcel.2019.11.017_PMC7060811.xml.txt','AssociatedAM176')
    if have(c,source):
        text=' '.join(source.path.read_text().split())
        for label,needle in [('Actualmotoradapter','sequencing adapter with preloaded motor protein'),('PhysicalRNAreadout','The resulting cDNA does not get sequenced'),('Acquisition48h','sequenced using MINKNOW software for 48 hours')]:c(label,needle in text,True)
    c.report()
if __name__=='__main__':selftest()

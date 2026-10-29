"""Source-local bounded partial model; no inherited kit recipe or finished library claim."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from dataclasses import dataclass,replace
from chemdraw import Construct,Scene,Segment,feature,revcomp,panel,workflow_from_protocol_graph,workflow_panel
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem
from pcr import PCRGeometry
from page import head,info
PAGE_FILE='chem.html'

def unknown(label):return Construct([Segment(label,'[unassembled]',placeholder=True)])
def cleanup(g,key,x,label,procedure,state=None):
    sel=g.input(key+'_selection',procedure,K.SELECTION)
    kept=g.product(key+'_kept',label,K.MOLECULAR,state=state)
    waste=g.product(key+'_waste','Physically separated '+key+' rejects',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(sel,R.SELECTION)),(kept,waste))
    g.discard(waste,'Source-described physical separation; purity and quantitative recovery are not assumed')
    return kept

def aliquot(g,key,x,label):
    chosen=g.product(key+'_chosen',label,K.MOLECULAR)
    rest=g.product(key+'_remainder','Untransferred '+key+' material',K.MOLECULAR)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(chosen,rest))
    g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Source aliquot remainder fate is not specified in the inspected section; no disposal inferred'))
    return chosen

def limit(g,x):
    g.mark_terminal(Terminal(x.id,Disposition.UNRESOLVED,LIMIT))
    g.declare_reviewed_scope(REVIEWED,SOURCE)
    g.add_source_item(SourceItem('unassembled_source','Remaining source library/input chemistry and full strand/read assembly',LIMIT))
    return g

def render_page():
    g=protocol_graph()
    return '\n'.join([head(TITLE),'<main class="wrap"><h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="'+NOTE+'">Research notes</a></p>',
        info('Private incomplete source model. '+LIMIT),detail(),
        ''.join(workflow_panel(workflow_from_protocol_graph(g,x)) for x in ENDPOINTS),'</main>'])

def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Unassembled source and endpoint remain unpublished',g.publishable,False)
    c('Private preview has no completion marker','data-protocol-complete=' in render_page(),False)
    c('Unreviewed geometry does not assert source absence','[adapter overhang not reported]' in render_page(),False)
    p=sources_root()/Path(__file__).resolve().parents[1].name/SOURCE_FILE
    if p.exists():
        source=''.join(p.read_text().split()).upper()
        for seq in TRANSCRIBED:c('Original transcription '+seq,seq.upper() in source,True)
    else:c.skip('Source transcription','External original unavailable')
    source_checks(c,g)
    c.report()
from adapter_ligation import AdapterLigationGeometry
def unreviewed_ligation_geometry(reason):
    sc=Scene()
    sc.strand('unassembled end',[Segment('source geometry not assembled','[unreviewed molecular geometry]',placeholder=True)],label='Local reading/assembly boundary')
    sc.note('unassembled end',reason+' This is not a claim that the original source fails to report an end.')
    return AdapterLigationGeometry('unavailable','unavailable',(sc,),reason)

TITLE='RAID — RNA reporter conjugation and fixed-cell first-strand prefix'
NOTE='01_raid.html'
SOURCE='Published PMC6365588 DOI10.1038/s41598-018-37977-7 Methods62–90; published SI FigS2 and DetailedProtocol2018-05-15 pp8–10'
LIMIT='Source RNA reporter and reversible fixation/conjugation are local. Full pCp/linker atom geometry, historical CEL-seq2 full templates, downstream second-strand/IVT/PCR/read landing remain unassembled. SI prints N8+cellbarcode8 but main reports Read2=14; full primer/run compatibility is unresolved, not silently changed to a usual CEL-seq2 primer.'
SOURCE_FILE='s41598-018-37977-7_PMC6365588.xml.txt'
COMPENSATION='GGGAGACCGACGAAACTGTTAACG'
ARC_MARKER='ATCAGTCAACAG'
TRANSCRIBED=(COMPENSATION,ARC_MARKER)
ENDPOINTS=('firststrand_clean_kept','arc_mixture')
@dataclass(frozen=True)
class Reporter:
    label3:str='OH'
    antibody_attached:bool=False
    capped5:bool=True
    rna_backbone_retained:bool=True
    def __post_init__(self):
        if self.label3 not in {'OH','source pCp-Azide','source pCp-linker-TCO','source antibody-linker conjugate','source reduced RNA-side linker'}:raise ValueError('Unknown source reporter stage')
        if not self.capped5 or not self.rna_backbone_retained:raise ValueError('Selected tag coupling/reduction preserves source RNA/cap, not a DNA or deleted-base substitute')
        if self.antibody_attached and self.label3!='source antibody-linker conjugate':raise ValueError('Antibody attachment requires source RNA-side covalent coupling')
    def pcp(self):
        if self.label3!='OH':raise ValueError('Source starts with reporter3primeOH')
        return replace(self,label3='source pCp-Azide')
    def tco(self):
        if self.label3!='source pCp-Azide':raise ValueError('Azide label precedes DBCO-PEG12-TCO')
        return replace(self,label3='source pCp-linker-TCO')
    def conjugate(self):
        if self.label3!='source pCp-linker-TCO':raise ValueError('TCO reporter precedes tetrazine-antibody coupling')
        return replace(self,label3='source antibody-linker conjugate',antibody_attached=True)
    def reduce(self):
        if not self.antibody_attached:raise ValueError('Selected antibody linker is already free/unassigned')
        return replace(self,label3='source reduced RNA-side linker',antibody_attached=False)
    def scene(self):
        sc=Scene();sc.strand('reporter',[Segment('unassembled5prime regions','[reporter5prime regions]',placeholder=True),
            Segment('reporter-template UMI','N'*15,'umi',placeholder=True,feature=feature('raid.arc_umi','umi','random',note='Source reporter-template identifier; IVT copies need not have distinct sequences')),
            Segment('antibody identity','N'*10,'cbc',placeholder=True,feature=feature('raid.antibody','feature_barcode','unknown')),
            Segment('identification marker',ARC_MARKER.replace('T','U')),Segment('unassembled extension','[extension]',placeholder=True),Segment('polyA tail','A'*28)],mod5='source cap',mod3=self.label3,label='RNA reporter; pCp/linker chemistry is a separate typed annotation, not an extra DNA sequence guess')
        return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    reporter=g.input('reporter','Source IVT-produced capped polyA28 RNA reporter; DNA-template/IVT details upstream',K.MOLECULAR,primary=True,state=Reporter().scene(),chemical_state=Reporter())
    pcp=g.input('pcp','Source pCp-Azide; detailed cytidyl/linker structure unassembled',K.MOLECULAR)
    ligase=g.input('rna_ligase','T4RNA ligase withATP/DMSO;16°Covernight',K.CATALYST)
    azide=g.product('azide','RNA bearing source3prime pCp-Azide label plus free labeling reagents',K.MOLECULAR,chemical_state=Reporter().pcp())
    g.react('label_azide','Source3prime RNA end labeling, not conversion of RNA to DNA',S.LIGATION,((reporter,R.TEMPLATE),(pcp,R.STRUCTURAL),(ligase,R.CATALYST)),(azide,),geometry=unreviewed_ligation_geometry('Source pCp-Azide labeling is explicit; complete donor cytidyl modification and phosphodiester/linker geometry have not been assembled.'))
    clean=cleanup(g,'pcp_cleanup',azide,'Column-recovered RNA-Azide','MicroprepRNA column removes excess pCp-Azide')
    linker=g.input('dbco_tco','20fold molar DBCO-PEG12-TCO;roomtemperatureovernight',K.MOLECULAR)
    tco_condition=g.input('tco_condition','Roomtemperature overnight per source',K.CONDITION)
    rna_tco=g.product('rna_tco','RNA-TCO plus excess linker molecules',K.MOLECULAR,chemical_state=Reporter().pcp().tco())
    g.react('install_tco','Source azide/DBCO coupling installs TCO-bearing linker',S.CHEMICAL_CONVERSION,((clean,R.TEMPLATE),(linker,R.STRUCTURAL),(tco_condition,R.CONDITION)),(rna_tco,))
    clean=cleanup(g,'tco_cleanup',rna_tco,'Column-recovered RNA-TCO','MicroprepRNA purification;borate-buffered salinepH8.4')
    antibody=g.input('antibody_tz','Source antibody previously functionalized with NHS-SS-Tetrazine; upstream activation unassembled',K.MOLECULAR)
    conjugation_condition=g.input('conjugation_condition','Roomtemperature1h in source boratebufferpH8.4',K.CONDITION)
    arc=g.product('arc_mixture','Antibody-RNA conjugate class plus unreacted RNA/antibody fractions',K.MOLECULAR,chemical_state=Reporter().pcp().tco().conjugate())
    g.react('arc_coupling','TCO/tetrazine conjugation1h roomtemperature',S.CHEMICAL_CONVERSION,((clean,R.TEMPLATE),(antibody,R.STRUCTURAL),(conjugation_condition,R.CONDITION)),(arc,))
    cells=g.input('cells','Keratinocyte suspension for fixed-cell source branch',K.BIOLOGICAL,primary=True)
    fix=g.input('fixatives','2.5mMDSP plus2.5mMSPDP;45min',K.CONDITION)
    fixed=g.product('fixed','Reversibly crosslinked cells with retained cellular material',K.MOLECULAR)
    g.react('fix','Source reducible cell fixation',S.CROSS_LINKING,((cells,R.TEMPLATE),(fix,R.CONDITION)),(fixed,))
    q=g.input('quench_buffer','100mMTris-HClpH7.5/150mMNaCl',K.CONDITION)
    quenched=g.product('quenched','Quenched fixed cells plus reagent species',K.MOLECULAR)
    g.react('quench','Quench fixation without assuming physical disposal',S.QUENCHING,((fixed,R.TEMPLATE),(q,R.CONDITION)),(quenched,))
    permeab=g.input('permeabilization','0.1%TritonX100/PFBB/tRNA/RNaseinhibitor source blocking buffer',K.CONDITION)
    accessible=g.product('permeabilized','Blocked/permeabilized fixed cells',K.MOLECULAR)
    g.react('permeabilize','Permit intracellular antibody access while retaining fixed-cell RNA',S.PERMEABILIZATION,((quenched,R.TEMPLATE),(permeab,R.CONDITION)),(accessible,))
    stained=g.product('stained','ARC-stained cells plus free ARC/RNA/reagent fractions',K.MOLECULAR,chemical_state=Reporter().pcp().tco().conjugate())
    g.react('stain','Overnight source ARC immunostaining',S.BINDING,((accessible,R.TEMPLATE),(arc,R.STRUCTURAL)),(stained,),changes_molecule=False)
    washed=cleanup(g,'cell_washes',stained,'Washed stained cells','Six gentle10ml washes;source0.1×PFBB/PBS')
    sort=g.input('sorter','Source FACS single-cell collection into primer/oil-preloaded384well plates',K.SELECTION)
    well=g.product('well','One source collected-cell well; mineral oil remains until aqueous transfer',K.MOLECULAR)
    other=g.product('other_cells','Other source cells/collection fractions',K.MOLECULAR)
    g.react('sort','Physically collect cells into separate barcode wells',S.SELECTION,((washed,R.TEMPLATE),(sort,R.SELECTION)),(well,other));g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Other wells/nonselected material not represented as discarded or a covalent pool'))
    dtt=g.input('dtt','Source DTT reverse-crosslinking additive;25°C45min then65°C5min',K.CONDITION)
    released=g.product('released','Released endogenousRNA and reporterRNA with reduced RNA-side linker; reduced protein/fixative material remains',K.MOLECULAR,chemical_state={'reporter':Reporter().pcp().tco().conjugate().reduce(),'RNA_classes':'endogenous and synthetic RNA remain physically distinct'})
    g.react('release','Reduce source crosslinks and SS-containing ARC connection',S.CHEMICAL_CONVERSION,((well,R.TEMPLATE),(dtt,R.CONDITION)),(released,))
    primers=g.input('rt_primers','Original CEL-seq2-compatible cell-barcode primers; N8/CBC8/run-width compatibility unresolved',K.OLIGO)
    rt=g.input('rt','MaximaHminus source reverse transcriptase',K.CATALYST)
    cdna=g.product('firststrand','Gene/ARC first strands; original RNA and reaction species still coexist',K.MOLECULAR)
    g.react('first_rt','Source first-strand synthesis links source cell identity without fusing RNA species',S.REVERSE_TRANSCRIPTION,((released,R.TEMPLATE),(primers,R.STRUCTURAL),(rt,R.CATALYST)),(cdna,))
    comp=g.input('compensation','Printed barcode-compensation primer',K.OLIGO,state=Construct([Segment('compensation',COMPENSATION)]))
    prepared=g.product('prepared','First-strand reaction with compensation primer added/annealed; exact substrate geometry unassembled',K.MOLECULAR)
    g.react('compensation_step','Add compensation primer;85°C5min then steppedcooling',S.HYBRIDIZATION,((cdna,R.TEMPLATE),(comp,R.STRUCTURAL)),(prepared,))
    others=g.input('other_wells','Other source first-strand wells to pool; distinct well identities',K.MOLECULAR)
    pooled=g.product('pooled','Pooled first-strand preparations with mineral oil and all reaction species',K.MOLECULAR)
    g.react('pool','Pool source wells after first-strand identity introduction',S.POOL,((prepared,R.TEMPLATE),(others,R.TEMPLATE)),(pooled,))
    aqueous=cleanup(g,'oil_separation',pooled,'Recovered lower aqueous first-strand preparation','Remove upper mineral-oil phase;transfer lower aqueous phase')
    exo=g.input('exoi','ExonucleaseI;37°C20min',K.CATALYST)
    digested=g.product('exoi_products','Source ExoI-treated first strands plus DNA-digestion products; no removal by digestion alone',K.MOLECULAR)
    g.react('digest_primers','Source exonuclease treatment; substrate-specific protection unassembled',S.DIGESTION,((aqueous,R.TEMPLATE),(exo,R.CATALYST)),(digested,))
    final=cleanup(g,'firststrand_clean',digested,'Recovered first-strand preparation; residual RNA/reagents unquantified','Source AMPures/PEGNaCl bead-buffer purification and ethanol washes')
    return limit(g,final)

def detail():return info('Published SI printed primer has N8 plus highlighted8nt CBC; main reports14cycleRead2. Exact primer landing and field coverage are unresolved; source ARC15nt UMI is separate. No ordinary CEL-seq2 N6 substitution is made.')+panel(Reporter().scene().rows(),caption='Capped RNA reporter identifiers and polyA tail; exact end-label polymer unassembled')
def source_checks(c,g):
    p=Reporter().pcp().tco().conjugate();q=p.reduce()
    c('Reduction releases reporter without deleting RNA/cap',(p.antibody_attached,q.antibody_attached,q.rna_backbone_retained,q.capped5),(True,False,True,True))
    c("Reporter UMI and antibody barcode have different structured roles",tuple(s.feature.role for s in Reporter().scene().strands['reporter'].segs if s.feature),('umi','feature_barcode'))
    c('Transcribed reporter marker is RNA, not DNA/readout spelling',next(s.top for s in Reporter().scene().strands['reporter'].segs if s.name=='identification marker'),'AUCAGUCAACAG')
    c('ExoI products reach actual bead cleanup',g.steps['firststrand_clean'].inputs[0].material,'exoi_products')

REVIEWED=('reporter', 'pcp', 'rna_ligase', 'azide', 'pcp_cleanup_selection', 'pcp_cleanup_kept', 'pcp_cleanup_waste', 'dbco_tco', 'tco_condition', 'rna_tco', 'tco_cleanup_selection', 'tco_cleanup_kept', 'tco_cleanup_waste', 'antibody_tz', 'conjugation_condition', 'arc_mixture', 'cells', 'fixatives', 'fixed', 'quench_buffer', 'quenched', 'permeabilization', 'permeabilized', 'stained', 'cell_washes_selection', 'cell_washes_kept', 'cell_washes_waste', 'sorter', 'well', 'other_cells', 'dtt', 'released', 'rt_primers', 'rt', 'firststrand', 'compensation', 'prepared', 'other_wells', 'pooled', 'oil_separation_selection', 'oil_separation_kept', 'oil_separation_waste', 'exoi', 'exoi_products', 'firststrand_clean_selection', 'firststrand_clean_kept', 'firststrand_clean_waste', 'label_azide', 'pcp_cleanup', 'install_tco', 'tco_cleanup', 'arc_coupling', 'fix', 'quench', 'permeabilize', 'stain', 'cell_washes', 'sort', 'release', 'first_rt', 'compensation_step', 'pool', 'oil_separation', 'digest_primers', 'firststrand_clean')

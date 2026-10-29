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
    p=sources_root()/Path(__file__).resolve().parents[1].name/SOURCE_FILE
    if p.exists():
        source=''.join(p.read_text().split()).upper()
        for seq in TRANSCRIBED:c('Original transcription '+seq,seq.upper() in source,True)
    else:c.skip('Source transcription','External original unavailable')
    source_checks(c,g)
    c.report()
TITLE='Xdrop — TP53 detection, sorting and droplet MDA'
NOTE='01_xdrop.html'
SOURCE='Published DOI10.1002/humu.24063; canonical PMCCloudPDF, Methods2.1,2.5–2.9,Results3.2/3.4; original renderedp1672 verified'
LIMIT='Source-local TP53 branch stops at amplified enriched DNA. Full dMDA formulation, target landing and original SMRTbell1.0 preparation/adapter topology remain unassembled; PacBioRSII output is source-reported, not a finished model endpoint.'
SOURCE_FILE='humu.24063_PMC7496172.pdf.txt'
FORWARD='GGTGTGATGGGATGGATAAA'
REVERSE='CCCTGCATTTCTTTTGTTTG'
TRANSCRIBED=('GGTGTGAT','GGGATGGATAAA','CCCTGCATTTCTT','TTGTTTG')
ENDPOINTS=('second_release_kept',)
@dataclass(frozen=True)
class DetectionMixture:
    long_genomic_population:bool=True
    detector_products:bool=False
    sorted_positive:bool=False
    mda_rounds:int=0
    def __post_init__(self):
        if not self.long_genomic_population:raise ValueError('This selected source route retains genomic-template material, not detector amplicons alone')
        if self.sorted_positive and not self.detector_products:raise ValueError('Positive gate requires the source detector-PCR population')
        if self.mda_rounds not in (0,1,2):raise ValueError('Only the two source droplet-MDA rounds are represented')
    def detect(self):return replace(self,detector_products=True)
    def sort_positive(self):return replace(self,sorted_positive=True)
    def amplify(self):
        if not self.sorted_positive or self.mda_rounds==2:raise ValueError('Amplify only sorted material within the two-round source route')
        return replace(self,mda_rounds=self.mda_rounds+1)
    def scene(self):
        sc=Scene();sc.strand('genomic',[Segment('long genomic cohort','[variable genomic DNA]',placeholder=True)],label='Long genomic template/product cohort; no guaranteed intact length')
        if self.detector_products:sc.strand('detector',[Segment('short detector amplicons','[short TP53 amplicons]',placeholder=True)],label='Physically separate short detector-PCR products; not genomic extensions')
        sc.note('genomic','MDA can also amplify co-transferred material; no assertion of genomic-only output, purity, loss or specificity')
        return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    dna=g.input('dna','4ng Jurkat genomic DNA in source double-emulsion dPCR mix',K.MOLECULAR,primary=True,state=DetectionMixture().scene())
    f=g.input('forward','TP53 detection forward primer',K.OLIGO,state=Construct([Segment('TP53 forward',FORWARD)]))
    r=g.input('reverse','TP53 detection reverse primer',K.OLIGO,state=Construct([Segment('TP53 reverse',REVERSE)]))
    enzyme=g.input('detector_mix','Samplix2×dPCR mix;40cycles94°C3s/60°C30s after94°C2min',K.CATALYST)
    detected_state=DetectionMixture().detect()
    detected=g.product('detected','Double-emulsion PCR mixture retains long template and short detector products',K.MOLECULAR,state=detected_state.scene(),chemical_state=detected_state)
    g.react('detect','Target-detection PCR in double-emulsion droplets',S.PCR,((dna,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(enzyme,R.CATALYST)),(detected,),geometry=PCRGeometry.unavailable(unknown('TP53 genomic template'),reason='Source intron4 detector primers are known, but complete genomic landing/product geometry is unassembled. Detection amplicon is not the long library insert.'))
    stain=g.input('stain','Source DropletStain and dPCR buffer; TP53 intercalating-dye branch',K.MOLECULAR)
    stained=g.product('stained','Stained double-emulsion reaction populations',K.MOLECULAR,state=detected_state.scene())
    g.react('stain_drops','Add DNA-intercalating droplet stain afterPCR',S.BINDING,((detected,R.TEMPLATE),(stain,R.STRUCTURAL)),(stained,),changes_molecule=False)
    gate=g.input('sorter','Source S3e FACS gates; enrich positive double-emulsion droplets',K.SELECTION)
    selected_state=detected_state.sort_positive()
    selected=g.product('selected','Positive-gated droplets and all their carried DNA species',K.MOLECULAR,state=selected_state.scene())
    reject=g.product('sort_reject','Uncollected gating fractions',K.WASTE)
    g.react('sort','Physically collect gated fluorescent droplets',S.SELECTION,((stained,R.TEMPLATE),(gate,R.SELECTION)),(selected,reject));g.discard(reject,'Physical FACS collection boundary; positive gating does not prove all selected molecules contain target')
    recovered=cleanup(g,'first_release',selected,'Recovered aqueous DNA after droplet coalescence','20µl PicoBreak;mix/centrifuge;aqueous transfer',selected_state.scene())
    prior=aliquot(g,'first_transfer',recovered,'3–10µl aqueous template for first dMDA')
    for n in (1,2):
        mix=g.input('mda_mix_'+str(n),'Source Samplix dMDA reaction mixture;Phi29 stated in Results; full formulation unresolved',K.CATALYST)
        state=selected_state.amplify() if n==1 else selected_state.amplify().amplify()
        product=g.product('mda_'+str(n),'Droplet-MDA reaction '+str(n)+'; genomic and detector-derived populations unquantified',K.MOLECULAR,state=state.scene(),chemical_state=state)
        g.react('amplify_'+str(n),'Emulsify on source dMDA cartridge;30°C16h then65°C10min',S.EXTENSION,((prior,R.TEMPLATE),(mix,R.CATALYST)),(product,))
        prior=cleanup(g,'first_mda_release' if n==1 else 'second_release',product,'Coalesced dMDA aqueous product','20µl PicoBreak and source phase recovery',state.scene())
        if n==1:prior=aliquot(g,'second_transfer',prior,'About6µl first dMDA template for second round')
    return limit(g,prior)

def detail():return panel(DetectionMixture().detect().scene().rows(),caption='Long genomic molecules and detector-PCR products are separate species')
def source_checks(c,g):
    p=DetectionMixture().detect().sort_positive().amplify().amplify()
    c('Two MDA rounds do not erase the detector population',(p.mda_rounds,p.detector_products,p.long_genomic_population),(2,True,True))
    c('Positive sorting does not convert amplicon into genomic DNA',len(p.scene().strands),2)
    c('Phase recovery follows detection and physical sorting',g.steps['first_release'].inputs[0].material,'selected')

REVIEWED=('dna', 'forward', 'reverse', 'detector_mix', 'detected', 'stain', 'stained', 'sorter', 'selected', 'sort_reject', 'first_release_selection', 'first_release_kept', 'first_release_waste', 'first_transfer_chosen', 'first_transfer_remainder', 'mda_mix_1', 'mda_1', 'first_mda_release_selection', 'first_mda_release_kept', 'first_mda_release_waste', 'second_transfer_chosen', 'second_transfer_remainder', 'mda_mix_2', 'mda_2', 'second_release_selection', 'second_release_kept', 'second_release_waste', 'detect', 'stain_drops', 'sort', 'first_release', 'first_transfer', 'amplify_1', 'first_mda_release', 'second_transfer', 'amplify_2', 'second_release')

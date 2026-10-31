from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Scene, Segment, feature, workflow_from_protocol_graph, workflow_panel
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan
from page import head,info
PAGE_FILE='chem.html'

def operations(g):
    def op(key,x,label,reagent=None,state=None,kind=S.MAGIC,*,additional=(),output_kind=K.MOLECULAR):
        inputs=[(x,R.TEMPLATE)]
        if reagent:
            z=g.input(key+'_reagent',reagent,K.MOLECULAR);inputs.append((z,R.STRUCTURAL))
        inputs.extend(additional)
        kw={}
        if kind is S.MAGIC:kw['reason']='Actual source operation; full molecular inputs/products, geometry and recovery unassembled.'
        if kind is S.GENERIC:kw['changes_molecule']=False
        if kind is S.CROSS_LINKING:
            condition=g.input(key+'_conditions',label+'; source exposure, not an extra chemical dose',K.CONDITION);inputs.append((condition,R.CONDITION))
        y=g.product(key+'_product',label+'; preparation includes unremoved reaction components',output_kind,state=state)
        g.react(key,label,kind,inputs,(y,),**kw);return y
    def clean(key,x,label):
        z=g.input(key+'_selection',label,K.SELECTION)
        y=g.product(key+'_kept','Physically retained source fraction: '+label,K.MOLECULAR)
        w=g.product(key+'_removed','Physically removed fractions: '+label,K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(z,R.SELECTION)),(y,w));g.discard(w,'Source physical workup only; no perfect specificity or purity assertion.');return y
    return op,clean

def end(g,x,label):g.sequence(x,SequencingPlan(acquisition='unspecified',limitation=label))
def incomplete(g):
    # Literal reviewed source scope; unread Visium capture/library internals
    # and whole products are not covered by iterating over constructed nodes.
    g.declare_reviewed_scope(('section','antibodies','antibody_buffer_exchange_selection','antibody_buffer_exchange','fix_reagent','fix_conditions','fix','block_reagent','block','saponin_reagent','saponin','polyt_block','stain','stain_wash_selection','stain_wash','image','permeabilize_reagent','permeabilize','permeabilization_wash_selection','permeabilization_wash','additive','second_strand','qpcr_allocation','qpcr_reagent','qpcr','amplification_additive','amplify_reagent','amplify','size_reagent','size_separate','adt_cleanup1_selection','adt_cleanup1','adt_cleanup2_selection','adt_cleanup2','adt_amp_mix','adt_pcr','adt_final_selection','adt_final'),SOURCE)
    g.add_source_item(SourceItem('unassembled','Unreviewed linked originals and full molecular geometry',LIMIT));return g

from dataclasses import dataclass
import illumina as il
TITLE='SPOTS — spatial RNA and antibody-derived tags (>20-antibody buffer-exchange route)'
SOURCE='Published SI MOESM1 Library Preparation Workflow pp.14–21 and MOESM3 antibody table; accepted PMC10272089 Methods separately qualified'
LIMIT='Private source-limited operational model. Published SI is primary here; accepted manuscript differs in additive-primer volumes and is not merged. Historical Visium guide, whole capture/tag polymers and exact product/read geometry remain unassembled.'
@dataclass(frozen=True)
class PhysicalOligo:
    bases:str
    notation:str
    modifications:tuple
    terminal3:str=''
    def scene(self):
        sc=Scene();sc.strand('oligo',[Segment(self.notation,self.bases)],label=self.notation,mod3=self.terminal3)
        sc.note('oligo','Physical source oligo: '+self.notation+'; '+', '.join(self.modifications))
        sc.source_oligo=self;return sc
BLOCK=PhysicalOligo('T'*25,'T23*T*T*/3InvdT/',('terminal phosphorothioate linkages as printed','3-prime inverted dT'),'inverted dT (blocked)')
ADDITIVE=PhysicalOligo('CCTTGGCACCCGAGAATTCCA','CCTTGGCACCCGAGAATT*C*C*A',('three printed phosphorothioate linkages',))
def tag_inventory():
    sc=Scene();sc.strand('tag',[Segment('CD4 antibody identity example','AACAAGACCCTTGAG',feature=feature('spots.cd4','feature_barcode','fixed'))]);sc.note('tag','Published MOESM3 CD4 RM4-5 identity only; whole antibody-conjugated oligo and conjugation geometry unassembled.');return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE);op,clean=operations(g)
    tissue=g.input('section','10um frozen OCT tissue on Visium spatial capture slide',K.BIOLOGICAL,primary=True)
    antibodies=g.input('antibodies','TotalSeq-A antibody panel; source CD4 identity example',K.MOLECULAR,state=tag_inventory())
    antibodies=clean('antibody_buffer_exchange',antibodies,'Optional source >20-antibody pool 50kDa ultrafiltration/3xSSC wash; recover retentate')
    x=op('fix',tissue,'1% formaldehyde fixation10min;5 SSC dips','Actual source formaldehyde/PBS and SSC contributions',kind=S.CROSS_LINKING)
    x=op('block',x,'Source Fc/RVC blocking15min4C','Actual source blocking mixture; molecular constituents not assigned here',kind=S.BINDING)
    x=op('saponin',x,'Saponin pre-permeabilization10min4C','Source saponin mixture',kind=S.GENERIC)
    block=g.input('polyt_block','Published protected dT25 blocker',K.OLIGO,state=BLOCK.scene())
    stained=g.product('stained','Antibody-bound tissue plus unbound panel/blocker',K.MOLECULAR)
    g.react('stain','Antibody staining90min4C;polyT blocker limits ADT capture during staining',S.BINDING,((x,R.TEMPLATE),(antibodies,R.STRUCTURAL),(block,R.STRUCTURAL)),(stained,))
    x=clean('stain_wash',stained,'Four antibody washes and20 SSC dips;retain slide/tissue')
    x=op('image',x,'Fluorescence imaging;optional coverslip removed in SSC',kind=S.GENERIC)
    x=op('permeabilize',x,'Tissue-removal enzyme/SDS permeabilization;time optimized per tissue','Source enzyme/SDS/SSC')
    x=clean('permeabilization_wash',x,'Three0.1xSSC washes;retain captured slide-associated species')
    x=op('capture_rt',x,'Historical Visium capture and reverse transcription;exact guide not reviewed','Historical RT reagents')
    add=g.input('additive','Published phosphorothioate FB additive primer',K.OLIGO,state=ADDITIVE.scene())
    y=g.product('second_strand_product','Second-strand reaction including ADT and RNA-derived species',K.MOLECULAR)
    g.react('second_strand','Second-strand synthesis;publishedSI2uL100uM additive/77uL total',S.MAGIC,((x,R.TEMPLATE),(add,R.STRUCTURAL)),(y,),reason='Published additive dose established; historical Visium synthesis reagents and full strand geometry unreviewed.')
    from protocol_graph import Terminal,Disposition
    q=g.product('qpcr_aliquot','Source qPCR aliquot;volume unassembled',K.MOLECULAR)
    remaining=g.product('remaining_cdna','Remaining cDNA preparation for amplification',K.MOLECULAR)
    g.react('qpcr_allocation','Allocate source qPCR material;not all cDNA consumed in cycle determination',S.SPLIT,((y,R.TEMPLATE),),(q,remaining))
    measured=op('qpcr',q,'Source qPCR determines cDNA amplification cycles','Historical kit qPCR mixture',output_kind=K.QC)
    g.mark_terminal(Terminal(measured.id,Disposition.QC,'Analytical qPCR contribution; not the library substrate.'))
    x=remaining
    amp_add=g.input('amplification_additive','Actual fresh source1uL0.2uM FB additive for cDNA amplification, distinct from second-strand dose',K.OLIGO,state=ADDITIVE.scene())
    x=op('amplify',x,'cDNA amplification;publishedSI1uL0.2uM additive/102uL total','Actual source51uL AmpMix and15uL cDNA primers; molecular stock contributions, not a new premix',additional=((amp_add,R.STRUCTURAL),))
    sel=g.input('size_reagent','0.6xSPRI selection',K.SELECTION)
    gene=g.product('gene_beads','Bead fraction enriched for >300bp RNA-derived cDNA',K.MOLECULAR)
    adt=g.product('adt_supernatant','Retained supernatant enriched for ~180bp antibody products;contains soluble reagents',K.MOLECULAR)
    g.react('size_separate','0.6xSPRI physically separates two retained branches',S.CLEANUP,((x,R.TEMPLATE),(sel,R.SELECTION)),(gene,adt))
    gene=op('gene_library',gene,'Historical Visium cDNA cleanup and spatial GEX library construction;operational recipe unreviewed','Actual source-delegated Visium GEX library molecular supplies; exact identities and packaging unassembled')
    end(g,gene,'Source GEX NextSeq500/Novaseq6000;spatial identifier and UMI classes known,full geometry/run cycles unassembled.')
    adt=clean('adt_cleanup1',adt,'First1.9xSPRI;discard supernatant,one80%ethanol wash;resuspend beads50uL')
    adt=clean('adt_cleanup2',adt,'Second1.9xSPRI95uL;two80%ethanol washes;elute45uL water and remove beads')
    f=g.input('si_pcr','Published SI-PCR',K.OLIGO,state=Construct([Segment('SI-PCR',il.P5+'TCTTTCCCTACACGACGCTC')]))
    r=g.input('rpi','Published RPI_B1 example;no assignment to experiment',K.OLIGO,state=Construct([Segment('P7',il.P7),Segment('index','AAGATTAC',feature=feature('spots.sample','sample_index','fixed')),Segment('sourceRPI arm','GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA')]))
    y=g.product('adt_pcr_product','ADT PCR population;primer modifications are not copied as phosphorothioates',K.MOLECULAR)
    mix=g.input('adt_amp_mix','Actual source50uL AmpMix in100uL ADT PCR; internal catalyst/nucleotide inventory unassembled',K.MOLECULAR)
    g.react('adt_pcr','SI-PCR/RPI indexing PCR6–10cycles in publishedSI',S.MAGIC,((adt,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(mix,R.STRUCTURAL)),(y,),reason='Source operational PCR and actual amplification mix known; exact kit components, primer landing and whole products remain unassembled.')
    adt=clean('adt_final',y,'1xSPRI;two80%ethanol washes;30uL water elution and remove beads')
    end(g,adt,'Source ADT NextSeq500/Novaseq6000;spatial address/UMI/antibody barcode not equated;full read-primer geometry unassembled.')
    return incomplete(g)
def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Private partial graph',g.publishable,False)
    c('ADT cleanup uses retained supernatant',g.steps['adt_cleanup1'].inputs[0].material,'adt_supernatant')
    c('Second cleanup follows first physical recovery',g.steps['adt_cleanup2'].inputs[0].material,'adt_cleanup1_kept')
    c('Blocker inversion retained',g.materials['polyt_block'].state.source_oligo.modifications[-1],'3-prime inverted dT')
    c('Actual blocker strand carries terminal block',g.materials['polyt_block'].state.strands['oligo'].mod3,'inverted dT (blocked)')
    c('Additive physical modifications retained',g.materials['additive'].state.source_oligo.modifications,ADDITIVE.modifications)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'s41587-022-01536-3_41587_2022_1536_MOESM1_ESM.pdf.txt'
    if p.exists():
        t=p.read_text();c('Original retained supernatant instruction','DO NOT DISCARD!' in t,True);c('Modified source oligos','TTTTTTTTTTTTTTTTTTTTTTT*T*T*/3InvdT/' in t and ADDITIVE.notation in t,True)
    else:c.skip('Source transcription','External source absent')
    c.report()

def render_page():
    g=protocol_graph()
    return head(TITLE)+'<main class="wrap"><h1>'+TITLE+'</h1>'+info(LIMIT)+'<p><a href="01_spots.html">Research notes</a></p>'+''.join(workflow_panel(workflow_from_protocol_graph(g,x)) for x in ('gene_library_product','adt_final_kept'))+'</main>'

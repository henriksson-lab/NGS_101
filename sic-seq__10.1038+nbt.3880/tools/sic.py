"""SiC-seq gel retention and cell-barcode workflow; no diffusion/yield simulation."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
import illumina as il
import nextera as nx
from chemdraw import Construct, Scene, Segment, feature, protocol_graph_panels
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SequencingPlan
from pcr import PCRGeometry, PCRPrimer
from page import head, info
PAGE_FILE='sic.html'
SOURCE='https://doi.org/10.1038/nbt.3880; Results39–43 and Online Methods107–131'
BARCODE_PREFIX='GCAGCTGGCGTAATAGCGAGTACAATCTGCTCTGATGCCGCATAG'
BARCODE_TAIL='TAAGCCAGCCCCGACACT'
FL127=il.P5+nx.S5
FL128=nx.ME_RC+nx.S7_RC+'GTGTCGGGGCTGGCTTA'
FL129=il.P7+'CAGCTGGCGTAATAGCG'
CUSTOM_INDEX='GCCCACGAGACGTGTCGGGGCTGGCTTA'
def barcode_template():
    return Construct([Segment('barcode5prime flank',BARCODE_PREFIX),
        Segment('cell barcode','N'*15,'cbc',placeholder=True,
                feature=feature('sic.cell','cell_barcode','random')),
        Segment('barcode3prime flank',BARCODE_TAIL)],name='Source single barcode template')
def barcode_geometry():
    return PCRGeometry.exact(barcode_template(),PCRPrimer('FL129',FL129,anneal_nt=len(FL129)-len(il.P7)),
        PCRPrimer('FL128',FL128,anneal_nt=len(FL128)-len(nx.ME_RC)-len(nx.S7_RC)),
        statement='Source primer3prime constant components located on source barcode template; this checks landing, not PCR yield or product population.')

@dataclass(frozen=True)
class GelRetainedGenome:
    stage:str
    basis:str=SOURCE
    def __post_init__(self):
        if self.stage not in {'native','transposase-associated fragments'}:
            raise ValueError('Only the source-declared native/associated-fragment states are represented')
    def scene(self):
        sc=Scene.duplex([Segment('genomic material','X'*40,placeholder=True)],label=self.stage)
        sc.note('top','Qualitative retention in agarose gel is source-declared; length, loss and diffusion are not calculated')
        return sc
NATIVE=GelRetainedGenome('native')
TAGMENTED=GelRetainedGenome('transposase-associated fragments')

def cleanup(g,key,template,reagent,label,kind=K.MOLECULAR):
    r=g.input(key+'_selection',reagent,K.SELECTION)
    out=g.product(key+'_kept',label,kind);w=g.product(key+'_waste','Unretained '+key+' material',K.WASTE)
    g.react(key,key.replace('_',' '),S.CLEANUP,((template,R.TEMPLATE),(r,R.SELECTION)),(out,w))
    g.discard(w,'Removed at the stated physical separation, not before it');return out

def protocol_graph():
    g=ProtocolGraph('SiC-seq — synthetic microbial community branch',authority=SOURCE)
    microgels=g.input('cell_microgels','Source-prepared microbial cells in agarose microgels, transferred to aqueous buffer',K.BIOLOGICAL,primary=True)
    lytic=g.input('lytic_enzymes','Zymolyase,lysostaphin,mutanolysin and lysozyme in source buffer',K.CATALYST)
    lysed=g.product('lysed','Cell-wall-digested microgels and soluble lysate',K.MOLECULAR)
    g.react('cell_wall_lysis','Digest microbial cell walls overnight at37°C',S.DIGESTION,
        ((microgels,R.TEMPLATE),(lytic,R.CATALYST)),(lysed,))
    retained=cleanup(g,'first_lysate_separation',lysed,'Centrifuge microgels and remove supernatant','Microgels after first lysate separation')
    protease=g.input('protease','ProteinaseK with0.5% lithium dodecyl sulfate/EDTA',K.CATALYST)
    proteolysed=g.product('proteolysed','Microgel genomes with digested cellular proteins and soluble material',K.MOLECULAR)
    g.react('proteolysis','Proteolyse/solubilize at50°C30min',S.DEPROTEINIZATION,
        ((retained,R.TEMPLATE),(protease,R.CATALYST)),(proteolysed,))
    washed=cleanup(g,'microgel_washes',proteolysed,'Centrifugal Tween/ethanol/Tween wash series','Washed genome-containing gels')
    # Attach a local physical state only after source-declared cleanup; it does
    # not claim perfect retention or predict how much DNA remains.
    genome=g.product('genome','Gel-retained native genomic material',K.MOLECULAR,state=NATIVE.scene(),chemical_state=NATIVE)
    g.react('retained_state','Carry washed microgels into source tagmentation droplets',S.GENERIC,
        ((washed,R.TEMPLATE),),(genome,),changes_molecule=False,
        reason='Physical re-encapsulation is source-described; droplet occupancy and any unencapsulated fraction are unassigned, not assumed quantitative.')
    tn5=g.input('tn5','Nextera DNA Library Prep FC-121-1030 transposase',K.CATALYST)
    adapters=g.input('loaded_adapters','Source kit loaded-adapter pool; full loaded geometry/outcome classes unresolved',K.ADAPTER)
    fragments=g.product('fragments','Tagmented genomic fragments retained as associated macromolecular material',K.MOLECULAR,
        state=TAGMENTED.scene(),chemical_state=TAGMENTED)
    g.react('tagment','Re-encapsulated Nextera tagmentation:50°C1h',S.TAGMENTATION,
        ((genome,R.TEMPLATE),(adapters,R.STRUCTURAL),(tn5,R.CATALYST)),(fragments,))
    template=g.input('barcode_template','Single random15-base cell barcode template',K.OLIGO,state=barcode_template())
    f128=g.input('fl128','FL128 barcode PCR primer',K.OLIGO,state=Construct([Segment('FL128',FL128)]))
    f129=g.input('fl129','FL129 P7-containing primer',K.OLIGO,state=Construct([Segment('FL129',FL129)]))
    phusion=g.input('phusion','Phusion Hot Start Flex PCR mix in barcode droplets',K.CATALYST)
    clonal=g.product('clonal_barcodes','Clonal barcode products isolated within barcode droplets',K.MOLECULAR)
    g.react('barcode_pcr','Digital droplet barcodePCR:40cycles',S.PCR,
        ((template,R.TEMPLATE),(f129,R.STRUCTURAL),(f128,R.STRUCTURAL),(phusion,R.CATALYST)),(clonal,),geometry=barcode_geometry())
    f127=g.input('fl127','FL127 P5-containing primer',K.OLIGO,state=Construct([Segment('FL127',FL127)]))
    mix=g.input('barcoding_mix','Platinum Multiplex PCR+Bst2.0; NexteraXT FC-131-1024 NT buffer, Tween/PEG',K.CATALYST)
    joined=g.product('joined','Cell-barcoded genomic fragments in PCR droplets',K.MOLECULAR)
    g.react('genome_barcoding','Merge gel/barcode/PCR droplets;65°C5min then30-cycle barcodingPCR',S.PCR,
        ((fragments,R.TEMPLATE),(clonal,R.TEMPLATE),(f127,R.STRUCTURAL),(f129,R.STRUCTURAL),(mix,R.CATALYST)),(joined,),
        geometry=PCRGeometry.unavailable(TAGMENTED.scene(),reason='Source describes handle-mediated splicing; full genomic substrates, adapter-end combinations and final product geometry are unassigned in this capped reading.'))
    pipette=g.input('pipette','Micropipette separation of large coalesced droplets',K.SELECTION)
    retained_drops=g.product('retained_drops','Remaining source barcoding droplets',K.MOLECULAR)
    coalesced=g.product('coalesced','Large coalesced droplets physically removed',K.WASTE)
    g.react('drop_selection','Remove large coalesced droplets afterPCR',S.SELECTION,
        ((joined,R.TEMPLATE),(pipette,R.SELECTION)),(retained_drops,coalesced))
    g.discard(coalesced,'Physically removed by micropipette before breaking remaining emulsion')
    recovered=cleanup(g,'emulsion_recovery',retained_drops,'Perfluorooctanol/centrifugation; aqueous transfer; Zymo DNA Clean&Concentrator','Recovered barcoded DNA')
    library=cleanup(g,'size_selection',recovered,'AMPure selection for200–600bp','SiC-seq genomic sequencing library',K.LIBRARY)
    g.sequence(library,SequencingPlan(acquired_roles=('Read 1','Custom barcode index'),
        read_lengths={'Read 1':150},
        limitation='Synthetic-community MiSeq single-end150 and custom index primer specified. Source barcode width15 is a molecular feature, not an assumed index-cycle schedule; final strand/primer coordinates remain unassigned.'))
    g.declare_reviewed_scope(('cell_microgels','lytic_enzymes','lysed','cell_wall_lysis','first_lysate_separation_selection','first_lysate_separation_kept','first_lysate_separation_waste','first_lysate_separation',
        'protease','proteolysed','proteolysis','microgel_washes_selection','microgel_washes_kept','microgel_washes_waste','microgel_washes','genome','retained_state',
        'tn5','loaded_adapters','fragments','tagment','barcode_template','fl128','fl129','phusion','clonal_barcodes','barcode_pcr','fl127','barcoding_mix','joined','genome_barcoding',
        'pipette','retained_drops','coalesced','drop_selection','emulsion_recovery_selection','emulsion_recovery_kept','emulsion_recovery_waste','emulsion_recovery',
        'size_selection_selection','size_selection_kept','size_selection_waste','size_selection'),SOURCE)
    return g

def render_page():
    return '\n'.join([head('SiC-seq — bounded gel/barcoding model'),'<div class="wrap"><h1>SiC-seq</h1>',
        '<p class="research-notes"><a href="01_sic_seq.html">Research notes</a></p>',
        info('Gel retention is a source-declared qualitative state, not a loss/diffusion prediction. Complete genomic splicing/read geometry remains unassigned. The cell barcode is not a UMI or sample index.'),
        info('The publisher supplement carries a21June2017 correction to Figure11a taxonomic display; this molecular model does not use that taxonomy panel.'),
        protocol_graph_panels(protocol_graph(),cls='long'),'</div>'])

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('Barcode role is cell identity, not UMI',barcode_template().get('cell barcode').feature.role,'cell_barcode')
    c('Coalesced droplets removed before emulsion recovery',g.steps['emulsion_recovery'].inputs[0].material,'retained_drops')
    c('Native and tagmented physical retention remain different states',NATIVE.stage==TAGMENTED.stage,False)
    from storage import sources_root
    p=sources_root()/Path(__file__).resolve().parents[1].name/'nbt.3880_PMC5531050.xml.txt'
    if p.exists():
        text=''.join(p.read_text().split())
        for seq in (barcode_template().top(),FL127,FL128,FL129,CUSTOM_INDEX):c('Source oligo transcription '+seq,seq in text)
    else:c.skip('Source oligos','External defining article unavailable')
    c.report(exit_on_fail=True)

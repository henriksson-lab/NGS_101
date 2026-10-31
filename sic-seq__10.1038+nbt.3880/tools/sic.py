"""SiC-seq gel retention and cell-barcode workflow; no diffusion/yield simulation."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
import illumina as il
import nextera as nx
from chemdraw import Construct, Scene, Segment, feature, protocol_graph_panels
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SequencingPlan, Terminal, Disposition
from pcr import PCRGeometry, PCRPrimer
from page import head, info, table
from protocol_model import role
from protocol_views import WorkflowSection, PlanLimitationSection, OligoSection, OligoEntry, assemble_page
PAGE_FILE='sic.html'
SOURCE='Associated accepted manuscript PMC5531050, Results workflow/synthetic-community validation and Online Methods; publisher corrected SI identity/Fig1caption only. No whole-VOR equivalence or complete linked-source review is asserted.'
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
        sc=Scene.duplex([role('genomic material; length and positions unassigned')],label=self.stage,unpaired=("genomic material; length and positions unassigned'",))
        sc.note('top','Qualitative retention in agarose gel is source-declared; length, loss and diffusion are not calculated')
        return sc
NATIVE=GelRetainedGenome('native')
TAGMENTED=GelRetainedGenome('transposase-associated fragments')

def cleanup(g,key,template,reagent,label,kind=K.MOLECULAR):
    r=g.input(key+'_selection',reagent,K.SELECTION)
    out=g.product(key+'_kept',label,kind);w=g.product(key+'_waste','Unretained '+key+' material',K.WASTE)
    g.react(key,key.replace('_',' '),S.CLEANUP,((template,R.TEMPLATE),(r,R.SELECTION)),(out,w))
    g.discard(w,'Removed at the stated physical separation, not before it');return out

def _qc(g,key,stock,label):
    analytical=g.product(key+'_analytical',label,K.QC,state=stock.state)
    observation=key=='barcode_qc'
    continuing=g.product(key+'_continuing','Continuing stock; barcode sampling/stain exposure or library aliquot volume unassigned',stock.kind,state=stock.state)
    g.react(key,'Source analytical observation; physical barcode sampling unspecified'if observation else'Reserve source analytical aliquot; volume unassigned',S.GENERIC if observation else S.SPLIT,((stock,R.TEMPLATE),),(analytical,continuing),changes_molecule=False)
    g.mark_terminal(Terminal(analytical.id,Disposition.QC,'Analytical boundary only. Barcode microscopy does not establish physical sampling or whether continuing droplets contacted stain; assay chemistry remains unexpanded.'))
    return continuing


def _parts(g,key,stock):
    enzyme=g.product(key+'_enzyme','Catalytic component of the same actual prepared stock',K.CATALYST)
    nt=g.product(key+'_nt','Nucleotide-bearing component of that same actual stock',K.MOLECULAR)
    buffer=g.product(key+'_buffer','Nonincorporated component of that same actual stock',K.CONDITION)
    g.react(key+'_components','Copresent inventory of one source prepared stock; no separation or extra additions',S.GENERIC,((stock,R.TEMPLATE),),(enzyme,nt,buffer),changes_molecule=False)
    return enzyme,nt,buffer


def protocol_graph():
    g=ProtocolGraph('SiC-seq — synthetic microbial community branch',authority=SOURCE,representation='published_instance')
    cells=g.input('cells','Source artificial microbial community; individual genomes and occupancy unassigned',K.BIOLOGICAL,primary=True)
    agarose=g.input('agarose','Source molten agarose stream; exact preparation unassembled',K.MOLECULAR)
    encapsulation=g.input('encapsulation_conditions','Source equal-volume cell/agarose coflow into oil droplets, cooling then aqueous transfer; settings unexpanded',K.CONDITION)
    microgels=g.product('cell_microgels','Source cell-containing microgel collection transferred into aqueous carrier; no perfect occupancy claim',K.BIOLOGICAL)
    carrier=g.product('initial_carrier','Removed initial oil/carrier material; exact composition and losses unassigned',K.MOLECULAR)
    g.react('prepare_microgels','Source aggregate cell encapsulation, cooling and transfer to aqueous phase',S.GENERIC,((cells,R.TEMPLATE),(agarose,R.STRUCTURAL),(encapsulation,R.CONDITION)),(microgels,carrier),changes_molecule=False)
    g.mark_terminal(Terminal(carrier.id,Disposition.UNRESOLVED,'Source oil-to-aqueous transfer; carrier disposal and any cell/DNA losses unassigned.'))
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
    g.react('retained_state','Record native genome retention after source microgel washes',S.GENERIC,
        ((washed,R.TEMPLATE),),(genome,),changes_molecule=False,
        reason='Qualitative reference-state annotation only; actual source re-encapsulation with reagents is included in the following tagmentation boundary. Occupancy and losses remain unassigned.')
    tag_stock=g.input('tagmentation_stock','One actual source8µl Nextera enzyme preparation; molecular internals unassembled',K.MOLECULAR)
    tn5=g.product('tn5','Transposase component of that same source enzyme preparation',K.CATALYST)
    adapters=g.product('loaded_adapters','Tag-donating component of that same preparation; loaded end geometry and packaging unassembled',K.ADAPTER)
    g.react('tagmentation_components','Inventory of the one actual source tagmentation supply; no fresh adapter addition or physical fractionation',S.GENERIC,((tag_stock,R.TEMPLATE),),(tn5,adapters),changes_molecule=False)
    tag_buffer=g.input('tagmentation_buffer','Source tagmentation buffer/water conditions; actual carrier and droplet membership unassigned',K.CONDITION)
    fragments=g.product('fragments','Tagmented genomic fragments retained as associated macromolecular material',K.MOLECULAR,
        state=TAGMENTED.scene(),chemical_state=TAGMENTED)
    g.react('tagment','Source re-encapsulation with Nextera reagents then50°C1h tagmentation',S.TAGMENTATION,
        ((genome,R.TEMPLATE),(adapters,R.STRUCTURAL),(tn5,R.CATALYST),(tag_buffer,R.CONDITION)),(fragments,))
    template=g.input('barcode_template','Single random15-base cell barcode template',K.OLIGO,state=barcode_template())
    f128=g.input('fl128','FL128 barcode PCR primer',K.OLIGO,state=Construct([Segment('FL128',FL128)]))
    f129=g.input('fl129','FL129 P7-containing primer',K.OLIGO,state=Construct([Segment('FL129',FL129)]))
    phusion=g.input('phusion','Actual Phusion Hot Start Flex Master Mix M0536L stock',K.MOLECULAR)
    phusion_enzyme,phusion_nt,phusion_buffer=_parts(g,'phusion',phusion)
    factory_conditions=g.input('barcode_droplet_conditions','Source limiting-dilution droplet formation; HFE7500 replaced by FC40 before cycling; Tween/PEG conditions',K.CONDITION)
    factory=g.product('barcode_droplets','Actual barcode-droplet collection with copresent PCR supplies; occupancy unassigned',K.MOLECULAR)
    factory_carrier=g.product('barcode_removed_carrier','Source pipetted-out HFE carrier; exact carryover unassigned',K.MOLECULAR)
    g.react('barcode_encapsulation','Source barcode mixture encapsulation and carrier exchange before PCR',S.GENERIC,((template,R.TEMPLATE),(f128,R.STRUCTURAL),(f129,R.STRUCTURAL),(phusion_enzyme,R.CATALYST),(phusion_nt,R.STRUCTURAL),(phusion_buffer,R.CONDITION),(factory_conditions,R.CONDITION)),(factory,factory_carrier),changes_molecule=False)
    g.mark_terminal(Terminal(factory_carrier.id,Disposition.UNRESOLVED,'Source carrier removed/replaced; no DNA-specific loss or downstream carrier use assigned.'))
    ft=g.product('factory_template','Same-droplet barcode template',K.OLIGO,state=template.state)
    ff=g.product('factory_fl129','Same-droplet FL129 contribution',K.OLIGO,state=f129.state)
    fr=g.product('factory_fl128','Same-droplet FL128 contribution',K.OLIGO,state=f128.state)
    fe=g.product('factory_enzyme','Same-droplet Phusion component',K.CATALYST)
    fn=g.product('factory_nt','Same-droplet nucleotide component',K.MOLECULAR)
    fb=g.product('factory_buffer','Same-droplet nonincorporated conditions',K.CONDITION)
    g.react('factory_components','Inventory of copresent barcode-droplet components; no physical fractionation or fresh doses',S.GENERIC,((factory,R.TEMPLATE),),(ft,ff,fr,fe,fn,fb),changes_molecule=False)
    clonal=g.product('clonal_barcodes','Source barcode-droplet product collection; occupancy and clonality not predicted',K.MOLECULAR)
    g.react('barcode_pcr','Digital droplet barcodePCR:40cycles',S.PCR,
        ((ft,R.TEMPLATE),(ff,R.STRUCTURAL),(fr,R.STRUCTURAL),(fe,R.CATALYST),(fn,R.STRUCTURAL),(fb,R.CONDITION)),(clonal,),geometry=barcode_geometry())
    clonal=_qc(g,'barcode_qc',clonal,'Source fluorescent encapsulation-rate observation reference; sampling and production stain exposure unassigned')
    f127=g.input('fl127','FL127 P5-containing primer',K.OLIGO,state=Construct([Segment('FL127',FL127)]))
    f129_join=g.input('fl129_join','Actual fresh FL129 contribution in the later barcoding mix; not factory primer stock reused',K.OLIGO,state=f129.state)
    mix=g.input('barcoding_mix','Actual Platinum Multiplex PCR Master Mix4464268 stock',K.MOLECULAR)
    enzyme,nt,buffer=_parts(g,'barcoding_mix',mix)
    bst=g.input('bst','Source separately added Bst2.0 WarmStart polymerase',K.CATALYST)
    nt_buffer=g.input('nt_buffer','Source NexteraXT NT buffer, Tween and PEG barcoding conditions',K.CONDITION)
    merged=g.product('merged_droplets','Actual sequentially merged gel/barcode/PCR droplet collection; internal membership unassigned',K.MOLECULAR)
    merger=g.input('merger_conditions','Source microfluidic merger then HFE-to-FC40 carrier replacement before thermal cycling',K.CONDITION)
    removed_phase=g.product('merger_removed_carrier','Source carrier removed before barcoding PCR; carryover unassigned',K.MOLECULAR)
    g.react('droplet_merger','Physically merge gel/barcode/PCR droplets before molecular barcoding',S.GENERIC,((fragments,R.TEMPLATE),(clonal,R.TEMPLATE),(f127,R.STRUCTURAL),(f129_join,R.STRUCTURAL),(enzyme,R.CATALYST),(nt,R.STRUCTURAL),(buffer,R.CONDITION),(bst,R.CATALYST),(nt_buffer,R.CONDITION),(merger,R.CONDITION)),(merged,removed_phase),changes_molecule=False)
    g.mark_terminal(Terminal(removed_phase.id,Disposition.UNRESOLVED,'Source carrier replacement; no quantitative recovery or perfectly isolated compartments implied.'))
    mt=g.product('merged_genome','Copresent tagmented genomic class in merged droplets',K.MOLECULAR,state=TAGMENTED.scene())
    mb=g.product('merged_barcode','Copresent barcode-factory product class and possible residuals',K.MOLECULAR)
    mf=g.product('merged_fl127','Copresent FL127 from the actual later mix',K.OLIGO,state=f127.state)
    mr=g.product('merged_fl129','Copresent fresh FL129 from the actual later mix',K.OLIGO,state=f129_join.state)
    me=g.product('merged_enzyme','Copresent Platinum component',K.CATALYST)
    mn=g.product('merged_nt','Copresent prepared-mix nucleotide component',K.MOLECULAR)
    mc=g.product('merged_conditions','Copresent nonincorporated conditions of the same droplets',K.CONDITION)
    ms=g.product('merged_bst','Copresent separately supplied Bst component',K.CATALYST)
    g.react('merged_components','Inventory view of the same merged droplets; not separated physical fractions',S.GENERIC,((merged,R.TEMPLATE),),(mt,mb,mf,mr,me,mn,mc,ms),changes_molecule=False)
    joined=g.product('joined','Cell-barcoded genomic fragments in PCR droplets',K.MOLECULAR)
    g.react('genome_barcoding','Source molecular barcoding after physical merger:65°C5min then30cycles',S.PCR,
        ((mt,R.TEMPLATE),(mb,R.TEMPLATE),(mf,R.STRUCTURAL),(mr,R.STRUCTURAL),(me,R.CATALYST),(mn,R.STRUCTURAL),(mc,R.CONDITION),(ms,R.CATALYST)),(joined,),
        geometry=PCRGeometry.unavailable(TAGMENTED.scene(),reason='Source describes handle-mediated splicing; full genomic substrates, adapter-end combinations and final product geometry are unassigned in this capped reading.'))
    pipette=g.input('pipette','Micropipette separation of large coalesced droplets',K.SELECTION)
    retained_drops=g.product('retained_drops','Remaining source barcoding droplets',K.MOLECULAR)
    coalesced=g.product('coalesced','Large coalesced droplets physically removed',K.WASTE)
    g.react('drop_selection','Remove large coalesced droplets afterPCR',S.SELECTION,
        ((joined,R.TEMPLATE),(pipette,R.SELECTION)),(retained_drops,coalesced))
    g.discard(coalesced,'Physically removed by micropipette before breaking remaining emulsion')
    recovered=cleanup(g,'emulsion_recovery',retained_drops,'Perfluorooctanol/centrifugation; aqueous transfer; Zymo DNA Clean&Concentrator','Recovered barcoded DNA')
    library=cleanup(g,'size_selection',recovered,'AMPure selection for200–600bp','SiC-seq genomic sequencing library',K.LIBRARY)
    library=_qc(g,'library_qc',library,'Source Bioanalyzer analytical allocation')
    loading=g.product('loading','Source sequencing-loading contribution; amount unassigned',K.LIBRARY)
    stock=g.product('library_stock','Untransferred library stock; amount unassigned',K.LIBRARY)
    g.react('loading_allocation','Select the source loading contribution',S.SPLIT,((library,R.TEMPLATE),),(loading,stock))
    g.mark_terminal(Terminal(stock.id,Disposition.UNRESOLVED,'No entire-stock consumption or additional sample pool invented.'))
    g.sequence(loading,SequencingPlan(acquired_roles=('Read 1','Custom barcode index'),
        read_lengths={'Read 1':150},
        limitation='Synthetic-community MiSeq single-end150 and custom index primer specified. Source barcode width15 is a molecular feature, not an assumed index-cycle schedule; final strand/primer coordinates remain unassigned.'))
    g.declare_reviewed_scope(('cell_microgels','lytic_enzymes','lysed','cell_wall_lysis','first_lysate_separation_selection','first_lysate_separation_kept','first_lysate_separation_waste','first_lysate_separation',
        'protease','proteolysed','proteolysis','microgel_washes_selection','microgel_washes_kept','microgel_washes_waste','microgel_washes','genome','retained_state',
        'tn5','loaded_adapters','fragments','tagment','barcode_template','fl128','fl129','phusion','clonal_barcodes','barcode_pcr','fl127','barcoding_mix','joined','genome_barcoding',
        'pipette','retained_drops','coalesced','drop_selection','emulsion_recovery_selection','emulsion_recovery_kept','emulsion_recovery_waste','emulsion_recovery',
        'size_selection_selection','size_selection_kept','size_selection_waste','size_selection')+ADDITIONAL_INVENTORY,SOURCE)
    return g.require_publishable()

# Curator-written obligations from the two inspected source scopes, independent
# of graph enumeration. Component views expose the same physical stocks/droplets.
ADDITIONAL_INVENTORY=('tagmentation_stock','tagmentation_components','tagmentation_buffer','cells','agarose','encapsulation_conditions','prepare_microgels','initial_carrier',
 'phusion_components','phusion_enzyme','phusion_nt','phusion_buffer','barcode_droplet_conditions','barcode_droplets','barcode_removed_carrier','barcode_encapsulation',
 'factory_template','factory_fl129','factory_fl128','factory_enzyme','factory_nt','factory_buffer','factory_components',
 'barcode_qc','barcode_qc_analytical','barcode_qc_continuing','fl129_join','barcoding_mix_components','barcoding_mix_enzyme','barcoding_mix_nt','barcoding_mix_buffer',
 'bst','nt_buffer','merged_droplets','merger_conditions','merger_removed_carrier','droplet_merger','merged_genome','merged_barcode','merged_fl127','merged_fl129','merged_enzyme','merged_nt','merged_conditions','merged_bst','merged_components',
 'library_qc','library_qc_analytical','library_qc_continuing','loading','library_stock','loading_allocation')


def render_page():
    g=protocol_graph()
    return assemble_page(g,(head('SiC-seq'),'<div class="wrap" data-protocol-complete="true"><h1>SiC-seq</h1>',
      '<p class="research-notes"><a href="01_sic_seq.html">Research notes</a></p>',
      info('<a href="https://doi.org/10.1038/nbt.3880">Lan et al. (2017)</a>. '+SOURCE),
      '<h2>Keep genomes together while adding a cell barcode</h2>',
      '<p>The source encapsulates microbial cells in agarose, transfers the gels into aqueous carrier, then lyses and washes them. High-molecular-weight DNA is retained in the gel. Re-encapsulation with Nextera reagents fragments and tags the genome; retention of the associated material is a qualitative source declaration, not a simulated diffusion or recovery result.</p>',
      '<p>A separate droplet factory amplifies random cell-barcode templates. The source merges gel droplets, barcode droplets and PCR solution before splicing the barcode onto genomic fragments. Coalesced droplets are removed; the remaining emulsion is broken, its aqueous DNA recovered, size-selected and sequenced.</p>',
      info('Source-qualified synthetic-community preparation, not a quantitative single-cell capture model. Droplet occupancy, clonality, cross-contamination, yields and complete occupied library geometry remain unassigned. Aggregate component views expose the same stocks and droplets; they do not introduce fresh doses or separated physical fractions. Native and tagmented gel retention are distinct qualitative states. Barcode microscopy sampling and stain exposure of continuing droplets are unassigned. No arbitrary genomic length is assigned.'),
      OligoSection('Known source oligos',tuple(OligoEntry(name,tuple(parts))for name,parts in [('Barcode template',barcode_template()),('FL127',[Segment('FL127',FL127)]),('FL128',[Segment('FL128',FL128)]),('FL129',[Segment('FL129',FL129)]),('Custom barcode-index primer',[Segment('custom index primer',CUSTOM_INDEX)])]),'The 15-base random window identifies a cell-barcode group. It is not a UMI or sequencer sample index; its molecular width does not set the index-cycle schedule.'),
      '<h2>Readout and source limits</h2>',
      '<p>The associated manuscript reports single-end 150-base MiSeq reads for the artificial microbial community. The custom barcode-index primer is known; exact final primer coordinates and index cycles remain unassembled. Read grouping and quality filtering occur computationally after sequencing.</p>',
      '<p>The corrected publisher supplement records a 21 June 2017 correction to the taxonomy display in Figure 11a. That display is outside this molecular model; no taxonomy result or whole-version equivalence is inferred from it.</p>',
      '<details><summary>Full preparation workflow and unresolved chemistry</summary>',WorkflowSection('Cell microgels → genome retention → tagmentation → droplet barcoding → library'),'</details>',PlanLimitationSection('loading'),'</div>'))

def selftest():
    import re
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
        read_length=re.search(r'single-endreadsof(\d+)bp',text)
        c('Reported synthetic-community read length remains source-derived',g.sequencing['loading'].read_lengths['Read 1'],int(read_length[1]) if read_length else None)
        c('Source reports fresh primer contributions in both factory and barcoding mixtures',text.count('400nMprimersFL')>=2,True)
    else:c.skip('Source oligos','External defining article unavailable')
    c('Later FL129 supply is not the consumed barcode-factory contribution',any(p.material=='fl129_join' for p in g.steps['droplet_merger'].inputs) and all(p.material!='fl129' for p in g.steps['droplet_merger'].inputs),True)
    c.report(exit_on_fail=True)

"""Source-owned CAP-seq microbial DNA workflow, explicit bioRxiv v3 scope."""
from dataclasses import dataclass,replace
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import nextera as nx
from chemdraw import Scene,feature,workflow_from_protocol_graph,workflow_panel
from protocol_model import role
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SequencingPlan,SourceItem
from pcr import PCRGeometry
from page import head,info
PAGE_FILE='cap-long.html'
SOURCE='Li et al.,bioRxiv10.1101/2024.09.10.612220v3,posted13October2025;Methodspp19–23,Fig1b/SuppFig1,originalp22rendered'
MIX_COMPONENTS_UL=(120,16,16,12,64,64,64,16)
PRINTED_MIX_TOTAL_UL=348
COMPONENT_TOTAL_UL=sum(MIX_COMPONENTS_UL)
ANALYSIS_MOTIF=nx.S7
CAP_BARCODE=feature('cap_microbe.capsule','cell_barcode','unknown')
ONT_BARCODE=feature('cap_microbe.native_index','sample_index','unknown')
DISCREPANCY=f'Preprintv3 p22 calls the linear-amplification mix {PRINTED_MIX_TOTAL_UL}µL, but its listed volumes sum to {COMPONENT_TOTAL_UL}µL. Both source values are retained; no normalized formulation/concentrations are selected.'

@dataclass(frozen=True)
class Preparation:
    """Selected preparation class; features/carryover flags are not purity or yield predictions."""
    compartment:str
    dna_class:str
    capsule_barcoded:bool=False
    capsule_remnants:bool=False
    native_barcoded:bool=False
    def __post_init__(self):
        if self.compartment not in {'CAP','droplet+CAP','bulk'}:raise ValueError('source compartment classes only')
        if self.capsule_remnants and self.compartment!='bulk':raise ValueError('dissolution remnants belong to the represented released bulk preparation')
    def pair(self):
        if self.compartment!='CAP':raise ValueError('source bead pairing starts from CAP preparation')
        return replace(self,compartment='droplet+CAP')
    def barcode_selected(self):
        if self.compartment!='droplet+CAP':raise ValueError('source barcode incorporation occurs in paired droplets')
        return replace(self,capsule_barcoded=True,dna_class='New barcoded DNA copies plus original/tagmented/unreacted DNA')
    def break_emulsion(self):
        if self.compartment!='droplet+CAP':raise ValueError('source barcoding emulsion required')
        return replace(self,compartment='CAP')
    def dissolve(self):
        if self.compartment!='CAP':raise ValueError('source capsule dissolution acts on recovered CAPs')
        return replace(self,compartment='bulk',capsule_remnants=True)
    def rows(self):
        sc=Scene();sc.strand('dna',[role(self.dna_class)],label=self.compartment+' preparation; heterogeneity/carryover unquantified')
        if self.capsule_barcoded:
            sc.strand('feature',[role('capsule identity barcode',molecular_feature=CAP_BARCODE)],label='Known copied feature inventory; strand/order/width unassembled')
        if self.native_barcoded:
            sc.strand('native_feature',[role('native sample barcode',molecular_feature=ONT_BARCODE)],label='Source native-barcode kit identity inventory; exact product/position/width unassembled')
        if self.capsule_remnants:
            sc.strand('polymer',[role('dissolved capsule polymers/reaction species')],label='Possible carryover;source recovery/purity unquantified')
        sc.note('dna','CAP retention is source size-dependent (~400bp assay cutoff), not an exact molecular sieve or pure single-cell guarantee.')
        return sc.rows()

def barcode_inventory(which):
    sc=Scene();f=CAP_BARCODE if which=='capsule' else ONT_BARCODE
    sc.strand('feature',[role(which+' barcode component',molecular_feature=f)],label='Component identity only; full oligo polymer/order/width not assembled')
    return sc

def protocol_graph():
    g=ProtocolGraph('CAP-seq microbial long-read route — preprintv3',authority=SOURCE,representation='published_instance')
    scope=[]
    def inp(key,label,kind,**kw):scope.append(key);return g.input(key,label,kind,**kw)
    def op(key,label,kind,material,reagents,state,**kw):
        out=g.product(key+'_mixture',label+' reaction preparation',K.MOLECULAR,state=state)
        g.react(key,label,kind,((material,R.TEMPLATE),*reagents),(out,),**kw);scope.append(key);return out
    def clean(key,material,label,state):
        reagent=inp(key+'_procedure',label,K.SELECTION)
        kept=g.product(key,'Source physically recovered preparation',K.MOLECULAR,state=state)
        waste=g.product(key+'_reject','Removed source workup fractions; DNA loss/carryover unquantified',K.WASTE)
        g.react(key+'_separation',label,S.CLEANUP,((material,R.TEMPLATE),(reagent,R.SELECTION)),(kept,waste));scope.append(key+'_separation');g.discard(waste,'Actual source extraction/washing/purification boundary;not perfect class-specific removal.');return kept
    start=Preparation('CAP','Genomic DNA in source cells; capsule occupancy and lysis efficiencies unassigned')
    cells=inp('cells','Source microbial-cell suspension;culture/stool processing upstream unassembled',K.BIOLOGICAL,primary=True)
    phases=inp('atps','Source dextran/PEGDA/LAP two-phase precursor and fluorinated oil/surfactant;full material chemistry unassembled',K.CONDITION)
    cur=op('encapsulate','Flow-focus source cell suspension into two-phase droplets',S.GENERIC,cells,((phases,R.CONDITION),),start)
    uv=inp('uv','Source200WUV4min;gel outerPEGDA phase',K.CONDITION)
    cur=op('gel','UV polymerize source capsule shell;DNA damage/occupancy outcomes not predicted',S.GENERIC,cur,((uv,R.CONDITION),),start)
    cur=clean('caps_recovered',cur,'Break initial oil with20%PFO/HFE7500;wash3xDPBS/0.1%Tween20',start)
    lyso=inp('wall_enzymes','Source lysozyme+lysostaphin/TE/NaCl/EDTA/Triton37°C2h',K.CATALYST)
    lysed=replace(start,dna_class='Source partially lysed cell/genomic-DNA and retained biomolecule population')
    cur=op('wall_lysis','Source enzymatic cell-wall lysis in capsules',S.DIGESTION,cur,((lyso,R.CATALYST),),lysed)
    cur=clean('wall_washed',cur,'Wash3xPBS/10mMEDTA',lysed)
    pk=inp('pk','SourceSDS/ProteinaseK55°C30min',K.CATALYST)
    lysate=replace(start,dna_class='Source capsule-confined genomic DNA/carryover after lysis')
    cur=op('protein_lysis','Source detergent/protease lysis;digested species persist until washing',S.DEPROTEINIZATION,cur,((pk,R.CATALYST),),lysate)
    cur=clean('lysis_washed',cur,'Wash3xsourceHEPES/Tween20/EDTA/NaCl;resuspendDPBS/Tween',lysate)
    mda=inp('mda_mix','SourcePhi29/pyrophosphatase/dNTPs/buffer/surfactants;full initiation/priming geometry unassembled',K.CATALYST)
    amplified=replace(lysate,dna_class='Original templates plus heterogeneous Phi29-amplified genomic DNA/carryover')
    cur=op('mda','Source MDA30°C2h then65°C10min;unassembled initiation chemistry',S.GENERIC,cur,((mda,R.CATALYST),),amplified)
    cur=clean('mda_washed',cur,'Wash3xsourcewashingbuffer;resuspendDPBS/Tween',amplified)
    tag=inp('tag_mix','SourceABclonalRK20239 TagBuffer/TagEnzyme50,55°C10min;historical ends/reaction geometry unassembled',K.CATALYST)
    tagmented=replace(amplified,dna_class='Original/untagged DNA and tagmented genomic fragments in source CAPs')
    adapter=inp('tag_adapter','Unassembled adapter component of source tagmentation kit;full source Table1 and loading/polymer geometry unreviewed',K.ADAPTER)
    cur=op('tagment','Source capsule-contained tagmentation',S.TAGMENTATION,cur,((tag,R.CATALYST),(adapter,R.STRUCTURAL)),tagmented)
    cur=clean('tag_washed',cur,'Wash3xsourcewashingbuffer;resuspendDPBS/Tween',tagmented)
    bead=inp('barcode_beads','DTT-dissolvable hydrogel barcode beads,PercentBioRRF0991;primers/adapters unassembled',K.OLIGO,state=barcode_inventory('capsule'))
    paired=tagmented.pair()
    cur=op('pair_bead','Co-encapsulate source CAPs and barcode beads in droplets;no universal one-cell-one-bead guarantee',S.GENERIC,cur,((bead,R.STRUCTURAL),),paired)
    mix=inp('linear_mix','PhusionHF+DeepVent+ETSSB+DTT/buffer/stabilizer/dNTPs;'+DISCREPANCY,K.CATALYST)
    dissolved=paired # DTT dissolves barcode beads, not the PEGDA capsule shell.
    cur=op('dissolve_beads','SourceDTT releases bead primers within paired droplets;capsules remain',S.GENERIC,cur,((mix,R.CATALYST),),dissolved)
    cur=clean('exchange_oil',cur,'Replace original emulsification oil withBioRadQX200oil before thermalcycling;capsules/droplets retained,removed oil fraction explicit',paired)
    copied=paired.barcode_selected()
    cur=op('linear_copy','Source72°C5min/98°C3min then15cycles98/60/72;linear barcoding copies DNA',S.EXTENSION,cur,((bead,R.STRUCTURAL),(mix,R.CATALYST)),copied)
    recovered=copied.break_emulsion()
    cur=clean('barcoded_caps',cur,'Break droplets20%PFO/HFE7500;washCAPs3xTE/0.1%Tween20',recovered)
    alkali=inp('capsule_alkali','Mix20µLCAPs with40µL0.8MNaOH;55°C10min',K.CONDITION)
    released=recovered.dissolve()
    cur=op('dissolve_capsules','Alkaline PEGDA capsule dissolution releases preparation to bulk',S.GENERIC,cur,((alkali,R.CONDITION),),released)
    acid=inp('acetic_acid','Add40µL0.8Maceticacid;source neutralization',K.CONDITION)
    cur=op('neutralize','Neutralize source released mixture;not physical purification',S.QUENCHING,cur,((acid,R.CONDITION),),released)
    # Cleanup has unquantified carryover; the retained state no longer asserts a polymer-free product.
    recovered_bulk=replace(released,dna_class='Recovered source DNA;residual dissolved polymers/reagents unquantified')
    cur=clean('bulk_recovered',cur,'SourceAMPureXP DNApurification',recovered_bulk)
    f=inp('pcr_f','Sourceforwardprimer,actualSupplementaryTable1 not yet located',K.OLIGO)
    r=inp('pcr_r','Sourcereverseprimer,actualSupplementaryTable1 not yet located',K.OLIGO)
    kapa=inp('kapa','KAPAhotstartHiFipolymerasemixture;source16cycles98/60/72',K.CATALYST)
    amplified_bulk=replace(recovered_bulk,dna_class='Amplified capsule-barcoded DNA plus original DNA/reagents and possible polymer carryover')
    cur=op('bulk_pcr','Source bulk16cycle amplification after capsule dissolution',S.PCR,cur,((f,R.STRUCTURAL),(r,R.STRUCTURAL),(kapa,R.CATALYST)),amplified_bulk,geometry=PCRGeometry.unavailable(Scene(),reason='Actual source primers/table and full tagmented/barcoded genomic strands not assembled;no prototype sequences supplied.'))
    cur=clean('pcr_recovered',cur,'Source magnetic-bead purification;BioanalyzerQC',amplified_bulk)
    ligation=inp('ont_ligation','Source SQK-LSK114 ligationkit;historical manual/end steps unassembled',K.CONDITION)
    indexing=inp('ont_indexing','Source SQK-NBD114.24 nativebarcodingkit;exact supplied barcode member unassigned',K.ADAPTER,state=barcode_inventory('native_sample'))
    lib=op('ont_handoff','Unassembled source ONTV14/nativebarcoding library handoff',S.MAGIC,cur,((ligation,R.CONDITION),(indexing,R.STRUCTURAL)),replace(amplified_bulk,native_barcoded=True),reason='Main names two ONTkits and a generic manufacturerURL;historical source-specific end/ligation/barcoding order must be traced,not imported from a current/default protocol.')
    g.sequence(lib,SequencingPlan(acquired_roles=('long-read DNA',),limitation='Source reports ONT long reads but mixes instrument/flow-cell labels:GridIONMK1/P2solo/PromethION/FLO-MIN114. Exact platform/run identity and molecule/adapter geometry remain unassembled.'))
    g.declare_reviewed_scope(scope,SOURCE)
    g.add_source_item(SourceItem('supplementary_table1','Actual source-cited SupplementaryTable1 primeroriginal,not located in capped pass','Mainv3Methods887 andsupplindex;no confirmed absence.'))
    g.add_source_item(SourceItem('historical_ont','Historical ONTkit/reagent/end geometry and instrumentidentity','Source generic manufacturerURL;do not substitute current vendor revision.'))
    g.add_source_item(SourceItem('source_recipe','Unreconciledsource mix totals/fullMDAinitiation/compartmentpolymers','Preprintv3p22 andrelatedMethods;no full molecular validation.'))
    return g

def render_page():
    g=protocol_graph()
    return head('CAP-seq microbial long-read DNA — preprintv3')+'<main class="wrap"><h1>CAP-seq (long-read microbiome)</h1><p><a href="01_cap-long.html">Research notes</a></p>'+info('Private incomplete model using bioRxivv3,13October2025. Exact primer table,full compartment/end chemistry and ONT handoff remain unresolved. This is a separate DNA method from RNA CapSeq and CXXC CAP-seq.')+info(DISCREPANCY)+workflow_panel(workflow_from_protocol_graph(g,'ont_handoff_mixture'))+'</main>'

def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();a=Preparation('CAP','genomicDNA');b=a.pair().barcode_selected().break_emulsion();d=b.dissolve()
    c('Bead dissolution/barcoding does not itself remove capsule',(b.compartment,b.capsule_barcoded),('CAP',True))
    c('Capsule release preserves copied identity and physical polymer remnants',(d.compartment,d.capsule_barcoded,d.capsule_remnants),('bulk',True,True))
    try:a.barcode_selected();bad=False
    except ValueError:bad=True
    c('Selected source barcoding cannot precede droplet pairing',bad,True)
    g=protocol_graph();c('Incomplete source/geometry stays unpublished',g.publishable,False)
    c('Private preview has no finished marker','data-protocol-complete=' in render_page(),False)
    p=sources_root()/HERE.parent.name/'2024.09.10.612220_v3.full.pdf.txt'
    if p.exists():
        text=' '.join(p.read_text().split());c('Exact v3 identity','this version posted October 13, 2025' in text,True)
        c('Published analysis motif,not whole sequencing primer',ANALYSIS_MOTIF in text,True)
    else:c.skip('Original v3 identity/motif','External source unavailable')
    c.report()

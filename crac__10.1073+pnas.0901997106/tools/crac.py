"""Original CRAC NGS operation order; original Sanger linkers are separate inventory."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
SOURCE='Granneman et al.2009 DOI10.1073/pnas.0901997106; published PMC2688437 Fig1, main library1143–1149 and Solexa1171–1179; SI/TableS5 unavailable'

@dataclass(frozen=True)
class NucleotideBlock:
    base: str
    sugar: str
    orientation: str
    def __post_init__(self):
        if (self.base,self.sugar,self.orientation) not in {('C','2,3-dideoxy','normal'),('T','2,3-dideoxy','reported_inverted')}:raise ValueError('only source-reported original linker blocks represented')

@dataclass(frozen=True)
class OriginalLinker:
    role: str
    body_polymer: str
    block_end: str
    block: NucleotideBlock
    def __post_init__(self):
        expected={'three_prime':('DNA','3-prime','C','normal'),'five_prime':('DNA–RNA hybrid','5-prime','T','reported_inverted')}
        if self.role not in expected or (self.body_polymer,self.block_end,self.block.base,self.block.orientation)!=expected[self.role]:raise ValueError('source linker role/polymer/block identity must agree')
    def ordinary_linear_dna(self):raise ValueError('source body and noncanonical terminal linkage unassembled; cannot flatten block into ordinary DNA')

def original_sanger_linkers():
    return (OriginalLinker('three_prime','DNA','3-prime',NucleotideBlock('C','2,3-dideoxy','normal')),
            OriginalLinker('five_prime','DNA–RNA hybrid','5-prime',NucleotideBlock('T','2,3-dideoxy','reported_inverted')))

def unreviewed_end(label):
    sc=Scene();sc.strand('adapter',[Segment('unassembled source end','[unreviewed '+label+' end]',placeholder=True)])
    return AdapterLigationGeometry('unavailable','unavailable',(sc,),'Original Solexa linker identity/end chemistry is in uninspected SI/TableS5; not a source-absence claim.')

def protocol_graph():
    g=ProtocolGraph('CRAC — original in-vivo/Solexa route',authority=SOURCE,representation='published_instance')
    prior=g.input('tev_eluate','In-vivo UV-crosslinked HTP-RNP TEV eluate; upstream exact dose/capture workup unassembled',K.MOLECULAR,primary=True)
    rnase=g.input('rnases','Source RNaseA/T1 partial RNA digestion',K.CATALYST)
    dig=g.product('partial_digest','RNP-protected RNA fragments plus other cleaved/uncut RNA/protein material; no physical discard',K.MOLECULAR)
    g.react('trim_rna','Source partial RNA digestion before denaturing nickel capture',S.DIGESTION,((prior,R.TEMPLATE),(rnase,R.CATALYST)),(dig,))
    denaturant=g.input('guanidine','Actual6Mguanidine-HCl molecular denaturant; covalentRNA-crosslinks not broken by this declaration',K.MOLECULAR)
    denatured=g.product('denatured','Denatured mixture retains covalent protein–RNA classes; noncovalent particles disrupted, purity unassigned',K.MOLECULAR)
    g.react('denature','Source denaturing condition; no covalent crosslink breaking inferred',S.GENERIC,((dig,R.TEMPLATE),(denaturant,R.STRUCTURAL)),(denatured,),changes_molecule=False,reason='Known denaturing operation, with uncomputed conformation/mixture composition')
    def select(key,x,label,procedure,kind=K.MOLECULAR):
        reagent=g.input(key+'_reagent',procedure,K.MOLECULAR if key=='tca_recovery' else K.SELECTION);keep=g.product(key+'_kept',label,kind);waste=g.product(key+'_waste','Physical '+key+' rejects; exact composition unassigned',K.WASTE)
        g.react(key,label,S.CLEANUP if key=='tca_recovery' else S.SELECTION,((x,R.TEMPLATE),(reagent,R.STRUCTURAL if key=='tca_recovery' else R.SELECTION)),(keep,waste));g.discard(waste,'Actual physical selection rejects only');return keep
    captured=select('nickel_capture',denatured,'Nickel-bound His-tagged RNP after source washing','Source nickel affinity resin and extensive guanidine-removing washes')
    ap=g.input('phosphatase','Source alkaline phosphatase for terminal2′/3′phosphates',K.CATALYST)
    repaired=g.product('repaired','Dephosphorylated RNA-containing RNP population; end fractions and exact structures unassigned',K.MOLECULAR)
    g.react('repair3','Source on-bead RNA-end dephosphorylation',S.DEPHOSPHORYLATION,((captured,R.TEMPLATE),(ap,R.CATALYST)),(repaired,))
    ligase=g.input('ligation_enzyme','Source on-bead RNA-linker ligation enzyme; exact original Solexa formulation unreviewed',K.CATALYST)
    a3=g.input('solexa_adapter3','Solexa-compatible3′linker; exact polymer/modifications uninspected in TableS5',K.ADAPTER)
    linked3=g.product('linked3','3′linker-bearing source RNA/RNP and residual substrate classes',K.MOLECULAR)
    g.react('ligate3','Source3′linker ligation; no Sanger-linker sequence inheritance',S.LIGATION,((repaired,R.TEMPLATE),(a3,R.STRUCTURAL),(ligase,R.CATALYST)),(linked3,),geometry=unreviewed_end('Solexa3′linker'))
    pnk=g.input('pnk','T4polynucleotide kinase',K.CATALYST);atp=g.input('gamma_atp','Source[γ-32P]ATP phosphate donor',K.MOLECULAR)
    labelled=g.product('labelled','5′phosphorylated/radiolabelled RNA-containing RNP; free nucleotide material retained until workup',K.MOLECULAR)
    g.react('label5','Source5′phosphorylation/radiolabelling between linker ligations',S.PHOSPHORYLATION,((linked3,R.TEMPLATE),(pnk,R.CATALYST),(atp,R.STRUCTURAL)),(labelled,))
    a5=g.input('solexa_adapter5','Solexa-compatible5′linker; exact TableS5 polymer/end topology unreviewed',K.ADAPTER)
    linked5=g.product('linked5','Both-linker-bearing RNP population; exact geometry/unreacted classes unassigned',K.MOLECULAR)
    g.react('ligate5','Source5′linker ligation on nickel beads',S.LIGATION,((labelled,R.TEMPLATE),(a5,R.STRUCTURAL),(ligase,R.CATALYST)),(linked5,),geometry=unreviewed_end('Solexa5′linker'))
    imidazole=g.input('imidazole','Actual molecular imidazole contribution forsource RNPelution; doseunreviewed',K.MOLECULAR)
    eluate=g.product('eluate','Eluted source RNP mixture, including unassigned soluble carryover',K.MOLECULAR)
    g.react('elute','Source nickel elution before TCA',S.GENERIC,((linked5,R.TEMPLATE),(imidazole,R.STRUCTURAL)),(eluate,),changes_molecule=False,reason='Reported affinity elution; complete recovery and binding changes not calculated')
    precip=select('tca_recovery',eluate,'TCA-precipitated protein/RNA complexes','Source TCA precipitation and recovery')
    selected=select('gel_membrane',precip,'Source radiolabel-selected nitrocellulose RNP region','BisTrisNuPAGE, nitrocellulose transfer, selected band excision; NGS-specific limits unreviewed')
    pk=g.input('proteinase','Source proteinaseK release from membrane slices',K.CATALYST)
    release=g.product('released','RNA-containing proteolysis mixture; peptide/crosslink remnants not silently removed',K.MOLECULAR)
    g.react('proteolyze','Source proteinaseK treatment',S.DEPROTEINIZATION,((selected,R.TEMPLATE),(pk,R.CATALYST)),(release,))
    ready=g.product('recovered_rna','Source recovered RNA for RT; exact SI extraction/workup unassembled',K.MOLECULAR)
    g.react('rna_recovery','Unassembled original RNA-recovery procedure',S.MAGIC,((release,R.TEMPLATE),),(ready,),reason='Main reports recovered RNA; gated SI determines exact extraction/cleanup and remnant fate.')
    primer=g.input('rt_primer','Original Solexa-linker RT primer; sequence/version unreviewed',K.OLIGO);rt=g.input('rt_enzyme','Source reverse transcriptase; exact SI identity/conditions unreviewed',K.CATALYST)
    rt_supplies=g.input('rt_supplies','Actual RT aggregate functional supplies; sourceSolexaidentity/amounts/stockpackaging remainunreviewed',K.MOLECULAR)
    rt_nt=g.product('rt_nts','Copresent DNA nucleotide functional component, not an extra dose',K.MOLECULAR)
    rt_buffer=g.product('rt_buffer','Copresent RT buffer functional component; exact recipe unreviewed',K.CONDITION)
    g.react('rt_inventory','Actual copying-supply bookkeeping; not physical separation or claimed premix',S.GENERIC,((rt_supplies,R.TEMPLATE),),(rt_nt,rt_buffer),changes_molecule=False)
    hybrid=g.product('hybrid','Source cDNA with RNA template/carryover retained; observed mutation/deletion signatures are copying outcomes',K.MOLECULAR,chemical_state={'RNA_retained':True,'cDNA_effect':'crosslink-associated errors inferred, no DNA mutation assigned to RNA'})
    g.react('rt','Source linker-primed cDNA synthesis',S.REVERSE_TRANSCRIPTION,((ready,R.TEMPLATE),(primer,R.STRUCTURAL),(rt,R.CATALYST),(rt_nt,R.STRUCTURAL),(rt_buffer,R.CONDITION)),(hybrid,))
    f=g.input('pcr_forward','Original Solexa-compatible PCR forward primer; unreviewed TableS5',K.OLIGO);r=g.input('pcr_reverse','Original Solexa-compatible PCR reverse primer; unreviewed TableS5',K.OLIGO);poly=g.input('polymerase','Source PCR enzyme/formulation; exact SI details unreviewed',K.CATALYST)
    pcr_supplies=g.input('pcr_supplies','Actual PCR aggregate functional supplies; sourceSolexaformulation/stockpackaging remainunreviewed',K.MOLECULAR)
    pcr_nt=g.product('pcr_nts','Copresent actual PCR nucleotide functional class',K.MOLECULAR)
    pcr_buffer=g.product('pcr_buffer','Copresent actual PCR buffer functional class',K.CONDITION)
    g.react('pcr_inventory','Actual amplification-supply bookkeeping; no independent fresh component doses',S.GENERIC,((pcr_supplies,R.TEMPLATE),),(pcr_nt,pcr_buffer),changes_molecule=False)
    amplified=g.product('amplified','Source amplified cDNA plus original templates and other reaction material',K.MOLECULAR)
    g.react('pcr','Source PCR; no inferred cycle count or full primer geometry',S.PCR,((hybrid,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(poly,R.CATALYST),(pcr_nt,R.STRUCTURAL),(pcr_buffer,R.CONDITION)),(amplified,),geometry=PCRGeometry.unavailable(Scene(),reason='Original Solexa template/primer polymers and strand assembly unreviewed.'))
    final=select('pcr_gel',amplified,'Source gel-purified PCR material for Solexa; applicable NGS size window unassembled','Source PCR-product gel purification; Sanger>60bp criterion not silently assigned to NGS',K.LIBRARY)
    g.sequence(final,SequencingPlan(acquired_roles=('Read 1',),limitation='Source direct Solexa single-end readout; actual run cycles and full library/run-primer geometry unreviewed.'))
    g.declare_reviewed_scope(['tev_eluate','rnases','trim_rna','guanidine','denature','nickel_capture_reagent','nickel_capture','phosphatase','repair3','ligation_enzyme','solexa_adapter3','ligate3','pnk','gamma_atp','label5','solexa_adapter5','ligate5','imidazole','elute','tca_recovery_reagent','tca_recovery','gel_membrane_reagent','gel_membrane','proteinase','proteolyze','rt','pcr','pcr_gel_reagent','pcr_gel'],SOURCE)
    g.declare_reviewed_scope(('rt_supplies','rt_inventory','rt_nts','rt_buffer','pcr_supplies','pcr_inventory','pcr_nts','pcr_buffer'),SOURCE)
    g.add_source_item(SourceItem('unreviewed_original_si','Original SI/TableS5 and exact Solexa-specific end/primer/workup/read geometry','PMC2688437 original SI index/PDF returned404'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();a,b=original_sanger_linkers()
    c('original inverted endpoint remains a thymidine nucleotide',(b.block.base,b.block.orientation,b.block_end),('T','reported_inverted','5-prime'))
    try:b.ordinary_linear_dna();bad=False
    except ValueError:bad=True
    c('inverted/body-unknown linker cannot silently flatten into ordinary DNA',bad,True)
    g=protocol_graph();c('RNA remains present after modeled copying',g.materials['hybrid'].chemical_state['RNA_retained'],True)
    c('NGS adapter material does not adopt separate Sanger inventory',g.materials['solexa_adapter5'].chemical_state,None)
    source=Source(sources_root()/'crac__10.1073+pnas.0901997106'/'pnas.0901997106_PMC2688437.html.txt','published main1147/1171')
    if have(c,source,label='source platform distinction'):
        text=' '.join(source.path.read_text().split());c('original explicitly changes linkers for Solexa','linkers that are compatible with Illumina Solexa' in text,True)
    c('Guanidine is a molecular contribution',g.materials['guanidine'].kind,K.MOLECULAR)
    c('TCA is a molecular precipitation reagent',g.materials['tca_recovery_reagent'].kind,K.MOLECULAR)
    c('Actual Solexa-specific SI item remains uncovered','unreviewed_original_si' in g.source_coverage,False)
    c('Copying NT components have actual supply provenance',g.producer_of('rt_nts').id,'rt_inventory')
    c.report(exit_on_fail=True)
if __name__=='__main__':selftest()

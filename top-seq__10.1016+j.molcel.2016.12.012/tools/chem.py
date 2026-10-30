"""Protocol-owned source-qualified incomplete graph; no publication marker."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry

def scene(label):
    sc=Scene();sc.strand('preparation',[Segment('population','['+label+']',placeholder=True)]);return sc

def op(g,key,x,label,kind=S.MAGIC,*,reagent='source reagents; geometry unassembled',product=None,state=None,catalyst=False):
    stock=g.input(key+'_stock',reagent,K.CATALYST if catalyst else K.CONDITION)
    out=g.product(key+'_product',product or label+'; molecular mixture/remnants unassembled',K.MOLECULAR,state=state)
    kw={'reason':'Source operation established; exact molecular ends, companion species and geometry unassembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
    g.react(key,label,kind,((x,R.TEMPLATE),(stock,R.CATALYST if catalyst else R.CONDITION)),(out,),**kw);return out

def cleanup(g,key,x,label):
    stock=g.input(key+'_condition',label,K.CONDITION)
    out=g.product(key+'_kept','Source recovered preparation; recovery and carryover unquantified',K.MOLECULAR)
    waste=g.product(key+'_waste','Actual source-separated workup fractions',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(stock,R.CONDITION)),(out,waste));g.discard(waste,'Source physical cleanup, no perfect partition implied.');return out

def uncertain_ligation():return AdapterLigationGeometry('unavailable','unavailable',(scene('unreviewed exact ends and junction geometry'),),'Known source ligation, complete molecular geometry not assembled.')

SOURCE='Associated accepted manuscript NIHMS838518/PMC5291818 mainDesign64, Methods153–177, KeyResources333–361; DOI10.1016/j.molcel.2016.12.012; final VOR applicability/chemistry SI unreviewed.'
A1='GATTGGAAGAGTGGTTCAGCAGGAATGCTGAG'
A2='ACACTCTTTCCCTACATGACACTCTTCCAATCT'
TO_METHODS='T(alkyneU)TTATATATTTATTGGAGACTGACTACCAGATGTAACA'
TO_TABLE='(alkyneU)TTATATATTTATTGGAGACTGACTACCAGATGTAACA'
@dataclass(frozen=True)
class SelectedTether:
    dna:str
    cytosine:int
    stage:str='native'
    def __post_init__(self):
        if not self.dna or set(self.dna)-set('ACGT') or type(self.cytosine)is not int or not 0<=self.cytosine<len(self.dna)-1 or self.dna[self.cytosine:self.cytosine+2]!='CG':raise ValueError('explicit selected unmodified CpG cytosine required')
        if self.stage not in {'native','azide','tethered'}:raise ValueError('source tag stage required')
    def tag(self):
        if self.stage!='native':raise ValueError('tag only native selected uCG class')
        return replace(self,stage='azide')
    def tether(self):
        if self.stage!='azide':raise ValueError('CuAAC tethering needs azide-tagged DNA')
        return replace(self,stage='tethered')
    def selected_copied_region(self,template_path):
        # Relation for explicitly chosen template coordinates, not full EP-derived
        # first strand or prediction of unresolved proximity-driven invasion.
        if self.stage!='tethered' or not template_path or any(type(i)is not int or not 0<=i<len(self.dna) for i in template_path) or any(a<=b for a,b in zip(template_path,template_path[1:])):raise ValueError('explicit3to5path on tethered DNA required')
        return ''.join({'A':'T','C':'G','G':'C','T':'A'}[self.dna[i]] for i in template_path)

def adapter(seq,name,phosphate=False):
    sc=Scene();sc.strand(name,[Segment(name,seq)],mod5='phosphate' if phosphate else '');return sc

def protocol_graph():
    g=ProtocolGraph('TOP-seq',authority=SOURCE,representation='published_instance')
    x=g.input('genomic','Source human or S.aureus genomicDNA; methylation state/site positions unassigned',K.MOLECULAR,primary=True)
    x=op(g,'sonicate',x,'Source Bioruptor peak150–200bp',S.FRAGMENTATION,reagent='source sonication/EB')
    x=op(g,'mtag',x,'Source eM.SssI0.5–1µM/Ado6azide200µM30°C1h',S.CHEMICAL_CONVERSION,product='Azide-derivatized selected unmodified CpG DNA plus other genomicDNA/enzyme/cofactor remnants; no C substitution')
    x=op(g,'pk',x,'Source ProteinaseK0.2mg/ml55°C30min',S.DEPROTEINIZATION,product='AzideDNA plus proteolytic debris and other material retained until column')
    x=cleanup(g,'tag_clean',x,'Source GeneJetPCRpurification')
    x=op(g,'end_repair',x,'Source DNAEndRepairKit end filling',S.END_REPAIR,reagent='historical DNAEndRepairKit enzymes',catalyst=True)
    x=cleanup(g,'end_clean',x,'Source GeneJetpurification')
    x=op(g,'atail',x,'Source Klenowexo-/dATP0.5mM37°C45min;75°C15min',S.TAILING,reagent='sourceKlenowexo-/dATP',catalyst=True)
    x=cleanup(g,'atail_clean',x,'Source GeneJetpurification')
    a=g.input('a1','Source5phosphate A1',K.ADAPTER,state=adapter(A1,'a1',True));b=g.input('a2','Source A2; no invented phosphate',K.ADAPTER,state=adapter(A2,'a2'));enzyme=g.input('t4ligase','SourceT4DNAligase15U',K.CATALYST)
    y=g.product('ligated','Source A1/A2-adapter-ligatedDNA classes plus unligatedDNA/adapters/reagents; exact fork/junction geometry unassembled',K.MOLECULAR)
    g.react('ligate','Source partially complementary A1/A2 annealing then22°Covernightligate;65°C10min',S.LIGATION,((x,R.TEMPLATE),(a,R.STRUCTURAL),(b,R.STRUCTURAL),(enzyme,R.CATALYST)),(y,),geometry=uncertain_ligation())
    x=cleanup(g,'lig_clean',y,'Source ZymoDNAcleanConcentrator5')
    to=g.input('to','Source alkyneU-bearingTO; Methods165 leadingT versus KeyResources333 unresolved; neither body silently selected',K.OLIGO)
    catalyst=g.input('cu','Source CuBr:TBTA3.3mM/60–65%DMSO',K.CATALYST);y=g.product('tether_mix','Covalently tethered genomicDNA/TO classes plus freeTO/adapters/reagent remnants',K.MOLECULAR)
    g.react('tether','Source CuAAC45°C2.5or6h',S.CHEMICAL_CONVERSION,((x,R.TEMPLATE),(to,R.STRUCTURAL),(catalyst,R.CATALYST)),(y,))
    x=cleanup(g,'tether_clean',y,'Source diluteDMSO<1% thenZymoDNAcleanConcentrator25')
    ep=g.input('ep','SourceEP bodyTGTTACATCTGGTAGTCAGTCTCCAATAAATATAT with customLNA positions unassembled',K.OLIGO);pfu=g.input('pfu','SourcePfu5U',K.CATALYST);y=g.product('top_extension','Conditional EP-primed genomicDNA copies plus original covalently tagged templates/TO/EP/reagents; physical substrate not erased',K.MOLECULAR)
    g.react('top','SourcePfu/nonhomologous proximitypriming95°C1min65°C1min72°C1min',S.EXTENSION,((x,R.TEMPLATE),(ep,R.STRUCTURAL),(pfu,R.CATALYST)),(y,))
    aliquot=g.product('pcr_aliquot','Source5µl TOPmix allocatedtoPCR',K.MOLECULAR);rest=g.product('top_remainder','Remaining TOPmixture; laterfate unassigned',K.MOLECULAR);g.react('allocate','Source5µl of50µlTOPreaction',S.SPLIT,((y,R.TEMPLATE),),(aliquot,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Unallocated source mixture retained.'))
    f=g.input('fusion_f','Source63nt barcode fusionprimer with phosphorothioates; exactPS positions unreviewed',K.OLIGO);r=g.input('fusion_r','Source45nt trP1-A2 primer with phosphorothioates; exactPS positions unreviewed',K.OLIGO);poly=g.input('maxima','SourceMaximaHotStartMasterMix',K.CATALYST);y=g.product('amp','Source15cycleindexedPCRproducts plus original tagged templates/primers/reagents; complete backbone unassembled',K.MOLECULAR)
    g.react('pcr','Source94°C4min then15cycles95/60/72°C1min',S.PCR,((aliquot,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(poly,R.CATALYST)),(y,),geometry=PCRGeometry.unavailable(scene('unreviewedLNA/PS/TOalternative and completePCRgeometry'),reason='Source sequences partially available, exact modifications/end/junction layout unresolved.'))
    x=cleanup(g,'size',y,'Source Zymocolumn thenMagJet~250bp sizeselection')
    g.sequence(x,SequencingPlan(limitation='SourceIonProtonPI single-end acquisition; raw cycle count/primer geometry unreviewed. Source trims70/165nt computationally; these are not physical run lengths.'))
    g.add_source_item(SourceItem('pending','ChemistrySI/table originals and source TO alternatives','Actual PMC gates in note; MethodsTO has extra leadingT versus table. LNA/PS and complete read geometry unresolved, no absence claim.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();native=SelectedTether('AACGTT',2);tag=native.tag().tether()
    c('covalent tethering preserves actual genomic nucleotide inventory',tag.dna,native.dna)
    c('explicit selected template path produces complementary copied region only',tag.selected_copied_region((5,4,3)), 'AAC')
    try:native.tether();bad=False
    except ValueError:bad=True
    c('tether cannot precede installed azide',bad,True)
    try:SelectedTether('AACTTT',2);bad=False
    except ValueError:bad=True
    c('target declaration rejects nonCpG site',bad,True)
    src=Source(sources_root()/'top-seq__10.1016+j.molcel.2016.12.012'/'j.molcel.2016.12.012_PMC5291818.xml.txt','Standard get_sources; actual Methods165 versus KeyResources333–336')
    if have(c,src,label='source TO alternative transcription'):
        text=src.path.read_text();table=''.join(text.split('TO (tethered oligonucleotide):',1)[1].split('A1 (adapter):',1)[0].split())
        c('actual source Methods and TO-table preserve the distinct leadingT alternatives',(TO_METHODS in text,TO_TABLE in table,TO_METHODS in table),(True,True,False))
    g=protocol_graph();c('unallocated TOPreaction remains independent',g.terminals['top_remainder'].disposition,Disposition.UNRESOLVED)
    c('source-limited graph cannot publish',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()

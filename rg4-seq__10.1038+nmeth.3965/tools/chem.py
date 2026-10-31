"""Protocol-owned source-qualified partial model; no completed-publication claim."""
from dataclasses import dataclass,replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
import illumina as il
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from pcr import PCRGeometry

def state(label):
    sc=Scene();sc.strand('population',[Segment(label,'['+label+']',placeholder=True)]);return sc

def clean(g,key,current,procedure):
    reagent=g.input(key+'_reagent',procedure,K.SELECTION)
    kept=g.product(key+'_kept','Recovered '+procedure+'; exact molecular recovery/carryover unquantified',K.MOLECULAR)
    waste=g.product(key+'_waste','Physical separation fractions; no complete target specificity assumed',K.WASTE)
    g.react(key,procedure,S.CLEANUP,((current,R.TEMPLATE),(reagent,R.SELECTION)),(kept,waste));g.discard(waste,'Actual source physical separation only.');return kept
from adapter_ligation import AdapterLigationGeometry
SOURCE='Kwok et al.2016 DOI10.1038/nmeth.3965; published SI Fig1/2,Notes2/3 and Table1; main Methods inaccessible. In-vitro folding assay, not proof of cellular G4.'
A3='AGATCGGAAGAGCACACGTCTG'
RT='CAGACGTGTGCTCTTCCGATCT'
A5='AGATCGGAAGAGCGTCGTGTAGCTCTTCCGATCT'
@dataclass(frozen=True)
class FoldedRNA:
    sequence:str
    condition:str
    selected_stall:int|None=None
    def __post_init__(self):
        if not self.sequence or set(self.sequence)-set('ACGU') or self.condition not in {'Li+','K+','K++PDS'}:raise ValueError('explicit RNA and source folding condition required')
        if self.selected_stall is not None and (type(self.selected_stall)is not int or not 0<=self.selected_stall<len(self.sequence)):raise ValueError('selected stall must identify a template position; no fold prediction')
    def first_dna(self):
        start=0 if self.selected_stall is None else self.selected_stall+1
        return revcomp(self.sequence[start:].replace('U','T'))

def adapter_scene(which):
    sc=Scene()
    if which=='3':sc.strand('a3',[Segment('source DNA3adapter',A3)],label='DNA3′adapter',mod5='rApp',mod3='C3 spacer')
    elif which=='5':sc.strand('a5',[Segment('source DNA5adapter',A5),Segment('source six random nucleotides','NNNNNN',placeholder=True)],label='DNA5′adapter; random sequence not assigned identifier meaning',mod5='phosphate',mod3='C3 spacer')
    else:raise ValueError('source adapter3 or5 only')
    return sc

PCR_INDEXES = {'I2':'ACATCG','I4':'TGGTCA','I5':'CACTGT','I6':'ATTGGC','I7':'GATCTG','I12':'TACAAG'}

def pcr_primer_scene(which, member=None):
    sc=Scene()
    if which=='forward':
        sc.strand('primer',[Segment('source forward primer',il.NEBNEXT_UNIVERSAL_PRIMER)])
    elif which=='reverse':
        if member is not None and member not in PCR_INDEXES:raise ValueError('actual published Table1 index member required')
        index='N'*6 if member is None else PCR_INDEXES[member]
        sc.strand('primer',[Segment('source P7 prefix',il.P7),Segment('unassigned source sample index' if member is None else 'source '+member+' index',index,placeholder=member is None,feature=feature('rg4.sample_index','sample_index','whitelist',whitelist='Published Table1 I2/I4/I5/I6/I7/I12',note='No index member is assigned to an experimental library.')),Segment('source reverse primer arm',il.TRUSEQ_READ2)])
    else:raise ValueError('forward or reverse primer required')
    return sc

def graph_for_condition(condition='K+'):
    if condition not in {'Li+','K+','K++PDS'}:raise ValueError('source conditions kept separate')
    g=ProtocolGraph('rG4-seq — source SI-qualified '+condition,authority=SOURCE,representation='published_instance')
    current=g.input('rna','Source RNA library input; exact fragmentation/end preparation blocked on main Methods',K.MOLECULAR,primary=True)
    a3=g.input('adapter3','Published Table1DNA3adapter5rApp/3C3',K.ADAPTER,state=adapter_scene('3'));lig=g.input('ligase3','Source3adapter ligation catalyst/recipe unreviewed',K.CATALYST)
    out=g.product('linked3','RNA/DNA3adapter-ligated population plus unligated substrates; actual workup unreviewed',K.MOLECULAR)
    g.react('ligate3','Published SI Fig1RNA3adapter ligation',S.LIGATION,((current,R.TEMPLATE),(a3,R.STRUCTURAL),(lig,R.CATALYST)),(out,),geometry=AdapterLigationGeometry('unavailable','unavailable',(adapter_scene('3'),),'Source adapter disclosed; exact substrate ends/workup/sealed junction assembly unreviewed.'))
    fold=g.input('fold_condition','Source separate '+condition+' in-vitro folding; concentrations/heating program blocked on main',K.CONDITION)
    folded=g.product('folded','Source condition-dependent RNA conformations; no automatic G-rich-sequence fold assignment',K.MOLECULAR)
    fold_inputs=((out,R.TEMPLATE),(fold,R.CONDITION))
    if condition=='K++PDS':
        pds=g.input('pds','Actual source pyridostatin molecular ligand; concentration/occupancy and whole interactions unassembled',K.MOLECULAR)
        fold_inputs+=((pds,R.STRUCTURAL),)
    g.react('fold','Published SI Fig1RNA folding after3adapter',S.GENERIC,fold_inputs,(folded,),changes_molecule=False,reason='Noncovalent RNA conformation changes; neither RNA sequence mutation nor guaranteed G4 formation.')
    primer_scene=Scene();primer_scene.strand('rt',[Segment('source RTprimer',RT)])
    primer=g.input('rt_primer','Published Table1reverse-transcription primer',K.OLIGO,state=primer_scene);rt=g.input('rt_enzyme','Source RTase; identity/conditions unreviewed',K.CATALYST)
    hybrid=g.product('hybrid','RNA retained alongside cDNA population with selected length/stop classes; exact purification unreviewed',K.MOLECULAR)
    g.react('rt','Published SI Fig1G4-associated RT stalling, not physical RNA cleavage',S.REVERSE_TRANSCRIPTION,((folded,R.TEMPLATE),(primer,R.STRUCTURAL),(rt,R.CATALYST)),(hybrid,))
    a5=g.input('adapter5','Published Table1DNA5adapter5phosphate/3C3/randomN6',K.ADAPTER,state=adapter_scene('5'));lig5=g.input('ligase5','Source cDNA5adapter ligation recipe unreviewed',K.CATALYST)
    linked=g.product('linked5','Source cDNA-linked5adapter products plus RNA/unligated/background material; cleanup unassembled',K.MOLECULAR)
    g.react('ligate5','Published SI Fig1cDNA5adapter ligation; name does not fix exact junction orientation',S.LIGATION,((hybrid,R.TEMPLATE),(a5,R.STRUCTURAL),(lig5,R.CATALYST)),(linked,),geometry=AdapterLigationGeometry('unavailable','unavailable',(adapter_scene('5'),),'Exact cDNA/substrate ends and source adapter junction processing unassembled.'))
    f=g.input('forward','Published Table1PCRforward primer; complete assembled landing unreviewed',K.OLIGO,state=pcr_primer_scene('forward'));r=g.input('reverse','Published Table1 indexed PCRreverse primer; unassigned source sample index, geometry unassembled',K.OLIGO,state=pcr_primer_scene('reverse'));poly=g.input('polymerase','Source PCR enzyme/conditions unreviewed',K.CATALYST)
    library=g.product('library','Source amplified library plus original/carryover template material; no assumed purification',K.LIBRARY)
    g.react('pcr','Published SI Fig1PCR',S.PCR,((linked,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(poly,R.CATALYST)),(library,),geometry=PCRGeometry.unavailable(Scene(),reason='Complete source primer/substrate/end/selection geometry unassembled.'))
    g.sequence(library,SequencingPlan(acquisition='unspecified',limitation='Published SI establishes NGS; exact platform, acquired roles/cycles and index/read primer geometry blocked on uninspected main. No single/paired-end assignment made.'))
    g.add_source_item(SourceItem('blocked_main','Original full Methods/workup/run layout','https://doi.org/10.1038/nmeth.3965'))
    return g

def protocol_graph():return graph_for_condition()
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();x=FoldedRNA('ACGGGUA','K+',3)
    c('RT selected stall copies downstream RNA suffix in complementary orientation',x.first_dna(),'TAC')
    c('RT stall does not truncate physical RNA',x.sequence,'ACGGGUA')
    c('source condition does not automatically guarantee a stall',FoldedRNA('ACGGGUA','K+').first_dna(),'TACCCGT')
    c('nonstabilizing Li is a comparison, not an enforced no-stall law',FoldedRNA('ACGGGUA','Li+',3).first_dna(),'TAC')
    c('source RTprimer complementary to complete3adapter',revcomp(RT),A3)
    c('PDS is an actual molecular contributor only on its source branch',('pds' in graph_for_condition('K++PDS').materials,'pds' in graph_for_condition('K+').materials),(True,False))
    c('Unassigned sample index is a placeholder, not synthesized N6',pcr_primer_scene('reverse').strands['primer'].segs[1].placeholder,True)
    c('all three condition graphs retain incomplete source boundary',all(graph_for_condition(v).publication_problems() for v in ('Li+','K+','K++PDS')),True)
    src=Source(sources_root()/'rg4-seq__10.1038+nmeth.3965'/'nmeth.3965_41592_2016_BFnmeth3965_MOESM356_ESM.xlsx.txt','Published Table1')
    if have(c,src,label='modified adapter source transcription'):
        t=src.path.read_text();c('Both PCR primer families transcribe actual source whole bodies',il.NEBNEXT_UNIVERSAL_PRIMER in t and all(il.P7+x+il.TRUSEQ_READ2 in t for x in PCR_INDEXES.values()),True);c('source adapters preserve5rApp/5phosphate and both3C3 blocks',all(v in t for v in ('/5rApp/'+A3+'/3SpC3/','/5Phos/'+A5+'NNNNNN/3SpC3/')),True)
    c.report()
if __name__=='__main__':selftest()

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



def qc(g,key,x,label):
    keep=g.product(key+'_retained','Sourceproductionremainder; analyticalaliquotamountnotinvented',K.MOLECULAR);test=g.product(key+'_assay','Separateanalyticalaliquotforthe sourceassay',K.QC)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(keep,test));g.mark_terminal(Terminal(test.id,Disposition.QC,'Analyticalmaterialdoesnotbecomethepreparativesequencingtemplate.'));return keep

def unknown_pcr(g,key,x,label,reagent):
    f=g.input(key+'_f','Sourceforwardindexingprimerfamily; memberandfullbodyunreviewed',K.OLIGO);r=g.input(key+'_r','Sourcereverseindexingprimerfamily; memberandfullbodyunreviewed',K.OLIGO);cat=g.input(key+'_enzyme',reagent,K.CATALYST);y=g.product(key+'_product','SourcePCRcopiesplusoriginaltemplates/unextendedprimersandreactionremnants; exactbackbonesunassembled',K.MOLECULAR)
    g.react(key,label,S.PCR,((x,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(cat,R.CATALYST)),(y,),geometry=PCRGeometry.unavailable(scene('Sourcefulladapter/primer/endgeometryunassembled'),reason='KnownPCRoperationdoesnotestablishcompletefirstcopyorfinalstrandgeometry.'));return y


from chemdraw import revcomp,feature
from illumina import P5,P7,TRUSEQ_READ1,TRUSEQ_READ2,INDEX1_PRIMER
SOURCE='Published2017PMC5696378Methods90–106;publishedMOESM1standardHaeIII/two-stepligationprotocol734–957,originaloligosp10–11/fixationp12rendered. Wholejunction/readgeometryunassembled.'
@dataclass(frozen=True)
class BiotinBridge:
    sequence:str='CGCGATATCTTATCTGACT'
    biotin_t:int=9
    five_prime:str='P'
    def __post_init__(self):
        if self.sequence!='CGCGATATCTTATCTGACT' or self.biotin_t!=9 or self.sequence[self.biotin_t]!='T' or self.five_prime!='P':raise ValueError('SelectedsourceS2-Factualinternalbiotin-dTand5Pmustbepreserved')
    def partner(self):return 'GTCAGATAAGATATCGCGT'
    def scene(self):
        sc=Scene();sc.strand('s2_f',[Segment('f_core',self.sequence[:-1]),Segment('f_t',self.sequence[-1])],mod5='P',label='S2-F; internalbiotin-dT10 retained as nucleotide')
        sc.anneal('s2_r',[Segment('r_core',self.partner()[:-1]),Segment('r_t',self.partner()[-1])],to='s2_f',pair=('r_core','f_core'),mod5='P');return sc

def oligo_state(name):
    seq={'y_f':INDEX1_PRIMER,'y_r':'T'+TRUSEQ_READ1,'universal':P5+TRUSEQ_READ1[4:20],'ntr_f':P5[:20],'ntr_r':P7[:20]}
    sc=Scene()
    if name=='index':sc.strand('index',[Segment('p7',P7),Segment('source_window','CGTGATGT',feature=feature('bl_source_index','sample_index','fixed')),Segment('tail','GACTGGAGTTCAGACGTGT')]);return sc
    if name not in seq:raise ValueError('Actualsourceoligorequired')
    sc.strand(name,[Segment(name,seq[name])],mod5='P' if name=='y_f' else '');return sc

def split(g,key,x,label,storage=False):
    a=g.product(key+'_used','Source selected production aliquot',K.MOLECULAR);b=g.product(key+'_remainder','Independent source remainder',K.MOLECULAR);g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(a,b));g.mark_terminal(Terminal(b.id,Disposition.STORED if storage else Disposition.UNRESOLVED,'Source explicitly −20Cstored bead aliquots.' if storage else 'Unassigned fraction fate, no implicit discard.'));return a

def protocol_graph():
    g=ProtocolGraph('BL-Hi-C — source standard HaeIII/two-stepS2 route',authority=SOURCE,representation='published_instance')
    x=g.input('cells','Source cells; native genomic sequence/unbounded contact partners not assigned',K.MOLECULAR,primary=True)
    x=op(g,'fix',x,'SourceFA10minRT/glycine10minRT+5minice; numerical formulation unresolved',reagent='SIp12 prints1mlmedium+278ul37%FA(final1%)+894ul2.5Mglycine(final0.2M); do not repair to10ml',product='Fixed cellular material and reaction companions; exact numeric fixation composition unassembled')
    x=cleanup(g,'cell_wash',x,'Source300g5min4C/PBSwash2')
    x=cleanup(g,'lysis',x,'Source0.1%SDSlysis/pellet800g twice,1%SDS/pellet3000g,0.1%SDSwash; actual pellet retained')
    x=op(g,'open',x,'Source0.5%SDS62C10min/Triton37C15min',product='Source nuclear preparation with companions, no DNA removal assumed')
    x=op(g,'haeiii',x,'SourceHaeIII100U37C2–16h',S.DIGESTION,reagent='SourceHaeIII/buffer',catalyst=True,product='Source blunt genomic fragments plus uncut/crosslinked DNA and enzyme retained')
    x=op(g,'a_tail',x,'SourcedATP/Klenowexo−37C40min',S.TAILING,reagent='SourcedATP/Klenowexo−',catalyst=True,product='Selected3A-tailed genomic ends plus unprocessed blunt DNA and companions')
    bridge=BiotinBridge();linker=g.input('s2_bridge','SourceS2 duplex;both5P/internalbiotin-dT;19ntperbackbone/18pairedcore/terminal3T;notbarcode',K.ADAPTER,state=bridge.scene(),chemical_state=bridge);cat=g.input('s2_ligase','SourceT4DNA ligase',K.CATALYST);y=g.product('ligated','Selectedgenomic–S2–genomic contact molecules plus one-sided/freeS2,uncut/otherDNA and ligase; completejunctionunassembled',K.MOLECULAR)
    g.react('s2_ligate','SourcecontinuousS2/T4ligation16C4h;kinetictwo-joins not twoinventedreactions',S.LIGATION,((x,R.TEMPLATE),(linker,R.STRUCTURAL),(cat,R.CATALYST)),(y,),geometry=uncertain_ligation())
    x=cleanup(g,'postligation_pellet',y,'Source3500g5min4C pellet retained/remove supernatant')
    x=op(g,'exo',x,'Sourcelambdaexonuclease/ExoI37C1h',S.DIGESTION,reagent='Sourcelambdaexonuclease/ExoI',catalyst=True,product='ContactDNA plus unprocessedDNA and free-linker digestion debris; no debris disappears before extraction')
    x=op(g,'proteinase',x,'SourcePK/SDS55Catleast2h;NaCl68C2h',S.DEPROTEINIZATION,reagent='SourceProteinaseK/SDS/NaCl',catalyst=True,product='DNA and proteolytic fragments/crosslink remnants until physical extraction')
    x=cleanup(g,'extract',x,'SourcePCI/MaXtractaqueous→GlycoBlue/NaOAc/isopropanol−80C1h→pellet/75%EtOHwash2→130ulEB')
    x=op(g,'shear',x,'SourceCovaris400bptarget',S.FRAGMENTATION,product='Biotin-linker-bearing and untagged genomic fragments; no immediate capture selection')
    beads=g.input('beads','SourceM280streptavidin beads',K.MOLECULAR);beads=cleanup(g,'bead_wash',beads,'SourceB&Wwash2/I-Block45minRT/B&Wwash2')
    block=g.input('salmon_dna','Source1ugSalmonSpermDNA blocking stock; independent DNA remains possible carryover',K.MOLECULAR);b=g.product('blocked','Blocked beads plus blocking DNA; no covalent genomic linkage asserted',K.MOLECULAR);g.react('block','SourceSalmonSpermDNA30minRT',S.BINDING,((beads,R.TEMPLATE),(block,R.STRUCTURAL)),(b,));b=cleanup(g,'block_wash',b,'SourceB&Wwash2 after DNA blocking; recovery/carryover not quantified')
    mix=g.product('capture_mix','BeadboundbiotinDNA and freeDNA/blocker/other companions',K.MOLECULAR);g.react('capture','Source45minRTstreptavidin capture',S.BINDING,((x,R.TEMPLATE),(b,R.STRUCTURAL)),(mix,));x=cleanup(g,'capture_wash',mix,'Source5xSSC/SDSwash,2xB&W,EB; beadassociatedDNA retained, blocker recovery unassigned')
    x=op(g,'repair',x,'SourcebeadboundPNK/T4pol/Klenow37C30min',S.END_REPAIR,reagent='SourcePNK/T4pol/Klenow',catalyst=True,product='Selected repaired genomic ends plus possible blockingDNA/enzymes/remnants')
    x=cleanup(g,'repair_wash',x,'Source2xTWB55C2min')
    x=op(g,'second_a',x,'SourceKlenowexo−Atailing37C30min',S.TAILING,reagent='SourceKlenowexo−/dATP',catalyst=True,product='Selected3Aends plus unprocessedDNA/possibleblockingDNA and remnants')
    x=cleanup(g,'a_wash',x,'Source2xTWB/Quickligationbufferwash')
    yf=g.input('y_f','Source5P Y-adapter-F',K.OLIGO,state=oligo_state('y_f'));yr=g.input('y_r','SourceY-adapter-R;5P not asserted',K.OLIGO,state=oligo_state('y_r'));ya=g.product('y_adapter','SourceannealedYadapter; complete ligating junction geometry unassembled',K.ADAPTER);temp=g.input('anneal_temperature','SourceTNE/95C2min then programmed cooling',K.CONDITION);g.react('anneal_y','SourceYadapter annealing',S.HYBRIDIZATION,((yf,R.TEMPLATE),(yr,R.STRUCTURAL),(temp,R.CONDITION)),(ya,))
    lig=g.input('quick_ligase','SourceQuickDNA ligase',K.CATALYST);y=g.product('adapter_ligated','SourceYadaptedselectedDNA plus possible blockerDNA/freeadapters/reagents untilwash',K.MOLECULAR);g.react('y_ligate','SourceQuickligationRT15min',S.LIGATION,((x,R.TEMPLATE),(ya,R.STRUCTURAL),(lig,R.CATALYST)),(y,),geometry=uncertain_ligation());x=cleanup(g,'ligate_wash',y,'Source2xTWB/EB')
    x=split(g,'pcr_allocate',x,'Source60ulEBdivided20ulaliquots;selectone20ul,othersstored−20C',True)
    inputs=[(x,R.TEMPLATE)]
    for name in ('universal','index','ntr_f','ntr_r'):
        primer=g.input('primer_'+name,'Sourcefour-primerQ5 '+name,K.OLIGO,state=oligo_state(name));inputs.append((primer,R.STRUCTURAL))
    cat=g.input('q5','SourceQ5HotStartHigh-FidelityDNApolymerase0.5ul; separateQ5buffer/dNTP sourceconditions',K.CATALYST);inputs.append((cat,R.CATALYST));y=g.product('pcr_product','DNAcopies plus originalDNA/possibleblockingDNA,freeprimers and beads; exact firstcopy/readgeometry unassembled',K.MOLECULAR);g.react('pcr','SourceQ5four-primerPCR9–12cycles',S.PCR,inputs,(y,),geometry=PCRGeometry.unavailable(scene('Fullsourcegenomic/Y/PCRgeometryunassembled'),reason='Source-known oligos do not establish complete library/backbone geometry.'))
    x=cleanup(g,'remove_dynabeads',y,'SourcePCRproduct transferred withoutDynabeads')
    x=cleanup(g,'ampure_large',x,'Source0.6xAMPure: retainSUPERNATANT/removebeads')
    x=cleanup(g,'ampure_small',x,'Sourceadditional7.5ul/cumulative0.75xAMPure: retainBEADS/removeSUPERNATANT;80%EtOHwash2/elute20ulwater')
    x=qc(g,'qc',x,'SourceQubit/Bioanalyzer separate analytical aliquots')
    lib=g.product('library','SourceBL-Hi-C DNAsequencinglibrary; fullstrand/primer/readgeometryunassembled',K.LIBRARY);g.react('identity','Sourcephysicalfinal library',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='SourceHiSeqXTenPE150 acquisition. Structured literal8nt /CGTGATGT/ sourcewindow retained, no assumed6nt correction or validated index-readgeometry. Contact mapping computational; no cellbarcode assumed.'))
    g.add_source_item(SourceItem('pending','Fixationnumericformulation/fulljunction/end/readgeometry and unreviewed sourcealternatives','PublishedSIp12values internally incompatible aswritten; no silentvolume repair. Sourceone-stepS1/overhangenzymevariants outside selectedHaeIII/two-stepS2 route; privateincomplete.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();a=BiotinBridge();c('InternalbiotinoccupiesactualT',(a.sequence[a.biotin_t],len(a.sequence),len(a.partner())),('T',19,19));a.scene()
    try:BiotinBridge(biotin_t=8);bad=False
    except ValueError:bad=True
    c('BiotinmarkcannotmigrateoffsourceT',bad,True)
    g=protocol_graph();c('Linkerdebrisretaineduntilphysicalextraction',g.steps['proteinase'].inputs[0].material,'exo_product')
    c('BlockerDNAhasindependentphysicalinput',g.steps['block'].inputs[1].material,'salmon_dna')
    c('PCRusesonlyallocatedbeadaliquot',g.steps['pcr'].inputs[0].material,'pcr_allocate_used')
    c('Actualaliquotremainderstored',g.terminals['pcr_allocate_remainder'].disposition,Disposition.STORED)
    c('FirstAMPurestageretainssupernatant','retainSUPERNATANT' in g.steps['ampure_large'].label,True)
    c('Analyticalmaterialdoesnotfeedlibrary',g.steps['identity'].inputs[0].material,'qc_retained')
    c('Unreviewedgeometrypreventspromotion',bool(g.publication_problems()),True)
    src=Source(sources_root()/'bl-hi-c__10.1038+s41467-017-01754-3'/'s41467-017-01754-3_supp_41467_2017_1754_MOESM1_ESM.pdf.txt','get_sources --slugs bl-hi-c__10.1038+s41467-017-01754-3')
    if have(c,src,label='publishedSIoligos'):
        text=src.path.read_text();c('SourceS2fullmodifiedstrand','CGCGATATC/iBIOdT/TATCTGACT' in text,True);c('Sourceindexwindowpreserved','/CGTGATGT/' in text,True)
    c.report()
if __name__=='__main__':selftest()

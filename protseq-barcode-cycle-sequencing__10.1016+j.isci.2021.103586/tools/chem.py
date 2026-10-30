from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Scene, Segment, feature, workflow_from_protocol_graph, workflow_panel
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan
from page import head,info
PAGE_FILE='chem.html'

def operations(g):
    def op(key,x,label,reagent=None,state=None,kind=S.GENERIC,reagent_kind=K.CONDITION,reagent_role=R.CONDITION):
        inputs=[(x,R.TEMPLATE)]
        if reagent:
            z=g.input(key+'_reagent',reagent,reagent_kind);inputs.append((z,reagent_role))
        y=g.product(key+'_product',label+'; preparation includes unremoved reaction components',K.MOLECULAR,state=state)
        g.react(key,label,kind,inputs,(y,));return y
    def clean(key,x,label):
        z=g.input(key+'_selection',label,K.SELECTION)
        y=g.product(key+'_kept','Physically retained source fraction: '+label,K.MOLECULAR)
        w=g.product(key+'_removed','Physically removed fractions: '+label,K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(z,R.SELECTION)),(y,w));g.discard(w,'Source physical workup only; no perfect specificity or purity assertion.');return y
    return op,clean

def end(g,x,label):g.sequence(x,SequencingPlan(acquisition='unspecified',limitation=label))
def incomplete(g):
    g.declare_reviewed_scope(tuple(g.materials)+tuple(g.steps),SOURCE)
    g.add_source_item(SourceItem('unassembled','Unreviewed linked originals and full molecular geometry',LIMIT));return g


from dataclasses import dataclass,replace
from protocol_graph import Terminal,Disposition

def join(g,key,inputs,label):
    y=g.product(key+'_product',label+';unremoved companion species remain',K.MOLECULAR)
    g.react(key,label,S.GENERIC,tuple((x,R.TEMPLATE if i==0 else R.STRUCTURAL) for i,x in enumerate(inputs)),(y,));return y

def render_page():
    return head(TITLE)+'<main class="wrap"><h1>'+TITLE+'</h1>'+info(LIMIT)+'<p><a href="01_chem.html">Research notes</a></p>'+workflow_panel(workflow_from_protocol_graph(protocol_graph(),ENDPOINT))+'</main>'

from illumina import P5_RC
TITLE='ProtSeq — demonstrated chip barcode-cycle sequencing'
SOURCE='Published PMC8717419 Results262–278,STARMethods588–609/643–657/669–673;publishedSI TableS2/S4 andFigS3'
LIMIT='Private incomplete model of demonstrated DNA-DNA barcode-cycle NGS. This is not complete amino-acid/protein sequencing. Peptide–nanobody capture and TargetSwitchSELEX are separate demonstrations. Source-selected populations,not perfect cognate binding or six successful ligations per molecule. Complete foundation/cololinker/binder junction geometry remains unassembled.'
ENDPOINT='ready_product'
FOUNDATION='CGACTGCGAGCTGATGGCCTTGATGATAACG'
NGS_PREFIX='AGATCGGAAGAGCGTCGTGTAGGGAAAGA'
NGS_ADAPTER=NGS_PREFIX+P5_RC
BRIDGE_PREFIX='CTGCGCCTATAGGAATTCGTTATC';BRIDGE_SUFFIX='GGACACGGCCGTTATC'
@dataclass(frozen=True)
class UniversalBridge:
    prefix:str=BRIDGE_PREFIX
    universal_residues:tuple=('5-Nitroindole',)*12
    suffix:str=BRIDGE_SUFFIX
    def notation(self):return self.prefix+'/i5NitInd/'*len(self.universal_residues)+self.suffix
    def __post_init__(self):
        if self.universal_residues!=('5-Nitroindole',)*12:raise ValueError('Source12actual5NitIndresidues,not12randomDNA or inosines')
@dataclass(frozen=True)
class ObservedChain:
    recorded:tuple=()
    def record_observed_event(self,binder,cycle):
        if not binder or type(cycle)is not int or cycle<1 or (self.recorded and cycle<=self.recorded[-1][1]):raise ValueError('Explicit increasing source observation needed')
        return replace(self,recorded=self.recorded+((binder,cycle),))
def foundation():
    sc=Scene();sc.strand('foundation',[Segment('common source foundationbody','CGACTGCGAGCTGATG'),Segment('sourceFd7identity','GCCTTGAT',feature=feature('protseq.fd7','feature_barcode','fixed',note='Construct/target foundation identity,not cell or UMI.')),Segment('source foundationtail','GATAACG')],mod5='phosphate',label='Source TableS2 Fd7 example;exact experimental target assignment unassembled');return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE);op,clean=operations(g)
    targets=g.input('targets','Source DNA target population for demonstrated1/6cycle experiment;source target–foundation pairing unassigned',K.OLIGO,primary=True)
    fd=g.input('foundation','PublishedFd7 example molecularinventory;assignmenttoselectedtarget unassembled',K.OLIGO,state=foundation())
    fc=g.input('forward_cololinker','SourceFC exactwhole polymer unassembled',K.OLIGO);rc=g.input('reverse_cololinker','SourceRC exactwhole polymer unassembled',K.OLIGO)
    x=join(g,'hybridize',(targets,fd,fc,rc),'Source sequencing-unit hybridization95→25°C stagedcool;FC10nM,othercomponents excess')
    chip=g.input('chip','SourceMiSeq flowcell P5/P7 surface-bound oligos;historical exactgeometryunreviewed',K.MOLECULAR)
    x=join(g,'seed',(x,chip),'Source foundationmix120pM/T4BluntTA;two30µL loads,28°C15min;selectedsource tethered classes')
    x=clean('seed_wash',x,'Source100%formamidewash40°C90s;retain chip-associated classes')
    x=op('block',x,'SourceBlockingBuffer/ChipBlockingSolution,37°C15min;surface occupancynotpredicted','SourceP5/P7complements/POC-tailblocking DNA andbuffer;whole exactinputs unassembled')
    bridge=g.input('universal_bridge','Source12actual5-Nitroindole universal-base residues between fixedDNAflanks;typed orderedchemical inventory,not N12/inosine;whole visual pairing unassembled',K.OLIGO,state=UniversalBridge())
    for cycle in range(1,7):
        name='cycle'+str(cycle)
        binder=g.input(name+'_binder','Source round-specific DNA-binder/barcode population;exactwhole polyermer/junction unassembled',K.OLIGO)
        mix=join(g,name+'_bridge',(binder,bridge),'Source binder:bridge1:2 hybridization95°C5min,RT1h;BSA100µg/mL')
        x=join(g,name+'_bind',(x,mix),'SourceDNA-DNA binding incubationRT30min;no perfect specificity or incorporation predicted')
        x=clean(name+'_wash',x,'Source washout ofunbound binders;chip-associated populationretained')
        x=op(name+'_ligate',x,'SourceT4ligase joins eligiblebinderbarcode tofoundation/chain;newbond retainsoriginalDNA bodies','Source2xT4ligase/bridge geometry unassembled')
        x=op(name+'_digest',x,'SourceEcoRI+complementaryrestrictionbridge cleavesbinder region frombarcode;cleavagefragments remainuntilwash','EcoRI/CutSmart+sourcebridge;no peptide/aminoacid degradation assigned')
        x=clean(name+'_release_wash',x,'Sourceformamide/wash removesreleasedDNA/binderclasses;retain chip-linkedbarcodechains without predicting yield')
    ad=g.input('ngs_adapter','Published5′phosphorylated universalNGSadapter',K.OLIGO,state=adapter_scene())
    abr=g.input('ngs_bridge','Published universalNGSadapterbridge9/5',K.OLIGO,state=Construct([Segment('source bridge','TTCCGATCTCGTTA')]))
    x=join(g,'ngs_ligate',(x,ad,abr),'Source2:1adapter:bridge,BluntTAligation;two30µLloads40°C2min45s;retainDNA bodies')
    x=clean('adapter_wash',x,'Two500µL nucleasefreewaterwashes,90s apart;retainchipDNA')
    x=op('ready',x,'Source chip goes directlytoMiSeq;adaptedrun omitsinitialwashing/libraryloading;PhiXcontrol cartridge input kept distinct')
    end(g,x,'PublishedMiSeqV2chemistry,45–600cycle read;experiment-specific cyclecount/primer geometryunassembled. Sourceonchip bridge-ligation route,not separateSELEX177bpPippinPCRprep.')
    return incomplete(g)
def adapter_scene():
    sc=Scene();sc.strand('adapter',[Segment('sourcebody',NGS_PREFIX),Segment('sourcecanonicalP5reversecomplement',P5_RC)],mod5='phosphate');return sc

def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();g=protocol_graph();n=ObservedChain().record_observed_event('sourcebinder',1).record_observed_event('sourceother',3)
    c('Observed missingcycle retainedwithout invented incorporation',tuple(z[1] for z in n.recorded),(1,3))
    c('5NitInd is actualmodifiedresidue,not DNA N/I',set(UniversalBridge().universal_residues),{'5-Nitroindole'})
    c('Foundation retains5phosphate',foundation().strands['foundation'].mod5,'phosphate')
    c('Foundation identity is not UMI',foundation().strands['foundation'].segs[1].feature.role,'feature_barcode')
    c('Nextcycle receiveswashedretainedchip material',g.steps['cycle2_bind'].inputs[0].material,'cycle1_release_wash_kept')
    c('NGSadapter aftersixsourcecycles',g.steps['ngs_ligate'].inputs[0].material,'cycle6_release_wash_kept')
    c('NGSadapter physical5phosphate',adapter_scene().strands['adapter'].mod5,'phosphate')
    c('No aminoaciddegradation step in demonstrated DNA branch','aminoacid_degrade' in g.steps,False)
    c('Incomplete sourcegeometry staysprivate',g.publishable,False)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'j.isci.2021.103586_PMC8717419.xml.txt'
    if p.exists():
        t=p.read_text();c('Wholemodifieduniversalbridge source transcription',UniversalBridge().notation() in t,True);c('WholeNGSadapter verbatim',NGS_ADAPTER in t,True)
    else:c.skip('Published originaltranscription','Externaloriginalabsent')
    c.report()

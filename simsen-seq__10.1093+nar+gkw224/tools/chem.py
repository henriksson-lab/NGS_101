from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from dataclasses import dataclass,replace
from chemdraw import Scene,Segment,Construct,feature,workflow_from_protocol_graph,workflow_panel
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from page import head,info
PAGE_FILE='chem.html'

class Operations:
    def __init__(self,g):self.g=g
    def op(self,key,x,label,reagent=None,state=None,chemical_state=None):
        g=self.g;inputs=[(x,R.TEMPLATE)]
        if reagent:
            y=g.input(key+'_condition',reagent,K.CONDITION);inputs.append((y,R.CONDITION))
        y=g.product(key+'_out',label+'; preparation includes unremoved companion species',K.MOLECULAR,state=state,chemical_state=chemical_state)
        g.react(key,label,S.GENERIC,inputs,(y,));return y
    def join(self,key,inputs,label,state=None,chemical_state=None):
        g=self.g;y=g.product(key+'_out',label+'; reaction components remain until an explicit separation',K.MOLECULAR,state=state,chemical_state=chemical_state)
        g.react(key,label,S.GENERIC,tuple((x,R.TEMPLATE if i==0 else R.STRUCTURAL) for i,x in enumerate(inputs)),(y,));return y
    def clean(self,key,x,label):
        g=self.g;z=g.input(key+'_selection',label,K.SELECTION)
        y=g.product(key+'_kept','Source physically retained fraction: '+label,K.MOLECULAR)
        w=g.product(key+'_removed','Source physically removed fractions; exact composition/purity unassigned',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(z,R.SELECTION)),(y,w));g.discard(w,'Only the source physical separation; no perfect selectivity asserted.');return y
    def split(self,key,x,label,other_label,other=Disposition.UNRESOLVED):
        g=self.g;y=g.product(key+'_production',label,K.MOLECULAR)
        z=g.product(key+'_other',other_label,K.MOLECULAR)
        g.react(key,'Source physical allocation: '+label,S.SPLIT,((x,R.TEMPLATE),),(y,z));g.terminals[z.id]=Terminal(z.id,other,other_label);return y,z
    def qc(self,key,x,label):
        g=self.g;y=g.product(key+'_production','Production preparation retained independently from analytical aliquot',K.MOLECULAR)
        z=g.product(key+'_assay',label,K.QC)
        g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,z));g.terminals[z.id]=Terminal(z.id,Disposition.QC,'Analytical aliquot does not become the sequencing substrate.');return y
    def opaque(self,key,x,label):
        g=self.g;y=g.product(key+'_out',label,K.MOLECULAR)
        g.react(key,label,S.MAGIC,((x,R.TEMPLATE),),(y,),reason='Selected originals establish this operational continuation; full molecular inventory/geometry remains unassembled.');return y

def finish(g,x,readout):
    g.sequence(x,SequencingPlan(acquisition='unspecified',limitation=readout))
    g.declare_reviewed_scope(tuple(g.materials)+tuple(g.steps),SOURCE)
    g.add_source_item(SourceItem('unassembled','Unassembled source inputs, molecular geometry and directly linked targets',LIMIT));return g

def render_page():
    return head(TITLE)+'<main class="wrap"><h1>'+TITLE+'</h1>'+info(LIMIT)+'<p><a href="01_chem.html">Research notes</a></p>'+workflow_panel(workflow_from_protocol_graph(protocol_graph(),ENDPOINT))+'</main>'


from chemdraw import revcomp
from hairpin import Hairpin
from checks import Source,have
from storage import sources_root
import illumina as il
TITLE='SiMSen-Seq — source hairpin-protected molecular barcoding'
SOURCE='PublishedPMC4914102 Methods65–75/Results79–93/Fig1; TableS1File008workbook universalprimers andTP2row,File006SIcaptions'
LIMIT='Private incomplete source model. Source N12 is a molecular identifier and six-basePCRindex identifies samples. Hairpin stem complementarity is computed,not thermal folding probability or perfectmisprimingprotection. Maintext74.0±0.3°C versusFig1caption70.01±0.24°C hairpinTm remains a visible source discrepancy; no value chosen for a calculation. StandardGGprimerinventoryisnotalltargetprimers:TableS1SMA1/SMA4lackGGandSMA2hasGGG;fullvariantorders/allocationunassembled. Fulltarget mixtures/PCR products/readlanding andoriginalkitdetails unassembled.'
ENDPOINT='pool_out'
STEM='GGACACTCTTTCCC'
ADAPTER='GG'+il.TRUSEQ_READ1
TAIL=revcomp(STEM)
UMI=feature('simsen.molecule','umi','random',note='SourceN12molecularbarcode;consensusfamiliesarecomputational,notphysicalone-to-oneyieldguarantee')
INDEX=feature('simsen.sample','sample_index','unknown',note='Source6ntPCRsampleindex;memberunassigned')
@dataclass(frozen=True)
class StandardProtectedPrimer:
    target:str
    barcode:str='N'*12
    def __post_init__(self):
        if not self.target or set(self.target)-set('ACGT') or len(self.barcode)!=12 or set(self.barcode)-set('ACGTN'):raise ValueError('SourceDNA targetarm and12base molecularbarcode required')
        Hairpin('sourcehairpinsequence',self.sequence,14,len(self.target))
    @property
    def sequence(self):return ADAPTER+self.barcode+'AT'+TAIL+self.target
    def hairpin(self):return Hairpin('sourceconditionalstem',self.sequence,14,len(self.target))

def primer_inventory():
    sc=Scene();sc.strand('source',[Segment('sourceGGstabilizer','GG'),Segment('sourceadapterbody',il.TRUSEQ_READ1),Segment('sourceN12moleculebarcode','N'*12,placeholder=True,feature=UMI),Segment('forcedmismatchdestabilizer','AT'),Segment('complementarystem',TAIL),Segment('unassembledtargetarm','[source target primer]',placeholder=True)],label='SourceTableS1standardGGforwardcomponentexample;otherstabilizervariants/wholemanufacturedtargetoligounassembled')
    return sc

def sample_primer():
    return Construct([Segment('P7',il.P7),Segment('source6ntsampleindex','N'*6,placeholder=True,feature=INDEX),Segment('sourcePCRarm','GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT')])

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE);o=Operations(g)
    x=g.input('dna','Source5–100nggDNAorcfDNA;DNAextraction/shearing/cell-versusplasmaapplicationnotforcedonepreparation',K.MOLECULAR,primary=True)
    f=g.input('forward','Sourcehairpin-protectedforwardprimerpool;shownstandardGGcomponentexampleonly,SMAnoGG/GGGvariantsandactualpanelallocationunassembled',K.OLIGO,state=primer_inventory())
    r=g.input('reverse','Sourceuniversalreverse-targetprimer;fullspecifictargetbodyunassembled',K.OLIGO)
    x=o.join('barcode_pcr',(x,f,r),'AccuPrime3cyclebarcoding,98°C10s/62°C6min/72°C30s after98°C3min;newDNAcopiescarrysourcebarcodewhileoriginaltemplates/reagentsremain')
    x=o.op('dilute_protease',x,'Add20µLTE/protease30ngperµL;65°C15mindigestsTaq,not physicalDNAcleanup','TE/protease')
    x=o.op('inactivate',x,'95°C15mininactivatesprotease; fragments/unusedprimers/barcodedcopiesremain')
    x,rest=o.split('pcr2_allocate',x,'Source10µLbarcodingreactiontransferinto40µLPCR2','Otherbarcodingreactionaftertransfer;laterfateunassigned')
    f=g.input('illumina_f','PublishedIlluminaforwardPCRprimer',K.OLIGO,state=Construct([Segment('sourceprimer',il.NEBNEXT_UNIVERSAL_PRIMER)]))
    r=g.input('illumina_r','PublishedP7/six-baseindex/universalreversearmprimer;memberunassigned',K.OLIGO,state=sample_primer())
    x=o.join('library_pcr',(x,f,r),'Q5HotStartPCR18–30cyclesaddscompleteIlluminalibraryarms;source80→72→76°Cslowramp,notidenticalbarcodingcycleprogram')
    x,rest=o.split('clean_allocate',x,'Source36µLPCR2productforpurification','OtherPCR2materialnotforcedintocleanup')
    x=o.clean('ampure',x,'Source0.83–1.0xAMPuredependingampliconlength;20µLTEelution')
    x=o.qc('fragment_analyzer',x,'SeparateFragmentAnalyzeranalyticalsizeassessment')
    x=o.op('normalize',x,'Sourceconcentrationnormalizationbeforemultiplexsequencing;exactworkupunassembled')
    peers=g.input('other_samples','Othersourceindependentlypreparednormalizedlibraries;batchmembershipunassigned',K.MOLECULAR)
    x=o.join('pool',(x,peers),'Sourcepoollibrariesforsequencing;noone-to-onerawread/moleculeguarantee')
    return finish(g,x,'Publishedlow-depthMiSeqNanoV2SE150versusfullMiSeq/HiSeq2500SEorPE150;actualcohortallocation/fullprimerlandingunassembled. N12molecularfamilies/consensuscallingnotphysicalDNAdeletionorcorrection.')

def selftest():
    from checks import Check
    c=Check();p=StandardProtectedPrimer('GACCTCAGGCGGCTCATAG');h=p.hairpin();g=protocol_graph()
    c('Sourcehairpinrightstemcomputedfromleft',h.right_stem,revcomp(h.left_stem))
    c('Sourcehairpintargetarmremainsoutsidefold',h.overhang_3,p.target)
    c('BarcodeandATretainedinhairpinloop',h.loop.endswith('N'*12+'AT'),True)
    c.raises('Wrongbarcodewidthrejected',lambda:StandardProtectedPrimer('ACGT','N'*11))
    c('ProteasedigestiondoesnotdiscardPCR1DNA',g.steps['inactivate'].inputs[0].material,'dilute_protease_out')
    c('PCR2receivesactualtransferredfraction',g.steps['library_pcr'].inputs[0].material,'pcr2_allocate_production')
    c('PCR2leftovernotdiscardedbyassumption',g.terminals['clean_allocate_other'].disposition,Disposition.UNRESOLVED)
    c('Pooldoesnotreceiveanalyticalsizealiquot',g.steps['normalize'].inputs[0].material,'fragment_analyzer_production')
    src=Source(sources_root()/'simsen-seq__10.1093+nar+gkw224'/'nar_gkw224_supp_gkw224_Supplementary_Data_nar-03291-met-g-2015-File008.xlsx.txt','PublishedTableS1universalforward/hairpin/targetTP2rows')
    if have(c,src,label='Sourcehairpincomponenttranscription'):
        t=src.path.read_text();c('Wholegenericforwardbodymatchesoriginal',ADAPTER+'N'*12+'AT'+TAIL+'-fwd target primer' in t,True)
        c('Selectedtargetarmoccursinsourcerow',p.target in next(x for x in t.splitlines() if x.startswith('TP_2\t')),True)
    c('Incompletegeometrynopublicpromotion',g.publishable,False);c.report()

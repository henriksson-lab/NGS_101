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


import re
from chemdraw import revcomp
from storage import sources_root
from checks import Source,have
TITLE='INSPIIRED — source integration-site library preparation'
SOURCE='PublishedPMC5363316 Methods177–185; publishedProtocolSI mmc1 steps1–10/TablesS9–12 and workbookmmc2 selectediSL1row'
LIMIT='Private incomplete source graph. This is the wet-preparation paper .002,not companionanalysis .003. iSL1 linker/blocker are source inventory examples,not one assigned clinical sample. N12random linker region is not automatically a UMI;12ntGolaybarcode identifiesreplicate. MainMiSeq cycle allocation179+12+143 differs from the300cyclekit label and is retained as printed. SI elutesligationin15µL but lists4PCR1 reactions eachusing15µL; actualscaling/allocationunassembled,not silentlyrepaired. Fullend/primer/read geometry and originalkitrevisions remain unreviewed.'
ENDPOINT='screen_qpcr_production'
COMMON='CTCCGCTTAAGGGACT'
PREFIX='GCCGAACATGAGCTATACGACTTAAACTTCCGCATGGCGT'
BLOCK='(bC)A(bG)(bG)GACTTGAAA(bG)(bC)GAAA(bG)(bG)(bG)AAA(bC)C'
SOURCE_OLIGO_EXCEPTIONS={
    'rep'+str(i)+'_block'+str(j):'PublishedProtocolSI j.omtm.2016.11.002_supp_mmc1.pdf originalp6TableS9,wholeHIVrow: '+BLOCK+'-3primeAmMO. Coordinatorandindependentreaderinspectedoriginal /tmp/inspiired19-blocker.png andconfirmedallbases/ninebridgedpositions;bC=iBNA-meC,bG=iBNA-G. Existingparserdoesnotrecognize(bN)sourceformat. ExactBLOCKsource-transcriptiontestremains;actualstrandmodifiedpositionsand3AmMOretained. Onlythisone-strandblockinginputexempt,notgeometry/absence/otheroligos;nosyntheticsourcetwin.'
    for i in range(1,5) for j in (1,2)
}
@dataclass(frozen=True)
class BlockingOligo:
    as_written:str=BLOCK
    terminal:str='3AmMO'
    def parts(self):
        tokens=re.findall(r'\(b[ACGT]\)|[ACGT]',self.as_written)
        if ''.join(tokens)!=self.as_written:raise ValueError('Preserve complete sourceBNA notation')
        return [Segment('iBNA-meC' if x=='(bC)' else ('iBNA-'+x[2] if x.startswith('(b') else 'DNA'),x[2] if x.startswith('(b') else x,note='PublishedS9 '+x+';bridgedsugar modification,not copied toPCRdaughters') for x in tokens]
    def scene(self):
        sc=Scene();sc.strand('blocker',self.parts(),label='PublishedS9HIVblockingexample;actualmodifiedbasespreserved',mod3=self.terminal);return sc
    @property
    def modified_count(self):return self.as_written.count('(b')

def linker():
    sc=Scene();sc.strand('positive',[Segment('iSL1sample-specificfixedregion',PREFIX,feature=feature('inspiired.linker.sample','sample_index','fixed',note='InventoryiSL1;actualsampleassignmentunreviewed')),Segment('SourceN12randomregion','N'*12,placeholder=True,note='Publishedcallsrandomsequence;identifierfunctionunassigned,notUMIbydefault'),Segment('annealedcommoncore',COMMON[:-1]),Segment('3primedToverhang',COMMON[-1])])
    sc.anneal('negative',[Segment('commoncorecomplement',revcomp(COMMON[:-1]))],to='positive',pair=('commoncorecomplement','annealedcommoncore'),mod5='P',mod3='AmMO');return sc

def protocol_graph():
    from adapter_ligation import AdapterLigationGeometry
    g=ProtocolGraph(TITLE,authority=SOURCE);o=Operations(g)
    dna=g.input('genomic','ConditionalsourceHIV-vectorintegrationgDNAentry;notassignedclinicalsample,normalgDNAcontrolseparate,notallDNAhasintegration',K.MOLECULAR,primary=True)
    x=o.op('sonicate',dna,'CovarisM220 sonicationtargets800–900bp; oldintegrationjunctionsremaininoriginalgDNA','Source50W/5percent/200cyclesburst/60s25°C')
    x=o.clean('fragment_clean',x,'Source0.7xAMPurepurification ofsonicatedDNA')
    x=o.op('endprep',x,'NEBNextUltraEndRepair/dA-Tailing20°C30min→65°C30minprepares5phosphate/3dAends; exactpopulationgeometryunassembled','HistoricalNEBNextUltramodule')
    a=g.input('linker','SourceiSL1positive/negativeannealedlinkerexample;actualsampleassignmentunassembled',K.OLIGO,state=linker())
    x=o.join('ligate',(x,a),'SourceBlunt/TA ligase16°CovernightTAjoinscustomlinker;3AmMOnegativearmretainedblocked')
    x=o.clean('lig_clean',x,'Source0.7xAMPureligatedDNApurification; published15µLelution')
    reps=[g.product('rep'+str(i),'IndependentallocatedsourcePCR1replicate'+str(i)+'; exactscalefrom15µLelutionunassembled',K.MOLECULAR) for i in range(1,5)]
    g.react('replicate_allocate','Mainatleast4replicates; sourcevolumeinconsistencykeptunresolved',S.SPLIT,((x,R.TEMPLATE),),tuple(reps))
    kept=[]
    for i,rep in enumerate(reps,1):
        k='rep'+str(i)
        f=g.input(k+'_ltr1','SourceHIVLTRprimer1; wholepolymerunassembled',K.OLIGO);r=g.input(k+'_linker1','Sourcesample-specificPCR1linkerprimer; fullpolymerunassembled',K.OLIGO)
        block=g.input(k+'_block1','SourceoptionalHIVblockingoligoforthisconditionalHIVentry;actualstudysampleunassigned',K.OLIGO,state=BlockingOligo().scene(),chemical_state=BlockingOligo())
        x=o.join(k+'_pcr1',(rep,f,r,block),'SourcePCR1Advantage2,5linear+20exponentialcycles,HIV70/67°Calternative; exactperstageprogramunassembled,notMLV80°Cdefault')
        x=o.qc(k+'_qc1',x,'Separateapproximately5µLPCR1agaroseassessment; productionremainderretained')
        x,rem=o.split(k+'_pcr2_allocate',x,'Source2µLPCR1templateintoPCR2','OtherPCR1DNAafterQC/PCR2allocation; laterfateunassigned')
        f=g.input(k+'_ltr2','SourceHIVLTRprimer2with12ntreplicateGolaybarcode;actualmember/bodyunassembled',K.OLIGO)
        r=g.input(k+'_linker2','Sourcesample-specificPCR2linkerprimeraddingP5arm;fullbodyunassembled',K.OLIGO)
        block=g.input(k+'_block2','SeparateoptionalblockingstockaliquotforPCR2; carriedPCR1companionsnotassumedremoved',K.OLIGO,state=BlockingOligo().scene(),chemical_state=BlockingOligo())
        x=o.join(k+'_pcr2',(x,f,r,block),'SourcePCR2Advantage2,5linear+15exponentialcycles,replicatebarcodeandIlluminaarms;newDNAcopiesnotBNAmanufacturedblocker')
        x=o.qc(k+'_qc2',x,'Separateapproximately5µLPCR2agaroseassessment')
        x,rem=o.split(k+'_pool_allocate',x,'Source10µLofPCR2replicateforwithin-samplepool','OtherPCR2materialretainedwithunresolvedlaterfate')
        kept.append(x)
    x=o.join('replicate_pool',tuple(kept),'Sourcecombine10µLfromeachPCR2replicateforonesample')
    x=o.clean('pool_clean',x,'Source0.7xAMPurepurification,20µLelution')
    x=o.qc('picogreen',x,'SeparatePicoGreenanalyticalmassmeasurement')
    others=g.input('other_samples','Otherindependentlypreparedsourcepurifiedsamplelibraries;actualbatchmembershipunassembled',K.MOLECULAR)
    x=o.join('sample_pool',(x,others),'Sourceequal-masssamplepoolingorcoverage-guidedallocation;exactbatchamountsunassigned')
    x=o.qc('screen_qpcr',x,'SeparateD1000ScreenTape/KAPASYBRFASTqPCRanalyticalaliquots; productionpoolretained')
    return finish(g,x,'PublishedMiSeqV2,source179R1/12I1/143R2cyclesasprinted;300cyclekitnamealsoasprinted. Fullprimerlanding/currentruncompatibilityunassembled. Linkerjunctioncountsproxyclonalabundancecomputationally,notphysicalUMIsbyassumption.')

def selftest():
    from checks import Check
    c=Check();g=protocol_graph();b=BlockingOligo()
    c('SourceHIVblockerretainsninebridgedsugarpositions',b.modified_count,9)
    c('Sourceblocker3primeisnonextendibleamino',b.terminal,'3AmMO')
    c('ReplicatePCR2receivesallocatedPCR1fraction',g.steps['rep1_pcr2'].inputs[0].material,'rep1_pcr2_allocate_production')
    c('PCR2copiedDNAdoesnotinheritblockersugarchemistry',g.materials['rep1_pcr2_out'].chemical_state,None)
    c('Pooldoesnotreceiveanalyticalgelaliquots',g.steps['replicate_pool'].inputs[0].material,'rep1_pool_allocate_production')
    c('SequencingsubstrateexcludesqPCRaliquot','screen_qpcr_production' in g.sequencing,True)
    src=Source(sources_root()/'inspiired__10.1016+j.omtm.2016.11.002'/'j.omtm.2016.11.002_supp_mmc1.pdf.txt','PublishedProtocolSI TablesS9/S10 renderedoriginal /tmp/inspiired19-blocker.png')
    if have(c,src,label='Modifiedsourceblockertranscription'):
        t=''.join(src.path.read_text().split());c('WholesourceHIVBNAstringmatchesoriginal',BLOCK in t,True)
        c('Negativearmwholephosphate/aminoorderingmatchesoriginal','5’Phos-'+revcomp(COMMON[:-1])+'-3’AmMO' in t,True)
    c('Fullgeometryunreviewednopromotion',g.publishable,False);c.report()

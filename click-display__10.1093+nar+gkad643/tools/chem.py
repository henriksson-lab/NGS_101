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

TITLE='Click display — corrected linker and mock-library NGS'
SOURCE='Published Zeng etal.2023 PMC10484664 Table1/Methods134–174; correction2024 PMC11662653 originalp1/table andp2text; originalSI TableS1'
LIMIT='Private incomplete source graph. Mandatory2024correction reverses original2023azide/DBCO labels and oligo names: correctedoligo1=5primeazide/internalTAMK/3primepuromycin;oligo2=5primeDBCO/internalCnvK/RTarm. OriginalVOR is not rewritten. Model uses correction. Its correctedFig1C labels BstI,whereas mainTable1/Methods specify BsaI-HF; model uses the explicit Table1/Methods enzyme and leaves the figure-label discrepancy unresolved. Modified groups are not normal DNA bases;5prime-to5primeclick bridge is not an ordinary phosphodiester. Exact full template/library geometry,IMAC workup variant/read cycles remain unassembled; no perfectclick/display/binding predicted.'
ENDPOINT='pool_out'
def sample_primer(n):
    sc=Scene();sc.strand('primer',[Segment('unassembledIlluminaadapterflank','[adapter flank]',placeholder=True),Segment('source2ntsample-roundindex',('AA','AC','AT','AG')[n],feature=feature('click.round','sample_index','fixed',note='SourceR0–R3sample/roundidentifier,notUMI')),Segment('unassembledlocusarm','[locus arm]',placeholder=True)],label='PartialsourceNGSforwardprimerinventory;wholepolymerunassembled');return sc
@dataclass(frozen=True)
class CorrectedOligo:
    number:int
    mod5:str
    left:str
    internal:str
    right:str
    mod3:str=''
    def __post_init__(self):
        expected={1:('5AzideN','AAAA','i5-TAMK','A'*21+'CC','3Puro'),2:('5DBCON','AAATTCCA','i3Cyan','GCCGCCCCCCGTCCT','')}
        if self.number not in expected or (self.mod5,self.left,self.internal,self.right,self.mod3)!=expected[self.number]:raise ValueError('Exact corrected source identity/modifications required')
    @property
    def as_written(self):return '/'+self.mod5+'/'+self.left+'/'+self.internal+'/'+self.right+('/'+self.mod3+'/' if self.mod3 else '')
    def scene(self,reacted=False):
        sc=Scene();sc.strand('oligo'+str(self.number),[Segment('source5primeDNAflank',self.left),Segment('sourceinternal/'+self.internal+'/',self.internal,placeholder=True,note='Actualsource noncanonical internal modifier;not ordinary DNAsequence'),Segment('source3primeDNAflank',self.right)],label='Correctedsourceoligo'+str(self.number),mod5='5prime-to5primetriazolejunction (consumed '+self.mod5+')' if reacted else self.mod5,mod3=self.mod3)
        return sc
O1=CorrectedOligo(1,'5AzideN','AAAA','i5-TAMK','A'*21+'CC','3Puro')
O2=CorrectedOligo(2,'5DBCON','AAATTCCA','i3Cyan','GCCGCCCCCCGTCCT')
@dataclass(frozen=True)
class DisplayClass:
    linker_bodies:tuple=(O1,O2)
    click_junction:str='nonphosphodiester5prime-to5primetriazole'
    reactive_click_groups:bool=False
    rna_crosslink:bool=False
    peptide_attached:bool=False
    cdna_copy:bool=False
    def crosslink(self):return replace(self,rna_crosslink=True)
    def transfer(self):
        if not self.rna_crosslink:raise ValueError('Selected route requires covalently associated RNA before peptide transfer')
        return replace(self,peptide_attached=True)
    def reverse_transcribe(self):
        if not self.peptide_attached:raise ValueError('Selected source RT follows transfer')
        return replace(self,cdna_copy=True)
    def scene(self):
        sc=Scene()
        for ol in self.linker_bodies:
            sc.strand('arm'+str(ol.number),[Segment('sourcearm'+str(ol.number)+'5primeflank',ol.left),Segment('source/'+ol.internal+'/',ol.internal,placeholder=True),Segment('sourcearm'+str(ol.number)+'3primeflank',ol.right)],label='Separate corrected covalentlinkerarm'+str(ol.number),mod5='shared5prime-to5primetriazolebridge',mod3='puromycin'+('–nascentpeptide' if self.peptide_attached else '') if ol.number==1 else ('newcDNAextension;sequenceunassembled' if self.cdna_copy else '3primeOH RTarm'))
        sc.note('arm2','Two arms connectedthroughsource5primeclickjunction;no ordinarylinearoligonucleotide concatenation. CnvK–RNA crosslink='+str(self.rna_crosslink)+'. RNA/template/reagentcarryover remains untilphysicalworkup.')
        return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE);o=Operations(g)
    a=g.input('oligo1','Corrected2024Oligo1',K.OLIGO,state=O1.scene(),chemical_state=O1)
    b=g.input('oligo2','Corrected2024Oligo2',K.OLIGO,state=O2.scene(),chemical_state=O2)
    ml=o.join('click',(a,b),'Equal0.5mMstockvolumes→250µMeach,roomtemperaturedark10min; successfulMLclass5prime-to5primetriazole,unreactedarmsremain',state=DisplayClass().scene(),chemical_state=DisplayClass())
    ml=o.op('ml_dilute',ml,'Source10folddilutioninDNAduplexbuffer and−80°Cstorage; no purification or100percentconversionassertion')
    ml_aliquots=[g.product('ml_round'+str(i),'Separate sourceML aliquotforround'+str(i)+';sourceTable1uses2µL,allocationfromstockexplicit',K.MOLECULAR) for i in range(1,4)]
    ml_rest=g.product('ml_remaining','OtherMLstockpreparation; laterfateunassigned',K.MOLECULAR)
    g.react('ml_allocate','Allocate separateMLstockaliquots tosourcepanningrounds',S.SPLIT,((ml,R.TEMPLATE),),(*ml_aliquots,ml_rest));g.terminals[ml_rest.id]=Terminal(ml_rest.id,Disposition.UNRESOLVED,'Stockremainingamount/futureuseunassigned.')
    x=g.input('input_template','SourcepreparationofnaiveN3CDARPinDNA+0.01percentD3G86 spike; primerassembly/templatefullgeometryunassembled',K.MOLECULAR,primary=True)
    production,reference=o.split('r0_allocate',x,'MockdisplayproductionDNA allocated to threepanningrounds','SourceR0referencePCRproduct forNGS; exact allocation/preceding amplificationunassembled')
    del g.terminals[reference.id];libs=[]
    for n in range(4):
        if n==0:sample=reference
        else:
            k='r'+str(n);d=DisplayClass()
            x=o.join(k+'_load',(production,ml_aliquots[n-1]),'SourceTable1DNA/ML/PURExpress/RNaseinhibitor mixture; actualcompanionsretained')
            x=o.op(k+'_transcribe',x,'37°C20minT7transcription makesnewRNA fromDNA; originalDNA/ML remain','PURExpresssolutionsA/B')
            d=d.crosslink();x=o.op(k+'_uv',x,'Source30secondUV crosslinksCnvK toassociatedRNA; no quantitativeconversioninference',state=d.scene(),chemical_state=d)
            x=o.op(k+'_bsai',x,'AddBsaI-HFv2;37°C10minwithoutshaking digests sourceDNAtemplate/terminatesfurthertranscription; fragmentsremain','BsaI-HFv2')
            x=o.op(k+'_translate',x,'37°C30minshaking translation; RNA/linker/peptide classesremain')
            d=d.transfer();x=o.op(k+'_transfer',x,'Roomtemperature20minpromotespuromycin-mediatednascentpeptidetransfer; untransferredclassesremain',state=d.scene(),chemical_state=d)
            x=o.op(k+'_edta',x,'EDTA30mM5mininactivatesribosomes,not physicalribosomeremoval','EDTA')
            x=o.op(k+'_mg',x,'AddMgCl2 thenRTbuffer/dNTP/DTTandwater; sourceexpandsreaction','MgCl2,Firststrandbuffer,dNTP,DTT')
            d=d.reverse_transcribe();x=o.op(k+'_rt',x,'SuperScriptRT42or50°C30min producesencodingcDNAextendedfromarm2,peptideremainsattachedthrougharm1; originalRNAcompanionsremain','SuperScriptIIorIV unresolvedselectedvariant',state=d.scene(),chemical_state=d)
            x=o.opaque(k+'_imac_prepare',x,'SourceIMACpreparation: ultrafiltration/PBSexchange OR10folddilutiontoPBS+10mMimidazole; actualmockexperimentvariantunassigned')
            x=o.clean(k+'_imac',x,'Sourceone-stepNi-NTAIMACretainsHis-taggedpreparation; standardprotocolworkup/eluategeometryunreviewed')
            target=g.input(k+'_target','Source200nMbiotinylatedeGFP mockselectionbait; attachedbiotinnotDNAbarcode',K.MOLECULAR)
            x=o.join(k+'_bind',(x,target),'Source1hroomtemperaturebinding')
            x=o.clean(k+'_pull',x,'SourceT1orneutravidinmagneticbeadpull-down and3–8washes; actualround-beadidentity/washcountunassigned,retainbeadassociatedpopulation')
            x=o.op(k+'_resuspend',x,'Resuspendbeads1xStandardTaqbuffer; noelution/separationinvented')
            f=g.input(k+'_f','Source1938T7promoterPCRprimer; wholebodyunassembled',K.OLIGO);r=g.input(k+'_r','Source2643T7terminatorPCRprimer; wholebodyunassembled',K.OLIGO)
            x=o.join(k+'_pcr',(x,f,r),'TaqPCRcopiesencodingDNAfrombeadpreparation; newDNAdoesnotinheritpuromycin/CnvK/TAMRA/triazole chemistry')
            x=o.clean(k+'_gel',x,'SourcegelpurificationofpanningPCRDNA')
            sample,production=o.split(k+'_sample',x,'Source roundPCRDNAaliquot forNGS','OtherroundPCRDNAforcontinuedpanningorunassignedremainder')
            if n<3:del g.terminals[production.id]
        f=g.input('ngs'+str(n)+'_f','SourceDARPin_NGS_F'+str(n)+'; actual2ntsample/roundbarcode aa/ac/at/ag,fullpolymerunassembled',K.OLIGO,state=sample_primer(n))
        r=g.input('ngs'+str(n)+'_r','SourceDARPin_NGS_R,fullpolymerunassembled',K.OLIGO,state=sample_primer(n))
        lib=o.join('ngs'+str(n),(sample,f,r),'Q5NGSPCRaddsIlluminaadapterandtwo-nucleotideroundbarcode,notUMI')
        lib=o.clean('ngs'+str(n)+'_gel',lib,'SourceZymoCleanGelDNARecovery')
        lib,remainder=o.split('ngs'+str(n)+'_allocate',lib,'Source150ngofthisroundNGSlibrary forpool','OtherpurifiedNGSlibraryremainder; laterfateunassigned')
        libs.append(lib)
    x=o.join('pool',tuple(libs),'Poolsource150ngofeachR0–R3NGSlibrary; no chemicaldeletionofreadsorbarcode')
    return finish(g,x,'PublishedAzentaAmplicon-EZMiSeqpaired-end,cycles/fullstrandlanding unreviewed. 2ntsamplebarcodes identifyrounds; FLASHmerging/translation are computational analyses,not molecularreactions.')

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();d=DisplayClass();post=d.crosslink().transfer().reverse_transcribe()
    c('Corrected reactivegroups areazide1/DBCO2',(O1.mod5,O2.mod5),('5AzideN','5DBCON'))
    c('Clickretainsbothoriginalarms',post.linker_bodies,(O1,O2))
    c('5to5junctionnotlinearphosphodiester',post.click_junction,'nonphosphodiester5prime-to5primetriazole')
    c('Consumedclickgroupsnotreactive',post.reactive_click_groups,False)
    c('RTretainsphysicalcrosslinkandpeptide',(post.rna_crosslink,post.peptide_attached,post.cdna_copy),(True,True,True))
    try:d.transfer();bad=False
    except ValueError:bad=True
    c('SelectedsourcecannottransferbeforeRNAassociation',bad,True)
    c('UVprecedesDNAtemplatedigestion',g.steps['r1_bsai'].inputs[0].material,'r1_uv_out')
    c('DNAdigestiondoesnotautomaticallydiscardfragments',g.steps['r1_translate'].inputs[0].material,'r1_bsai_out')
    c('OriginalproteinlinkerdoesnotbecomeallPCRdaughterchemistry',g.materials['r1_pcr_out'].chemical_state,None)
    c('EachpanningroundusesseparateMLaliquot',g.steps['r2_load'].inputs[1].material,'ml_round2')
    c('NGSpoolusesphysical150ngallocation',g.steps['pool'].inputs[0].material,'ngs0_allocate_production')
    c('Postselectiontemplatecontinuesafterphysicalgel',g.steps['r2_load'].inputs[0].material,'r1_sample_other')
    c('Thirdroundremainderretainedunresolved',g.terminals['r3_sample_other'].disposition,Disposition.UNRESOLVED)
    src=Source(sources_root()/'click-display__10.1093+nar+gkad643'/'nar_gkae1172_PMC11662653.xml.txt','Published2024correction doi10.1093/nar/gkae1172,Tablep1,originalrendered /tmp/click18-correction.png')
    if have(c,src,label='Publishedcorrectedwholemodifiedoligotranscription'):
        txt=src.path.read_text();table=txt[txt.index('Name	Sequence	Vendor'):txt.index('The description of Oligo')];collapsed=''.join(table.split())
        c('WholecorrectedO1modifiersandbasesmatchoriginal',O1.as_written in collapsed,True)
        c('WholecorrectedO2modifiersandbasesmatchoriginal',O2.as_written in collapsed,True)
    c('No finishedclaimfrompartialgeometry',g.publishable,False);c.report()

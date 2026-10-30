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
from checks import Source,have
from storage import sources_root
TITLE='LAM-HTGTS — source-qualified linear amplification and bridge capture'
SOURCE='AssociatedauthorMS PMC4320661 Methods107 LAM-specificparagraph; publishedMOESM10Table8 originalp30 (selectedmodifiedoligos); definingVORapplicabilityunresolved'
LIMIT='Private incomplete source graph. Operational2015LAMrecipe is explicitlyassociatedauthor-manuscript scope; publishedSI oligos independentlycontrolinventory. OriginalemulsionHTGTSend-repair/A-tailrecipe is not imported. Bothsourcebridgevariantsareinventory;actualend-conditionedallocationisexplicitlyunassembled,notallusedonebatch. EcoNI cuts only recognition-site-containingmolecules; digestiondoesnotphysicallydiscardthem. Fullbait/rearrangedDNAgeometry,predecessorbridgeinstructionsandkitrevisionsremainunreviewed.'
ENDPOINT='loading_out'
CORE='GCGACTATAGGGCACGCGTGG'
@dataclass(frozen=True)
class DigestClass:
    has_econi:bool
    cut:bool=False
    def __post_init__(self):
        if type(self.has_econi)is not bool or type(self.cut)is not bool or (self.cut and not self.has_econi):raise ValueError('A selectedEcoNIcut requires an explicitrecognitionsite')
    def digest(self):return replace(self,cut=self.has_econi)

def upper():
    sc=Scene();sc.strand('upper',[Segment('bridgecore',CORE),Segment('sourceT','T'),Segment('fiveNrandomsplint','N'*5,placeholder=True,note='SourceTNNNNNbridgeannealingwindow,notUMI')],mod3='AmMO');return sc

def lower(extra_a=False):
    sc=Scene();sc.strand('lower',[Segment('source5primeA' if extra_a else 'sourceempty','A' if extra_a else ''),Segment('bridgecorecomplement',revcomp(CORE))],mod5='P',mod3='AmMO');return sc

def bridge_inventory():
    sc=Scene();sc.strand('upper',[Segment('bridgecore',CORE),Segment('sourceTNNNNN','T'+'N'*5,placeholder=True,note='Six-basepartiallydegeneratesplint;notidentifier')],mod3='AmMO')
    sc.anneal('lower',[Segment('bridgecorecomplement',revcomp(CORE))],to='upper',pair=('bridgecorecomplement','bridgecore'),mod5='P',mod3='AmMO')
    sc.note('lower','PublishedssA-lowerintendedforA-tailedDNA;5AprimealternativeisaseparateinventoryforG/C/T-tailedorbluntDNA. Thisinventorydoesnotassigngenomicendclassesorcompletejunctions.');return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE);o=Operations(g)
    x=g.input('genomic','SourceisolatedgenomicDNAaftercellularbait/preybreakjoining; cellularjoinnotaninvitronewligationstep',K.MOLECULAR,primary=True)
    x=o.op('sonicate',x,'SourcegenomicDNAsonication; LAM-specificsizeallocationunassembled','Sonication')
    p=g.input('bio_primer','Sourcesingle5biotinLAMprimer; chosenlocus/wholepolymerunassembled',K.OLIGO)
    x=o.join('lam50',(x,p),'Taqsingle-primerLAM50cyclesmakesnew5biotinDNAcopies; originalgDNA/unusedprimerstillpresent')
    x=o.op('lam100',x,'AddmoreTaqandperformadditional50cycles,notreplacebyone-cyclePEMextension','AdditionalTaqpolymerase')
    x=o.clean('capture',x,'SourceMyOneC1streptavidinbindingretainsbiotin-bearingpreparation; exactwash/stranddenaturationunreviewed')
    u=g.input('upper','PublishedssA-upper-TN',K.OLIGO,state=upper())
    lo=g.input('lower','PublishedssA-lower; intendedA-tailedDNAinventory',K.OLIGO,state=lower())
    other=g.input('lower_a','PublishedssA-lower-A; intendedG/C/T-tailedorbluntDNAinventory;actualsourceallocationunassembled',K.OLIGO,state=lower(True))
    assigned=g.product('assigned_bridge','Sourceend-conditionedbridgepreparation;lowerforA-tailedDNAorlower-Aforotherends,actualallocationunassembled',K.MOLECULAR)
    remainder=g.product('bridge_inventory_remainder','Otheradapterstock/preparationsnotassignedthisselectedreaction',K.MOLECULAR)
    g.react('bridge_allocate','Unreviewedend-conditionedallocationofsourcebridgeinventory',S.MAGIC,((u,R.TEMPLATE),(lo,R.STRUCTURAL),(other,R.STRUCTURAL)),(assigned,remainder),reason='PublishedTable8distinguishesA-tailedversusotherendclasses;actualsubstrate/variantamountsnotassignedinselectedMethods. NoDNAreceivesbothlowerstrandsbyassumption.');g.terminals[remainder.id]=Terminal(remainder.id,Disposition.UNRESOLVED,'Unassignedstockfateremainsunresolved.')
    x=o.join('bridge',(x,assigned),'Sourceon-beadovernightbridge-adapterligation15percentPEG8000/1mMhexamminecobalt; alternativeadapterallocation/wholejunctionunassembled,notforcedbothvariantsoneDNA')
    x=o.op('lig_temp',x,'PublishedMSligationtemperatures25°C60min→22°C120min→16°C8–12h; beads/reagentsremain')
    f=g.input('nested_bait','Sourcesite-specificnestedprimerwithI5-relatedtail;wholepolymer/memberunassembled',K.OLIGO)
    r=g.input('nested_adapter','SourceadapterprimerwithI7-relatedtail;wholepolymerunassembled',K.OLIGO)
    x=o.join('nested',(x,f,r),'Source15cyclenested-locusPCR;newDNAcopiedfrombiotin/linker-bearingtemplates,biotin/blockersnotautomaticallycopied')
    x=o.op('econi',x,'EcoNIblockingdigestcleavessite-containinggermlineclasses; site-freeclassesremainintact,allfragmentsremainuntilphysicalworkup','EcoNI',chemical_state=(DigestClass(True).digest(),DigestClass(False).digest()))
    f=g.input('p5','SourceTable8P5PCRIIIprimer;fullbodyunassembled',K.OLIGO);r=g.input('p7','SourceTable8P7PCRIIIprimer;fullbodyunassembled',K.OLIGO)
    x=o.join('final_pcr',(x,f,r),'Source10–15cyclesP5/P7step-outPCRreconstructsMiSeqsequencearms; intact/cuttemplates/reagentsremainwithoutinventedcleanup')
    x=o.opaque('loading',x,'KnownMiSeqlibraryacquisitionhandoff; productioncleanup/QC/sizeallocationunreviewedinselectedparagraph')
    return finish(g,x,'PublishedassociatedsourceMiSeqpaired-end;runlength/fullprimerlandingunreviewed. Molecularrearrangementjunctionsareexistinggenomictemplates;mapping/duplicatefiltersnotchemicalseparations.')

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('Recognition-site-freeclassnotartificiallycut',DigestClass(False).digest().cut,False)
    c('Site-containingclassselectedcut',DigestClass(True).digest().cut,True)
    c('Originalsiteclassnotmutatedbyconditionaldigest',DigestClass(True).cut,False)
    c('Additional50cyclesretainfirstreaction',g.steps['lam100'].inputs[0].material,'lam50_out')
    c('Bridgeoperatesonsourcecapturedfraction',g.steps['bridge'].inputs[0].material,'capture_kept')
    c('FinalPCRreceivesunremoveddigestmixture',g.steps['final_pcr'].inputs[0].material,'econi_out')
    c('Noend-repairrecipeimportedfromolderHTGTS','end_repair' in g.steps,False)
    src=Source(sources_root()/'lam-htgts__10.1038+nbt.3101'/'nbt.3101_41587_2015_BFnbt3101_MOESM10_ESM.pdf.txt','PublishedTable8p30render /tmp/lam19-bridge.png')
    if have(c,src,label='Wholemodifiedbridgeoligotranscription'):
        t=''.join(src.path.read_text().split())
        c('Sourceupperwholeblockedbody',CORE+'TNNNNN-/3AmMO/' in t,True)
        c('Sourcelowerwhole5phosphate3aminobody','/5Phos/'+revcomp(CORE)+'-/3AmMO/' in t,True)
        c('Sourcealternatewhole5Aphosphate3aminobody','/5Phos/A'+revcomp(CORE)+'/3AmMO/' in t,True)
    c('Incompletegeometry/sourceversionpreventpromotion',g.publishable,False);c.report()

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

from chemdraw import revcomp
SOURCE='AssociatedAMNIHMS338311/PMC3243364 Methods132,Results40/Fig1;DOI10.1016/j.cell.2011.11.013. DefiningpublishedCellversionapplicability andchemistrySI unreviewed/gated;originalSOLiDscope,notlaterIlluminaadaptation.'
@dataclass(frozen=True)
class SelectedResection:
    remaining:str
    removed:str
    original:str
    inferred:bool=True
    def __post_init__(self):
        if self.removed+self.remaining!=self.original:raise ValueError('ResectionmustretainwholeoriginalDNAinventory')
    def selected_complement_copy(self):return revcomp(self.remaining)
def selected_5prime_resection(sequence,stop):
    if not sequence or set(sequence)-set('ACGT') or type(stop)is not int or not 0<=stop<len(sequence):raise ValueError('Explicitsingle-strandDNAandcallerselectedphysicalstoprequired')
    return SelectedResection(sequence[stop:],sequence[:stop],sequence)
def protocol_graph():
    g=ProtocolGraph('ChIP-exo — associated2011SOLiDmain-methodscope',authority=SOURCE,representation='published_instance')
    x=g.input('chip','SourceimmunoprecipitatedcrosslinkedchromatinonSepharose;priorstandardChIPVenters2009methodunreviewed',K.MOLECULAR,primary=True)
    x=op(g,'polish',x,'Sourceon-resinDNApolishing',S.END_REPAIR,reagent='ExactchemistrySIenzyme/endsunreviewed',catalyst=True)
    a=g.input('p2_adapter','SourceoriginalSOLiDP2adapterfamily;exactsequence/modifications/endgeometryunreviewed',K.OLIGO);cat=g.input('p2_ligase','Sourceadapterligationcatalyticconditionsunreviewed',K.CATALYST);y=g.product('p2_product','P2-bearingchromatinDNAplusunligatedDNA/freeP2/resin/reactioncompanions',K.MOLECULAR)
    g.react('p2_ligate','Sourceon-resinP2adapterligation',S.LIGATION,((x,R.TEMPLATE),(a,R.STRUCTURAL),(cat,R.CATALYST)),(y,),geometry=uncertain_ligation())
    x=op(g,'lambda',y,'Sourcelambda5to3exonucleaseresectiontowardcrosslinkbarriers',S.DIGESTION,reagent='Sourcelambdaexonuclease;stopdistribution/yieldnotpredicted',catalyst=True,product='SelectedresectedDNAclassesplusremovedDNAproducts,unresectedDNA,P2/protein/resin/reagents;digestionnotphysicaldiscard')
    x=cleanup(g,'elute',x,'Sourceelutionoflambda-treatedDNAfromresin;exactworkupunreviewed')
    x=op(g,'pk',x,'Sourcecrosslinkreversal/proteinaseK65C',S.DEPROTEINIZATION,reagent='SourcePK/65C;detailedconditionsunreviewed',catalyst=True,product='DNAincludingselectedresectionproductsplusprotein/peptideremnantsandreagentsuntilprecipitation')
    x=cleanup(g,'ethanol',x,'SourceethanolprecipitationafterPK')
    x=op(g,'denature',x,'Source95CDNAdenaturation',S.THERMAL_TREATMENT,reagent='Source95Cconditions',product='SeparatephysicalDNAstrandsplusretainedcompanions;denaturationdoesnotdeleteoriginalDNA')
    primer=g.input('p2_primer','SourceoriginalSOLiDP2primer;exactbody/endunreviewed',K.OLIGO);enzyme=g.input('extension_enzyme','Sourceprimer-extensionpolymeraseunreviewed',K.CATALYST);y=g.product('extended','SelectednewcomplementaryP2-primedcopiesplusoriginalresectedDNA,unextendedprimerandcompanions;fullprimer/productgeometryunassembled',K.MOLECULAR)
    g.react('p2_extend','SourceP2-primerextensionafter95Cdenaturation',S.EXTENSION,((x,R.TEMPLATE),(primer,R.STRUCTURAL),(enzyme,R.CATALYST)),(y,))
    a=g.input('p1_adapter','SourceoriginalSOLiDP1adapterfamily;exactbody/endgeometryunreviewed',K.OLIGO);cat=g.input('p1_ligase','SourceP1ligationcatalyticconditionsunreviewed',K.CATALYST);x=g.product('p1_product','SourceP1/P2libraryclassesplusunligatedtemplates/adapters/remnants;fulllibraryunassembled',K.MOLECULAR)
    g.react('p1_ligate','SourceP1adapterligationatlambda-definedend',S.LIGATION,((y,R.TEMPLATE),(a,R.STRUCTURAL),(cat,R.CATALYST)),(x,),geometry=uncertain_ligation())
    x=unknown_pcr(g,'pcr',x,'SourcePCRamplification;cycles/completegeometryunreviewed','SourceoriginalSOLiDlibraryPCR')
    x=cleanup(g,'gel',x,'Sourcegelpurificationofamplifiedlibrary;windowunreviewed')
    lib=g.product('library','Sourcephysicalsequencinglibrary;fullmolecular/readgeometryunassembled',K.LIBRARY);g.react('library_identity','Sourcefinalpreparationidentity',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='SourceactualSOLiDsequencer,notIllumina. Exactreadcycles/colourspaceadapterlandingunreviewed. Shiftedtagcoordinatesandpeakpairsarecomputational;notgenomicDNAedits orproofeverymoleculestopsperfectlyatonebarrier.'))
    g.add_source_item(SourceItem('pending','PublishedCelloriginal,chemistrySupplementDOCandcitedVenters2009ChIP/vendororiginals','AssociatedAMmainonlysourcechemistryread;actualSIbarrierrecorded. Noacceptedabsence/nofamilyrecipeimport.'))
    return g
def selftest():
    from checks import Check
    c=Check();a=selected_5prime_resection('ACGTGCC',3)
    c('resectionretainsremovedandremainingDNA',a.removed+a.remaining,a.original)
    c('selectednewcopyiscomplementnotnativeidentity',a.selected_complement_copy(),'GGCA')
    c('callerselectedreferenceinferred',a.inferred,True)
    try:selected_5prime_resection('ACGT',4);bad=False
    except ValueError:bad=True
    c('stopbeyondretainedDNArejected',bad,True)
    g=protocol_graph();c('digestiondebrisreachessourceelution',g.steps['elute'].inputs[0].material,'lambda_product')
    c('primerextensionfollowsphysical95Cstrandseparation',g.steps['p2_extend'].inputs[0].material,'denature_product')
    c('P1ligationfollowsP2copying',g.steps['p1_ligate'].inputs[0].material,'extended')
    c('sourceSI/geometryboundariespreventpromotion',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()

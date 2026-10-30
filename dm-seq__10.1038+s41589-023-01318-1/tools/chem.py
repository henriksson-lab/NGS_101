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

SOURCE='AssociatedPMC10763687 optimizedMethods199/GBMacquisition227; publishedMOESM1 Table1renderedp2 anddetailedprocedure572–800; DOI10.1038/s41589-023-01318-1. Main50µlA3A andSIcommercial100µlformamidebranchesdistinct.'
COPY_WRITTEN='GA/5mC/TGGAGTT/5mC/AGA/5mC/GTGTG/5mC/T/5mC/TT/5mC/GAT/5mC/'
COPY_BODY=COPY_WRITTEN.replace('/5mC/','C')
@dataclass(frozen=True)
class ModifiedDNA:
    sequence:str
    marks:tuple=()
    def __post_init__(self):
        if not self.sequence or set(self.sequence)-set('ACGTU'):raise ValueError('explicitDNA nucleobasesrequired')
        if len(dict(self.marks))!=len(self.marks) or any(type(i)is not int or not 0<=i<len(self.sequence) or self.sequence[i]!='C' or mark not in {'5mC','5hmC','5gmC','5cxmC','5pyC'} for i,mark in self.marks):raise ValueError('namedmodifiedcytosinepositionsrequired')
    def glucosylate_selected(self):return replace(self,marks=tuple((i,'5gmC' if m=='5hmC' else m) for i,m in self.marks))
    def carboxymethylate_selected(self,positions):
        marks=dict(self.marks)
        for i in positions:
            if type(i)is not int or i<0 or self.sequence[i:i+2]!='CG' or i in marks:raise ValueError('selectednativeunmodifiedCpGrequired; notsynthetic5mC orotherCcontext')
            marks[i]='5cxmC'
        return replace(self,marks=tuple(sorted(marks.items())))
    def deaminate_selected(self):
        marks=dict(self.marks);seq=list(self.sequence);retained=[]
        for i,b in enumerate(seq):
            if b!='C':continue
            m=marks.get(i)
            if m=='5hmC':raise ValueError('unprotected5hmC fate not assigned bythisselectedclass')
            if m=='5mC':seq[i]='T'
            elif m is None:seq[i]='U'
            else:retained.append((i,m))
        return ModifiedDNA(''.join(seq),tuple(retained))
    def copied_insert_reference(self):
        # Sequence relation only: fullBst product also includesactualmethylatedprimer.
        from chemdraw import revcomp
        seq=revcomp(self.sequence.replace('U','T'));return ModifiedDNA(seq,tuple((i,'5mC') for i,b in enumerate(seq) if b=='C'))

def copy_primer():
    sc=Scene();sc.strand('copy_primer',[Segment('v2',COPY_BODY,note='EveryCexplicitlywritten/5mC/inpublishedTable1; DNAprimer, notnativegenomicmethylation')]);return sc

def graph_for_branch(branch='optimized_gbm'):
    if branch not in {'optimized_gbm','si_commercial'}:raise ValueError('sourceoptimizedGBMorSIcommercialalternativebranchrequired')
    g=ProtocolGraph('DM-Seq — '+branch,authority=SOURCE,representation='published_instance')
    x=g.input('genomic','SourceextractedgenomicDNA; originalcytosinemodificationsnotinstantiated',K.MOLECULAR,primary=True)
    spike=g.input('spikes','SourceextractedunmethylatedpUC19/CpGmethylatedlambda/T4hmC DNAcontrolmix; targetoccupancymeasurednotassumedperfect',K.MOLECULAR);y=g.product('spiked','SourceDNA/controlmixture; independentgenomicDNA speciesnotonecovalentstrand',K.MOLECULAR)
    g.react('spike','Source1:100controlmix:genomicDNA massratio',S.GENERIC,((x,R.TEMPLATE),(spike,R.SPIKE_IN)),(y,));x=y
    x=op(g,'shear',x,'SourceCovarisshearing; detailedSI~550bp example',S.FRAGMENTATION)
    x=cleanup(g,'shear_clean',x,'Source1.2xSPRI')
    x=qc(g,'input_qubit',x,'SourceanalyticalQubitquantification')
    x=op(g,'endprep',x,'SourceUltraIIendprep20°C30min65°C30min',S.END_REPAIR,reagent='SourcehistoricalUltraIIEndPrepenzymes',catalyst=True)
    adapter=g.input('adapter','SourcecustomIDTxGENYadapterswithallC=5propynylC; fullbodies/endgeometryunreviewed',K.ADAPTER);cat=g.input('ligase','SourceUltraIILigationMasterMix/enhancer',K.CATALYST);y=g.product('ligated','Source5pyCadapterligatedDNAplusunligatedDNA/adapter/reagents; sourcePCRscaffoldnotmCTruSeqdefault',K.MOLECULAR)
    g.react('ligate','Sourceannealed5pyCadapterligation20°C15min',S.LIGATION,((x,R.TEMPLATE),(adapter,R.STRUCTURAL),(cat,R.CATALYST)),(y,),geometry=uncertain_ligation())
    x=cleanup(g,'adapter_clean',y,'Source1.2xSPRI')
    x=qc(g,'ligated_qubit',x,'SourceanalyticalQubitofligatedDNA')
    used=g.product('copy_input','Source'+('10ng' if branch=='optimized_gbm' else '~1–10ng')+'ligatedDNAallocatedtocopyreaction',K.MOLECULAR);rest=g.product('copy_remainder','UnallocatedligatedDNA; laterfateunassigned',K.MOLECULAR);g.react('copy_allocate','Sourcecopyinputallocation',S.SPLIT,((x,R.TEMPLATE),),(used,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Sourceallocatedmassknown; noassumptionofperfectligationyieldorzero remainder.'))
    primer=g.input('copy_primer','Publishedfullymethylatedcopyprimer(v2)',K.OLIGO,state=copy_primer(),chemical_state=ModifiedDNA(COPY_BODY,tuple((i,'5mC') for i,b in enumerate(COPY_BODY) if b=='C')));cat=g.input('bst','SourceBstlargefragment/dATP/dGTP/dTTP/dmCTP',K.CATALYST);y=g.product('copied','OriginalDNAplusselectedmethylatedBstcopy-strands retainingactualv2primer andpossibleuncopied/displacedDNA/enzymes; sequencecopydoesnotpreservenativemethylation',K.MOLECULAR)
    g.react('copy','Sourcecopyprimeranneal thenBst65°C30min; notunoptimizedKlenow/v1',S.EXTENSION,((used,R.TEMPLATE),(primer,R.STRUCTURAL),(cat,R.CATALYST)),(y,))
    x=op(g,'glucosylate',y,'SourcebetaGT/UDPglucose37°C1hprotectoriginal5hmC',S.CHEMICAL_CONVERSION,product='Selected5hmC→5gmCDNA withsyntheticmCcopy/originalbiologicalmCandothercytosines/remnantsretained')
    x=op(g,'mungbean',x,'SourceMungBeanNuclease30°C30min',S.DIGESTION,reagent='SourceMungBeanNuclease; exactefficiency/enddistributionunreviewed',catalyst=True,product='SelectedcopiedDNAclassespluscleavageproducts/uncleavedDNA/reagents; noinstantremovalorperfectcopyselection')
    x=cleanup(g,'mung_clean',x,'Source1.2xSPRI')
    x=op(g,'cx',x,'SourceMBP-M.MpeI-N374K/CxSAM37°Covernight;95°C5min',S.CHEMICAL_CONVERSION,product='SelectedoriginalunmodifiedCpGs→5cxmC; original5mCand5gmC/syntheticmCcopy/5pyCadaptersremainotherclasses; carryoverretained')
    x=op(g,'pk',x,'SourceProteinaseK37°C15min',S.DEPROTEINIZATION,reagent='SourceProteinaseK',catalyst=True,product='ModifiedDNAplusproteolysisproducts/buffersuntilSPRI')
    x=cleanup(g,'pk_clean',x,'Source1.2xSPRI')
    if branch=='si_commercial':
        x=op(g,'denature',x,'SourceSI20%formamide85°C10min thencoldblock',S.THERMAL_TREATMENT,reagent='Source16µlDNA+4µlformamide; coolingblocknotationreportedinsource')
        x=op(g,'deaminate',x,'SourcecommercialAPOBEC100µltotal;4°C→50°Cslowramp/hold',S.CHEMICAL_CONVERSION,product='Selectedoriginal5mC→T, unprotectedC→U;5cxmC/5gmC/5pyCretained; deaminatedsyntheticcopyandotherDNA/reagentsremainphysical')
    else:
        x=op(g,'denature',x,'Sourceoptimizedmain snapcoolingafterDNAdenaturation',S.THERMAL_TREATMENT,reagent='SourceoptimizedA3AconditionsperdirectACEreference; nocommercial100µlformamidemerge')
        x=op(g,'deaminate',x,'SourceoptimizedA3A50µlreactionperreference25',S.CHEMICAL_CONVERSION,product='Selectedoriginal5mC→T, unprotectedC→U;5cxmC/5gmC/5pyCretained; deaminatedsyntheticcopyandotherDNA/reagentsremainphysical')
    x=cleanup(g,'deamination_clean',x,'Source1.2xSPRI')
    x=qc(g,'qpcr',x,'SourceSIanalyticalqPCRaliquotsetsindexcycles; materialnotreusedforlibrary')
    x=unknown_pcr(g,'index_pcr',x,'SourceindexingPCR; detailedSIadaptivecycles/Table2example10cycles','SourceKAPAHiFiHotStartUracil+ mix; originaltemplatesandcopy-strandsnotassumedphysicallygone')
    x=cleanup(g,'pcr_clean',x,'Source0.8xSPRI')
    x=qc(g,'final_qc',x,'SourceQubit/Bioanalyzeranalyticalmaterialseparatefromproductionlibrary')
    g.sequence(x,SequencingPlan(limitation=('ActualGBMNextSeqHighOutputv2.5single-end75bp; branchsourceoptimizedmain50µlA3A. ' if branch=='optimized_gbm' else 'SourceSIcommercialAPOBECalternative usedExtendedData7/8; kitmateriallistMiSeqnotassignedasactualrun. ')+'Exactadapter/copy/PCR/readgeometryandSI10µlelutionversus19µlindexDNA scalingunassembled. Copystrandpooramplificationandcomputationalfiltersarenotphysicaldiscard.'))
    g.add_source_item(SourceItem('pending','Historical5pyCadapterfulloligos/ACEA3Aconditions/indexprimeroriginals','Exactlinkedreagents/primer geometryunreviewed; sourcevolumeallocation/copystrandfateslimitedasnote. NogeneralEM-seqorTAPSrecipeimported.'))
    return g

def protocol_graph():return graph_for_branch()
def selftest():
    from checks import Check
    c=Check();original=ModifiedDNA('ACGCCGC',((4,'5mC'),(6,'5hmC')));protected=original.glucosylate_selected().carboxymethylate_selected((1,));product=protected.deaminate_selected()
    c('deaminationdistinguishesphysicalnative5mCTfromunprotectedCU',(product.sequence[4],product.sequence[3]),('T','U'))
    c('selectedprotectedcytosinesretainactualadducts',dict(product.marks),{1:'5cxmC',6:'5gmC'})
    c('nativeoriginalisnotrewrittenbyconstructingchemicalproduct',original.sequence,'ACGCCGC')
    copy=original.copied_insert_reference();c('syntheticcopycontains5mCindependentofnativepattern',set(dict(copy.marks)),{i for i,b in enumerate(copy.sequence) if b=='C'})
    try:original.carboxymethylate_selected((4,));bad=False
    except ValueError:bad=True
    c('native5mCcannotbepassedasunmodifiedCxMTasesubstrate',bad,True)
    g=protocol_graph();c('qPCRreactionnotpreparativePCRinput',g.steps['index_pcr'].inputs[0].material,'qpcr_retained')
    c('uncopiedDNAdebrisretaineduntilSPRI',g.steps['mung_clean'].inputs[0].material,'mungbean_product')
    c('commercialalternativekeepsseparatedenaturationcondition','formamide' in graph_for_branch('si_commercial').steps['denature'].label and 'formamide' not in g.steps['denature'].label,True)
    c('incompletegeometrypreventspromotion',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()

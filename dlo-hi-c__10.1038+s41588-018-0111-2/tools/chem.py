"""Protocol-owned source-qualified incomplete graph; no publication marker."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,feature,revcomp
import illumina
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry

def scene(label):
    sc=Scene();sc.strand('preparation',[Segment('population','['+label+']',placeholder=True)]);return sc

def op(g,key,x,label,kind=S.MAGIC,*,reagent='source reagents; geometry unassembled',product=None,state=None,catalyst=False):
    stock=g.input(key+'_stock',reagent,K.CATALYST if catalyst else K.CONDITION)
    out=g.product(key+'_product',product or label+'; molecular mixture/remnants unassembled',K.MOLECULAR,state=state)
    kw={'reason':'Source operation established; exact molecular ends, companion species and geometry unassembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
    inputs=((x,R.TEMPLATE),(stock,R.CATALYST if catalyst else R.CONDITION))
    if kind==S.LIGATION:kw['geometry']=uncertain_ligation();inputs+=((x,R.STRUCTURAL),)
    g.react(key,label,kind,inputs,(out,),**kw);return out

def cleanup(g,key,x,label):
    stock=g.input(key+'_condition',label,K.CONDITION)
    out=g.product(key+'_kept','Source recovered preparation; recovery and carryover unquantified',K.MOLECULAR)
    waste=g.product(key+'_waste','Actual source-separated workup fractions',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(stock,R.CONDITION)),(out,waste));g.discard(waste,'Source physical cleanup, no perfect partition implied.');return out

def uncertain_ligation():return AdapterLigationGeometry('unavailable','unavailable',(scene('unreviewed exact ends and junction geometry'),),'Known source ligation, complete molecular geometry not assembled.')



def qc(g,key,x,label):
    keep=g.product(key+'_retained','Sourceproductionremainder; analyticalaliquotamountnotinvented',K.MOLECULAR);test=g.product(key+'_assay','Separateanalyticalaliquotforthe sourceassay',K.QC)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(keep,test));g.mark_terminal(Terminal(test.id,Disposition.QC,'Analyticalmaterialdoesnotbecomethepreparativesequencingtemplate.'));return keep

def source_adapter(which=1):
    if which not in (1,2):raise ValueError('Published FigureS2 adapter1 or adapter2 required')
    body=illumina.TRUSEQ_READ1 if which==1 else illumina.INDEX1_PRIMER
    sc=Scene()
    if which==1:
        sc.strand('top',[Segment('known arm',body),Segment('random cohesive end','NN',placeholder=True)],mod5='OH',mod3='OH')
        sc.anneal('bottom',[Segment('paired arm',revcomp(body))],to='top',pair=('paired arm','known arm'),mod5='OH',mod3='OH')
        sc.ligation_end('top','random cohesive end',kind='cohesive',polarity='3-prime',owner='adapter')
    else:
        sc.strand('top',[Segment('known arm',body)],mod5='OH',mod3='OH')
        sc.anneal('bottom',[Segment('paired arm',revcomp(body)),Segment('random cohesive end','NN',placeholder=True)],to='top',pair=('paired arm','known arm'),mod5='OH',mod3='OH')
        sc.ligation_end('bottom','random cohesive end',kind='cohesive',polarity='3-prime',owner='adapter')
    return sc

def source_pcr_primer(which=1):
    if which not in (1,2):raise ValueError('Published FigureS2 primer1 or primer2 required')
    sc=Scene()
    if which==1:sc.strand('primer',[Segment('FigureS2 primer1',illumina.TRUSEQ_P5_FULL)])
    else:sc.strand('primer',[Segment('source P7',illumina.P7),Segment('unassigned sample index','[index]',placeholder=True,feature=feature('dlo.sample','sample_index','unknown',note='FigureS2 index window; width and selected member unassigned')),Segment('FigureS2 exact tail',illumina.TRUSEQ_READ2[:-1],note='Source printed primer2 ends GATC; no extra terminal T imported from canonical arm')])
    return sc

def unknown_pcr(g,key,x,label,reagent):
    stock=g.input(key+'_supplies','Actual source PCR inventory; packaging/count/recipe unassembled',K.MOLECULAR)
    f=g.product(key+'_f','Published FigS2 primer1 whole ordering',K.OLIGO,state=source_pcr_primer(1))
    r=g.product(key+'_r','Published FigS2 primer2 fixed flanks and unassigned index window',K.OLIGO,state=source_pcr_primer(2))
    cat=g.product(key+'_enzyme',reagent,K.CATALYST)
    nts=g.product(key+'_nts','Source PCR nucleotide substrates; actual identities/doses unreviewed',K.MOLECULAR)
    buffer=g.product(key+'_buffer','Source PCR buffer components; actual composition unreviewed',K.CONDITION)
    g.react(key+'_inventory','Bookkeeping of actual PCR supplies; no prepared master mix asserted',S.GENERIC,((stock,R.STRUCTURAL),),(f,r,cat,nts,buffer),changes_molecule=False)
    y=g.product(key+'_product','Source PCR copies and original templates/unextended primers/remnants; whole library geometry unassembled',K.MOLECULAR)
    g.react(key,label,S.PCR,((x,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(cat,R.CATALYST),(nts,R.STRUCTURAL),(buffer,R.CONDITION)),(y,),geometry=PCRGeometry.unavailable(scene('Source full template/primer/end junction geometry unassembled'),reason='Known primer bodies do not establish full occupied library or all first-copy strands.'))
    return y

SOURCE='Published DLOHiC2018 DOI10.1038/s41588-018-0111-2;MOESM1 FigS1–5/Table1–2(originalp3,p11rendered). Mainpublisherpageispreview,notMethods.'
HALVES={'A':'GATAGCTCAAGAACTCCGAC','B':'CTTAGCTCAAGAACTCCGAC'}
def half_linker(which='A'):
    from chemdraw import revcomp
    if which not in HALVES:raise ValueError('sourceK562HindIII linkerAorBrequired')
    body=HALVES[which];sc=Scene();sc.strand('top',[Segment('half linker',body)],mod5='P')
    sc.anneal('bottom',[Segment('HindIII cohesive end','AGCT'),Segment('half complementary body',revcomp(body))],to='top',pair=('half complementary body','half linker'),mod5='P')
    sc.ligation_end('bottom','HindIII cohesive end',kind='cohesive',polarity='5-prime',owner='adapter');return sc
@dataclass(frozen=True)
class ContactTagClass:
    left: str
    right: str
    def __post_init__(self):
        if self.left not in HALVES or self.right not in HALVES:raise ValueError('Explicitsourcehalflinkerpairrequired')
    @property
    def linker_top(self):
        from chemdraw import revcomp
        return revcomp(HALVES[self.left])+HALVES[self.right]
    @property
    def heterotypic(self):return self.left!=self.right
    @property
    def genomic_tag_lengths(self):return (20,20)

def protocol_graph():
    g=ProtocolGraph('DLOHiC — selectedK562HindIII publishedSIroute',authority=SOURCE,representation='published_instance')
    cells=g.input('cells','SourceK562cells;notTHP1/BglIIorin-situMseIbranch',K.MOLECULAR,primary=True)
    fix_supplies=g.input('fix_supplies','Actual source double-crosslinking reagent inventory; packaging and doses unassigned',K.MOLECULAR)
    egs=g.product('fix_egs','Source EGS crosslinking reagent',K.MOLECULAR)
    fa=g.product('fix_fa','Source formaldehyde crosslinking reagent',K.MOLECULAR)
    g.react('fix_inventory','Bookkeeping of source EGS/formaldehyde contributions; not a master mix factory',S.GENERIC,((fix_supplies,R.STRUCTURAL),),(egs,fa),changes_molecule=False)
    fix_environment=g.input('fix_environment','Actual main crosslinking/lysis conditions unreviewed; selected SI establishes operations',K.CONDITION)
    x=g.product('fix_product','Crosslinked chromatin plus cellular/reagent companions',K.MOLECULAR)
    g.react('fix','Source EGS/formaldehyde double crosslinking and lysis interface',S.MAGIC,((cells,R.TEMPLATE),(egs,R.STRUCTURAL),(fa,R.STRUCTURAL),(fix_environment,R.CONDITION)),(x,),reason='Selected SI establishes molecular crosslinkers; full main dose/order/lysis recipe remains unreviewed.')
    tops=[]
    for which in ('A','B'):
        from chemdraw import revcomp
        a=g.input('linker_'+which.lower()+'_top','SourceK562HindIII'+which+'top;T4PNKphosphorylated',K.OLIGO,state=Scene());a.state.strand('oligo',[Segment('source half',HALVES[which])],mod5='P')
        b=g.input('linker_'+which.lower()+'_bottom','SourceK562HindIII'+which+'bottom;T4PNKphosphorylated',K.OLIGO,state=Scene());b.state.strand('oligo',[Segment('source bottom','AGCT'+revcomp(HALVES[which]))],mod5='P')
        h=g.product('half_'+which.lower(),'Sourceindependentlypreparedphosphorylatedhalflinker'+which,K.MOLECULAR,state=half_linker(which))
        g.react('anneal_'+which.lower(),'SourceT4PNKphosphorylation/95°C5mincool5°Cminto25°C;inputstatesarepostPNK',S.HYBRIDIZATION,((a,R.TEMPLATE),(b,R.STRUCTURAL)),(h,));tops.append(h)
    enzyme=g.input('hind_t7','SourceHindIII/T7DNAligase;fullconditionsunreviewed',K.CATALYST)
    linked=g.product('half_linked','Sourcehalf-linker-bearingchromatinpopulationplusunligatedDNA/freehalflinkers/enzymes;sourceA/Ballocation/mixinggeometryunreviewed',K.MOLECULAR)
    g.react('digest_half_ligate','SourceHindIII/T7simultaneousdigestionandcohesivehalf-linkerligation;A/Bbatchallocationunassembled',S.MAGIC,((x,R.TEMPLATE),(tops[0],R.STRUCTURAL),(tops[1],R.STRUCTURAL),(enzyme,R.CATALYST)),(linked,),reason='PublishedSIestablisheshalf-linkerligationandfourresultinglinkerclasses;preciseA/Ballocation/mixingtimingisoutsideinspectedmainMethods,notinvented.')
    mixed=g.product('pellet','Sourcecentrifugedlinker-bearingchromatinpellet;purityandA/Bcompositionunquantified',K.MOLECULAR);sup=g.product('supernatant','Actualsourceexcludedexcess-linkersupernatantwithpossiblecompanions',K.WASTE)
    g.react('spin','Sourcecentrifugationremovesexcesslinkersbeforenextligation',S.CLEANUP,((linked,R.TEMPLATE),),(mixed,sup));g.discard(sup,'Sourcephysicallyexcludedfraction;notallunligatedDNAassumedremoved.')
    x=op(g,'proximity',mixed,'Sourcein-gelT4DNAligaseproximityligation',S.LIGATION,reagent='SourceT4DNAligase;genomicintermolecularproductsnotpredicted',product='Sourcecontact-junctionclassesAA/AB/BA/BBplusunjoinedDNA/linkers/proteinsandreagents',catalyst=True)
    x=op(g,'reverse_crosslink',x,'SourceTable2reversecrosslinking;exactmainconditionsunread',product='DNAcontactproductsplusreleasedproteins/cellular/reagentcompanions;noinstantpurification')
    x=cleanup(g,'dna_purification',x,'SourceTable2DNApurificationafterreversecrosslinking')
    x=op(g,'mmei',x,'SourceMmeIdigestionreleasesapproximately80bpcontacttags',S.DIGESTION,reagent='SourceMmeI;FigS2two-base3primegenomiccohesiveends',catalyst=True,product='Selected20bpgenomic+40bplinkercore+20bpgenomiccontactclassesplusothercleavedDNA/remnants/undigestedmaterial')
    x=cleanup(g,'page',x,'SourceFigS3nativePAGErecoveryofapproximately80bpcontactDNA;notperfectpurity')
    adapter_stock=g.input('sequencing_adapter_inventory','Actual published FigureS2 PE-adapter inventory; packaging/amounts unassigned',K.MOLECULAR)
    adapters=[g.product('sequencing_adapter_'+str(i),'Published FigureS2 PE-adapter '+str(i)+' with 3prime NN cohesive end and 5prime OH',K.ADAPTER,state=source_adapter(i)) for i in (1,2)]
    g.react('sequencing_adapter_components','Bookkeeping of source adapter inventory; not a new annealing recipe or fresh dose',S.GENERIC,((adapter_stock,R.STRUCTURAL),),adapters,changes_molecule=False)
    cat=g.input('adapter_ligase','Sourcesequencingadapterligationreagentsunreviewed',K.CATALYST);y=g.product('adapter_product','Sourceadapter-ligatedcontactDNAplusunligatedcontactDNA/adapters/remnants;noA-tailorUSERopeninginvented',K.MOLECULAR)
    g.react('adapter_ligate','SourceNN3primecohesivesequencingadapterligation',S.LIGATION,((x,R.TEMPLATE),(adapters[0],R.STRUCTURAL),(adapters[1],R.STRUCTURAL),(cat,R.CATALYST)),(y,),geometry=AdapterLigationGeometry('known','cohesive',(source_adapter(1),source_adapter(2)),'Known source adapter local 3prime NN ends only; whole genomic tag/junction geometry and allocation unassembled.'))
    x=unknown_pcr(g,'pcr',y,'SourceTable2libraryamplification;FigS2primer2sample-indexmember/lengthunassigned','SourcePCRreagents/cyclecountunreviewed')
    x=qc(g,'early_qc',x,'Sourceearlyrestriction-basednoiseQCindependentlyallocated;mainassaytimingunassembled')
    g.sequence(x,SequencingPlan(limitation='SourceactualNGScontactreadoutestablished;exactexperimentalrun/readprimerlandingunreviewed. SIcost-tableSE100/150arenotassignedasanactualrun. Linkerclassificationandgenomicmappingtailadjustmentsarecomputational,notphysicalDNAremoval.'))
    g.add_source_item(SourceItem('pending','DefiningfullmainMethodsandhistoricalsequencingadapter/PCR/readgeometry','ExactNaturepublisherHTMLisabstract/figurepreview,notMethods;publishedchemistrySIreadincludingrenderedoligos. Nofamilyinheritanceoracceptedabsence.'))
    return g

def selftest():
    from checks import Check
    from chemdraw import revcomp
    c=Check();g=protocol_graph();classes=[ContactTagClass(a,b) for a in HALVES for b in HALVES]
    c('sourcefourlinkerclassesretained',len(classes),4)
    c('heterotypicQCclassificationdoesnotrewriteDNA',sum(x.heterotypic for x in classes),2)
    c('selectedcontactcorederivesfromtwoactualhalves',ContactTagClass('A','B').linker_top,revcomp(HALVES['A'])+HALVES['B'])
    c('allselectedcoresremain40bp',{len(x.linker_top) for x in classes},{40})
    c('sourcehalvesbothphosphorylated',tuple(x.mod5 for x in half_linker().strands.values()),('P','P'))
    c('MmeIremaindersreachphysicalPAGE',g.steps['page'].inputs[0].material,'mmei_product')
    c('restrictionQCnotsequencingsubstrate',g.terminals['early_qc_assay'].disposition,Disposition.QC)
    c('Both actual adapters derive from same source inventory',g.producer_of('sequencing_adapter_1').id==g.producer_of('sequencing_adapter_2').id=='sequencing_adapter_components',True)
    c('Both source adapter classes consumed in ligation',{'sequencing_adapter_1','sequencing_adapter_2'} <= {z.material for z in g.steps['adapter_ligate'].inputs},True)
    c('Exact FigureS2 primer1 fixture',source_pcr_primer(1).strands['primer'].segs[0].top,'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT')
    c('Exact FigureS2 primer2 tail fixture',source_pcr_primer(2).strands['primer'].segs[-1].top,'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC')
    c('PCR consumes same-inventory nucleotide class','pcr_nts' in {z.material for z in g.steps['pcr'].inputs},True)
    c('Adapter1 source whole top including random3primeNN', ''.join(x.top for x in source_adapter(1).strands['top'].segs),'ACACTCTTTCCCTACACGACGCTCTTCCGATCTNN')
    c('Adapter2 source whole top ordering', ''.join(x.top for x in source_adapter(2).strands['top'].segs),'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC')
    c('Crosslinkers are molecular source inventory products',g.materials['fix_egs'].kind==K.MOLECULAR and g.producer_of('fix_fa').id=='fix_inventory',True)
    c('gatedmainandgeometrypreventpromotion',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()

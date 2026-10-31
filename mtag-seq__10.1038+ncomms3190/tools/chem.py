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

def op(g,key,x,label,kind=S.MAGIC,*,reagent='source reagents; geometry unassembled',product=None,state=None,catalyst=False,extra_inputs=()):
    stock=g.input(key+'_stock',reagent,K.CATALYST if catalyst else K.CONDITION)
    out=g.product(key+'_product',product or label+'; molecular mixture/remnants unassembled',K.MOLECULAR,state=state)
    kw={'reason':'Source operation established; exact molecular ends, companion species and geometry unassembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
    g.react(key,label,kind,((x,R.TEMPLATE),(stock,R.CATALYST if catalyst else R.CONDITION))+tuple(extra_inputs),(out,),**kw);return out

def cleanup(g,key,x,label):
    stock=g.input(key+'_condition',label,K.CONDITION)
    out=g.product(key+'_kept','Source recovered preparation; recovery and carryover unquantified',K.MOLECULAR)
    waste=g.product(key+'_waste','Actual source-separated workup fractions',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(stock,R.CONDITION)),(out,waste));g.discard(waste,'Source physical cleanup, no perfect partition implied.');return out

def uncertain_ligation():return AdapterLigationGeometry('unavailable','unavailable',(scene('unreviewed exact ends and junction geometry'),),'Known source ligation, complete molecular geometry not assembled.')


SOURCE='Published Nature ncomms3190.html fullmainMethods883–904, publishedSI MOESM598 S19lines333–366; DOI10.1038/ncomms3190. SelectedNGShighthroughputclick route, notmicroarraylibrary.'
@dataclass(frozen=True)
class SelectedCpG:
    dna:str
    cytosine:int
    stage:str='unmodified'
    def __post_init__(self):
        if not self.dna or set(self.dna)-set('ACGT') or type(self.cytosine)is not int or self.cytosine<0 or self.dna[self.cytosine:self.cytosine+2]!='CG':raise ValueError('explicit selected unmodifiedCpG required, not an arbitrary genomicC')
        if self.stage not in {'unmodified','azide','SSbiotin','thiol-remnant'}:raise ValueError('source local chemicalclass required')
    def tag(self):
        if self.stage!='unmodified':raise ValueError('eM.SssI tag requires declaredunmodifiedsubstrateclass')
        return replace(self,stage='azide')
    def click(self):
        if self.stage!='azide':raise ValueError('DBCO needs installedazide')
        return replace(self,stage='SSbiotin')
    def reduce(self):
        if self.stage!='SSbiotin':raise ValueError('DTT needsdisulfidecapturetag')
        return replace(self,stage='thiol-remnant')

def protocol_graph():
    g=ProtocolGraph('mTAG-seq — selectedhighthroughputclickNGS',authority=SOURCE,representation='published_instance')
    x=g.input('genomic','SourcehumanbrainDNA; unmodifiedCpG positionsunassigned',K.MOLECULAR,primary=True)
    x=op(g,'sonicate',x,'SourceNGSsamples shearedaverage300bp',S.FRAGMENTATION,reagent='sourceNGSmain902; SI200bpCovarisoptimizationnot silentlysubstituted')
    inputdna=g.product('input_control','Actualsource shearedunenrichedinput aliquot',K.MOLECULAR);selected=g.product('label_input','ShearedDNA allocatedtoenrichment',K.MOLECULAR)
    g.react('input_allocate','Sourcepairedshearedinput versus mTAG-enrichedlibrary',S.SPLIT,((x,R.TEMPLATE),),(selected,inputdna));g.mark_terminal(Terminal(inputdna.id,Disposition.UNRESOLVED,'Separateactualinputlibrary followsownpreparation; no duplicated labelledDNA or missingcontrolfate.'));x=selected
    x=op(g,'end_repair',x,'SourcehistoricalFastDNAEndRepair20°C10min;75°C20min/stepcool',S.END_REPAIR,reagent='sourceFastDNAEndRepairenzymeMix',catalyst=True)
    x=cleanup(g,'repair_clean',x,'SourceZR96DNAcleanConcentrator5')
    em=g.input('tag_enzyme','Actual source Mediumclick eM.SssI90nM contribution, not amine or optimization conditions',K.CATALYST)
    donor=g.input('tag_donor','Actual source Ado6azide25uM donor contribution',K.MOLECULAR)
    x=op(g,'tag',x,'SourceMediumIntensityeM.SssI90nM/Ado6azide25µM37°C30or60min;80°C10min',S.CHEMICAL_CONVERSION,reagent='Source Mediumclick37C30or60min/80C10min and buffer conditions; no inferred identity for unnamed mainbuffer component',product='SelectedazidoalkylatedCpG DNA plus otherDNA/enzyme/cofactorremnants; tagnotglucosyl5hmC',extra_inputs=((em,R.CATALYST),(donor,R.STRUCTURAL)))
    addition=g.input('click_addition','Actual source5ul highthroughput addition containing4.4ulTris and0.6ul20mM DBCO-S-S-PEG3biotin inDMF',K.MOLECULAR)
    linker=g.product('click_linker','DBCO-S-S-PEG3biotin contribution from same actual5ul addition',K.MOLECULAR)
    solvent=g.product('click_solvent','Tris/DMF conditions from same actual5ul addition',K.CONDITION)
    g.react('click_components','Expose actual addition components, not fresh stocks or physical separation',S.GENERIC,((addition,R.STRUCTURAL),),(linker,solvent),changes_molecule=False)
    x=op(g,'click',x,'SourcehighthroughputDBCO-S-S-PEG3biotin37°C2h',S.CHEMICAL_CONVERSION,reagent='Sourcehighthroughput37C2h programme after actual5ul addition',product='CovalentSSbiotin-taggedunmodifiedCpGDNA plus unreactedDNA/freeprobe/reagents',extra_inputs=((linker,R.STRUCTURAL),(solvent,R.CONDITION)))
    x=cleanup(g,'click_clean',x,'SourceZR96DNAcleanConcentrator5/elute32µl')
    beads=g.input('beads','SourceMyOneStreptavidinT1 concentratedstock8µl/32µlDNA; mainC1route separate',K.SELECTION);bound=g.product('bound','SourcecapturedDNA/beadcomplex plus possible nonspecificspecies',K.MOLECULAR);removed=g.product('capture_wash','Physicallyremovedsupernatant and3highsalt/2lowsaltwashfractions',K.WASTE)
    g.react('capture','SourceT1capture2hRT900rpm andwash',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(beads,R.SELECTION)),(bound,removed));g.discard(removed,'Source actualremoved binding/wash fractions, not perfectcapture.')
    dtt=g.input('dtt_reductant','Actual source50mMDTT in20ul highthroughput release solution; not standaloneC1 200mM route',K.MOLECULAR)
    x=op(g,'dtt',bound,'Sourcehighthroughput50mMDTT/60mMTrispH7.8,1hmixing',S.CHEMICAL_CONVERSION,reagent='Source60mMTrispH7.8/1h900rpm programme',product='ReleasedDNA withthiol-sidechainremnant plusbeads/biotinlinkerproducts/DTT; allremainuntilmagnet',extra_inputs=((dtt,R.STRUCTURAL),))
    stock=g.input('magnet','Sourcemagneticplate',K.CONDITION);released=g.product('released','SourcecollectedDNA retainingcovalenttagremnant; not nativeunmodifiedC or reportedmethylationfraction',K.MOLECULAR);waste=g.product('bead_fraction','Source separatedbeads/linkerremnants and possibleretainedDNA',K.WASTE)
    g.react('collect','SourcephysicalmagneticseparationafterDTT',S.CLEANUP,((x,R.TEMPLATE),(stock,R.CONDITION)),(released,waste));g.discard(waste,'Sourceactualbeadremoval; noquantitativereleaseassumption.');x=released
    x=op(g,'ngs_transfer',x,'SourcepurifiedmTAG-enrichedDNA submittedtoBGI',reagent='ExactsourceNGStransfer/workupandhistoricalkitboundary; nopropanolormicroarrayPCRdefault')
    x=op(g,'atail',x,'SourceNGS3primeadenylation perIlluminaMultiplexingkit',S.TAILING,reagent='historicalkitadenylationenzymesunreviewed',catalyst=True)
    adapter=g.input('adapter','HistoricalIlluminaMultiplexingadapter; exactbody/endgeometryunreviewed',K.ADAPTER);ligase=g.input('ligase','SourceNGSadapterligationenzymesunreviewed',K.CATALYST);y=g.product('ligated','SourceadapterligatedDNA plusunligatedDNA/adapter/reagents; tagremnantsnoterased',K.MOLECULAR)
    g.react('ligate','SourceNGSadapterligation',S.LIGATION,((x,R.TEMPLATE),(adapter,R.STRUCTURAL),(ligase,R.CATALYST)),(y,),geometry=uncertain_ligation())
    beta=g.input('beta_me_reductant','Actual sourceNGS2mercaptoethanol10mM contribution; not microarray32or100mM',K.MOLECULAR)
    x=op(g,'beta_me',y,'SourceNGS10mM2mercaptoethanolRT10min BEFOREPCR',S.CHEMICAL_CONVERSION,reagent='SourceNGSroomtemperature10minbeforePCR programme',product='Source reducinglibrarymixture; originalmodifiedDNA/adapters/reagentsretained',extra_inputs=((beta,R.STRUCTURAL),))
    f=g.input('forward','Sourceuniqueindexprimerfamilyforward; bodyunreviewed',K.OLIGO);r=g.input('reverse','Sourceuniqueindexprimerfamilyreverse; bodyunreviewed',K.OLIGO);poly=g.input('polymerase','HistoricalNGSkitpolymeraseunreviewed',K.CATALYST);y=g.product('amp','Source14cycleNGSPCRcopies plus originalmodifiedDNA/primers/reagents; newDNAdoesnotuniversallyinheritcovalentmarks',K.MOLECULAR)
    g.react('pcr','Source14cycleNGSPCR; notmicroarraytwo-roundamplification',S.PCR,((x,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(poly,R.CATALYST)),(y,),geometry=PCRGeometry.unavailable(scene('Exacthistoricalkitindex/primer/fullstrandgeometryunassembled'),reason='Do notsubstituteavailablemicroarrayA25/Ad1/Ad2oligosforNGSkit.'))
    g.sequence(y,SequencingPlan(limitation='Publishedmain904actualIlluminaHiSeqpairedend90bp. NativeCpGtagging/captureisnotdirectbaseconversionorabsoluteCpGmethylationyield; mappedread-density/PCRdedupcomputational. Exactprimer/index/readgeometryunassembled.'))
    g.add_source_item(SourceItem('pending','HistoricalIlluminaMultiplexingSamplePreparationOligonucleotidekit','ActualsourceNGSkit/versionnotreviewed; noadapter/UMI/defaultpolymeraseinvented. Inputcontrolandmicroarray/aminebranchesoutofselectedgraph.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();native=SelectedCpG('AACGTT',2);released=native.tag().click().reduce()
    c('labelandrelease preservegenomicDNAsequenceandsite',(released.dna,released.cytosine),('AACGTT',2))
    c('DTTdoesnotrestorenativeunmodifiedC',released.stage,'thiol-remnant')
    try:native.click();bad=False
    except ValueError:bad=True
    c('capturetagcannotprecedeazide',bad,True)
    g=protocol_graph();c('physicalbeadseparationfollowsDTT',g.steps['collect'].inputs[0].material,'dtt_product')
    c('sourcepairedinputisindependent',g.terminals['input_control'].disposition,Disposition.UNRESOLVED)
    src=Source(sources_root()/'mtag-seq__10.1038+ncomms3190'/'ncomms3190.html.txt','get_sourcescanonicalNaturearticle')
    if have(c,src,label='actualNGSversusmicroarraylibrary'):
        t=src.path.read_text().split('mTAG-seq\n')[-1];c('publishedNGSparagraphspecifies14PCRcyclesandPE90','only 14 cycles' in t and '90' in t and 'Paired-end' in t,True)
    c('Source click reagent derives from actual5ul addition',g.steps['click_components'].inputs[0].material,'click_addition')
    c('NGS PCR predecessor consumes actual10mM reducing contributor',g.steps['beta_me'].inputs[2].material,'beta_me_reductant')
    c('Historical NGS kit remains uncovered','pending' in g.source_coverage,False)
    c('unreviewedhistoricalgeometrypreventspromotion',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()

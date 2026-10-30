"""Protocol-owned selected chemical classes; incomplete source-qualified graph."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Scene, Segment, revcomp, feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry

def scene(name,text):
    sc=Scene();sc.strand(name,[Segment(name,text,placeholder=True)]);return sc

def clean(g,key,x,label):
    condition=g.input(key+'_condition',label,K.CONDITION)
    kept=g.product(key+'_kept','Source recovered material; no perfect purity or quantitative partition assumed',K.MOLECULAR)
    waste=g.product(key+'_waste','Actually separated cleanup fractions with unquantified molecular carryover',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(condition,R.CONDITION)),(kept,waste));g.discard(waste,'Source physical separation, not perfect molecular classification.');return kept

def op(g,key,x,label,kind=S.MAGIC,*,reagent='source mixture',catalyst=False,product=None,state=None,chemical_state=None):
    stock=g.input(key+'_reagent',reagent,K.CATALYST if catalyst else K.CONDITION)
    y=g.product(key+'_product',product or label+'; all physical reaction companions retained until separation',K.MOLECULAR,state=state,chemical_state=chemical_state)
    kwargs={}
    if kind==S.MAGIC:kwargs['reason']='Source operation established; exact internal chemistry/geometry not assembled from unreviewed historical originals.'
    if kind==S.PCR:kwargs['geometry']=PCRGeometry.unavailable(Scene(),reason='Exact source template/primer endpoint geometry unassembled, no ordinary complete library asserted.')
    if kind in (S.LIGATION,S.PROXIMITY_LIGATION):kwargs['geometry']=AdapterLigationGeometry('unavailable','unavailable',(scene('unassembled','[source-known ligation; full substrate/end geometry unassembled]'),),'Unreviewed geometry; not a claim of unavailable disclosure.')
    inputs=[(x,R.TEMPLATE),(stock,R.CATALYST if catalyst else R.CONDITION)]
    if kind in {S.DIGESTION,S.END_REPAIR,S.EXTENSION,S.TAILING,S.PROXIMITY_LIGATION,S.REVERSE_TRANSCRIPTION} and not catalyst:
        enzyme=g.input(key+'_enzyme','Explicit catalytic component of source mixture: '+reagent,K.CATALYST);inputs.append((enzyme,R.CATALYST))
    if kind==S.REVERSE_TRANSCRIPTION:
        primer=g.input(key+'_primer','SourceRTprimer;exacthistoricalgeometryunassembled',K.OLIGO,state=scene(key+'_primer','[source RT primer]'));inputs.append((primer,R.STRUCTURAL))
    if kind==S.LIGATION:
        adapter=g.input(key+'_adapter','Source ligating molecular adapter; exact historical oligo/end inventory unassembled',K.ADAPTER,state=scene(key+'_adapter','[source adapter geometry unassembled]'));inputs.append((adapter,R.STRUCTURAL))
    if kind==S.PCR:
        for arm in ('f','r'):
            primer=g.input(key+'_primer_'+arm,'Source PCR primer family; exact oligo/endpoint unassembled',K.OLIGO,state=scene(key+'_primer_'+arm,'[source primer body/index geometry unassembled]'));inputs.append((primer,R.STRUCTURAL))
    g.react(key,label,kind,inputs,(y,),**kwargs);return y

def qc_split(g,key,x,label):
    assay=g.product(key+'_assay','Physically separate source analytical aliquot',K.QC)
    kept=g.product(key+'_retained','Unassayed retained preparation for downstream library work',K.MOLECULAR)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(assay,kept));g.mark_terminal(Terminal(assay.id,Disposition.QC,'Source analytical aliquot never feeds retained library.'));return kept
from illumina import TRUSEQ_SMALL_RNA_RA3,TRUSEQ_SMALL_RNA_RA5,P5,P7
from chemdraw import Construct
SOURCE='PublishedMOESM157SupplementaryTable1(originalPDFp10rendered),associatedAMPMC3560312 overview73–79/Methods129–141. Finalmainapplicabilityunreviewed;publishedSIoligoinventorydistinctfromunassembledwholelibrary.'
SOURCE_OLIGO_EXCEPTIONS={'capture_oligo':'Generic guard does not expand source CU5T45 plus prose into actual DNA T45/RNA U5. Coordinator independently read associated original PMC3560312 Results73 and Methods131 on2026-10-29:45thymidines at5primeportion,5uridines at3primeportion,5primebiotin. One capture-oligo macro exemption only; all actual sugars/width/orientation retained. Whole descriptor Source/have check present. Not an absence audit.'}
UMI=feature('reads3.random4','umi','random',note='PublishedTable1footnote:4randomnucleotidesfacilitateclustersanddistinguishRNAfragments/removePCRredundancy.')
INDEX=feature('reads3.i7','sample_index','fixed',note='PublishedbracketCGTGATmultiplexingindexexample;notstudysampleassignment.')
@dataclass(frozen=True)
class CaptureOligo:
    def scene(self):
        sc=Scene();sc.strand('capture',[Segment('DNA45','T'*45),Segment('RNA5','U'*5)],mod5='biotin');return sc
    @property
    def ordered_sugars(self):return (('DNA','T'*45),('RNA','U'*5))
    def selected_rna_release(self,body,tail,remain):
        if not body or set(body)-set('ACGU') or type(tail)is not int or type(remain)is not int or not 0<=remain<=tail:raise ValueError('Explicitselectedtailcleavageboundaryrequired,notperfectRNaseHprediction')
        return (body+'A'*remain,'A'*(tail-remain))
@dataclass(frozen=True)
class BlockedAdapter:
    direction:str='reverse'
    identifier:str='NNNN'
    def __post_init__(self):
        if self.direction not in {'forward','reverse'} or len(self.identifier)!=4 or set(self.identifier)-set('ACGTN'):raise ValueError('Actualsourceorientationand4ntidentifierrequired')
    def construct(self):
        return Construct(([Segment('identifier',self.identifier,placeholder='N' in self.identifier,feature=UMI),Segment('body',revcomp(TRUSEQ_SMALL_RNA_RA5))] if self.direction=='reverse' else [Segment('body',TRUSEQ_SMALL_RNA_RA3)]))
    def scene(self):
        sc=Scene();sc.strand('adapter',self.construct().segments+[Segment('terminal_ddC','C')],mod5='rApp',mod3='2prime3primedideoxy (blocked)');return sc
    @property
    def terminal_block(self):return ('2prime3primedideoxycytidine','C')
    def complementary_reference(self):
        # Abstract DNA-complement reference only, not the operational RT product: primer landing remains unassembled.
        return Construct([Segment('ddC_complement','G')]+[Segment('complement_'+s.name,revcomp(s.top),placeholder=s.placeholder,feature=s.feature) for s in reversed(self.construct().segments)])

def rna_adapter(direction):return revcomp(TRUSEQ_SMALL_RNA_RA3).replace('T','U') if direction=='forward' else TRUSEQ_SMALL_RNA_RA5.replace('T','U')
def rt_primer(direction):return 'CTAGCAGCCTGACATCTTGAGACTTG' if direction=='forward' else 'GCCTTGGCACCCGAGAATTCCA'
def reverse_pcr_primer():return Construct([Segment('P7',P7),Segment('i7','CGTGAT',feature=INDEX),Segment('fixed','GTGACTGGAGTT'+revcomp(TRUSEQ_SMALL_RNA_RA3))])
def graph_for_direction(direction='reverse'):
    if direction not in {'forward','reverse'}:raise ValueError('ActualpublishedTable1orientationrequired')
    g=ProtocolGraph('3′READS — '+direction+'source-tablebranch',authority=SOURCE,representation='published_instance')
    x=g.input('total_rna','SourceextractedtotalRNA;study-specificsamplecontextunassigned',K.MOLECULAR,primary=True)
    x=op(g,'poly_a',x,'SourceonePolyAPuristMAGselection',reagent='Historicalsourcekitinternalsunreviewed',product='SourcepolyA-enrichedpreparationwithkitremnantfatesunassembled')
    x=op(g,'fragment',x,'SourceAmbionRNAfragmentation70C5min',S.FRAGMENTATION,reagent='SourceRNAfragmentationbuffer',product='AllRNAfragments/shorttails/remnantsuntilsourcecaptureselection')
    oligo=CaptureOligo();a=g.input('capture_oligo','Source5biotin-DNA45T/RNA5Uchimera',K.OLIGO,state=oligo.scene(),chemical_state=oligo);beads=g.input('beads','SourceMyOnestreptavidinC1beads',K.SELECTION);coated=g.product('coated_beads','Sourcebiotin-oligoassociatedbeads;freeoligocarryoverunquantified',K.SELECTION)
    g.react('coat_beads','Sourcecaptureoligoboundthrough5biotin',S.BINDING,((a,R.TEMPLATE),(beads,R.SELECTION)),(coated,))
    bound=g.product('bound','SourcepolyAtail/captureoligo/beadcomplexesplusnonidealcarryover',K.MOLECULAR);unbound=g.product('capture_unbound','Sourcephysicallyseparatedunboundcapturefraction',K.WASTE);g.react('capture','SourceRT1hRNAcapture/low-saltformamidewash',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(coated,R.SELECTION)),(bound,unbound));g.discard(unbound,'Sourcecapturedfractionselected,notperfectpolyAtailclassifier.')
    x=op(g,'rnase_h',bound,'SourceRNaseH37C1htailcleavage/release',S.DIGESTION,reagent='SourceRNaseH5U',catalyst=True,product='ReleasedRNAwithselectedterminalAspluscleavedRNAresidues/beads/chimericcaptureoligos/reagents;exactcleavageboundarynotfixedtofiveAs')
    released=g.product('released','SourcecollectedelutedRNAfraction',K.MOLECULAR);remains=g.product('capture_remains','Remainingbead/captureoligo/tailfragmentcompanions;laterfateunassigned',K.MOLECULAR);g.react('collect','SourceelutedRNAphysicalrecovery',S.CLEANUP,((x,R.TEMPLATE),),(released,remains));g.mark_terminal(Terminal(remains.id,Disposition.UNRESOLVED,'Nofabricatedbead/remnantstorageordiscard.'))
    x=clean(g,'organic',released,'Sourcephenol/chloroformandEtOH')
    x=op(g,'pnk',x,'SourceT4kinase5endphosphorylation',S.END_REPAIR,reagent='SourceT4kinase',catalyst=True,product='5phosphorylatedRNAclassesplusenzymaticremnants;notadjacentRNAseqSAP3enddephosphorylation')
    x=clean(g,'rneasy',x,'SourceRNeasycleanup')
    adapter=BlockedAdapter(direction);a=g.input('three_adapter','PublishedTable1'+direction+'3adapterrApp/terminalddC',K.ADAPTER,state=adapter.scene(),chemical_state=adapter);cat=g.input('rnl2','SourcetruncatedT4RNLII',K.CATALYST);y=g.product('three_ligated','RNA/3adapterclassesretainactualterminalddCplusunligatedRNA/freeadapter/remnants',K.MOLECULAR)
    g.react('three_ligate','Sourcepreadenylated3adapterligation',S.LIGATION,((x,R.TEMPLATE),(a,R.STRUCTURAL),(cat,R.CATALYST)),(y,),geometry=AdapterLigationGeometry('unavailable','unavailable',(a.state,),'Actualpublishedadapterinventoryknown;completeRNAend/junctiongeometryunassembled.'))
    a=g.input('five_adapter','PublishedTable1'+direction+'RNA5adapter',K.ADAPTER,state=Scene());a.state.strand('five',[Segment('sourceRNA',rna_adapter(direction))]);cat=g.input('rnl1','SourceT4RNLI',K.CATALYST);z=g.product('five_ligated','DoubleadapterRNAclassesplusoldRNA/freeadapters/enzymaticremnants',K.MOLECULAR)
    g.react('five_ligate','Source5adapterligationafter3adapter',S.LIGATION,((y,R.TEMPLATE),(a,R.STRUCTURAL),(cat,R.CATALYST)),(z,),geometry=AdapterLigationGeometry('unavailable','unavailable',(a.state,),'ActualRNAadapterinventoryknown;wholeligatedorientationunassembled.'))
    p=g.input('rt_primer','PublishedTable1'+direction+'RTprimer',K.OLIGO,state=Scene());p.state.strand('primer',[Segment('sourceRT',rt_primer(direction))]);cat=g.input('ssiii','SourceSuperScriptIII',K.CATALYST);x=g.product('rt_product','SourcecDNAplusoriginalRNA/adapters/templates/unextendedprimer/remnants;wholefirstcopyunassembled',K.MOLECULAR)
    g.react('rt','SourceRTwithSSIII',S.REVERSE_TRANSCRIPTION,((z,R.TEMPLATE),(p,R.STRUCTURAL),(cat,R.CATALYST)),(x,))
    f=g.input('pcr_f','PublishedTable1forwardPCRprimer',K.OLIGO,state=Scene());f.state.strand('f',[Segment('sourcePCR',P5+TRUSEQ_SMALL_RNA_RA5[:-5])]);r=g.input('pcr_r','PublishedTable1reversePCRprimerexampleCGTGAT',K.OLIGO,state=reverse_pcr_primer());cat=g.input('phusion','SourcePhusion',K.CATALYST);y=g.product('pcr_product','Source12cycleDNAcopiesplusoriginalDNA/RNAandadapter/remnantcompanions',K.MOLECULAR)
    g.react('pcr','SourcePhusion12cyclePCR',S.PCR,((x,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(cat,R.CATALYST)),(y,),geometry=PCRGeometry.unavailable(scene('geometry','[Table1oligosknown;wholeorientedlibrarygeometryunassembled]'),reason='Sourcebranchlabelsdonotjustifyinventingprimerlandingorasilentlyrepairedlibrary.'))
    lib=g.product('library','Sourcephysical3READSDNAsequencingpreparation;fulltemplategeometryunassembled',K.LIBRARY);g.react('identity','Sourcefinalpreparationidentity',S.GENERIC,((y,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='AssociatedmainGAIIx1x72nt;publishedTable1forward/reverseinventoriesretainedbutcompleteRT/PCRlandingunassembled. SourceN4israndomfragmentidentifier/clusteringaid;PASSclassificationdoesnotphysicallyselectperfectpolyAsites.'))
    g.add_source_item(SourceItem('pending','Finalmainapplicability/wholeTable1orientation/linkedhistoricalsources','AssociatedmainandpublishedSIreadseparately;noassumedfullversionequivalence. AdjacentCLIP/RNAseqrecipesnotimported.'));return g
def protocol_graph():return graph_for_direction()
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();a=CaptureOligo();c('SourcecaptureoligoretainsDNAandRNAsugars',a.ordered_sugars,(('DNA','T'*45),('RNA','U'*5)))
    kept,lost=a.selected_rna_release('ACGU',20,3);c('SelectedRNaseHconservestailresidues',kept+lost,'ACGU'+'A'*20)
    b=BlockedAdapter('reverse','ACGT');c('ActualddCremainsnucleotideplusblock',b.terminal_block,('2prime3primedideoxycytidine','C'))
    c('ComplementaryreferenceincludesddCbasewithoutcopyingblock',b.complementary_reference().top()[0],'G')
    c('Randomidentifieridentitypreservedincomplementaryreference',b.complementary_reference().segments[-1].feature.id,b.construct().segments[0].feature.id)
    g=protocol_graph();c('5adapterligationfollows3adapterligation',g.steps['five_ligate'].inputs[0].material,'three_ligated')
    c('CaptureRNAreleasedbeforephosphorylation',g.steps['pnk'].inputs[0].material,'organic_kept')
    c('AdjacentstandardRNAseqSAPnotimported','sap' not in g.steps,True)
    c('Unassembledprimergeometrypreventsfinishedclaim',bool(g.publication_problems()),True)
    src=Source(sources_root()/'3-reads__10.1038+nmeth.2288'/'nmeth.2288_41592_2013_BFnmeth2288_MOESM157_ESM.pdf.txt','get_sources --slugs3-reads__10.1038+nmeth.2288;originalp10rendered')
    if have(c,src,label='publishedTable1'):
        t=''.join(src.path.read_text().split());c('Wholepublishedreverseadapterwithterminalmods','5’-rApp/NNNN'+revcomp(TRUSEQ_SMALL_RNA_RA5)+'/3ddC' in t,True)
        c('Publishedidentifierfunctionandindexexplicit','eliminateredundantreadscausedbyPCR' in t and 'nucleotidesinbracketaretheindex' in t,True)
    main=Source(sources_root()/'3-reads__10.1038+nmeth.2288'/'nmeth.2288_PMC3560312.xml.txt','get_sources --slugs3-reads__10.1038+nmeth.2288;associatedAMoverview73/Methods131')
    if have(c,main,label='sourcecapturecomposition'):
        text=main.path.read_text();c('Wholecapturedescriptorestablishes45DNA5RNAand5biotin','45 thymidines (Ts) at the 5′ portion and 5 uridines (Us) at the 3′ portion' in text and 'through biotin at its 5′ end' in text,True)
    c.report()
if __name__=='__main__':selftest()

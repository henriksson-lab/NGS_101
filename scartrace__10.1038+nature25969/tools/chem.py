from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Scene, Segment, feature, workflow_from_protocol_graph, workflow_panel
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan
from page import head,info
PAGE_FILE='chem.html'

def operations(g):
    def op(key,x,label,reagent=None,state=None,kind=S.GENERIC,reagent_kind=K.CATALYST,reagent_role=R.CATALYST):
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
def incomplete(g,v):
    # Literal independently reviewed source scope, not graph enumeration. Linked
    # transcriptome library manufacture remains unreviewed and unaccounted.
    g.declare_reviewed_scope(('cells','sort_reagent','sort','freeze_reagent','freeze','proteinase_reagent','proteinase','scar_pcr1_reagent','scar_pcr1','scar_primer','scar_pcr2','pool_wells','oil_separation_selection','oil_separation','scar_adapter_pcr_reagent','scar_adapter_pcr','scar_adapter_pcr_product'),SOURCE)
    if v=='joint':
        g.declare_reviewed_scope(('rt_primer','rt_catalyst','rt','second_strand_reagent','second_strand','allocate','ivt_reagent','ivt','rna_library_product'),SOURCE)
    g.add_source_item(SourceItem('unassembled','Unreviewed linked originals and full molecular geometry',LIMIT));return g

TITLE='ScarTrace — genomic GFP scars and same-cell transcriptomes'
SOURCE='Published institutional mirror ofNature25969 Methods1385–1402,Fig1d,andpublishedSI MOESM6/7 originalworkbooks'
LIMIT='Private partial graph. Published main operational order and two source primer tables inspected. Direct detailed ProtocolExchange2018.017 and SORT-seq predecessor remain unreviewed;full nested PCR,scar primer assembly,IVT continuation,library polymers and read-primer geometry unassembled. No LINNAEUS/scGESTALT recipe inherited.'
RT_PREFIX='GCCGGTAATACGACTCACTATAGGCCTTGGCACCCGAGAATTCCA'
SCAR_ARM='AGGACGACGGCAACTACAAG'
def rt_primer():return Construct([Segment('source RT handle',RT_PREFIX),Segment('UMI','N'*8,placeholder=True,feature=feature('scartrace.rna_umi','umi','random')),Segment('RNA cell barcode example','CGTCTAAT',feature=feature('scartrace.rna_cell','cell_barcode','fixed')),Segment('polyT','T'*25),Segment('anchor','VN',placeholder=True)])
def recorder():
    sc=Scene();sc.strand('scar',[Segment('GFP genomic scar','[source-dependent GFP indel allele]',placeholder=True,feature=feature('scartrace.scar','lineage_barcode','unknown'))]);sc.note('scar','Histone-GFP transgene copies;multiple scar alleles percell possible;not a unique cell barcode or predicted indel.');return sc

from dataclasses import dataclass
@dataclass(frozen=True)
class ScarPrimerInventory:
    barcode:str='GTCACCGG'
    arm:str=SCAR_ARM
    def rows(self):
        sc=Scene()
        sc.strand('barcode_component',[Segment('Scar well barcode inventory',self.barcode,feature=feature('scartrace.scar_cell','cell_barcode','fixed'))],label='Component of ONE unassembled source primer;not a free oligo')
        sc.strand('arm_component',[Segment('GFP arm inventory',self.arm)],label='Same source primer inventory;order/junction not reconstructed')
        sc.note('barcode_component','Separate inventory rows do not declare independent supplied molecular species or a validated covalent junction.')
        return sc.rows()

def variant_graph(v):
    if v not in ('joint','dna_only'):raise ValueError(v)
    g=ProtocolGraph(TITLE+' / '+v,authority=SOURCE,representation='illustrative_example',example_description='One symbolic well and selected source primer inventory illustrate the joint/gDNA-only source workflow. No actual384-cell pool, sample identity or edited allele is assigned.');op,clean=operations(g)
    cells=g.input('cells','Source DAPI-negative single cells after recorder editing;one symbolic well',K.BIOLOGICAL,primary=True,state=recorder())
    x=op('sort',cells,'FACS into384well plate,mineral oil,dNTPs,spike-ins andRNase inhibitor','Actual source mineral-oil/dNTP/spike-in/RNase-inhibitor allocations; composition/stockcounts unassembled',reagent_kind=K.MOLECULAR,reagent_role=R.STRUCTURAL)
    x=op('freeze',x,'Spin plate,storeminus80C;source65C5minlysis','Source65C5min lysis and storage/spin programme',reagent_kind=K.CONDITION,reagent_role=R.CONDITION)
    if v=='joint':
        primer=g.input('rt_primer','PublishedCELSeq2-BC-001 member;wellassignmentnotexperimental sampleidentity',K.OLIGO,state=rt_primer())
        out=g.product('rt_product','Well mixture:RNA/cDNAhybrids,genomicDNA,freeprimers/reagents',K.MOLECULAR)
        rt=g.input('rt_catalyst','Actual source reverse-transcription catalyst allocation at42C; identity/dose require detailed original, not an imported SORT-seq recipe',K.CATALYST)
        g.react('rt','Source reverse transcription42C',S.REVERSE_TRANSCRIPTION,((x,R.TEMPLATE),(primer,R.STRUCTURAL),(rt,R.CATALYST)),(out,));x=op('second_strand',out,'Source second-strand synthesis16C;genomicDNAremainscopresent','Source second-strand mixture unassembled')
    x=op('proteinase',x,'ProteinaseK makesgenomicDNAaccessible;retainwellcDNAifpresent','Source proteinaseK')
    x=op('scar_pcr1',x,'FirstsourceGFP-targetedPCR;gDNAandotherwellmaterialremaincopresent','Source firstPCRprimers/enzymes unassembled')
    primer=g.input('scar_primer','One source barcode-primer reagent;Table2 component inventory,whole junction unassembled',K.OLIGO,state=ScarPrimerInventory())
    out=g.product('scar_pcr2_product','NestedGFPproductsbearingwellscarbarcodesplusotherwellmaterial',K.MOLECULAR)
    g.react('scar_pcr2','SecondnestedPCRintroducesuniquescarbarcodeperwell',S.GENERIC,((x,R.TEMPLATE),(primer,R.STRUCTURAL)),(out,))
    # One representativewell is displayed;populationpooling doesnot inventadditionalindependentmolecularinputs.
    x=op('pool_wells',out,'Poolallsourcewells;symbolicwellrepresentspopulation,notanewbarcodeassociationreaction')
    x=clean('oil_separation',x,'Separateaqueouspoolfrommineraloil;retainaqueousmaterial')
    if v=='joint':
        scar=g.product('scar_aliquot','Scar-libraryaliquotofpooledaqueousmixture',K.MOLECULAR)
        gene=g.product('rna_aliquot','Transcriptome-libraryaliquotofsamepooledmixture',K.MOLECULAR)
        g.react('allocate','Splitpooledmaterialintotwolibrarypreparations;notDNA/RNApurification',S.SPLIT,((x,R.TEMPLATE),),(scar,gene))
        gene=op('ivt',gene,'SourceFig1dIVTonpooledcDNA;SORT-seqdownstreammolecularstepsunreviewed','SourceSORT-seqIVTcomponents')
        gene=op('rna_library',gene,'SourceSORT-seqtranscriptomelibrarycontinuation;exactlinkedrecipeunreviewed','HistoricalSORT-seqcomponents')
        end(g,gene,'PublishedPE75NextSeq;RNARead1cellbarcode/UMI,Read2transcript. Wholelibrary/readprimergeometryunassembled;computationalassociationnotphysicalfusion.')
    else:scar=x
    scar=op('scar_adapter_pcr',scar,'PCRintroducingIlluminaTruSeqadapters;noinventedintermediatepurifications','SourceadaptorPCRprimers/enzymesunassembled')
    end(g,scar,'PublishedPE75NextSeq;scarwellbarcodeandGFPscarreadout. Wholelibrary,indexcycles/readprimergeometryunassembled;scarbarcode≠inheritedscar.')
    return incomplete(g,v)
def protocol_graph():return variant_graph('joint')
def render_page():return head(TITLE)+'<main class="wrap"><h1>'+TITLE+'</h1>'+info(LIMIT)+info('ILLUSTRATIVE EXAMPLE — one symbolic well represents protocol stages, not a verified384-cell pool or an edited sample allele. SourceFigure1d schematizes80percenttranscriptome/20percentscar; mainMethods does not specify transfer quantities. The mixture split is not DNA/RNA purification.')+'<p><a href="01_scartrace.html">Research notes</a></p>'+''.join(workflow_panel(workflow_from_protocol_graph(protocol_graph(),x)) for x in ('scar_adapter_pcr_product','rna_library_product'))+'</main>'
def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();g=protocol_graph();d=variant_graph('dna_only')
    c('Scar barcode and arm are one input reagent',tuple(i.material for i in g.steps['scar_pcr2'].inputs),('scar_pcr1_product','scar_primer'))
    c('Symbolic well is explicitly illustrative',g.representation,'illustrative_example')
    c('Original detailed and SORT-seq obligation stays uncovered','unassembled' in g.source_coverage,False)
    c('Bothroutesprivate',g.publishable or d.publishable,False)
    c('DNAgDNArouteomitsRT','rt' in d.steps,False)
    c('Proteinaseaftersecondstrand',g.steps['proteinase'].inputs[0].material,'second_strand_product')
    c('Poolbeforeoilseparation',g.steps['oil_separation'].inputs[0].material,'pool_wells_product')
    c('Splitafterphysicaloilremoval',g.steps['allocate'].inputs[0].material,'oil_separation_kept')
    p=sources_root()/Path(__file__).resolve().parents[1].name/'nature25969_41586_2018_BFnature25969_MOESM6_ESM.xlsx.txt'
    if p.exists():c('SourceRTprimerwholetranscription',''.join(seg.top for seg in rt_primer().segments) in p.read_text(),True)
    else:c.skip('RTprimertranscription','Externaloriginalabsent')
    table=sources_root()/Path(__file__).resolve().parents[1].name/'nature25969_41586_2018_BFnature25969_MOESM7_ESM.xlsx.txt'
    if table.exists():c('Published scar primer components,without joining their junction','\t'+ScarPrimerInventory().barcode+'\t'+ScarPrimerInventory().arm in table.read_text(),True)
    else:c.skip('Scar primer component transcription','External original absent')
    c.report()

"""Published easyCLIP workbook case; source alternatives and full mixtures stay explicit."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
from chemdraw import Construct,Segment,Scene,feature
import illumina as il
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SequencingPlan,SourceItem
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry
SOURCE='PublishedPMC7946914 Methods229–231; published SupplementaryData1/MOESM3 workbook LysisToRT/RTtolibraries/L5linkers/Otheroligos; SI/MOESM1 read-processing paragraph. Three originals only.'
L5_PRECURSOR_HANDLE='CTTGTGTGTGTGAAGCGAAGGGTA'
P6_HANDLE='CTTGTGTGTGAAGCGAAGGGTA'
P6=il.P5+'TCTTTCCC'+P6_HANDLE
P3=il.PE_PCR_PRIMER_2
RT=il.PE_READ2_PRIMER[il.PE_READ2_PRIMER.index('GGCATTCCT'):]
L5_MACRO='/5AzideN/'+L5_PRECURSOR_HANDLE+'TTAGGC'+'rN'*7
L3_MACRO='CTGAGATCGGAAGAGCGGTTCAG'+'A'*12+'/iAzideN/'+'A'*12+'/3Bio/-3′'

def precursor_l5():
    """As-written precursor inventory, not the complete dye-conjugated ligation reagent."""
    s=Scene();s.strand('workbook precursor',[Segment('source DNA handle',L5_PRECURSOR_HANDLE),Segment('source L5 sample tag','TTAGGC',feature=feature('easy_l5_3','sample_index','fixed')),Segment('source RNA random window','N'*7,placeholder=True,feature=feature('easy_l5_random','umi','random'),note='Workbook7rN versus SI read-processing hexamer; RNA synthesis pattern, not seven ordinary DNA bases')],mod5='source azide precursor; final fluorophore coupling unassembled');return s

def unknown(label):
    s=Scene();s.strand('mixture',[Segment('unassembled','['+label+']',placeholder=True)]);return s

def op(g,key,x,label,kind=S.MAGIC,reagent='Source conditions; full mixture unassembled'):
    role=R.CATALYST if kind==S.DEPROTEINIZATION else R.CONDITION
    stock=g.input(key+'_stock',reagent,K.CATALYST if role==R.CATALYST else K.CONDITION)
    y=g.product(key+'_product',label+'; physical substrate, products and reagent remnants retained until source separation',K.MOLECULAR)
    kw={'reason':'Source operation known; complete substrate/product/condition inventory unassembled.'}if kind==S.MAGIC else {}
    g.react(key,label,kind,((x,R.TEMPLATE),(stock,role)),(y,),**kw);return y

def clean(g,key,x,label):
    r=g.input(key+'_reagent',label,K.SELECTION);y=g.product(key+'_kept','Source physically retained fraction; no perfect recovery/purity assertion',K.MOLECULAR);w=g.product(key+'_removed','Actual separated source wash/unbound fractions',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(r,R.SELECTION)),(y,w));g.discard(w,'Physical separation only; no molecular classification or perfect selection prediction.');return y

def ligate(g,key,x,label,adapter_label):
    a=g.input(key+'_adapter',adapter_label,K.ADAPTER);e=g.input(key+'_enzyme','Source T4RNA ligaseI/ATP/PEG400/DTT/RNaseinhibitor; exact reaction mixture/end geometry unassembled',K.CATALYST)
    y=g.product(key+'_product','Ligated/partly ligated/unligated RNA-protein classes, unused adapters and all reaction remnants',K.MOLECULAR)
    geom=AdapterLigationGeometry('unavailable','unavailable',(unknown('Source RNA/DNA/dye/end inventory not fully assembled'),),'Actual ligation known; incomplete model geometry is not a publication-absence claim.')
    g.react(key,label,S.LIGATION,((x,R.TEMPLATE),(a,R.STRUCTURAL),(e,R.CATALYST)),(y,),geometry=geom);return y

def protocol_graph():
    g=ProtocolGraph('easyCLIP — selected published workbook3prime-first path',authority=SOURCE,representation='published_instance')
    x=g.input('cells','Source RBP-expressing cultured cells; exact study sample/amount unassigned',K.BIOLOGICAL,primary=True)
    x=op(g,'uv',x,'Workbook ice/PBS/UV0.3Jcm2; source crosslinked and uncrosslinked populations',S.CROSS_LINKING)
    x=op(g,'lysis',x,'Workbook CLIPlysis/freeze/thaw/sonicate20–30s in7sbursts; main6x5s alternative separate')
    x=clean(g,'clarify',x,'14krcf10min4C; transfer source supernatant')
    x=clean(g,'ip',x,'AntiHA orantibodyProteinG beadIP1h4C;highstringency/highsalt/lowsalt washes')
    e=g.input('rnase','WorkbookRNaseONE0.01–0.1Uul,30C3–10min;main0.025Uul8–12min separate',K.CATALYST);y=g.product('digested','Cleaved RNAs, protein complexes and reagent remnants; no predicted exact fragments',K.MOLECULAR);g.react('rnase_digest','Selected workbook on-bead RNA digestion',S.DIGESTION,((x,R.TEMPLATE),(e,R.CATALYST)),(y,));x=clean(g,'rnase_wash',y,'RemoveRNasemixture/addcoldHighStr then2NT2washes')
    x=op(g,'pnk3',x,'Source PNK3prime dephosphorylation37C30min; ATP-freepH6.8mix')
    x=clean(g,'pnk3_remove',x,'Magnetic removalPNKreaction; no extra wash before3ligation')
    x=ligate(g,'lig3',x,'Workbook3prime ligation16C2h-to-overnight','Source L3barcoded/IR800labelled reagent; precursor azide/3biotin inventory known, final conjugate/end state unassembled')
    x=clean(g,'lig3_wash',x,'OneHighSaltorHighStr and2NT2 washes')
    x=op(g,'pnk5',x,'Source5primePNK/ATP37C15min')
    x=clean(g,'pnk5_wash',x,'Source2NT2 washes')
    x=ligate(g,'lig5',x,'Workbook5prime ligation16C2h','Source L5/IR680labelled chimeric reagent;7rNprecursor window versus hexamerprocessing retained; no precursor-dye equivalence assumed')
    x=clean(g,'postlig_wash',x,'Source1–2NT2 washes;stepdownNT1. Optional10percentWB branch not selected')
    others=g.input('other_samples','Source independent indexed preparations typically combined before gel; identities/amounts unassigned',K.MOLECULAR,primary=True);y=g.product('pooled','Physical indexed-sample pool; not a covalent assembly',K.MOLECULAR);g.react('sample_pool','Main both-ended libraries typically pooled before single gel lane',S.POOL,((x,R.TEMPLATE),(others,R.TEMPLATE)),(y,))
    x=op(g,'lds',y,'Workbook1.6xLDS75C15min')
    x=clean(g,'membrane',x,'NuPAGE180V30–50min/transfernitrocellulose400mA27min;imageandcutselectedminimalbandandhigher. Main200V45min/400–500mA25min alternatives separate')
    x=op(g,'pk',x,'Workbook375ulPKbuffer+25ulProteinaseK50–55C45min',S.DEPROTEINIZATION)
    x=clean(g,'dt_capture',x,'TransferPKextract to20uloligodTbeads/600ulBIB,4C20min;BIB/NT2/PBSwashes. Mainroomtemperature orKCl alternatives unselected')
    primer=g.input('rt_primer','Source P3.tall.V4 RTprimer',K.OLIGO,state=Construct([Segment('source RTprimer',RT)]))
    y=g.product('eluted','Source13.3ul transferred RNA/primer supernatant',K.MOLECULAR);rem=g.product('elution_remainder','Beads and untransferred source liquid; exact composition/fate unassigned',K.MOLECULAR);c=g.input('elution_condition','Workbook12.4ulwater+3ul1uMprimer;95C3min,13.3ulsupernatant transfer/4C5min',K.CONDITION)
    g.react('elute','Source primer-assisted heatelution/physical transfer',S.STRAND_SEPARATION,((x,R.TEMPLATE),(primer,R.STRUCTURAL),(c,R.CONDITION)),(y,rem));g.mark_terminal(Terminal(rem.id,Disposition.UNRESOLVED,'No complete consumption or perfect capture assumption.'))
    rt=g.input('ssiv','SourceSSIV/dNTP/DTT/RNaseinhibitor mixture; workbook20ul total, subsection calls6.7ul addition',K.CATALYST);x=g.product('rt_product','cDNA plus original RNA/primer/RT reagents; directly toPCR',K.MOLECULAR);g.react('rt','Selected workbookparagraph53C40min/55C10min/80C10min; adjacent50/10/5 table unresolved',S.REVERSE_TRANSCRIPTION,((y,R.TEMPLATE),(primer,R.STRUCTURAL),(rt,R.CATALYST)),(x,))
    branches=[g.product('pcr_input_'+str(i),'Source PCR aliquot≤3.3ul to47ul mix; exact volumes unassigned',K.MOLECULAR)for i in range(7)];rest=g.product('rt_remainder','Source remainingRT stock; allocation not stoichiometrically assigned',K.MOLECULAR);g.react('pcr_allocate','Workbook7PCRtubes from8reactionmastermix; exact RT allocation unassigned',S.SPLIT,((x,R.TEMPLATE),),(*branches,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Seven upper bounds do not assert seven maximal aliquots from20ul.'))
    spare=g.input('unused_mastermix','Source8-reaction PCRmastermix prepared/7tubes distributed; remaining stock volume/composition not assembled',K.MOLECULAR)
    g.mark_terminal(Terminal(spare.id,Disposition.UNRESOLVED,'Unconsumed source reagent stock kept distinct; individual PCRinputs are source stock aliquots, not extra prepared reactions.'))
    products=[]
    for i,t in enumerate(branches):
        f=g.input('p6_'+str(i),'SourceP6solexa aliquot',K.OLIGO,state=Construct([Segment('source P6',P6)]));r=g.input('p3_'+str(i),'SourceP3solexa aliquot',K.OLIGO,state=Construct([Segment('source P3',P3)]));e=g.input('pcr_enzyme_'+str(i),'SourcePrimeSTARMax/SYBR mix;Phusionalternative unselected;8reactionstockremainder unassembled',K.CATALYST);y=g.product('pcr_'+str(i),'PCRcopies plusoriginalDNA/primers/reagents',K.MOLECULAR)
        g.react('amplify_'+str(i),'Workbook12–18cycles/exponentialstop;optional9+5–8 programme unselected',S.PCR,((t,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(e,R.CATALYST)),(y,),geometry=PCRGeometry.unavailable(unknown('Full RNA-derived/dye/primer-site geometry unassembled; source L5/P6handle mismatch unresolved'),reason='Whole source primer inventories do not establish an exact full-product assembly.'));products.append(y)
    # Unknown aggregation is a declared boundary, never seven invented cleanup tubes.
    lib=g.product('library','Collection of source purified PCR recoveries; physical aggregation/run assignment unassembled',K.LIBRARY)
    removed=g.product('pcr_cleanup_removed','Source physically removed cleanup fractions; individual tube/pool assignment unassembled',K.WASTE)
    wash=g.input('pcr_cleanup_stock','Source NucleoSpinGel/PCRcleanup;85percentethanol replaceskitwash',K.SELECTION)
    g.react('library_boundary','Source purification/quantification/acquisition boundary; pool timing and tube assignments unassembled',S.MAGIC,tuple((x,R.TEMPLATE)for x in products)+((wash,R.SELECTION),),(lib,removed),reason='The workbook says purify PCR reactions and quantify/sequence; the physical tube/aggregation scheme is unassigned rather than invented.')
    g.discard(removed,'Source cleanup removes physical fractions; no perfect purity or per-tube allocation asserted.')
    g.sequence(lib,SequencingPlan(limitation='SourceIlluminaHiSeq mentioned; exact run cycles/primer landing unassigned. BothL5/L3barcodes andPCRduplicate filtering computational; sevenRNA precursor positions versushexamer processing unresolved.'))
    g.add_source_item(SourceItem('pending','Fullfinaldye-conjugate/end geometry,exactPCRpooling/run/coverage anddirectcitedprotocol','Three-originalpassclosed;ProtocolExchange2021rs.3.pex-1333/v1/vendororiginals andunselectedexperimentsunread. No confirmedabsence.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();c('Source precursor inventory builds',bool(precursor_l5()),True);c('No complete publication contract claimed',bool(g.publication_problems()),True);c('RTremnant explicitly unresolved','rt_remainder'in g.terminals,True);c('Seven distinct physical PCR inputs',len(g.steps['pcr_allocate'].outputs),8);c('RTproduct feeds split without invented cleanup',g.steps['pcr_allocate'].inputs[0].material,'rt_product');c('Source elution residue retained','elution_remainder'in g.terminals,True)
    src=Source(sources_root()/'easyclip__10.1038+s41467-021-21623-4'/'s41467-021-21623-4_supp_41467_2021_21623_MOESM3_ESM.xlsx.txt','get_sources --slugs easyclip__10.1038+s41467-021-21623-4;publishedSupplementaryData1')
    if have(c,src,label='PublishedWorkbook'):
        text=src.path.read_text();c('FullsourceP3P6RTprimertranscription',all(seq in text for seq in (P3,P6,RT)),True);c('WholemodifiedL5precursorrow',L5_MACRO in text,True);c('WholemodifiedL3precursorrow',L3_MACRO in text,True)
    c.report()
if __name__=='__main__':selftest()

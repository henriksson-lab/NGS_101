"""Source-owned private partial LABS preparation; no historical kit inheritance."""
TITLE='LABS selected extracted-cfDNA route'
AUTHORITY='Published PMC11177480 main129/211–237,SI FigS1caption and TableS1cfDNAcomparison rows;3originals. Historicalkit details unreviewed.'
INITIAL='Source extracted cfDNA; plasma extraction/cfDNA nativefragment lengths not equated with150–400bp fragmentedgDNA comparison'
CELL_INPUT=False
READOUT='PublishedTableS1cfDNAcomparison rows specify SE75. No NextSeq500 SR80 assumption inherited from separately described MethylC/EpiGnome comparators. DNA uracil chemistry afterbisulfite is distinct from RNA U and sequencedDNA T; complete kit/index/primer landing geometry remains unassembled.'
DEFERRED='KAPA HyperPrepKK8504/MethylCodeMECOV50/EpiMarkM0490S/HiScribeE2040S/KAPARNAHyperPrepKK8541 historical guide originals and exact source-specific extension/kitfinish geometry unassembled; PCRcycles/indexoligos not assigned by this selected pass. No claimed absence.'
T7_RAW='/5Phos/iMe-dC/iMe-dC/iMe-dC/TATAGTGAGT/iMe-dC/GTATTAATTT/iMe-dC/G/iMe-dC/GGGG/iMe-dC/T'
T7='CCCTATAGTGAGTCGTATTAATTTCGCGGGGCT'
HELPER='CGACTCACTATAGGGT'
T7_PRIMER='AGCCCCGCGAAATTAATACGACTCACTATAGGG'
OLIGOS={
 'Source methylatedT7 long strand':{'segments':lambda:[Segment('source methylated T7 adapter',T7)],'mod5':'P','chemistry':{'polymer':'DNA','five_prime':'P','modified_bases':tuple((i,'5mC')for i,b in enumerate(T7)if b=='C')}},
 'Source3P short helper':{'segments':lambda:[Segment('source helper',HELPER)],'mod3':'P','chemistry':{'polymer':'DNA','three_prime':'P'}},
 'Source T7 extension primer':{'segments':lambda:[Segment('source T7 primer',T7_PRIMER)],'chemistry':{'polymer':'DNA'}}}

from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from chemdraw import Scene,Segment,feature
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry

def scene(label):
    item=OLIGOS.get(label)
    if item is None:return None
    sc=Scene();sc.strand('oligo',item['segments'](),mod5=item.get('mod5',''),mod3=item.get('mod3',''));return sc

def chemistry(label):return OLIGOS.get(label,{}).get('chemistry')
def unknown_scene():
    sc=Scene();sc.strand('unassembled',[Segment('unassembled geometry','[source geometry unassembled]',placeholder=True)]);return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance')
    x=g.input('input',INITIAL,K.BIOLOGICAL if CELL_INPUT else K.MOLECULAR,primary=True)
    def stock(key,label):return g.input(key,label,K.OLIGO,state=scene(label),chemical_state=chemistry(label))
    def op(key,label,kind,reagent,product,oligos=(),symmetric=False):
        nonlocal x
        enzyme=kind in {S.DIGESTION,S.END_REPAIR,S.EXTENSION,S.TAILING,S.LIGATION,S.PCR,S.DEPROTEINIZATION,S.TAGMENTATION,S.TRANSCRIPTION}
        c=g.input(key+'_stock',reagent,K.CATALYST if enzyme else K.CONDITION)
        ins=[(x,R.TEMPLATE),(c,R.CATALYST if enzyme else R.CONDITION)]
        for i,o in enumerate(oligos):ins.append((stock(key+'_oligo_'+str(i),o),R.STRUCTURAL))
        kw={}
        if kind==S.MAGIC:kw['reason']='Source-specific compound transition or molecular geometry remains unassembled in this selected model.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        if kind==S.LIGATION:kw['geometry']=AdapterLigationGeometry('unavailable','unavailable',(unknown_scene(),),'Actual source ligation; full source-specific substrate/junction geometry unassembled, not a source absence claim.')
        if kind==S.PCR:kw['geometry']=PCRGeometry.unavailable(unknown_scene(),reason='Source amplification known; whole source-specific primer landing and product geometry unassembled.',symmetric_species=symmetric)
        y=g.product(key+'_product',product,K.MOLECULAR);g.react(key,label,kind,ins,(y,),**kw);x=y
    def clean(key,label):
        nonlocal x
        c=g.input(key+'_condition',label,K.CONDITION);y=g.product(key+'_kept','Physically recovered source preparation; exact recovery/carryover unassigned',K.MOLECULAR);w=g.product(key+'_waste','Physically separated source workup fractions',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(c,R.CONDITION)),(y,w));g.discard(w,'Source physical workup, not perfect molecular exclusion.');x=y
    def split(key,label,rest_label,disposition=Disposition.UNRESOLVED):
        nonlocal x
        y=g.product(key+'_used',label,K.MOLECULAR);r=g.product(key+'_rest',rest_label,K.QC if disposition==Disposition.QC else K.MOLECULAR);g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,r));g.mark_terminal(Terminal(r.id,disposition,rest_label));x=y
    def qc(key,label):split(key,'Production preparation retained after separate analytical aliquot',label,Disposition.QC)
    def external_anneal(key,primer,label,enzyme,programme,product,rt=False):
        nonlocal x
        p=stock(key+'_primer',primer);c=g.input(key+'_anneal_condition',label,K.CONDITION)
        t=g.product(key+'_template','Same-tube template component after annealing; external primer tracked separately, not an isolated fraction',K.MOLECULAR)
        retained=g.product(key+'_retained_primer','Same-tube co-retained external primer component, not a second stock addition or intramolecular geometry',K.OLIGO,state=scene(primer),chemical_state=chemistry(primer))
        g.react(key+'_anneal',label,S.THERMAL_TREATMENT,((x,R.TEMPLATE),(p,R.STRUCTURAL),(c,R.CONDITION)),(t,retained))
        e=g.input(key+'_enzyme',enzyme,K.CATALYST);y=g.product(key+'_product',product,K.MOLECULAR)
        g.react(key,programme,S.REVERSE_TRANSCRIPTION if rt else S.EXTENSION,((t,R.TEMPLATE),(retained,R.STRUCTURAL),(e,R.CATALYST)),(y,));x=y

    long=stock('t7_long','Source methylatedT7 long strand');short=stock('t7_helper','Source3P short helper');cond=g.input('adapter_anneal_condition','Source10mMTrispH8/0.1mMEDTA/50mMNaCl;95C5min -0.25Cmin to4C;50uMequimolarthen15uM/-20C',K.CONDITION)
    adapter=g.product('t7_adapter','Source preannealed partialduplex T7adapter; ends known by oligo inventory, complete ligation substrate/product geometry unassembled',K.ADAPTER)
    g.react('adapter_anneal','Source independent T7adapter stock preparation',S.HYBRIDIZATION,((long,R.TEMPLATE),(short,R.STRUCTURAL),(cond,R.CONDITION)),(adapter,))
    op('endprep','SourceKAPAHyperPrep20C30min65C30min endrepair/Atail',S.END_REPAIR,'Source KAPA endrepair/Atailenzymemix/buffer','Source endconditioned cfDNA classes; historical exact end distributions unassembled')
    ligase=g.input('ligase','Source KAPA ligase/buffer/water',K.CATALYST);lig=g.product('ligate_product','Source T7adapter-ligated DNA plus unligated/partial products; no yield assumption',K.MOLECULAR)
    g.react('ligate','Source T7ligation20C4h or4Covernight; selected20C4h',S.LIGATION,((x,R.TEMPLATE),(adapter,R.STRUCTURAL),(ligase,R.CATALYST)),(lig,),geometry=AdapterLigationGeometry('unavailable','unavailable',(unknown_scene(),),'Exact source-specific whole cfDNA ends/junction paths unassembled; not undisclosed.'));x=lig
    clean('ligate_clean','Source1.4xAMPure postligationcleanup')
    op('bisulfite','Source MethylCodeMECOV50 conversion;9–20uleluate',S.CHEMICAL_CONVERSION,'Source MethylCodekit formulation/thermal programme unreviewed','Source converted DNA with unmodifiedC→U and protected5mC; damaged/nicked fragments and incompleteconversion remain possible')
    op('extend','Source EpiMarkHotStartTaq/T7primer extension; thermal programme notassigned',S.EXTENSION,'Source EpiMarkHotStartTaqM0490S','Source T7promoter-duplex-bearing DNA templates; source-specific extent/end geometry unassembled',('Source T7 extension primer',))
    op('protease','SourceQiagenprotease50C2h75C30min',S.DEPROTEINIZATION,'Source Qiagenprotease19155','Source protease-treated extendedDNA/reagents; no unreported intermediate cleanup')
    op('ivt','SourceHiScribeT7HighYield+SUPERaseIn37C12–16h',S.TRANSCRIPTION,'Source HiScribeT7RNApolymerase/nucleotidekit andRNaseinhibitor','Source linear RNA copies and templateDNA/proteolysis products; no exact amplification yield')
    op('dnase','SourceDNaseI/buffer37C10min',S.DIGESTION,'Source NEBDNaseIM0303S','Source RNA and DNA digestion products; digestion is notphysicalDNAremoval')
    clean('rna_clean','SourceRNACleanandConcentratorR1013')
    qc('rna_quant','Separate source QubitRNAHS analytical aliquot')
    split('rna_aliquot','Source maximum100ngRNA allocated to librarykit','Remaining RNA mass/fate unassigned; no allRNA assumption')
    op('rna_library','SourceKAPARNAHyperPrepKK8541 RT/secondstrand/adaptorligation/amplification compound boundary',S.MAGIC,'Source historicalKAPA kit internals/oligos unreviewed; no extra enzyme stocks invented','Source amplified DNA library derived from IVTRNA; exact adapter/index/PCRgeometry unassembled')
    library=g.product('library','Source sequencing preparation; full source-specific library geometry unassembled',K.LIBRARY)
    g.react('loading','Physical source loading preparation',S.GENERIC,((x,R.TEMPLATE),),(library,),changes_molecule=False)
    g.sequence(library,SequencingPlan(limitation=READOUT))
    g.add_source_item(SourceItem('pending',DEFERRED,'Selected source scope only; incomplete linked-source/whole molecular review, no accepted absence claim.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Source methylated full T7 ordering',T7,'CCCTATAGTGAGTCGTATTAATTTCGCGGGGCT')
    c('Source helper full ordering',HELPER,'CGACTCACTATAGGGT')
    c('Source T7primer whole ordering',T7_PRIMER,'AGCCCCGCGAAATTAATACGACTCACTATAGGG')
    c('Helper blocking retained',chemistry('Source3P short helper')['three_prime'],'P')
    c('Ligation uses prepared adapter not newstock',g.steps['ligate'].inputs[1].material,'t7_adapter')
    c('Bisulfite follows postligation cleanup',g.steps['bisulfite'].inputs[0].material,'ligate_clean_kept')
    c('Extension before protease',g.steps['protease'].inputs[0].material,'extend_product')
    c('Protease before T7 IVT',g.steps['ivt'].inputs[0].material,'protease_product')
    c('RNAphysicalcleanup after DNase',g.steps['rna_clean'].inputs[0].material,'dnase_product')
    c('Kit receives allocated RNA',g.steps['rna_library'].inputs[0].material,'rna_aliquot_used')
    c('Untransferred RNA explicit',g.terminals['rna_aliquot_rest'].disposition,Disposition.UNRESOLVED)
    src=Source(sources_root()/Path(__file__).resolve().parents[1].name/'s13059-024-03262-2_PMC11177480.xml.txt','Publishedmain219 originalT7ordering doi10.1186/s13059-024-03262-2')
    if have(c,src,label='source whole modified T7 ordering'):
        import re
        text=src.path.read_text();c('Source complete modified ordering present',T7_RAW in text,True)
        normalized=re.sub(r'/5Phos/','',T7_RAW).replace('iMe-dC','C').replace('/','')
        c('Normalized whole modified T7 ordering',T7,normalized)
        c('All source methylated Cs retained',chemistry('Source methylatedT7 long strand')['modified_bases'],tuple((i,'5mC')for i,b in enumerate(normalized)if b=='C'))
    c('Private partial graph remains incomplete',bool(g.publication_problems()),True)
    c.report()
if __name__=='__main__':selftest()

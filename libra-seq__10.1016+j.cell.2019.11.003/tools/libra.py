"""Published LIBRA-seq13-base antigen branch; source-limited kit handoffs."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
import illumina as il
from chemdraw import Construct, Scene, Segment, feature, revcomp, protocol_graph_panels
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SequencingPlan
from pcr import PCRGeometry
from page import head, info
PAGE_FILE='libra.html'
SOURCE='https://doi.org/10.1016/j.cell.2019.11.003; published-version Methods198–220'
PREFIX=revcomp(il.TRUSEQ_SMALL_RNA_RA3)
CAPTURE_TAIL='CCCATATAAGAAA'
TAG_WRITTEN=PREFIX+'N'*13+'CCCATATAAGA*A*A'
ADDITIVE_WRITTEN='CCTTGGCACCCGAGAATT*C*C'
ADDITIVE=ADDITIVE_WRITTEN.replace('*','')
ANTIGEN_BARCODES={'BG505':'CATGATTGGCTCA','CZA97':'TGTCCGGCAATAA','H1 New Caledonia':'GATCGTAATACCA'}
def antigen_tag():
    return Construct([Segment('smallRNA-derived handle',PREFIX),
        Segment('antigen identity','N'*13,'cbc',placeholder=True,feature=feature('libra.antigen','feature_barcode','whitelist',whitelist=SOURCE+';13-base antigen dictionary')),
        Segment('beadTSO annealing sequence',CAPTURE_TAIL,note='Three PS bonds as written:CCCATATAAGA*A*A')],name='Source13-base antigen-tag pattern')
def cleanup(g,key,template,reagent,label,kind=K.MOLECULAR):
    r=g.input(key+'_reagent',reagent,K.SELECTION)
    kept=g.product(key+'_kept',label,kind);waste=g.product(key+'_waste','Unretained '+key+' material',K.WASTE)
    g.react(key,key.replace('_',' '),S.CLEANUP,((template,R.TEMPLATE),(r,R.SELECTION)),(kept,waste))
    g.discard(waste,'Removed during this physical recovery/selection');return kept

def protocol_graph():
    g=ProtocolGraph('LIBRA-seq —13-base antigen/BCR branch',authority=SOURCE)
    cells=g.input('cells','B cells for antigen-specific receptor profiling',K.BIOLOGICAL,primary=True)
    oligo=g.input('antigen_dna','5prime-amino antigen DNA pattern',K.OLIGO,state=antigen_tag())
    antigen=g.input('antigen','Biotinylated antigen proteins',K.BIOLOGICAL)
    linker=g.input('linkers','Solulink4FB/HyNic conjugation reagents',K.CONDITION)
    conjugated=g.product('conjugated','Antigen–DNA conjugates',K.MOLECULAR)
    g.react('conjugate','Conjugate5prime-amino oligo to antigen',S.CHEMICAL_CONVERSION,
        ((oligo,R.TEMPLATE),(antigen,R.ANALYTE),(linker,R.CONDITION)),(conjugated,))
    purified=cleanup(g,'fplc',conjugated,'AKTA FPLC','Purified antigen–DNA conjugates')
    pe=g.input('streptavidin_pe','Streptavidin–phycoerythrin',K.SELECTION)
    fluorescent=g.product('fluorescent','Fluorescent antigen–DNA complexes',K.MOLECULAR)
    g.react('fluorescent_complex','Combine biotinylated antigen with streptavidin-PE',S.BINDING,
        ((purified,R.TEMPLATE),(pe,R.SELECTION)),(fluorescent,))
    stained=g.product('stained','B cells stained with barcoded antigens',K.MOLECULAR)
    g.react('stain','Stain B cells with antigen complexes and markers',S.BINDING,
        ((cells,R.TEMPLATE),(fluorescent,R.STRUCTURAL)),(stained,))
    washed=cleanup(g,'wash',stained,'Three PBS/BSA washes','Washed antigen-stained cells')
    facs=g.input('facs','Antigen-positive FACS gate',K.SELECTION)
    selected=g.product('selected','Antigen-positive sorted B cells',K.MOLECULAR)
    rejected=g.product('facs_rejected','Cells not retained by antigen-positive selection',K.WASTE)
    g.react('sort','Select antigen-positive B cells',S.SELECTION,((washed,R.TEMPLATE),(facs,R.SELECTION)),(selected,rejected))
    g.discard(rejected,'Not collected in antigen-positive sorted input; no perfect binding-specificity claim')
    bead=g.input('bead','Cited Chromium B-cellV(D)J bead-delivered primers; historical revision unresolved',K.OLIGO)
    rt=g.input('rt','RT reagents of cited B-cellV(D)J solution; exact formulation unassigned',K.CATALYST)
    barcoded=g.product('barcoded','Cell-indexed BCR cDNA and antigen reporter products',K.MOLECULAR)
    g.react('barcode','Co-encapsulate B cells and index BCR/antigen material',S.REVERSE_TRANSCRIPTION,
        ((selected,R.TEMPLATE),(bead,R.STRUCTURAL),(rt,R.CATALYST)),(barcoded,))
    add=g.input('additive','Source PS-protected additive antigen primer',K.OLIGO,state=Construct([Segment('additive',ADDITIVE,note=ADDITIVE_WRITTEN)]))
    fw=g.input('cdna_fw','Cited cDNA amplification forward primer; design unresolved',K.OLIGO)
    rv=g.input('cdna_rv','Cited cDNA amplification reverse primer; design unresolved',K.OLIGO)
    pol=g.input('polymerase','Cited cDNA/PCR reagents; source-specific formulation unassigned',K.CATALYST)
    amplified=g.product('amplified','Amplified cellular/BCR and antigen products',K.MOLECULAR)
    g.react('cdna_pcr','Amplify cDNA with additive antigen primer',S.PCR,
        ((barcoded,R.TEMPLATE),(fw,R.STRUCTURAL),(rv,R.STRUCTURAL),(add,R.STRUCTURAL),(pol,R.CATALYST)),(amplified,),
        geometry=PCRGeometry.unavailable(Scene(),reason='Complete historical bead/cDNA strand geometry not assigned; no poly(A) antigen-capture recipe imported.'))
    spri=g.input('spri','Source SPRI size separation; initial ratio unassigned',K.SELECTION)
    bcr=g.product('bcr_fraction','Cellular/BCR cDNA fraction',K.MOLECULAR)
    tags=g.product('antigen_fraction','Antigen reporter fraction',K.MOLECULAR)
    g.react('size_split','Retain both BCR and antigen fractions',S.SPLIT,((amplified,R.TEMPLATE),(spri,R.SELECTION)),(bcr,tags))
    bcrlib=g.product('bcr_library','BCR sequencing library',K.LIBRARY)
    g.react('bcr_handoff','Cited10xV(D)J library preparation',S.MAGIC,((bcr,R.TEMPLATE),),(bcrlib,),
        reason='Source delegates to historical B-cellV(D)J solution; exact enrichment/adapter/primer version remains unreviewed.')
    tag_fw=g.input('tag_fw','Antigen library forward primer, exact design unassigned',K.OLIGO)
    tag_rv=g.input('tag_rv','Antigen library reverse primer, exact design unassigned',K.OLIGO)
    tag_pcr=g.product('tag_pcr','10–12-cycle antigen PCR products',K.MOLECULAR)
    g.react('tag_pcr_step','Amplify antigen library10–12cycles',S.PCR,
        ((tags,R.TEMPLATE),(tag_fw,R.STRUCTURAL),(tag_rv,R.STRUCTURAL),(pol,R.CATALYST)),(tag_pcr,),
        geometry=PCRGeometry.unavailable(Scene(),reason='Full antigen PCR product/library-adapter geometry remains unassigned while cited protocol/supplements are unread.'))
    taglib=cleanup(g,'tag_cleanup',tag_pcr,'1.6× purification as reported','Antigen sequencing library',K.LIBRARY)
    for lib in (bcrlib,taglib):
        g.sequence(lib,SequencingPlan(acquired_roles=('Read 1','Read 2'),
            limitation='NovaSeq6000 BCR/antigen outputs. Antigen analysis uses reverse-complementedR2 to recover cell barcode,followingUMI10 and antigenID13 in this branch; complete primer/cycle geometry and kit revision unassigned.'))
    g.declare_reviewed_scope(('cells','antigen_dna','antigen','linkers','conjugated','conjugate','fplc_reagent','fplc_kept','fplc_waste','fplc',
        'streptavidin_pe','fluorescent','fluorescent_complex','stained','stain','wash_reagent','wash_kept','wash_waste','wash','facs','selected','facs_rejected','sort',
        'bead','rt','barcoded','barcode','additive','cdna_fw','cdna_rv','polymerase','amplified','cdna_pcr','spri','bcr_fraction','antigen_fraction','size_split',
        'bcr_library','bcr_handoff','tag_fw','tag_rv','tag_pcr','tag_pcr_step','tag_cleanup_reagent','tag_cleanup_kept','tag_cleanup_waste','tag_cleanup'),SOURCE)
    return g

def render_page():
    return '\n'.join([head('LIBRA-seq — bounded antigen/BCR model'),'<div class="wrap"><h1>LIBRA-seq</h1>',
        '<p class="research-notes"><a href="01_libra_seq.html">Research notes</a></p>',
        info('13-base antigen branch. BCR and antigen libraries remain separate. Historical bead/TSO and final primer geometry are unassigned; no poly(A) antigen capture or physical VH:VL junction is invented.'),protocol_graph_panels(protocol_graph(),cls='long'),'</div>'])
def selftest():
    from checks import Check
    c=Check();g=protocol_graph();c('Separate BCR and antigen outputs',len(g.sequencing),2)
    c('No polyA antigen tag substituted','A'*20 in antigen_tag().top(),False)
    from storage import sources_root
    p=sources_root()/Path(__file__).resolve().parents[1].name/'j.cell.2019.11.003_PMC7158953.xml.txt'
    if p.exists():
        text=''.join(p.read_text().split())
        for seq in (TAG_WRITTEN,ADDITIVE_WRITTEN,*ANTIGEN_BARCODES.values()):c('Published Methods transcription '+seq,seq in text)
    else:c.skip('Source transcription','External published article unavailable')
    c.report(exit_on_fail=True)

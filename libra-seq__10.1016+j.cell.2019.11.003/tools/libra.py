"""Published LIBRA-seq13-base antigen branch; source-limited kit handoffs."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
import illumina as il
from chemdraw import Construct, Scene, Segment, feature, revcomp, panel
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SequencingPlan, Terminal, Disposition
from pcr import PCRGeometry
from page import head, info
from protocol_views import assemble_page, WorkflowSection
PAGE_FILE='libra.html'
SOURCE='https://doi.org/10.1016/j.cell.2019.11.003; final-version author manuscript Methods198–220; whole VOR equivalence not established'
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

def supplies(g,key,label,enzyme_label):
    inventory=g.input(key,label,K.MOLECULAR)
    enzyme=g.product(key+'_enzyme',enzyme_label,K.CATALYST)
    nts=g.product(key+'_nts','Copresent nucleotide class from the same actual supplies; no extra dose',K.MOLECULAR)
    buffer=g.product(key+'_buffer','Copresent buffer class; exact formulation unassembled',K.CONDITION)
    g.react(key+'_inventory','Declare copresent functional components; not physical separation or a claimed premix factory',S.GENERIC,((inventory,R.TEMPLATE),),(enzyme,nts,buffer),changes_molecule=False)
    return enzyme,nts,buffer

def conjugate_qc(g,stock,key,label):
    qc=g.product(key+'_qc',label,K.QC)
    remaining=g.product(key+'_stock','Actual conjugate stock after analytical allocations; amounts unassigned',K.MOLECULAR)
    g.react(key+'_allocate','Actual analytical aliquot allocation; no amount inferred',S.SPLIT,((stock,R.TEMPLATE),),(qc,remaining))
    g.mark_terminal(Terminal(qc.id,Disposition.QC,'Source analytical assay; not production material'))
    return remaining


def protocol_graph():
    g=ProtocolGraph('LIBRA-seq —13-base antigen/BCR branch',authority=SOURCE,representation='published_instance')
    cells=g.input('cells','B cells for antigen-specific receptor profiling',K.BIOLOGICAL,primary=True)
    oligo=g.input('antigen_dna','5prime-amino antigen DNA pattern',K.OLIGO,state=antigen_tag(),chemical_state={'source_ends':'5prime amino modification; three source-written phosphorothioate bonds; 3prime tail retained, no unreported terminal chemistry assigned'})
    antigen=g.input('antigen','Biotinylated antigen proteins',K.BIOLOGICAL)
    linker=g.input('linkers','Actual Solulink 4FB/S-HyNic molecular reactants; kit packaging/quantities unassigned',K.MOLECULAR)
    conjugated=g.product('conjugated','Antigen–DNA conjugates',K.MOLECULAR)
    g.react('conjugate','Desalt and modify amino oligo with 4FB and antigen with S-HyNic, then combine; factory intermediates unassembled',S.CHEMICAL_CONVERSION,
        ((oligo,R.TEMPLATE),(antigen,R.ANALYTE),(linker,R.STRUCTURAL),(g.input('conjugation_environment','Source kit-directed conjugation environment; exact buffer/condition values unassembled',K.CONDITION),R.CONDITION)),(conjugated,))
    assessed=conjugate_qc(g,conjugated,'conjugate_assay','Actual BCA and NanoDrop analytical allocations; concentration/HyNic substitution assays')
    fplc=g.input('fplc_reagent','AKTA FPLC operation; actual column/mobile-phase molecular inventory unassembled, apparatus is not a reagent',K.CONDITION)
    purified=g.product('fplc_kept','FPLC-retained antigen–DNA conjugates; whole conjugation structures unassembled',K.MOLECULAR)
    excess=g.product('fplc_waste','Excess oligonucleotides in unretained FPLC fractions',K.WASTE)
    g.react('fplc','FPLC removes excess oligo',S.CLEANUP,((assessed,R.TEMPLATE),(fplc,R.CONDITION)),(purified,excess));g.discard(excess,'Unretained chromatography material')
    purified=conjugate_qc(g,purified,'silver_stain','Actual conjugate SDS-PAGE/silver-stain analytical allocation')
    pe=g.input('streptavidin_pe','Actual streptavidin–phycoerythrin protein contribution; whole source stock allocation',K.MOLECULAR)
    fluorescent=g.product('fluorescent','Fluorescent antigen–DNA complexes',K.MOLECULAR)
    g.react('fluorescent_complex','Combine at antigen:streptavidin 5:1 molar ratio; add 1/5 of the same PE stock every 20 min at 4°C',S.BINDING,
        ((purified,R.TEMPLATE),(pe,R.STRUCTURAL)),(fluorescent,))
    stained=g.product('stained','B cells stained with barcoded antigens',K.MOLECULAR)
    g.react('stain','Stain B cells with antigen complexes and markers',S.BINDING,
        ((cells,R.TEMPLATE),(fluorescent,R.STRUCTURAL),(g.input('staining_markers','Actual source phenotyping antibodies/viability dye; selected 13-base cohort marker identities qualified in Methods',K.MOLECULAR),R.STRUCTURAL)),(stained,))
    wash_buffer=g.input('wash_reagent','Actual PBS/1% BSA wash-buffer protein/carrier contribution; three source washes',K.MOLECULAR)
    washed=g.product('wash_kept','Washed antigen-stained cells',K.MOLECULAR)
    wash_waste=g.product('wash_waste','Actual wash supernatants with unbound staining material',K.WASTE)
    g.react('wash','Three physical PBS/BSA washes',S.CLEANUP,((stained,R.TEMPLATE),(wash_buffer,R.STRUCTURAL)),(washed,wash_waste));g.discard(wash_waste,'Removed wash supernatants')
    facs=g.input('facs','Antigen-positive FACS gate',K.SELECTION)
    selected=g.product('selected','Antigen-positive sorted B cells',K.MOLECULAR)
    rejected=g.product('facs_rejected','Cells not retained by antigen-positive selection',K.WASTE)
    g.react('sort','Select antigen-positive B cells',S.SELECTION,((washed,R.TEMPLATE),(facs,R.SELECTION)),(selected,rejected))
    g.discard(rejected,'Not collected in antigen-positive sorted input; no perfect binding-specificity claim')
    bead=g.input('bead','Cited Chromium B-cellV(D)J bead-delivered primers; historical revision unresolved',K.OLIGO)
    rt,rt_nt,rt_buffer=supplies(g,'rt','Actual historical Chromium RT aggregate supplies; stock packaging and formulation unassigned','Actual reverse-transcriptase functional class')
    barcoded=g.product('barcoded','Cell-indexed BCR cDNA and antigen reporter products',K.MOLECULAR)
    g.react('barcode','Co-encapsulate B cells and index BCR/antigen material',S.REVERSE_TRANSCRIPTION,
        ((selected,R.TEMPLATE),(bead,R.STRUCTURAL),(rt,R.CATALYST),(rt_nt,R.STRUCTURAL),(rt_buffer,R.CONDITION)),(barcoded,))
    add=g.input('additive','Source PS-protected additive antigen primer',K.OLIGO,state=Construct([Segment('additive',ADDITIVE,note=ADDITIVE_WRITTEN)]))
    fw=g.input('cdna_fw','Cited cDNA amplification forward primer; design unresolved',K.OLIGO)
    rv=g.input('cdna_rv','Cited cDNA amplification reverse primer; design unresolved',K.OLIGO)
    pol,pcr_nt,pcr_buffer=supplies(g,'polymerase','Actual cDNA-PCR aggregate supplies; stock packaging/formulation unassigned','Actual cDNA-PCR polymerase functional class')
    amplified=g.product('amplified','Amplified cellular/BCR and antigen products',K.MOLECULAR)
    g.react('cdna_pcr','Amplify cDNA with additive antigen primer',S.PCR,
        ((barcoded,R.TEMPLATE),(fw,R.STRUCTURAL),(rv,R.STRUCTURAL),(add,R.STRUCTURAL),(pol,R.CATALYST),(pcr_nt,R.STRUCTURAL),(pcr_buffer,R.CONDITION)),(amplified,),
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
    tag_pol,tag_nt,tag_buffer=supplies(g,'tag_supplies','Actual separate antigen-PCR aggregate supplies; historical recipe/stock packaging unassembled','Actual antigen-PCR polymerase functional class')
    tag_pcr=g.product('tag_pcr','10–12-cycle antigen PCR products',K.MOLECULAR)
    g.react('tag_pcr_step','Amplify antigen library10–12cycles',S.PCR,
        ((tags,R.TEMPLATE),(tag_fw,R.STRUCTURAL),(tag_rv,R.STRUCTURAL),(tag_pol,R.CATALYST),(tag_nt,R.STRUCTURAL),(tag_buffer,R.CONDITION)),(tag_pcr,),
        geometry=PCRGeometry.unavailable(Scene(),reason='Full antigen PCR product/library-adapter geometry remains unassigned while cited protocol/supplements are unread.'))
    taglib=cleanup(g,'tag_cleanup',tag_pcr,'1.6× purification as reported','Antigen sequencing library',K.LIBRARY)
    for key,lib in (('bcr',bcrlib),('antigen',taglib)):
        loading=g.product(key+'_loading','Actual '+key+' library sequencing allocation; fraction/absolute amount unassigned',K.LIBRARY)
        stock=g.product(key+'_stock','Untransferred '+key+' prepared library stock; amount/fate unassigned',K.LIBRARY)
        g.react(key+'_load_allocate','Allocate actual NovaSeq loading contribution; no cross-branch physical pooling asserted',S.SPLIT,((lib,R.TEMPLATE),),(loading,stock))
        g.mark_terminal(Terminal(stock.id,Disposition.UNRESOLVED,'Untransferred library stock; no disposal inferred'))
        g.sequence(loading,SequencingPlan(acquired_roles=('Read 1','Read 2'),
            limitation='NovaSeq6000 BCR/antigen outputs. Antigen analysis uses reverse-complementedR2 to recover cell barcode,followingUMI10 and antigenID13 in this branch; complete primer/cycle geometry and kit revision unassigned.'))
    g.declare_reviewed_scope(('cells','antigen_dna','antigen','linkers','conjugated','conjugate','fplc_reagent','fplc_kept','fplc_waste','fplc',
        'streptavidin_pe','fluorescent','fluorescent_complex','stained','stain','wash_reagent','wash_kept','wash_waste','wash','facs','selected','facs_rejected','sort',
        'bead','rt','barcoded','barcode','additive','cdna_fw','cdna_rv','polymerase','amplified','cdna_pcr','spri','bcr_fraction','antigen_fraction','size_split',
        'bcr_library','bcr_handoff','tag_fw','tag_rv','tag_pcr','tag_pcr_step','tag_cleanup_reagent','tag_cleanup_kept','tag_cleanup_waste','tag_cleanup'),SOURCE)
    g.declare_reviewed_scope(('conjugation_environment','conjugate_assay_qc','conjugate_assay_stock','conjugate_assay_allocate','silver_stain_qc','silver_stain_stock','silver_stain_allocate','staining_markers',
        'rt_inventory','rt_enzyme','rt_nts','rt_buffer','polymerase_inventory','polymerase_enzyme','polymerase_nts','polymerase_buffer',
        'tag_supplies','tag_supplies_inventory','tag_supplies_enzyme','tag_supplies_nts','tag_supplies_buffer',
        'bcr_loading','bcr_stock','bcr_load_allocate','antigen_loading','antigen_stock','antigen_load_allocate'),SOURCE)
    return g.require_publishable()

def render_page():
    g=protocol_graph()
    sc=Scene();sc.strand('tag',list(antigen_tag()),mod5='amino',label='13-base antigen-tag pattern; PS bonds are source-written in the tail')
    return assemble_page(g,(head('LIBRA-seq — antigen specificity and B-cell receptors'),
        '<main class="wrap" data-protocol-complete="true"><h1>LIBRA-seq</h1><p>DNA-tagged antigens identify binding B cells. Droplet processing co-indexes antigen reporters and B-cell receptor transcripts; the two physical library branches are retained separately.</p><p><a href="https://doi.org/10.1016/j.cell.2019.11.003">Setliff et al. (2019)</a> · <a href="01_libra_seq.html">Research note</a></p>',
        info('Source-qualified 13-base antigen branch from the final-version author manuscript, Methods 198–220. Complete version-of-record equivalence is not established. Historical Chromium/CITE finishing, whole occupied cellular molecules and read-primer geometry remain unassembled; no modern kit recipe is substituted.'),
        '<h2>Known antigen reporter</h2>',panel(sc.rows(),caption='Source 5′ amino tag: small-RNA-derived handle, variable 13-base antigen identity and TSO-annealing tail. Source tail is CCCATATAAGA*A*A; stars retain phosphorothioate bond notation. This is a free tag pattern, not a reconstructed complete sequencing library.'),
        '<h2>Preparation and read setup</h2><p>4FB/HyNic conjugation and FPLC precede streptavidin-PE labeling, staining, washing and antigen-positive sorting. Chromium co-indexing and cDNA amplification with a source additive primer feed SPRI separation. The antigen fraction receives 10–12 additional PCR cycles and 1.6× purification; BCR finishing follows the cited historical solution.</p><p>Both branches were sequenced on NovaSeq6000. Paired-end FASTQ analysis reverse-complements antigen R2, identifies a cell barcode, then a ten-base UMI and antigen identity. This processing convention does not establish exact acquired cycles or physical primer landing. N90 uses a different 15-base antigen panel and is outside this 13-base model.</p>',
        '<details><summary>Full source workflow, analytical allocations and unresolved chemistry</summary>',WorkflowSection('Source-qualified preparation and separate library outputs'),'</details></main>'))

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
    c('Conjugation reactants are actual molecular supplies',g.materials['linkers'].kind,K.MOLECULAR)
    c('FPLC apparatus is not a molecular reagent',g.materials['fplc_reagent'].kind,K.CONDITION)
    c('Wash BSA remains molecular contribution',g.materials['wash_reagent'].kind,K.MOLECULAR)
    c('Actual conjugate QC is separate',g.terminals['conjugate_assay_qc'].disposition,Disposition.QC)
    c('Distinct stage-PCR supply provenance',g.producer_of('tag_supplies_enzyme').id,'tag_supplies_inventory')
    c('RT consumes actual copresent nucleotides','rt_nts' in tuple(x.material for x in g.steps['barcode'].inputs),True)
    c('Both library stock remainders retained',tuple(g.terminals[x].disposition for x in ('bcr_stock','antigen_stock')),(Disposition.UNRESOLVED,Disposition.UNRESOLVED))
    c('Publication gate',g.publication_problems(),[])
    c('Public complete marker','data-protocol-complete="true"' in render_page(),True)
    c.report(exit_on_fail=True)

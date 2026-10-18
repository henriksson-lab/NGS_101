"""mim-tRNAseq owns primer-dependent TGIRT, chemical ends and ss-circle input."""
from chemdraw import Construct,Scene,Segment,feature,revcomp
from circular import circularize_ssdna,CircularizationGeometry
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_model import molecule,role
from protocol_graph import InputRole,MaterialKind,ProtocolGraph,SequencingPlan,StepKind,Disposition,Terminal
SOURCE_OLIGO_EXCEPTIONS={"adapter_precursor":"Published TableS4 and Methods394–400 write terminal /ddC/. The source scraper strips that modification without adding its physical C base, while this state correctly retains terminal dideoxy-C as a nucleotide. Local selftest directly verifies the complete source-written adapter including /ddC/ against original TableS4; no altered sequence is substituted."}
TITLE='mim-tRNAseq — modification-induced misincorporation tRNA sequencing'
DOI='10.1016/j.molcel.2021.01.028'
SOURCE='Behrens et al., Molecular Cell(2021), <a href="https://doi.org/10.1016/j.molcel.2021.01.028">doi:10.1016/j.molcel.2021.01.028</a>; published Methods and TableS4.'
NOTES='01_mim_trnaseq.html'
SUMMARY='Barcoded3′adapter ligation, primer-dependent low-salt TGIRT and single-strand cDNA circularization retain native modification-dependent substitutions.'
CAVEAT='Selected route is primer-dependent TGIRT42°C16h after representative I1 adapter ligation; I2–I4 are source alternatives. The tRNA input includes source CC/CCA synthetic standards as independent RNA species. Template-switching/high-salt/short-RT comparators and charging oxidation/β-elimination are not substituted. ddC is a physical terminal nucleotide; iSp18 is a nonnucleotide spacer retained in physical cDNA/circles. Native first-copy PCR and final internal geometry remain unimplemented, not absent publication data.'
ADAPTER_BASES='GATATCGTCAAGATCGGAAGAGCACACGTCTGAA'
RT_LEFT='AGATCGGAAGAGCGTCGTGTAGGGAAAGAG'
RT_BIND='GTGACTGGAGTTCAGACGTGTGCTC'
MIM_FORWARD='AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTC'
MIM_REVERSE='CAAGCAGAAGACGGCATACGAGAT'+'N'*6+'GTGACTGGAGTTCAGACGTGTG'
BARCODE=feature('mim_adapter_barcode','sample_index','whitelist',whitelist='mim-tRNAseq published I1–I4 first10nt')
INDEX=feature('mim_pcr_index','sample_index','unknown')

def adapter_parts():
    return [Segment('I1 demultiplexing prefix',ADAPTER_BASES[:10],'cbc',feature=BARCODE),Segment('adapter remainder',ADAPTER_BASES[10:]),Segment('terminal ddC nucleotide','C',note='3′dideoxy chain terminator')]

def rt_parts():
    return [Segment('circularization bias-mitigation RN','RN',placeholder=True,length_bp=2,note='DNA IUPAC R=A/G and N; not a molecular identifier'),Segment('RT left PCR handle',RT_LEFT),role('non-nucleotide iSp18; not18bases'),Segment('RT annealing/right handle',RT_BIND)]

def _single(label,parts,**kw):
    sc=Scene();sc.strand('strand',parts,label=label,**kw);return sc

# Compute only the source-compatible local primer/adapter relation. The ddC
# beyond this binding region is retained in the physical adapter, not copied.
_BIND=max(n for n in range(1,min(len(RT_BIND),len(ADAPTER_BASES))+1) if revcomp(RT_BIND[-n:])==ADAPTER_BASES[-n:])
CDNA=Construct([*rt_parts(),Segment('copied adapter between RT site and barcode',revcomp(ADAPTER_BASES[10:-_BIND])),Segment('copied I1 demultiplexing prefix',revcomp(ADAPTER_BASES[:10]),'cbc',feature=BARCODE),role('full-length tRNA cDNA; modifications recorded as substitutions')],name='mim physical single-strand cDNA')
CLOSED=circularize_ssdna(CDNA,five_prime_phosphate=True)
class CircleScene(Scene):
    """Retain actual ss topology for both drawing and unavailable-PCR boundary."""
    def __init__(self):
        super().__init__();self.strand('ss circle',list(CDNA),label='actual single-strand circle; iSp18 retained')
    def rows(self):return CLOSED.rows('covalently closed ss cDNA; physical iSp18 retained')

FINAL_LIBRARY=Construct([Segment('source forward-primer end',MIM_FORWARD),role('circle-derived tRNA/adapter/RN interval; full internal strand geometry unassigned'),Segment('source reverse-primer annealing complement',revcomp('GTGACTGGAGTTCAGACGTGTG')),Segment('indexed-primer index complement','N'*6,placeholder=True,length_bp=6,feature=INDEX),Segment('source P7-primer end complement',revcomp('CAAGCAGAAGACGGCATACGAGAT'))],name='mim primer-derived endpoint reference')
SEQ_PRIMERS=()
READ_LENGTHS=None
FINAL_CAPTION='Source PCR primer-derived ends and six-base physical sample index; complete internal read geometry remains unassigned.'
SEQUENCING_INTRO='Source sequences150cycles on NextSeq; run read/index allocation and exact sequencing primer are not assigned here.'
SEQUENCING_ENDING=SEQUENCING_INTRO
SEQUENCING_UNAVAILABLE='150total reported cycles are not assigned to Read1 or index cycles without a source recipe.'

def oligos():
    return ('I1 ordered precursor:5′-p'+ADAPTER_BASES+'/ddC/-3′ (Mth preadenylation precedes RNA ligation)',
        'Primer-dependent RT:5′-pRN'+RT_LEFT+'/iSp18/'+RT_BIND+'-3′',
        'Forward:5′-'+MIM_FORWARD[:-1]+'*C-3′',
        'Indexed reverse:5′-'+MIM_REVERSE[:-1]+'*G-3′')

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    trna=g.input('native_trna','gel-purified60–100nt tRNA pool;CC/CCAstandards3:1 at0.06pmol/µg totalRNA; baseline deacylated/T4PNK-pretreated',MaterialKind.BIOLOGICAL,state=_single('native modified tRNA',[role('mature tRNA; sequence/CCA heterogeneity unassigned')],mod3='OH'),primary=True)
    adapter=g.input('adapter_precursor','ordered5′phosphorylated I1/ddC adapter',MaterialKind.ADAPTER,state=_single('I1 ordered precursor',adapter_parts(),mod5='P',mod3='ddC'),source='Published Methods394–400 and TableS4')
    mth=g.input('mth','Mth RNA ligase adapter preadenylation',MaterialKind.CATALYST)
    activated=g.product('activated_adapter','preadenylated I1/ddC adapter',MaterialKind.ADAPTER,state=_single('ligation-ready I1',adapter_parts(),mod5='App',mod3='ddC'))
    g.react('preadenylate','Preadenylate ordered5′P adapter before ligation',StepKind.CHEMICAL_CONVERSION,((adapter,InputRole.TEMPLATE),(mth,InputRole.CATALYST)),(activated,))
    ligase=g.input('adapter_ligase','T4RNA ligase2trKQ200U;3h25°C,25%PEG8000,4:1adapter:tRNA',MaterialKind.CATALYST)
    adapted=g.product('adapted_trna','I1-adapter-ligated tRNA',MaterialKind.MOLECULAR,state=_single('adapter-ligated tRNA',[role('native tRNA'),*adapter_parts()],mod3='ddC'))
    ast=activated.state;ast.ligation_end('strand','I1 demultiplexing prefix',kind='single-strand',polarity='5-prime',owner='adapter',text='5′App ligation end')
    g.react('adapter_ligation','Ligate preadenylated3′adapter to tRNA3′OH',StepKind.LIGATION,((trna,InputRole.TEMPLATE),(activated,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(adapted,),geometry=AdapterLigationGeometry.known('single-strand',ast,statement='Published5′P precursor is preadenylated by Mth; ddC prevents adapter concatemerization and remains the physical3′terminal base.'))
    gel=g.input('ligation_gel','10%denaturingPAGE;excise95–125nt',MaterialKind.SELECTION)
    clean=g.product('clean_trna','gel-purified adapter-ligated tRNA',MaterialKind.MOLECULAR,state=adapted.state)
    gw=g.product('ligation_gel_waste','excess adapters/out-of-window gel material',MaterialKind.WASTE)
    g.react('ligation_cleanup','Separate ligated tRNA from excess adapters',StepKind.SELECTION,((adapted,InputRole.TEMPLATE),(gel,InputRole.SELECTION)),(clean,gw),changes_molecule=False);g.discard(gw,'physically outside excised95–125nt band')
    primer=g.input('rt_oligo','primer-dependent5′pRN/iSp18 RT primer',MaterialKind.OLIGO,state=_single('physical RT primer',rt_parts(),mod5='P'),source='Published Methods404/TableS4')
    rt=g.input('tgirt','TGIRT500nM;low-salt75mMKCl,42°C16h;125nMprimer/template',MaterialKind.CATALYST)
    cdna=g.product('full_length_cdna','extended single-strand cDNA target',MaterialKind.MOLECULAR,state=_single('physical cDNA target',list(CDNA),mod5='P',mod3='OH'))
    rna=g.product('rt_template_rna','RNA templates remain after RT',MaterialKind.MOLECULAR,state=clean.state)
    g.react('reverse_transcription','Primer-dependent prolonged low-salt TGIRT readthrough',StepKind.REVERSE_TRANSCRIPTION,((clean,InputRole.TEMPLATE),(primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,rna))
    naoh=g.input('rna_hydrolysis','1µL5MNaOH;95°C3min',MaterialKind.CONDITION)
    ctarget=g.product('clean_cdna','cDNA target after RNA hydrolysis',MaterialKind.MOLECULAR,state=cdna.state)
    hydro=g.product('degraded_rna','RNA hydrolysis products; not removed by hydrolysis',MaterialKind.MOLECULAR)
    g.react('hydrolyse_rna','Hydrolyse RNA templates',StepKind.CHEMICAL_CONVERSION,((cdna,InputRole.TEMPLATE),(rna,InputRole.TEMPLATE),(naoh,InputRole.CONDITION)),(ctarget,hydro))
    size=g.input('cdna_size_selection','10%denaturingPAGE60–150nt;elution/ethanol precipitation',MaterialKind.SELECTION)
    sized=g.product('sized_cdna','size-selected physical cDNA',MaterialKind.MOLECULAR,state=ctarget.state)
    waste=g.product('size_waste','unextended primer/out-of-window gel material',MaterialKind.WASTE)
    g.react('select_cdna','Select60–150nt cDNA and recover by ethanol precipitation',StepKind.SELECTION,((ctarget,InputRole.TEMPLATE),(size,InputRole.SELECTION)),(sized,waste),changes_molecule=False);g.discard(waste,'physically excluded gel/precipitation fractions')
    g.mark_terminal(Terminal(hydro.id,Disposition.UNRESOLVED,'Hydrolysis creates fragments/nucleotides. Their individual recovery through source PAGE/ethanol purification is unquantified, not an asserted discard.'))
    enzyme=g.input('circligase','CircLigase;ATP/MgCl2/betaine,60°C3h;80°C10min',MaterialKind.CATALYST)
    circle=g.product('circular_cdna','ss full-length tRNA cDNA circle',MaterialKind.MOLECULAR,state=CircleScene(),chemical_state=CLOSED)
    g.react('circularize','Close5′phosphorylated cDNA intramolecularly',StepKind.CIRCULARIZATION,((sized,InputRole.TEMPLATE),(enzyme,InputRole.CATALYST)),(circle,),geometry=CircularizationGeometry.known('single-strand intramolecular ligation',statement='Published CircLigase closes RT-primer5′P to cDNA3′OH; physical iSp18 remains in the circle.'))
    used=g.product('pcr_aliquot','one-fifth circle sample forPCR',MaterialKind.MOLECULAR,state=circle.state)
    rest=g.product('unused_circles','remaining four-fifths circle sample',MaterialKind.MOLECULAR,state=circle.state)
    aliquot_rule=g.input('circle_aliquot_rule','one-fifth circularization sample',MaterialKind.CONDITION)
    g.react('sample_circles','Use one-fifth circularization reaction',StepKind.SAMPLING,((circle,InputRole.TEMPLATE),(aliquot_rule,InputRole.CONDITION)),(used,rest),changes_molecule=False);g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Source uses one-fifth; remaining sample handling unassigned.'))
    f=g.input('library_forward','source forward;terminal phosphorothioate',MaterialKind.OLIGO,state=MIM_FORWARD,source='Published Methods410/TableS4')
    r=g.input('library_reverse','six-base indexed reverse;terminal phosphorothioate',MaterialKind.OLIGO,state=Construct([Segment('P7 primer prefix',MIM_REVERSE[:24]),Segment('six-base sample index','N'*6,placeholder=True,length_bp=6,feature=INDEX),Segment('reverse annealing end',MIM_REVERSE[30:],note='terminal phosphorothioate')]),source='Published Methods410/TableS4')
    pol=g.input('library_polymerase','KAPAHiFi GCbuffer;5–6cycles98/62/72°C',MaterialKind.CATALYST)
    amp=g.product('amplified','primer-derived PCR product',MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
    g.react('library_pcr','Amplify source circle aliquot with printed primers',StepKind.PCR,((used,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(used.state,reason='Actual template is the closed single strand retaining iSp18. Native first copying, recopy and any spacer boundary are unimplemented; a completed duplex is not used as input.'))
    cleanup=g.input('library_cleanup','DNAclean&concentrator5;8%PAGE130–220bp;ethanol;10µLTris',MaterialKind.SELECTION)
    lib=g.product('library_library','purified mim-tRNAseq library',MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    w=g.product('library_cleanup_waste','column/gel/precipitation excluded material',MaterialKind.WASTE)
    g.react('library_purify','Purify and gel-select PCR library',StepKind.CLEANUP,((amp,InputRole.TEMPLATE),(cleanup,InputRole.SELECTION)),(lib,w),changes_molecule=False);g.discard(w,'source physical column/gel/precipitation transfers')
    g.sequence(lib,SequencingPlan(library=FINAL_LIBRARY,acquisition='unspecified',limitation=SEQUENCING_INTRO+' '+SEQUENCING_UNAVAILABLE))
    g.declare_reviewed_scope(('native_trna','adapter_precursor','mth','preadenylate','adapter_ligase','adapter_ligation','ligation_gel','ligation_cleanup','rt_oligo','tgirt','reverse_transcription','rna_hydrolysis','hydrolyse_rna','cdna_size_selection','select_cdna','circligase','circularize','circle_aliquot_rule','sample_circles','library_forward','library_reverse','library_polymerase','library_pcr','library_cleanup','library_purify','library_library'),SOURCE)
    return g.require_publishable()

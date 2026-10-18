"""Protocol-owned molecular declarations; source scope is deliberately bounded."""
from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from protocol_model import molecule, role
from protocol_graph import ProtocolGraph, MaterialKind as M, InputRole as I, StepKind as K, SequencingPlan, Terminal, Disposition
from pcr import PCRGeometry
import seqprimers as sp
import nextera as nx
TITLE = 'SMART-seq — original Tn5 branch'
NOTES = '01_smart-seq.html'
SOURCE = 'Ramskold et al. 2012, doi:10.1038/nbt.2282, PMC3467340 Methods and published Supplementary Table 1.'
SUMMARY = 'Original SMART cDNA synthesis followed by the modified Epicentre Nextera preparation; the alternative PE ligation route is outside this selected branch.'
CAVEAT = "Native unassigned oligos and adapter geometry remain role tokens; a later protocol's sequence does not establish the original reagent."
SEQ_PRIMERS = ()
ACQUIRED_ROLES = ()
SEQUENCING_UNAVAILABLE = "Native sequencing primer and index assignments remain unverified in this bounded source pass."
SEQUENCING_INTRO = ""
SEQUENCING_ENDING = "Illumina-compatible native cDNA library. Exact native sequencing geometry remains unassigned."
READ_SPECS = (sp.read_spec('Read 1','50 or 40','Single-end cDNA sequence; Supplementary Table 1 specifies single-end acquisition, with Tn5 rows at 50 nt (GAIIx) or 40 nt (MiSeq).','include'),)
FINAL_CAPTION = "Source-qualified endpoint; unknown ends do not imply canonical Nextera index architecture."
DT_SEQUENCE='AAGCAGTGGTATCAACGCAGAGTAC'+'T'*30+'VN'
DT=(Segment('CDS primer',DT_SEQUENCE,placeholder=True),)
TSO=(role('native SMARTer II A oligo'),)
COPIED_TSO=(role('copied native SMARTer II A boundary'),)
PREAMP_PRIMERS=Construct((role('native kit IS PCR primer'),))
AMPLIFIED=(role('native SMARTer-derived end'),role('amplified cDNA'),role('CDS-derived opposite end'))
FINAL_LIBRARY=Construct((role('native Epicentre Nextera left end'),role('cDNA fragment'),role('native Epicentre Nextera opposite end')),name='original SMART-seq Tn5 library')
def oligos(): return ('CDS primer — 5′-'+DT_SEQUENCE+'-3′','SMARTer II A and kit IS PCR primer: exact native bases remain unassigned in this pass.')

def single(label, parts):
    return molecule(label, tuple(parts))
def duplex(label, parts):
    return Scene.duplex(parts,unpaired=tuple(p.name+"'" for p in parts if p.is_role_token()))
def cleanup(g, current, ident, label):
    out=g.product(ident+'_recovered',label,M.MOLECULAR,state=current.state)
    waste=g.product(ident+'_waste','removed supernatant and washes',M.WASTE)
    reagent=g.input(ident+'_selection',label,M.SELECTION)
    g.react(ident,label,K.CLEANUP,((current,I.TEMPLATE),(reagent,I.SELECTION)),(out,waste),changes_molecule=False)
    g.discard(waste,'Physically removed during the stated purification')
    return out

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    rna=g.input('rna','poly(A) RNA in cell lysate',M.MOLECULAR,state=single('RNA',(role('RNA transcript'),role('native poly(A) tract'))),primary=True)
    dt=g.input('dt','Source CDS anchored oligo(dT)30VN',M.OLIGO,state=Construct(DT))
    rt=g.input('rt','MMLV-derived SMARTScribe reverse transcriptase',M.CATALYST)
    first=g.product('first','first-strand cDNA',M.MOLECULAR,state=single('first cDNA',(*DT,role('reverse-transcribed RNA'),role('native untemplated cDNA 3-prime boundary'))))
    first.state.strand('retained_RNA',(role('independent source RNA template'),),label='RNA retained after RT; extent and hybrid geometry unassigned')
    g.react('reverse_transcription','Reverse-transcribe poly(A) RNA using the printed CDS primer',K.REVERSE_TRANSCRIPTION,((rna,I.TEMPLATE),(dt,I.STRUCTURAL),(rt,I.CATALYST)),(first,))
    tso=g.input('tso','Native SMARTer II A commercial reagent',M.OLIGO,state=Construct(TSO))
    switched_state=single('template-switched single strand',(*DT,role('reverse-transcribed RNA'),*COPIED_TSO))
    switched_state.strand('retained_RNA',(role('independent source RNA template'),),label='RNA beside switched cDNA; no covalent RNA-cDNA join')
    switched=g.product('switched','single-stranded template-switched cDNA',M.MOLECULAR,state=switched_state)
    g.react('template_switch','Copy the native template-switch oligo in the same RT',K.TEMPLATE_SWITCHING,((first,I.TEMPLATE),(tso,I.STRUCTURAL),(rt,I.CATALYST)),(switched,))
    primers=g.input('preamp_primers','Native kit IS PCR primer',M.OLIGO,state=PREAMP_PRIMERS)
    pol=g.input('preamp_polymerase','Native kit cDNA amplification polymerase',M.CATALYST)
    amp=g.product('amplified','amplified full-length double-stranded cDNA',M.MOLECULAR,state=duplex('amplified',AMPLIFIED))
    residual_rna=g.product('residual_rna','RNA-template-derived material after PCR thermal cycling',M.MOLECULAR,state=single('RNA-derived population',(role('native RNA and thermal breakdown products; proportions unassigned'),)))
    g.react('preamp','PCR preamplification: 12 cycles at 1 ng, 15 at 100 pg, 18 at 10 pg or single-cell input',K.PCR,((switched,I.TEMPLATE),(primers,I.STRUCTURAL),(pol,I.CATALYST)),(amp,residual_rna),geometry=PCRGeometry.unavailable(switched_state,reason='Actual single-strand input is shown; native first-cycle conversion and complete primer-site geometry remain unvalidated.',symmetric_species=True))
    g.mark_terminal(Terminal(residual_rna.id,Disposition.UNRESOLVED,'RNA is not converted into cDNA or physically discarded by PCR. Native thermal fate and carryover through subsequent bead cleanups remain unassigned.'))
    clean=amp
    tn=g.input('transposome','Original Epicentre Nextera transposase and 5x HMW buffer; native oligos unassigned',M.ADAPTER)
    tn_enzyme=g.input('tn_enzyme','native Tn5 transposase',M.CATALYST)
    tagged_state=single('native tagged-fragment population',(role('native transferred adapter boundary'),role('cDNA fragment'),role('native opposite fragment end')))
    tagged=g.product('tagged','native tagmented cDNA populations',M.MOLECULAR,state=tagged_state)
    other=g.product('other_tagged','other native end combinations and retained terminal fragments',M.MOLECULAR,state=single('other tagged populations',(role('other native tagmented molecules'),)))
    g.react('tagment','Tagment approximately 5 ng cDNA: 55 C for 5 min',K.TAGMENTATION,((clean,I.TEMPLATE),(tn,I.STRUCTURAL),(tn_enzyme,I.CATALYST)),(tagged,other))
    g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Nonamplification or library incompatibility does not establish physical removal; native frequencies and fates are unassigned.'))
    released=cleanup(g,tagged,'strip_cleanup','35 ul PB strip followed by 88 ul SPRI XP, sample:beads 1:1.6')
    pcr_input=released
    ip=g.input('index_primers','Original standard Nextera PCR primers; native identities and index widths unassigned',M.OLIGO)
    ip_reverse=g.input('index_reverse','native opposite index PCR primer',M.OLIGO)
    ipol=g.input('index_polymerase','Original standard Nextera PCR polymerase',M.CATALYST)
    lib=g.product('indexed','indexed cDNA library',M.LIBRARY,state=FINAL_LIBRARY)
    g.react('index_pcr','Nine cycles of original standard Nextera PCR',K.PCR,((pcr_input,I.TEMPLATE),(ip,I.STRUCTURAL),(ip_reverse,I.STRUCTURAL),(ipol,I.CATALYST)),(lib,),geometry=PCRGeometry.unavailable(pcr_input.state,reason='Native tagged input, not completed indexed library, is the PCR substrate; native adapter/index primer geometry remains unvalidated.'))
    final=lib
    g.sequence(final,SequencingPlan(library=FINAL_LIBRARY,limitation=SEQUENCING_UNAVAILABLE))
    g.declare_reviewed_scope(REVIEWED,SOURCE)
    return g.require_publishable()

REVIEWED=('rna', 'dt', 'rt', 'first', 'reverse_transcription', 'tso', 'switched', 'template_switch', 'preamp_primers', 'preamp_polymerase', 'amplified', 'residual_rna', 'preamp', 'transposome', 'tn_enzyme', 'tagged', 'other_tagged', 'tagment', 'index_primers', 'index_reverse', 'index_polymerase', 'indexed', 'index_pcr', 'strip_cleanup', 'strip_cleanup_recovered', 'strip_cleanup_waste', 'strip_cleanup_selection')

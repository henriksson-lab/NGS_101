"""Protocol-owned molecular declarations; source scope is deliberately bounded."""
from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from protocol_model import molecule, role
from protocol_graph import ProtocolGraph, MaterialKind as M, InputRole as I, StepKind as K, SequencingPlan, Terminal, Disposition
from pcr import PCRGeometry
import seqprimers as sp
import nextera as nx
TITLE = 'SMART-seq3 — author protocol V.3'
NOTES = '01_smart-seq3.html'
SOURCE = 'Hagemann-Jensen et al. 2020 doi:10.1038/s41587-020-0497-0, published Supplementary Table 1; author protocol V.3 doi:10.17504/protocols.io.bcq4ivyw.'
SUMMARY = 'Biotinylated asymmetric RT oligos install the eight-base UMI; native tagmented populations enter custom-index PCR.'
CAVEAT = "Native unassigned oligos and adapter geometry remain role tokens; a later protocol's sequence does not establish the original reagent."
SEQ_PRIMERS = ()
ACQUIRED_ROLES = ()
SEQUENCING_UNAVAILABLE = "Native sequencing primer and index assignments remain unverified in this bounded source pass."
SEQUENCING_INTRO = ""
SEQUENCING_ENDING = "Illumina-compatible native cDNA library. Exact native sequencing geometry remains unassigned."
READ_SPECS = (sp.read_spec('Sequencing output',None,'Author V.3 permits single-end or paired-end Illumina acquisition; read lengths and native run index assignment are unassigned here. The retained 5-prime class carries the TSO tag and eight-base UMI; internal fragments do not.','include'),)
FINAL_CAPTION = "Source-qualified endpoint; unknown ends do not imply canonical Nextera index architecture."
DT_HANDLE='ACGAGCATCAGCAGCATACGA'
DT_SEQUENCE=DT_HANDLE+'T'*30+'VN'
TSO_HANDLE='AGAGACAGATTGCGCAATG'
FWD=nx.S5+nx.ME+'ATTGCGCAATG'
REV=DT_HANDLE
UMI=feature('transcript_umi','umi','random')
DT=(Segment('biotinylated dT handle',DT_HANDLE,note='5-prime /5BiosG/'),Segment('dT30','T'*30),Segment('VN','VN',placeholder=True))
TSO=(Segment('biotinylated TSO handle',TSO_HANDLE,note='5-prime /5Biosg/'),Segment('UMI','N'*8,placeholder=True,feature=UMI),Segment('rGrGrG','GGG',note='three riboguanosines'))
# complement_segments returns reversed segment order with reverse-complemented bases
# (5-prime CCC, copied UMI, copied TSO handle), preserving molecular features.
COPIED_TSO=tuple(complement_segments(TSO))
PREAMP_PRIMERS=Construct((Segment('Fwd PCR',FWD,note='3-prime AA*T*G: two terminal phosphorothioate bonds'),))
AMPLIFIED=(Segment('forward PCR-derived end',FWD),Segment('copied UMI','N'*8,placeholder=True,feature=UMI),Segment('template-switch boundary','GGG'),role('cDNA copied from RNA'),role('native opposite dT-derived end'))
I5=role('native custom indexed left end',molecular_feature=feature('sample_i5','sample_index','unknown'))
I7=role('native opposite custom indexed end',molecular_feature=feature('sample_i7','sample_index','unknown'))
FINAL_LIBRARY=Construct((I5,Segment('TSO tag','ATTGCGCAATG'),Segment('transcript UMI','N'*8,placeholder=True,feature=UMI),Segment('template-switch boundary','GGG'),role('cDNA fragment'),I7),name='SMART-seq3 retained-TSO class')
INTERNAL_LIBRARY=Construct((I5,role('internal cDNA fragment'),I7),name='SMART-seq3 internal class without original UMI')

def oligos(): return ('Smartseq3_OligodT30VN — 5′-/5BiosG/'+DT_SEQUENCE+'-3′','Smartseq3_N8_TSO — 5′-/5Biosg/'+TSO_HANDLE+'NNNNNNNNrGrGrG-3′','Fwd_PCR_primer — 5′-'+FWD[:-3]+FWD[-3]+'*'+FWD[-2]+'*'+FWD[-1]+'-3′','Rev_PCR_primer — 5′-'+REV[:-3]+REV[-3]+'*'+REV[-2]+'*'+REV[-1]+'-3′')

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
    dt_state=Scene();dt_state.strand('dT primer',DT,mod5='/5Biosg/')
    dt=g.input('dt','5-prime biotinylated Smartseq3_OligodT30VN',M.OLIGO,state=dt_state)
    rt=g.input('rt','Maxima H Minus RT',M.CATALYST)
    first=g.product('first','first-strand cDNA',M.MOLECULAR,state=single('first cDNA',(*DT,role('reverse-transcribed RNA'),role('native untemplated cDNA 3-prime boundary'))))
    first.state.strands['molecule'].mod5='/5Biosg/'
    first.state.strand('retained_RNA',(role('independent source RNA template'),),label='RNA retained after RT; extent and hybrid geometry unassigned')
    g.react('reverse_transcription','Lyse at 72 C 10 min; RT 42 C 90 min, 10x(50 C 2 min/42 C 2 min), 85 C 5 min',K.REVERSE_TRANSCRIPTION,((rna,I.TEMPLATE),(dt,I.STRUCTURAL),(rt,I.CATALYST)),(first,))
    tso_state=Scene();tso_state.strand('TSO',TSO,mod5='/5Biosg/')
    tso=g.input('tso','5-prime biotinylated Smartseq3_N8_TSO with rGrGrG',M.OLIGO,state=tso_state)
    switched_state=single('template-switched single strand',(*DT,role('reverse-transcribed RNA'),*COPIED_TSO))
    switched_state.strands['molecule'].mod5='/5Biosg/'
    switched_state.strand('retained_RNA',(role('independent source RNA template'),),label='RNA beside switched cDNA; no covalent RNA-cDNA join')
    switched=g.product('switched','single-stranded template-switched cDNA',M.MOLECULAR,state=switched_state)
    g.react('template_switch','Copy the native template-switch oligo in the same RT',K.TEMPLATE_SWITCHING,((first,I.TEMPLATE),(tso,I.STRUCTURAL),(rt,I.CATALYST)),(switched,))
    primers=g.input('preamp_primers','Phosphorothioate-protected Fwd and Rev PCR primers',M.OLIGO,state=PREAMP_PRIMERS)
    reverse_primer=g.input('preamp_reverse','Rev PCR primer',M.OLIGO,state=Construct((Segment('Rev PCR',REV),)))
    pol=g.input('preamp_polymerase','KAPA HiFi HotStart polymerase',M.CATALYST)
    amp=g.product('amplified','amplified full-length double-stranded cDNA',M.MOLECULAR,state=duplex('amplified',AMPLIFIED))
    residual_rna=g.product('residual_rna','RNA-template-derived material after PCR thermal cycling',M.MOLECULAR,state=single('RNA-derived population',(role('native RNA and thermal breakdown products; proportions unassigned'),)))
    g.react('preamp','Asymmetric preamplification with 0.5 uM Fwd and 0.1 uM Rev; input-dependent cycle count',K.PCR,((switched,I.TEMPLATE),(primers,I.STRUCTURAL),(reverse_primer,I.STRUCTURAL),(pol,I.CATALYST)),(amp,residual_rna),geometry=PCRGeometry.unavailable(switched_state,reason='Actual single-strand input is shown; native first-cycle conversion and complete primer-site geometry remain unvalidated.',symmetric_species=False))
    g.mark_terminal(Terminal(residual_rna.id,Disposition.UNRESOLVED,'RNA is not converted into cDNA or physically discarded by PCR. Native thermal fate and carryover through subsequent bead cleanups remain unassigned.'))
    clean=cleanup(g,amp,'cdna_cleanup','Selected optional cDNA cleanup: 0.6:1 PEG/Sera-Mag beads')
    tn=g.input('transposome','Nextera XT Amplicon Tagment Mix; native transferred-oligo geometry unvalidated',M.ADAPTER)
    tn_enzyme=g.input('tn_enzyme','native Tn5 transposase',M.CATALYST)
    tagged_state=Scene()
    tagged_state.strand('retained_TSO',(Segment('PCR-derived s5 ME tag',FWD),Segment('original transcript UMI','N'*8,placeholder=True,feature=UMI),Segment('template-switch boundary','GGG'),role('cDNA fragment'),role('native Tn5-transferred opposite boundary')))
    tagged_state.strand('internal',(role('native Tn5-transferred left boundary'),role('internal cDNA fragment without original UMI'),role('native Tn5-transferred right boundary')))
    tagged=g.product('tagged','native tagmented cDNA populations',M.MOLECULAR,state=tagged_state)
    other=g.product('other_tagged','other native end combinations and retained terminal fragments',M.MOLECULAR,state=single('other tagged populations',(role('other native tagmented molecules'),)))
    g.react('tagment','Tagment 100 pg normalized cDNA at 55 C for 10 min',K.TAGMENTATION,((clean,I.TEMPLATE),(tn,I.STRUCTURAL),(tn_enzyme,I.CATALYST)),(tagged,other))
    g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Nonamplification or library incompatibility does not establish physical removal; native frequencies and fates are unassigned.'))
    sds=g.input('sds','0.2% SDS addition, 5 min',M.CONDITION)
    released=g.product('released','Tn5-released fragments',M.MOLECULAR,state=tagged_state)
    g.react('strip','Strip Tn5 with SDS',K.DEPROTEINIZATION,((tagged,I.TEMPLATE),(sds,I.CONDITION)),(released,),changes_molecule=False)
    fill=g.input('fill_polymerase','Phusion HF; 72 C 3 min',M.CATALYST)
    filled_state=Scene()
    filled_state.strand('retained_TSO',(*tagged_state.strands['retained_TSO'].segs,role('native extended opposite adapter end')))
    filled_state.strand('internal',(*tagged_state.strands['internal'].segs,role('native extended adapter boundaries')))
    filled=g.product('filled','native extension products before index amplification',M.MOLECULAR,state=filled_state)
    g.react('gap_fill','First 72 C 3 min extension',K.EXTENSION,((released,I.TEMPLATE),(fill,I.CATALYST)),(filled,))
    pcr_input=filled
    ip=g.input('index_primers','Custom S50X/N70X index primers; V.3 recommends 10 or 8 bp, actual member and width unassigned',M.OLIGO)
    ip_reverse=g.input('index_reverse','native opposite index PCR primer',M.OLIGO)
    ipol=g.input('index_polymerase','Phusion HF polymerase',M.CATALYST)
    lib=g.product('indexed','indexed cDNA library',M.LIBRARY,state=FINAL_LIBRARY)
    internal=g.product('indexed_internal','indexed internal cDNA without original transcript UMI',M.LIBRARY,state=INTERNAL_LIBRARY)
    g.react('index_pcr','Twelve index PCR cycles after 72 C extension: 98/55/72 C',K.PCR,((pcr_input,I.TEMPLATE),(ip,I.STRUCTURAL),(ip_reverse,I.STRUCTURAL),(ipol,I.CATALYST)),(lib,internal),geometry=PCRGeometry.unavailable(pcr_input.state,reason='Native tagged input, not completed indexed library, is the PCR substrate; native adapter/index primer geometry remains unvalidated.'))
    final=g.product('library_cleanup_recovered','pooled and purified retained-TSO class',M.LIBRARY,state=FINAL_LIBRARY)
    final_internal=g.product('library_cleanup_internal','pooled and purified internal class',M.LIBRARY,state=INTERNAL_LIBRARY)
    beads=g.input('library_cleanup_selection','0.6:1 PEG/Sera-Mag beads',M.SELECTION)
    waste=g.product('library_cleanup_waste','removed supernatant and washes',M.WASTE)
    g.react('library_cleanup','Pool wells and purify both molecular classes using 0.6:1 PEG/Sera-Mag beads',K.CLEANUP,((lib,I.TEMPLATE),(internal,I.TEMPLATE),(beads,I.SELECTION)),(final,final_internal,waste),changes_molecule=False)
    g.discard(waste,'Physically removed bead supernatant and ethanol washes')
    g.sequence(final_internal,SequencingPlan(library=INTERNAL_LIBRARY,limitation=SEQUENCING_UNAVAILABLE))
    g.sequence(final,SequencingPlan(library=FINAL_LIBRARY,limitation=SEQUENCING_UNAVAILABLE))
    g.declare_reviewed_scope(REVIEWED,SOURCE)
    return g.require_publishable()

REVIEWED=('rna', 'dt', 'rt', 'first', 'reverse_transcription', 'tso', 'switched', 'template_switch', 'preamp_primers', 'preamp_reverse', 'preamp_polymerase', 'amplified', 'residual_rna', 'preamp', 'transposome', 'tn_enzyme', 'tagged', 'other_tagged', 'tagment', 'index_primers', 'index_reverse', 'index_polymerase', 'indexed', 'index_pcr', 'cdna_cleanup', 'cdna_cleanup_recovered', 'cdna_cleanup_waste', 'cdna_cleanup_selection', 'sds', 'released', 'strip', 'fill_polymerase', 'filled', 'gap_fill', 'library_cleanup', 'library_cleanup_recovered', 'library_cleanup_waste', 'library_cleanup_selection', 'indexed_internal', 'library_cleanup_internal')

"""Di et al. 2020; selected 200-ng purified-RNA SHERRY branch."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]; sys.path[:0]=[str(ROOT/'lib')]
from chemdraw import Construct, Scene, Segment, oligo, feature, revcomp
from protocol_model import role
from pcr import PCRGeometry
from protocol_graph import Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
import seqprimers as sp
import nextera as nx
TITLE='SHERRY — direct tagmentation of RNA/DNA hybrids'
NOTES='01_sherry.html'
SOURCE='Defining source: <a href="https://doi.org/10.1073/pnas.1919800117">Di et al., PNAS (2020)</a>; published SI pp.2–4.'
SUMMARY='Selected purified 200-ng RNA branch: poly(T)30VN RT, hybrid tagmentation, SuperScript II/Q5 gap extension, 12-cycle indexed PCR and 1:1 bead purification.'
CAVEAT='Native gap-fill topology and TD202 index-primer member/width are unassigned. Known transferred S5/S7 sequences are retained; no complete Nextera PCR shell or fixed index length is supplied. RNA-transfer test products and the single-cell branch are distinct experiments.'
SOURCE_OLIGO_EXCEPTIONS={'oligo_dt':'Published SI p.3 prints poly(T)30VN shorthand, verified directly; guard expands placeholder VN to NN and does not match this source shorthand. The published variable VN class is retained in the displayed ordering oligo.'}
SEQ_PRIMERS=()
FINAL_LIBRARY=Construct([role('PCR-completed end A; native geometry unassigned'),Segment('S5–ME',nx.ADAPTOR_S5),role('transcript-derived insert'),Segment('opposite S7–ME',revcomp(nx.ADAPTOR_S7)),role('PCR-completed end B; native geometry unassigned'),role('sample index acquisition/placement unassigned',molecular_feature=feature('sherry_sample','sample_index','unknown'))],name='SHERRY source-limited amplifiable heterologous endpoint')
FINAL_CAPTION='Role-level heterologous endpoint; complete TD202 primers and molecular PCR derivation remain unassigned.'
SEQUENCING_INTRO='Published instruments: NextSeq500 or HiSeq4000. Read and index acquisition lengths are unassigned in this bounded source scope.'
SEQUENCING_ENDING=SEQUENCING_INTRO
SHOW_READ_CONTENT=False
READ_SPECS=tuple(sp.read_spec(n,None,'Acquisition unassigned in selected source scope','undetermined') for n in ('Read 1','Read 2','Index 1 (i7)','Index 2 (i5)'))
def seg(*args,**kw): return Segment(*args,**kw)
def rna_scene():
 sc=Scene(); sc.strand('RNA',[role('poly(A) RNA')],label='independent RNA molecules'); return sc
def hybrid_scene():
 sc=Scene(); sc.strand('RNA',[role('RNA template')],label='RNA'); sc.strand('cDNA',[Segment('poly(T)30','T'*30),Segment('VN','VN',placeholder=True),role('first-strand transcript copy')],label='cDNA 5′→3′; pairing extent unassigned'); return sc
def oligos():
    yield oligo("poly(T)30VN RT primer", [Segment("poly-T", "T" * 30),
                                        Segment("anchored VN", "VN", placeholder=True)])
    yield oligo("ME loading oligo", [Segment("mosaic-end complement", nx.ME_RC)])
    yield oligo("Adaptor A", [Segment("S5–ME", nx.ADAPTOR_S5)])
    yield oligo("Adaptor B", [Segment("S7–ME", nx.ADAPTOR_S7)])


def tagged_scene():
 sc=Scene(); sc.strand('RNA',[role('adapter-transferred RNA fragments; strand-test result')],label='RNA population, not covalently joined to cDNA'); sc.strand('cDNA',[role('tagmented cDNA with source A/B ends; individual outcome unassigned')],label='cDNA population'); return sc

def extended_scene():
 sc=Scene(); sc.strand('DNA',[role('gap-extended tagmented cDNA; actual precursor topology unassigned')],label='DNA PCR substrate'); sc.strand('RNA',[role('residual RNA and RNA-transfer products')],label='RNA retained in reaction; fate unassigned'); return sc

def loaded_scene():
 sc=Scene()
 for name,prefix in (('A',nx.S5),('B',nx.S7)):
  sc.strand(name,[Segment('prefix '+name,prefix),Segment('ME '+name,nx.ME)],label='published adaptor '+name)
  sc.anneal('loading '+name,[Segment('ME complement '+name,nx.ME_RC)],to=name,pair=('ME complement '+name,'ME '+name),label='published loading oligo')
 return sc

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 rna=g.input('rna','200 ng purified total RNA',MaterialKind.MOLECULAR,state=rna_scene(),primary=True)
 primer=g.input('oligo_dt','poly(T)30VN; no ISPCR or TSO',MaterialKind.OLIGO,state=Construct([Segment('T30','T'*30),Segment('VN','VN',placeholder=True)]))
 rt=g.input('reverse_transcriptase','SuperScript II',MaterialKind.CATALYST)
 hybrid=g.product('hybrid','RNA/cDNA reaction mixture',MaterialKind.MOLECULAR,state=hybrid_scene())
 g.react('reverse_transcribe','RT at 42°C for 1.5 h without cycling',StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(hybrid,))
 loaded=g.input('loaded_transposome','0.05 µl mixed Tn5 loaded separately with ME+A or ME+B',MaterialKind.ADAPTER,state=loaded_scene())
 tagged=g.product('tagged','tagmented RNA/cDNA mixture',MaterialKind.MOLECULAR,state=tagged_scene())
 tn5=g.input('tn5','loaded Tn5 catalytic protein',MaterialKind.CATALYST)
 g.react('tagment','Tagment hybrid at 55°C for 30 min',StepKind.TAGMENTATION,((hybrid,InputRole.TEMPLATE),(loaded,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,))
 q5=g.input('q5','Q5 High-Fidelity Master Mix',MaterialKind.CATALYST)
 extended=g.product('extended','gap-extension mixture',MaterialKind.MOLECULAR,state=extended_scene())
 g.react('gap_extend','100 U SuperScript II plus Q5: 42°C 15 min, then 70°C 15 min',StepKind.EXTENSION,((tagged,InputRole.TEMPLATE),(rt,InputRole.CATALYST),(q5,InputRole.CATALYST)),(extended,))
 primers=g.input('indexed_primers','Vazyme TD202 first indexed common primer; member unassigned',MaterialKind.OLIGO)
 opposite=g.input('opposite_primer','Vazyme TD202 opposite indexed common primer; member unassigned',MaterialKind.OLIGO)
 amp=g.product('amplified','PCR library before purification',MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 residual=g.product('residual','RNA and non-amplified molecular products; physical partition unassigned',MaterialKind.MOLECULAR)
 g.react('index_pcr','12 cycles for selected 200-ng input',StepKind.PCR,((extended,InputRole.TEMPLATE),(primers,InputRole.STRUCTURAL),(opposite,InputRole.STRUCTURAL),(q5,InputRole.CATALYST)),(amp,residual),geometry=PCRGeometry.unavailable(extended.state,reason='Actual gap-extension mixture; exact selected TD202 primer geometry is unassigned.'))
 g.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,reason='Non-amplification does not demonstrate physical removal.'))
 beads=g.input('clean_beads','VAHTS DNA Clean Beads',MaterialKind.CONDITION)
 library=g.product('library',FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 waste=g.product('purify_waste','bead-purification non-retained fraction; composition unassigned',MaterialKind.WASTE)
 g.react('purify','Purify resulting library at 1:1 bead ratio',StepKind.CLEANUP,((amp,InputRole.TEMPLATE),(beads,InputRole.CONDITION)),(library,waste))
 g.mark_terminal(Terminal(waste.id,Disposition.DISCARD,reason='Non-retained physical fraction at this source-described purification; no species-specific depletion asserted.'))
 g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=(),limitation='Study platforms known; exact cycle allocation and selected indexed primer synthesis unassigned.'))
 g.declare_reviewed_scope(('tn5','opposite_primer','rna','oligo_dt','reverse_transcriptase','reverse_transcribe','hybrid','loaded_transposome','tagged','tagment','q5','extended','gap_extend','indexed_primers','amplified','index_pcr','clean_beads','purify','library'),SOURCE)
 return g.require_publishable()

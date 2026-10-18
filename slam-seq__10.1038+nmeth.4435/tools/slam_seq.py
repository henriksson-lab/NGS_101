"""Herzog et al. 2017; source-owned metabolic conversion and QuantSeq boundary."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from base_conversion import BasePath
from chemdraw import Construct, Scene
from protocol_model import role
from protocol_graph import Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
import seqprimers as sp
def slam_path(*, labelled: bool) -> BasePath:
    """SLAM readout path: terminal C/T is RNA-sense sequencing, not first-strand cDNA G/A."""
    return BasePath(("s4U", "IAA-s4U", "C")) if labelled else BasePath(("U", "T"))


TITLE='SLAM-seq — thiol-linked metabolic RNA sequencing'
NOTES='01_slam-seq.html'
SOURCE='Defining source: <a href="https://doi.org/10.1038/nmeth.4435">Herzog et al., Nature Methods (2017)</a>, published Online Methods.'
SUMMARY='Pulse-labelled total RNA is isolated, IAA-alkylated, quenched and ethanol-precipitated before QuantSeq 3′-end library preparation.'
CAVEAT='Native QuantSeq strand ends, internal second-strand/PCR order and kit oligos remain unassigned in this source scope. No TruSeq adapter shell, index count or fixed index length is inferred. Physical IAA-s4U is an RNA modification; the T-to-C sequencing signal is a subsequent readout.'
ALKYLATION={'IAA_mM':10,'DMSO_percent':50,'phosphate_mM':50,'pH':8,'temperature_C':50,'minutes':15}
FINAL_LIBRARY=Construct([role('QuantSeq native end A; sequence unassigned'),role('3′-end DNA insert; labelled-site readout C or unconverted T'),role('QuantSeq native end B; sequence unassigned')],name='SLAM-seq source-limited QuantSeq endpoint')
SEQ_PRIMERS=()
FINAL_CAPTION='Native kit-end roles; no standard TruSeq sites assigned. T>C detection is probabilistic and does not mean every labelled position converts.'
SEQUENCING_INTRO='HiSeq2500 SR50; Fig.3 transcriptional-output experiments SR100. Index acquisition is unassigned.'
SEQUENCING_ENDING=SEQUENCING_INTRO
SHOW_READ_CONTENT=False
READ_SPECS=(sp.read_spec('Read 1','50; 100 for Fig.3','3′ transcript read; reference-U positions can report C','include'),sp.read_spec('Index 1 (i7)',None,'Index geometry and cycle allocation unassigned','undetermined'))
def mixture(label):
 sc=Scene(); sc.strand('labelled',[role(label)],label='newly synthesized RNA population'); sc.strand('other',[role('unlabelled RNA containing ordinary U')],label='pre-existing and unlabelled RNA retained'); return sc

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 cells=g.input('cells','mESC population',MaterialKind.BIOLOGICAL,primary=True)
 s4u=g.input('s4u','100 µM s4U in medium; exchange every 3 h during pulse',MaterialKind.CONDITION)
 pulse=g.product('pulse','cells containing labelled and unlabelled RNA',MaterialKind.BIOLOGICAL,state=mixture('new RNA with stochastic s4U incorporation'))
 rnap=g.input('rnap','cellular RNA polymerase',MaterialKind.CATALYST)
 g.react('pulse_label','Metabolic pulse; selected duration is experiment-specific',StepKind.TRANSCRIPTION,((cells,InputRole.TEMPLATE),(s4u,InputRole.CONDITION),(rnap,InputRole.CATALYST)),(pulse,))
 trizol=g.input('trizol','TRIzol; 0.1 mM DTT during isopropanol precipitation',MaterialKind.CONDITION)
 total=g.product('total_rna','total RNA resuspended in 1 mM DTT',MaterialKind.MOLECULAR,state=mixture('s4U-bearing RNA'))
 waste=g.product('isolate_rna_waste','TRIzol extraction/precipitation non-retained material; composition unassigned',MaterialKind.WASTE)
 g.react('isolate_rna','Harvest and isolate total RNA; retain labelled and unlabelled species',StepKind.CLEANUP,((pulse,InputRole.TEMPLATE),(trizol,InputRole.CONDITION)),(total,waste))
 g.mark_terminal(Terminal(waste.id,Disposition.DISCARD,reason='Non-retained physical fraction at this source-described purification; no species-specific depletion asserted.'))
 iaa=g.input('iaa','10 mM IAA, 50% DMSO, 50 mM sodium phosphate pH8',MaterialKind.CONDITION)
 alk=g.product('alkylated','IAA-reacted total RNA',MaterialKind.MOLECULAR,state=mixture('RNA bearing carboxyamidomethyl-s4U; conversion yield unassigned'))
 g.react('alkylate','Typical 5 µg total RNA: 50°C for 15 min',StepKind.CHEMICAL_CONVERSION,((total,InputRole.TEMPLATE),(iaa,InputRole.CONDITION)),(alk,))
 dtt=g.input('dtt','excess DTT',MaterialKind.CONDITION)
 quenched=g.product('quenched','quenched total RNA reaction',MaterialKind.MOLECULAR,state=alk.state)
 g.react('quench','Quench excess IAA with DTT',StepKind.CHEMICAL_CONVERSION,((alk,InputRole.TEMPLATE),(dtt,InputRole.CONDITION)),(quenched,))
 ethanol=g.input('ethanol','ethanol precipitation reagents',MaterialKind.CONDITION)
 recovered=g.product('recovered','precipitated total RNA for QuantSeq',MaterialKind.MOLECULAR,state=alk.state)
 waste=g.product('precipitate_waste','ethanol precipitation supernatant; molecular composition unassigned',MaterialKind.WASTE)
 g.react('precipitate','Ethanol precipitate RNA after alkylation',StepKind.CLEANUP,((quenched,InputRole.TEMPLATE),(ethanol,InputRole.CONDITION)),(recovered,waste))
 g.mark_terminal(Terminal(waste.id,Disposition.DISCARD,reason='Non-retained physical fraction at this source-described purification; no species-specific depletion asserted.'))
 kit=g.input('quantseq','Lexogen QuantSeq 3′ mRNA library reagents',MaterialKind.CONDITION)
 cdna_sc=Scene(); cdna_sc.strand('RNA',[role('alkylated poly(A) RNA template')],label='RNA'); cdna_sc.strand('cDNA',[role('oligo(dT)-initiated first-strand cDNA; G can be incorporated opposite IAA-s4U, A opposite U')],label='cDNA; conversion is probabilistic')
 cdna=g.product('cdna','first-strand RNA/cDNA reaction',MaterialKind.MOLECULAR,state=cdna_sc)
 rtp=g.input('rt_primer','QuantSeq oligo(dT) RT primer; native ends unassigned',MaterialKind.OLIGO)
 rtenz=g.input('rt_enzyme','QuantSeq RT enzyme',MaterialKind.CATALYST)
 g.react('reverse_transcribe','QuantSeq oligo(dT) reverse transcription',StepKind.REVERSE_TRANSCRIPTION,((recovered,InputRole.TEMPLATE),(rtp,InputRole.STRUCTURAL),(rtenz,InputRole.CATALYST)),(cdna,))
 library=g.product('library',FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 residual=g.product('kit_residual','RNA templates and kit byproducts; physical fate unassigned',MaterialKind.MOLECULAR)
 g.react('kit_construct','Native QuantSeq library-construction boundary after verified RT',StepKind.MAGIC,((cdna,InputRole.TEMPLATE),(kit,InputRole.CONDITION)),(library,residual),reason='Inspected published Methods names QuantSeq and establishes RT readout, but does not establish its native second-strand/PCR order or intermediate strands. Version-matched kit guide remains deferred; no speculative internal route is substituted.')
 g.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,reason='The inspected scope does not assign RNA destruction or kit purification fractions.'))
 g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=(),limitation='SR50 or Fig.3 SR100, exact native primer and index acquisition unassigned.'))
 g.declare_reviewed_scope(('rnap','rt_primer','rt_enzyme','cells','s4u','pulse_label','pulse','trizol','isolate_rna','total_rna','iaa','alkylate','alkylated','dtt','quench','quenched','ethanol','precipitate','recovered','quantseq','reverse_transcribe','cdna','kit_construct','kit_residual','library'),SOURCE)
 return g.require_publishable()

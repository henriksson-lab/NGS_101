"""Protocol-owned source, oligos, molecular workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, Workflow, feature, revcomp
from topology import HairpinTopology
import illumina as il
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, Terminal, Disposition
from protocol_model import chain, molecule, role
import seqprimers as sp


SS2_HAIRPIN = Segment("Structure-seq2 hairpin donor", "TGAAGAGCCTAGTCGCTGTTCANNNNNNCTGCCCATAGAG", "r1", placeholder=True)

SS2_RT_SEQUENCE = "CAGACGTGTGCTCTTCCGATC" + "N" * 6

SS2_FORWARD_SEQUENCE = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTTGAACAGCGACTAGGCTCTTCA"

SS2_REVERSE_SEQUENCE = "CAAGCAGAAGACGGCATACGAGAT" + "N" * 6 + "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC"

# Reagent states for the graph inputs.  The RT primer is drawn only in Supplementary
# Figure S1 (step 1A, 3'-NNNNNNCTAGCCTTCTCGTGTGCAGAC-5', 27 nt as stated in Figure S4);
# the hairpin donor and both PCR primers are printed in the Methods.
SS2_RT_PRIMER = Construct((
    Segment("RT-primer fixed handle", SS2_RT_SEQUENCE[:-6], "r2"),
    Segment("random hexamer", SS2_RT_SEQUENCE[-6:],
            note="degenerate priming bases; they anneal at random RNA positions"),
), name="Structure-seq2 RT primer")

SS2_HAIRPIN_DONOR = Construct((
    Segment("hairpin 5' arm and loop", "TGAAGAGCCTAGTCGCTGTTCA", "r1",
            note="5' phosphate (Methods: '5'-p')"),
    Segment("degenerate six-base region", "NNNNNN",
            note="pairs with diverse cDNA 3' ends"),
    Segment("hairpin 3' arm", "CTGCCCATAGAG", "r1",
            note="3' C3 spacer (Methods: '3'-Spacer')"),
), name="Structure-seq2 hairpin donor")

SS2_REVERSE_PRIMER = Construct((
    Segment("P7", il.P7, "p7"),
    Segment("sample barcode", "NNNNNN", "cbc", placeholder=True,
            feature=feature("sample_index", "sample_index", "unknown",
                            note="Methods state a 6-8-nt barcode; the six-base member is shown")),
    Segment("TruSeq Read 2 / Index 1 arm", SS2_REVERSE_SEQUENCE[-33:], "r2"),
), name="Structure-seq2 indexed reverse primer (six-base member)")

if SS2_REVERSE_PRIMER.top() != SS2_REVERSE_SEQUENCE:
    raise AssertionError("indexed reverse primer state differs from its printed sequence")

# Read by tools/selftest_source_oligos.py.
SOURCE_OLIGO_EXCEPTIONS = {
    "rt_primer": ("Printed only as artwork in Supplementary Figure S1 (File006 p. S2, "
                  "step 1A/1B: 3'-NNNNNNCTAGCCTTCTCGTGTGCAGAC-5'); the extracted text "
                  "has no bases, so it was transcribed from the rendered figure."),
}

SS2_FORWARD = Segment("Structure-seq2 TruSeq forward", SS2_FORWARD_SEQUENCE, "p5")

SS2_REVERSE = Segment("Structure-seq2 indexed reverse", SS2_REVERSE_SEQUENCE, "p7", placeholder=True,
                      feature=feature("sample_index", "sample_index", "unknown"))

SS2_CUSTOM = Segment("Structure-seq2 custom sequencing primer", "TCTTCCGATCTTGAACAGCGACTAGGCTCTTCA", "r1")

SS2_FINAL_LIBRARY = Construct((
    Segment("forward PCR primer", SS2_FORWARD_SEQUENCE, "p5"),
    Segment("stop-defined cDNA", "X" * 42, placeholder=True),
    Segment("copied RT random-hexamer complement", "N" * 6, placeholder=True, note="Retained primer-derived N6, distinct from the biological stop-defined extension"),
    Segment("Index 1 primer site", revcomp(SS2_REVERSE_SEQUENCE[-33:]), "r2"),
    Segment("sample index complement", "N" * 6, "cbc", placeholder=True,
            feature=feature("sample_index", "sample_index", "unknown")),
    Segment("P7 reverse complement", il.P7_RC, "p7"),
), name="Structure-seq2 six-base-index library")

SS2_SEQ_PRIMERS = (
    sp.custom("Read 1", "Structure-seq2 custom Read 1", SS2_CUSTOM.top,
              "Ritchey et al. Supplementary Figure S1"),
    sp.TRUSEQ["I1"],
)


TITLE = 'Structure-seq2 — low-bias in-vivo RNA structure probing'
NOTES = '01_structure_seq2.html'
SOURCE = 'Ritchey et al., Nucleic Acids Research (2017); <a href="https://doi.org/10.1093/nar/gkx533">doi:10.1093/nar/gkx533</a>.'
SUMMARY = 'Random-primed DMS-stop cDNA receives a degenerate hairpin adapter by T4 DNA ligase; PAGE or biotin–streptavidin purification removes the primer–adapter by-product.'
CAVEAT = 'This graph follows the PAGE purification branch. The biotin/streptavidin variation and the original Structure-seq comparison libraries are separate experimental variants.'
FINAL_LIBRARY = SS2_FINAL_LIBRARY
FINAL_CAPTION = 'Representative six-base-index PCR product. The custom primer abuts the unique insert, bypassing the 22-nt low-diversity prefix.'
SEQ_PRIMERS = SS2_SEQ_PRIMERS
ACQUIRED_ROLES = ('Read 1', 'Index 1 (i7)')
READ_LENGTHS = {'Read 1': 150}
SHOW_READ_CONTENT = False
READ_SPECS = (sp.read_spec('Read 1',150,'Custom-primer cDNA sequence','include'),sp.read_spec('Index 1 (i7)','unassigned','Sample barcode; actual index acquisition length is unassigned','demultiplex'))
READ_SPECS_NOTE = 'Figure S1 shows an index primer. The six-base barcode in the drawn library is a representative physical member of the stated six-to-eight-base family; it does not establish an acquired six-cycle index read.'
READ_USES = {'Read 1': 'include', 'Index 1 (i7)': 'demultiplex'}
SEQUENCING_INTRO = 'The custom primer and standard multiplexing-index primer both land on source-exact PCR-derived sites.'
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''

def oligos():
    return ('RT primer — 5′-CAGACGTGTGCTCTTCCGATCNNNNNN-3′', 'hairpin donor — 5′-/5Phos/TGAAGAGCCTAGTCGCTGTTCANNNNNNCTGCCCATAGAG/3SpC3/-3′', 'TruSeq forward — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTTGAACAGCGACTAGGCTCTTCA-3′', 'indexed reverse — 5′-CAAGCAGAAGACGGCATACGAGAT[6–8-nt barcode]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-3′', 'custom sequencing primer — 5′-TCTTCCGATCTTGAACAGCGACTAGGCTCTTCA-3′')


DONOR_FOLD=HairpinTopology('Structure-seq2 donor',SS2_HAIRPIN_DONOR.top(),0,18,4)
TOPOLOGY_PANELS=(('Source-declared donor fold',DONOR_FOLD.folded_panel('Four-base TGAA/TTCA stem from Supplementary Figure S1. The donor has a 5-prime phosphate, six degenerate pairing bases after the stem and a blocked 3-prime C3 end. A compatible fold is not a thermodynamic prediction.')),)

# Exact source-written substrates are separate from PCR output orientation.
CDNA_PARTS=(*SS2_RT_PRIMER.segments,role('stop-defined cDNA extension'))
CDNA=molecule('single-strand cDNA',CDNA_PARTS)
LIGATED=molecule('ligated single-strand cDNA',(*CDNA_PARTS,*SS2_HAIRPIN_DONOR.segments),junctions=(('stop-defined cDNA extension',"hairpin 5' arm and loop",'T4 DNA ligation'),))
LIGATED.strands['molecule'].mod3='C3 spacer'
# The six-base index is a selected representative, not the study-specific value.
PCR_GEOMETRY=PCRGeometry.unavailable(LIGATED,reason='Actual hairpin-ligated ssDNA enters PCR. The four-base donor stem is separately checked; its degenerate substrate pairing and first-cycle conversion are not reconstructed in this pass; the completed library is not used as a template.')

def _purify(g,current,identifier,label,kind=MaterialKind.MOLECULAR):
    sel=g.input(identifier+'_selection',label,MaterialKind.SELECTION)
    out=g.product(identifier+'_recovered',label+' recovered molecules',kind,state=current.state)
    waste=g.product(identifier+'_waste','excluded gel fractions and recovery waste',MaterialKind.WASTE)
    g.react(identifier,label,StepKind.CLEANUP,((current,InputRole.TEMPLATE),(sel,InputRole.SELECTION)),(out,waste),changes_molecule=False)
    g.discard(waste,'Physically excluded gel fractions, filtration waste and precipitation supernatants')
    return out

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    start=g.input('input','native rice-shoot RNA in tissue',MaterialKind.BIOLOGICAL,state=molecule('RNA',(role('native cellular RNA'),)),primary=True)
    dms=g.input('dms','DMS treatment with matched minus-DMS control; DTT quench',MaterialKind.CONDITION)
    probed=g.product('probed','RNA carrying DMS adducts',MaterialKind.MOLECULAR,state=molecule('probed RNA',(role('RNA with physical DMS adducts on susceptible A/C'),)))
    g.react('probe','Probe in vivo with DMS; quench with DTT and wash',StepKind.CHEMICAL_CONVERSION,((start,InputRole.TEMPLATE),(dms,InputRole.CONDITION)),(probed,))
    rna=_purify(g,probed,'rna_purify','Plant RNA extraction and one round poly(A) selection')
    primer=g.input('rt_primer','Source Figure S1 random-hexamer RT primer',MaterialKind.OLIGO,state=SS2_RT_PRIMER)
    rt=g.input('rt_enzyme','SuperScript III; 55 C 50 min after 90 C denaturation',MaterialKind.CATALYST)
    hybrid=Scene();hybrid.strand('RNA',(role('DMS-probed RNA'),));hybrid.strand('cDNA',CDNA_PARTS)
    first=g.product('first','RNA/cDNA population after stop-defined RT',MaterialKind.MOLECULAR,state=hybrid)
    g.react('rt','Random-prime stop-defined cDNA',StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(first,))
    rnase=g.input('rnase','RNase H 37 C 20 min',MaterialKind.CATALYST)
    cdna=g.product('cdna','single-strand cDNA and residual free primer',MaterialKind.MOLECULAR,state=CDNA)
    hydro=g.product('rna_hydrolysate','RNA cleavage products',MaterialKind.MOLECULAR,state=molecule('RNA fragments',(role('RNase H cleavage products'),)))
    g.react('digest_rna','Cleave residual RNA with RNase H',StepKind.DIGESTION,((first,InputRole.TEMPLATE),(rnase,InputRole.CATALYST)),(cdna,hydro))
    g.mark_terminal(Terminal(hydro.id,Disposition.UNRESOLVED,'Hydrolysis is not physical separation; detailed recovery of RNA fragments during PAGE is unassigned.'))
    recovered=_purify(g,cdna,'rt_page','Pool two RT reactions; denaturing PAGE above 50 nt, crush/soak and ethanol recovery')
    donor_state=Scene();donor_state.strand('hairpin donor',SS2_HAIRPIN_DONOR.segments,mod5='p',mod3='C3 spacer')
    donor=g.input('hairpin','5-prime-phosphorylated, 3-prime-C3-blocked hairpin donor',MaterialKind.ADAPTER,state=donor_state)
    ligase=g.input('ligase','T4 DNA ligase, betaine and PEG 8000',MaterialKind.CATALYST)
    ligated=g.product('ligated','hairpin-ligated ssDNA',MaterialKind.MOLECULAR,state=LIGATED)
    other=g.product('ligate_other','unreacted donor and residual primer-donor by-product',MaterialKind.MOLECULAR,state=molecule('residual reaction classes',(role('unassigned residual donor/primer-donor population'),)))
    g.react('ligate','Ligate the cDNA 3-prime end to the donor 5-prime phosphate',StepKind.LIGATION,((recovered,InputRole.TEMPLATE),(donor,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(ligated,other),geometry=AdapterLigationGeometry.unavailable('Source-known hairpin donor sequence and modifications are shown. The four-base donor stem is separately checked; native degenerate substrate pairing/end geometry has not been reconstructed in this pass.'))
    g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'The following above-90-nt PAGE targets donor/by-product exclusion; exact residual carryover is unassigned.'))
    selected=_purify(g,ligated,'ligate_page','Denaturing PAGE above 90 nt; recover ligated cDNA')
    fwd=g.input('forward','Illumina TruSeq forward primer',MaterialKind.OLIGO,state=SS2_FORWARD)
    rev=g.input('reverse','Indexed reverse primer, selected six-base member',MaterialKind.OLIGO,state=SS2_REVERSE_PRIMER)
    pol=g.input('polymerase','Q5 High Fidelity DNA polymerase; PCR cycle test selects cycle count',MaterialKind.CATALYST)
    library=g.product('amplified','PCR-derived library before final PAGE',MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    g.react('pcr','Amplify with published forward and indexed reverse primers',StepKind.PCR,((selected,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),(rev,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(library,),geometry=PCRGeometry.unavailable(selected.state,reason=PCR_GEOMETRY.statement))
    final=_purify(g,library,'final_page','Denaturing PAGE 220-600 nt; exclude residual amplified by-product and primer dimers',MaterialKind.LIBRARY)
    g.sequence(final,SequencingPlan(library=FINAL_LIBRARY,primers=SEQ_PRIMERS,acquired_roles=ACQUIRED_ROLES,read_lengths=READ_LENGTHS,read_uses=READ_USES))
    g.declare_reviewed_scope(REVIEWED,SOURCE)
    return g.require_publishable()

REVIEWED=('input','dms','probed','probe','rna_purify_selection','rna_purify_recovered','rna_purify_waste','rna_purify','rt_primer','rt_enzyme','first','rt','rnase','cdna','rna_hydrolysate','digest_rna','rt_page_selection','rt_page_recovered','rt_page_waste','rt_page','hairpin','ligase','ligated','ligate_other','ligate','ligate_page_selection','ligate_page_recovered','ligate_page_waste','ligate_page','forward','reverse','polymerase','amplified','pcr','final_page_selection','final_page_recovered','final_page_waste','final_page')

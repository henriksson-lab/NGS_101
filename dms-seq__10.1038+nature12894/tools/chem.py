"""Protocol-owned source, oligos, molecular workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Segment, Workflow, feature, revcomp
from circular import CircularizationGeometry
import illumina as il
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role
import seqprimers as sp
from rna_end_readout import (StopReadout, BaseConversionReadout, SelectiveCleavageReadout,
                             _finish_pcr, _apply_linear_step)

DMS_LINKER = Segment("microRNA cloning linker-1", "CTGTAGGCACCATCAAT", "r2")


TITLE = 'DMS-seq — in-vivo RNA structure probing'
NOTES = '01_dms_seq.html'
SOURCE = 'Rouskin et al., Nature (2014); <a href="https://doi.org/10.1038/nature12894">doi:10.1038/nature12894</a>.'
SUMMARY = 'DMS marks solvent-accessible A and C in living cells; reverse-transcriptase termination sites report nucleotide accessibility after adapter ligation and cDNA circularization.'
CAVEAT = ''
FINAL_LIBRARY = Construct((role('P5'), role('DMS-stop cDNA'), role('sample index', tag='cbc', molecular_feature=feature('sample_index', 'sample_index', 'unknown')), role('P7')), name='DMS-seq library')
FINAL_CAPTION = 'Role-level PCR product; the paper publishes linker-1 but not the adapter PCR or sequencing-primer oligos.'
SEQ_PRIMERS = ()
READ_LENGTHS = None
SEQUENCING_INTRO = ''
SEQUENCING_ENDING = 'The defining experiment sequenced RT-termination libraries on an Illumina HiSeq 2000.'
SEQUENCING_UNAVAILABLE = 'Exact PCR and sequencing-primer oligos are not reproduced in the defining paper.'

def oligos():
    return ('microRNA cloning linker-1 — 5′-/5rApp/CTGTAGGCACCATCAAT/3ddC/-3′',)

def workflow() -> Workflow:
    return chain(molecule("poly(A)-selected RNA", (role("native cellular RNA structure"),)), (
        ("Probe accessible A and C with DMS",
         molecule("DMS-modified RNA", (role("RNA"), role("methylated exposed A/C"))),
         "Parallel in-vivo, refolded in-vitro and denatured samples separate cellular structure from intrinsic reactivity."),
        ("Fragment with Zn²⁺, select 60–80 nt, and repair 3′ phosphates",
         molecule("end-repaired DMS RNA fragment", (role("60–80-nt modified RNA"), role("3-prime hydroxyl"))),
         "T4 PNK resolves the 3′ phosphates left by chemical fragmentation."),
        ("Ligate preadenylated cloning linker-1",
         molecule("linker-ligated DMS RNA", (role("modified RNA fragment"), DMS_LINKER),
                  junctions=(("modified RNA fragment", "microRNA cloning linker-1", "Rnl2tr K227Q ligation"),)),
         "The linker is 5′ preadenylated and 3′ ddC-blocked, so ligation proceeds without ATP."),
        ("Reverse-transcribe and isolate prematurely terminated cDNA",
         molecule("DMS-stop cDNA", (role("adapter-primed cDNA ending before modified A/C"),)),
         "A 25–50-nt cDNA size selection enriches Superscript III termination products."),
        ("Circularize cDNA and add Illumina adapters by PCR",
         molecule("DMS-seq library", (role("P5"), role("DMS-stop cDNA"), role("sample index"), role("P7"))),
         "CircLigase precedes eight to ten PCR cycles and HiSeq 2000 sequencing."),
    ))


PCR_GEOMETRY = None
ACQUIRED_ROLES = ()
READ_USES = None

RECIPE = (('probe', 'Probe accessible A and C with DMS', StepKind.CHEMICAL_CONVERSION), ('fragment', 'Fragment DMS-modified RNA', StepKind.FRAGMENTATION), ('repair', 'Repair fragmentation-derived RNA ends', StepKind.PHOSPHORYLATION), ('ligate', 'Ligate the preadenylated cloning linker', StepKind.LIGATION), ('rt', 'Reverse-transcribe to DMS-dependent stops', StepKind.REVERSE_TRANSCRIPTION), ('circle', 'Circularize the stop-defined cDNA', StepKind.CIRCULARIZATION))

def protocol_graph() -> ProtocolGraph:
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    start = graph.input("input", TITLE + " source RNA", MaterialKind.BIOLOGICAL,
                        state=molecule("source RNA", (role("input RNA"),)), primary=True)
    reviewed = ["input"]
    current = start
    for identifier, label, kind in RECIPE:
        current = _apply_linear_step(graph, current, identifier, label, kind, reviewed)
    _finish_pcr(
        graph, "library", current, FINAL_LIBRARY,
        SEQUENCING_ENDING or FINAL_CAPTION, reviewed, label=TITLE,
        geometry=PCR_GEOMETRY,
        primers=SEQ_PRIMERS,
        acquired_roles=ACQUIRED_ROLES,
        read_lengths=READ_LENGTHS,
        read_uses=READ_USES)
    graph.declare_reviewed_scope(reviewed, SOURCE)
    return graph.require_publishable()

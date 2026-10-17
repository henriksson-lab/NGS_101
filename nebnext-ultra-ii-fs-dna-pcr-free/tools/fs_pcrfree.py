"""NEBNext Ultra II FS DNA PCR-free, E7430/E7435, manual Version 2.0_6/25 (inspected 2026-10-17).

Adaptor architecture: NEB E7395/E7874/E7876/E7878 manual Version 4.0_6/26.
A previously cited E7430 "v3.0 9/25" was not located; NEB serves Version 2.0_6/25,
whose revision history lists only 1.0 (9/21) and 2.0 (6/25)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
from adapter_ligation import reported_adapter_end
from chemdraw import Construct, Segment
from endprep import dA_tailed_scene, repair_and_dA_tail
from nebnext import STANDARD_PRIMERS, indexed_library
from protocol_graph import Highlight, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind

TITLE="NEBNext Ultra II FS DNA PCR-free"
NOTES="01_nebnext-ultra-ii-fs-pcr-free.html"
SOURCE='Commercial protocol: <a href="https://www.neb.com/en-us/products/e7435-nebnext-ultra-ii-fs-dna-pcr-free-library-prep-with-sample-purification-beads">NEB E7430/E7435 manual Version 2.0_6/25</a>; adaptor read layout from <a href="https://www.neb.com/en/-/media/nebus/files/manuals/manuale7395_e7874_e7876_e7878.pdf">NEB E7395/E7874/E7876/E7878 manual Version 4.0_6/26</a>.'
SUMMARY="Intact DNA is enzymatically fragmented, end-repaired and dA-tailed in one reaction, then receives full-length dual-index UMI adaptors. There is no PCR."
UMI_NT=12
I7_READ_NT=20
# High input (100-500 ng), manual section 1.1.5: (min at 37 C, 37, min at 65 C, 65); heated lid 75 C.
FS_PROGRAMS={350:(10,37,30,65),450:(8,37,30,65)}
# Second size-selection bead volume added to the 100 ul eluate, section 1.3.14.
SIZE_SELECT_BEADS_UL={350:65,450:58}
# Low input (50-99 ng), section 2.1.5 / 2.2: 300-bp program, adaptor diluted 1:10 to 2 uM, one 50 ul bead cleanup only (2.3).
LOW_INPUT_FS_PROGRAM={300:(20,37,30,65)}
DNA=Construct([Segment("genomic DNA", "X"*50, placeholder=True)], name="intact DNA")
END_PREP=repair_and_dA_tail(Construct([Segment("FS fragment", "X"*42, placeholder=True)], name="FS fragment"))
ADAPTOR_END=reported_adapter_end(
    "full-length UDI-UMI adaptor", kind="cohesive", polarity="3-prime",
    terminal_sequence="T", substrate_scenes=(dA_tailed_scene(END_PREP),),
    statement="The E7430/E7435 workflow specifies complementary dA/dT ligation ends; the indexed UMI adaptor body is proprietary.")
ADAPTOR=Construct([Segment("full-length UDI-UMI adaptor", "X"*18, "cbc", placeholder=True, inferred=True,
    note="current sequence proprietary; 8-base i5, 8-base i7 and 12-base UMI are documented")], name="pre-annealed adaptor")
STEPS=(("Fragment, repair ends and add dA", dA_tailed_scene(END_PREP), "The single FS reaction couples time-dependent fragmentation to end repair, 5′ phosphorylation and one-base 3′ dA tails."),
       ("Ligate full-length UDI-UMI adaptors", ADAPTOR, "The adaptor is already flow-cell competent: it carries both unique indices and a 12-base UMI. Dotted bases remain proprietary.", ADAPTOR_END),
       ("Clean up and size-select without amplification", Construct([Segment("adapter-ligated molecules only", "X"*50, placeholder=True)], name="PCR-free pool"), "A 50 µl bead cleanup is followed by a second bead binding (65 µl for 350 bp, 58 µl for 450 bp); supernatants are discarded. No polymerase amplification completes or copies this library."))
SEQ_PRIMERS=STANDARD_PRIMERS; REQUIRED_ROLES=tuple(p.role for p in SEQ_PRIMERS)
READOUT="Read 1 and Read 2 sequence the insert (user-defined lengths); Index 2 reads the 8-base i5. NEB specifies a 20-cycle Index 1 (i7) read when the UMI is used, i.e. the 8-base i7 plus the 12-base UMI; the drawn order of i7 before UMI is inferred, because the inspected text gives only the cycle total."
FINAL_CAPTION="Full-length PCR-free molecule drawn as a fully paired duplex for primer landing. NEB states that PCR-free libraries keep single-stranded adaptor regions (forked ends), which this drawing does not show. The adapter sequence was withdrawn as proprietary by NEB; documented functional blocks and computed primer sites are retained without inventing bases."
def final_library(): return indexed_library(Construct([Segment("DNA insert", "X"*42, placeholder=True)]), "Ultra II FS PCR-free library", umi_nt=UMI_NT, inferred_adaptors=True)

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    dna=g.input("dna","intact genomic DNA",MaterialKind.BIOLOGICAL,state=DNA,primary=True)
    fs=g.input("fs_mix","Ultra II FS fragmentation/end-prep enzyme mix",MaterialKind.CATALYST)
    fragmented=g.product("fragmented","FS-fragmented DNA",MaterialKind.MOLECULAR,state=END_PREP.insert)
    g.react("fragment","Fragment DNA",StepKind.FRAGMENTATION,((dna,InputRole.TEMPLATE),(fs,InputRole.CATALYST)),(fragmented,))
    repaired=g.product("repaired","end-repaired FS fragments",MaterialKind.MOLECULAR,state=END_PREP.insert)
    g.react("end_repair","Repair fragmented DNA ends",StepKind.END_REPAIR,((fragmented,InputRole.TEMPLATE),(fs,InputRole.CATALYST)),(repaired,))
    tailed=g.product("tailed","fragmented, repaired and dA-tailed DNA",MaterialKind.MOLECULAR,state=dA_tailed_scene(END_PREP))
    g.react("da_tail","Add one 3-prime dA to each repaired end",StepKind.TAILING,((repaired,InputRole.TEMPLATE),(fs,InputRole.CATALYST)),(tailed,))
    adapter=g.input("adapter","full-length UDI-UMI adapter",MaterialKind.ADAPTER,state=ADAPTOR,inferred=True,source=SOURCE+" (adapter body proprietary)")
    ligase=g.input("ligase","NEBNext Ultra II Ligation Master Mix + Ligation Enhancer (20 °C, 15 min)",MaterialKind.CATALYST)
    ligated=g.product("ligated","adapter-ligated PCR-free molecules",MaterialKind.LIBRARY,state=final_library(),inferred=True,source=SOURCE+" (adapter body proprietary)")
    g.react("ligate","Ligate full-length UDI-UMI adapters",StepKind.LIGATION,((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(ligated,),geometry=ADAPTOR_END)
    beads=g.input("beads","50 µl SPRIselect or NEBNext Sample Purification Beads (after adding 40 µl 0.1X TE)",MaterialKind.SELECTION)
    cleaned=g.product("cleaned","bead-bound adapter-ligated DNA, eluted in 102 µl 0.1X TE",MaterialKind.LIBRARY,state=final_library(),inferred=True,source=SOURCE+" (adapter body proprietary)")
    cleanup_waste=g.product("cleanup_waste","first-cleanup supernatant (unwanted DNA)",MaterialKind.WASTE)
    g.react("cleanup","Bead cleanup of the ligation reaction",StepKind.SELECTION,((ligated,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(cleaned,cleanup_waste),changes_molecule=False)
    g.discard(cleanup_waste,"supernatant is discarded after magnetic separation (manual 1.3.6)")
    size_beads=g.input("size_beads","second bead addition to 100 µl eluate: 65 µl (350 bp) or 58 µl (450 bp)",MaterialKind.SELECTION)
    library=g.product("library","size-selected Ultra II FS PCR-free library (eluted in 22 µl 0.1X TE)",MaterialKind.LIBRARY,state=final_library(),inferred=True,source=SOURCE+" (adapter body proprietary)")
    waste=g.product("size_waste","size-selection supernatant (unwanted DNA)",MaterialKind.WASTE)
    g.react("size_select","Bead size selection without amplification",StepKind.SELECTION,((cleaned,InputRole.TEMPLATE),(size_beads,InputRole.SELECTION)),(library,waste),changes_molecule=False)
    g.discard(waste,"supernatant is discarded after magnetic separation (manual 1.3.17)")
    g.sequence(library,SequencingPlan(library=final_library(),primers=tuple(SEQ_PRIMERS),acquired_roles=REQUIRED_ROLES,read_lengths={"Index 1 (i7)":I7_READ_NT}))
    g.annotate(Highlight("size_select","Forked adaptor ends remain single-stranded",
        "E7430 Version 2.0_6/25 Figures 1 and 2 state that PCR-free libraries migrate slower because of single-strand regions of the adaptors; the UDI-UMI manual Figure 1 draws the PCR-free product with unpaired adaptor arms. No PCR converts these ends to full duplex. The final drawing below is a fully paired simplification used only to show primer landing."))
    g.annotate(Highlight("fragment","Low-input branch not drawn",
        "For 50–99 ng (manual Section 2) the FS program is 20 min at 37 °C then 30 min at 65 °C for ~300 bp inserts, the adaptor is diluted 1:10 to 2 µM, and only the single 50 µl bead cleanup is performed (elution in 17 µl). The diagram follows the high-input 350/450 bp path. Whether the 65 °C step adds dA tails or inactivates enzymes is not stated; the end-repair/dA-tail split is inferred from the FS kit's stated end-prep function."))
    g.declare_reviewed_scope(["dna","fs_mix","fragment","end_repair","da_tail","adapter","ligase","ligate","beads","cleanup","size_beads","size_select","library"],SOURCE)
    return g.require_publishable()

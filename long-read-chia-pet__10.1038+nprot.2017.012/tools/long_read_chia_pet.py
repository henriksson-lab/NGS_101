"""Li et al. 2017 long-read ChIA-PET molecular model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature
from proximity_linker import ThreePrimeBridgeLinker
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role

def seg(name,top,tag=None,**kw):return Segment(name,top,tag,**kw)

TITLE = "Long-read ChIA-PET — bridge-linked chromatin contacts"
NOTES = "01_long-read-chia-pet.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/nprot.2017.012">'
          'Li et al., <i>Nature Protocols</i> (2017)</a>.')
SUMMARY = ("Protein-associated chromatin fragments are repaired and dA-tailed, joined in "
           "proximity through a double-ended biotin bridge, then tagmented; streptavidin "
           "retains only fragments that still contain the contact junction.")
CAVEAT = ("Original Methods (PMC5537732, workflow overview) place crosslink reversal and DNA release before Tn5 tagmentation; the graph now includes that transition. "
          "The bridge-linker bases and modification are published. INFERRED — canonical Nextera components are a conditional reference for the named FC-121-1031 kit, not exact historical physical primers or index widths.")

BRIDGE_F = "CGCGATATCTTATCTGACT"
BRIDGE_R = "GTCAGATAAGATATCGCGT"
BRIDGE = ThreePrimeBridgeLinker("biotin bridge linker", BRIDGE_F, BRIDGE_R,
                                overhang="T", modified_top_position=9)

def contact_segments() -> list[Segment]:
    core = BRIDGE.top_core
    return [
        seg("left contact tag", "X" * 28, placeholder=True),
        seg("left A:T ligation junction", "A", bottom="T"),
        seg("bridge core left", core[:9], "cbc"),
        seg("internal biotin-dT", core[9], "umi"),
        seg("bridge core right", core[10:], "cbc"),
        seg("right T:A ligation junction", "T", bottom="A"),
        seg("right contact tag", "X" * 28, placeholder=True),
    ]

CONTACT = Construct(contact_segments(), name="DNA–bridge–DNA proximity product")
# Explicit conditional reference for the named FC-121-1031 Nextera DNA kit.
# Source gives i5/i7 primer roles, not physical members/widths or loaded oligos.
FINAL_LIBRARY=Construct([
 seg("conditional P5",il.P5,"p5",inferred=True),
 role("i5 index; member/width unassigned",tag="cbc",inferred=True,molecular_feature=feature("chia_i5","sample_index","unknown")),
 seg("conditional S5",nx.S5,"s5",inferred=True),seg("conditional left mosaic end",nx.ME,"me",inferred=True),
 *contact_segments(),
 seg("conditional right mosaic end",nx.ME_RC,"me",inferred=True),seg("conditional S7 reverse complement",nx.S7_RC,"s7",inferred=True),
 role("i7 reverse complement; member/width unassigned",tag="cbc",inferred=True,molecular_feature=feature("chia_i7","sample_index","unknown")),
 seg("conditional P7 reverse complement",il.P7_RC,"p7",inferred=True)],name="Conditional source-kit-compatible long-read ChIA-PET reference")
SEQ_PRIMERS=(sp.NEXTERA["R1"],sp.NEXTERA["R2"])
if sp.verify(FINAL_LIBRARY,SEQ_PRIMERS,required_roles=("Read 1","Read 2")):
 raise ValueError("Conditional Nextera reference primers do not land")
SHOW_READ_CONTENT=False
READ_SPECS=(sp.read_spec("Read 1","150 or250;source alternatives","contact-tag side A; bridge may be reached","include"),
 sp.read_spec("Read 2","150 or250;source alternatives","contact-tag side B; bridge may be reached","include"),
 sp.read_spec("Index 1 (i7)",None,"acquisition cycles and index member/width unassigned","undetermined"),
 sp.read_spec("Index 2 (i5)",None,"acquisition cycles and index member/width unassigned","undetermined"))
FINAL_CAPTION = ("INFERRED — conditional source-kit-compatible ends with unassigned index widths. A bridge-containing proximity product after Nextera tagmentation, "
                 "streptavidin selection and indexed PCR.")
SEQUENCING_INTRO = ("Paired reads enter from the two Tn5 ends. The source describes 2 × 150 or 2 × 250 base reads; the exact run and index mask remain unassigned. A bridge sequence identifies a linker-bearing fragment; it does not by itself establish a biological contact.")
TAGMENTATION = (55, 5)
MAX_PCR_CYCLES = 13


def contact_scene() -> Scene:
    sc = Scene.duplex(list(CONTACT), label="proximity product")
    sc.junction("top", "left contact tag", "left A:T ligation junction", "ligation")
    sc.junction("top", "right T:A ligation junction", "right contact tag", "ligation")
    sc.mark("top", "internal biotin-dT", "streptavidin handle")
    return sc


def independent_loci():
 sc=Scene()
 sc.strand("A",[role("protein-associated genomic locus A")],label="independent genomic locus A")
 sc.strand("B",[role("protein-associated genomic locus B")],label="independent genomic locus B")
 return sc


def prepared_fragments():
 sc=Scene()
 for name in ("A","B"):
  key="prepared"+name
  sc.strand(key,[role("independent repaired chromatin fragment "+name),seg("source3primeA"+name,"A")],label="dA-prepared fragment "+name)
  sc.overhang(key,"source3primeA"+name,polarity="3-prime",owner="substrate",text="source3′dA")
  sc.note(key,"Representative end strand of source-described duplex; full opposite-end geometry unassigned.")
 return sc


def tagmented_precursor():
 sc=molecule("source Tn5-tagged contact before PCR fill-in",(
  role("transferred Nextera end; strand/gap geometry unassigned"),
  role("bridge-containing contact DNA; variable Tn5 insertion boundaries"),
  role("opposite transferred Nextera end; geometry unassigned")))
 sc.note("molecule","Actual captured substrate is not the final indexed duplex. Source initial72°C3min precedes amplification; native copying/gap geometry remains unassigned.")
 return sc


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    chromatin = graph.input(
        "crosslinked_chromatin", "dual-cross-linked chromatin", MaterialKind.BIOLOGICAL,
        state=independent_loci(), primary=True)
    sonication = graph.input("sonication", "chromatin-sonication condition",
                             MaterialKind.CONDITION)
    fragments = graph.product("chromatin_fragments", "sonicated cross-linked chromatin",
                              MaterialKind.MOLECULAR, state=chromatin.state)
    graph.react("fragment", "Sonicate dual-cross-linked chromatin",
                StepKind.FRAGMENTATION,
                ((chromatin, InputRole.TEMPLATE), (sonication, InputRole.CONDITION)),
                (fragments,))
    antibody = graph.input("target_antibody", "target-specific antibody and capture beads",
                           MaterialKind.SELECTION)
    enriched = graph.product("enriched_chromatin", "protein-associated chromatin fragments",
                             MaterialKind.MOLECULAR, state=fragments.state)
    flowthrough = graph.product("ip_flowthrough", "unbound chromatin fragments",
                                MaterialKind.WASTE)
    graph.react("immunoprecipitate", "Immunoprecipitate protein-associated chromatin",
                StepKind.AFFINITY_CAPTURE,
                ((fragments, InputRole.TEMPLATE), (antibody, InputRole.SELECTION)),
                (enriched, flowthrough))
    graph.discard(flowthrough, "is not retained by the protein-directed immunoprecipitation")
    repair_mix = graph.input("end_repair_mix", "DNA end-repair enzyme mix",
                             MaterialKind.CATALYST)
    repaired = graph.product("end_repaired", "end-repaired chromatin fragments",
                             MaterialKind.MOLECULAR, state=enriched.state)
    graph.react("end_repair", "Repair immunoprecipitated DNA ends", StepKind.END_REPAIR,
                ((enriched, InputRole.TEMPLATE), (repair_mix, InputRole.CATALYST)),
                (repaired,))
    tailing_polymerase = graph.input("tailing_polymerase", "dA-tailing polymerase",
                                     MaterialKind.CATALYST)
    dntp = graph.input("datp", "dATP tailing nucleotide", MaterialKind.MOLECULAR)
    tailed = graph.product("da_tailed", "3′-dA-tailed chromatin fragments",
                           MaterialKind.MOLECULAR, state=prepared_fragments())
    graph.react("da_tail", "Add one 3′ dA to repaired chromatin ends", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (dntp, InputRole.STRUCTURAL),
                 (tailing_polymerase, InputRole.CATALYST)), (tailed,))
    bridge = graph.input("bridge_linker", "published double-ended biotin bridge linker",
                         MaterialKind.ADAPTER, state=BRIDGE.scene())
    ligase = graph.input("bridge_ligase", "bridge-linker DNA ligase",
                         MaterialKind.CATALYST)
    contact = graph.product("bridge_contact", CONTACT.name, MaterialKind.MOLECULAR,
                            state=contact_scene())
    graph.react("bridge_ligation", "Join two dA-tailed fragments through the bridge linker",
                StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (bridge, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (contact,),
                geometry=AdapterLigationGeometry.known(
                    "cohesive", tailed.state, BRIDGE.scene(),
                    statement="The published phosphorylated bridge duplex exposes one 3′ dT at each end for ligation to genomic 3′ dA ends."))
    decrosslinking = graph.input("decrosslinking", "crosslink-reversal and DNA-release conditions", MaterialKind.CONDITION)
    released = graph.product("released", "released bridge-ligated DNA", MaterialKind.MOLECULAR, state=contact.state)
    graph.react("release_dna", "Reverse crosslinks and release bridge-ligated DNA before tagmentation", StepKind.DEPROTEINIZATION,
                ((contact, InputRole.TEMPLATE), (decrosslinking, InputRole.CONDITION)), (released,))
    loaded_adapter = graph.input("loaded_nextera_adapter", "Nextera FC-121-1031 adapter-loaded transposome; physical strands unassigned",
                                 MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "source Nextera Tn5,55°C5min", MaterialKind.CATALYST)
    tagged = graph.product("tagmented_contact", "tagmented bridge-containing contact fragments",
                           MaterialKind.MOLECULAR, state=tagmented_precursor())
    tn5_waste = graph.product("tn5_waste", "candidate non-library tagmentation products",
                              MaterialKind.MOLECULAR)
    graph.react("tagment", "Tagment the bridge-linked contact products",
                StepKind.TAGMENTATION,
                ((released, InputRole.TEMPLATE), (loaded_adapter, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, tn5_waste))
    graph.mark_terminal(Terminal(tn5_waste.id, Disposition.UNRESOLVED, "Candidate insertion products and physical recovery have not been resolved; capture is represented separately."))
    zymo=graph.input("post_tn5_cleanup","Steps51–56:Zymo GenomicDNA Clean&Concentrator column",MaterialKind.SELECTION)
    cleaned_tn5=graph.product("cleaned_tagmented","column-recovered tagmented contact products",MaterialKind.MOLECULAR,state=tagged.state)
    zymo_waste=graph.product("post_tn5_waste","column flowthrough and wash fractions",MaterialKind.WASTE)
    graph.react("purify_tagmented","Purify tagmented DNA before bridge capture",StepKind.CLEANUP,((tagged,InputRole.TEMPLATE),(zymo,InputRole.SELECTION)),(cleaned_tn5,zymo_waste))
    graph.discard(zymo_waste,"Physically separated by the source column transfers; molecular partition remains unassigned.")
    streptavidin = graph.input("streptavidin", "M280 streptavidin beads, iBlock and sheared-genomic-DNA surface blocking; source B&W washes",
                               MaterialKind.SELECTION)
    captured = graph.product("captured_contact", "biotin-bridge-bearing tagmented fragments",
                             MaterialKind.MOLECULAR, state=tagged.state)
    nonbridge = graph.product("nonbridge_fragments", "tagmented fragments lacking the biotin bridge",
                              MaterialKind.WASTE)
    graph.react("capture_bridge", "Capture fragments retaining the internal biotin bridge",
                StepKind.AFFINITY_CAPTURE,
                ((cleaned_tn5, InputRole.TEMPLATE), (streptavidin, InputRole.SELECTION)),
                (captured, nonbridge))
    graph.discard(nonbridge, "does not contain the biotinylated contact bridge")
    ppc=graph.input("ppc_primers","source PPC PCR primer cocktail; native strands unassigned",MaterialKind.OLIGO)
    forward = graph.input("pcr_forward", "Nextera indexed PCR forward primer",
                          MaterialKind.OLIGO)
    reverse = graph.input("pcr_reverse", "Nextera indexed PCR reverse primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "Nextera NPM;initial72°C3min,then minimal cycles,generally≤13",
                             MaterialKind.CATALYST)
    library = graph.product("sequencing_library", FINAL_LIBRARY.name,
                            MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    graph.react("library_pcr", "INFERRED reference product — fill/ index-amplify captured Tn5 contacts on bead",
                StepKind.PCR,
                ((captured, InputRole.TEMPLATE), (ppc,InputRole.STRUCTURAL), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    captured.state,
                    reason="Actual captured Tn5-tagged contacts precede the source72°C first-copy/fill phase. Physical PPC/index species, index widths and full precursor geometry remain unassigned; the conditional final library is not the template."))
    cleanup=graph.input("post_pcr_cleanup","Step69 AmpureXP cleanup;Box3 instructions",MaterialKind.SELECTION)
    cleaned=graph.product("cleaned_library","purified PCR library",MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
    cleanup_waste=graph.product("pcr_cleanup_waste","Ampure supernatant and washes",MaterialKind.WASTE)
    graph.react("cleanup_pcr","Purify PCR products",StepKind.CLEANUP,((library,InputRole.TEMPLATE),(cleanup,InputRole.SELECTION)),(cleaned,cleanup_waste))
    graph.discard(cleanup_waste,"Physically excluded during source bead cleanup.")
    size=graph.input("blue_pippin","Steps70–71:2%agarose BluePippin,300–600bp",MaterialKind.SELECTION)
    selected=graph.product("selected_library","300–600bp long-read ChIA-PET library",MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    offsize=graph.product("off_size","non-selected gel fractions",MaterialKind.WASTE)
    graph.react("size_select","Select300–600bp production library",StepKind.SELECTION,((cleaned,InputRole.TEMPLATE),(size,InputRole.SELECTION)),(selected,offsize))
    graph.discard(offsize,"Not included in the selected BluePippin fraction.")
    graph.sequence(selected, SequencingPlan(
        library=FINAL_LIBRARY, primers=tuple(SEQ_PRIMERS),
        acquired_roles=("Read 1","Read 2"),limitation="Source biological reads2×150 or2×250 are alternatives. Exact acquired index cycles/member/width remain unassigned; canonical priming sites are a conditional reference."))
    graph.declare_reviewed_scope((
        "crosslinked_chromatin", "sonication", "fragment", "target_antibody",
        "immunoprecipitate", "end_repair_mix", "end_repair", "tailing_polymerase",
        "datp", "da_tail", "bridge_linker", "bridge_ligase", "bridge_ligation",
        "decrosslinking", "release_dna", "loaded_nextera_adapter", "tn5", "tagment", "streptavidin",
        "post_tn5_cleanup", "purify_tagmented", "capture_bridge", "ppc_primers", "pcr_forward", "pcr_reverse", "pcr_polymerase",
        "library_pcr", "sequencing_library", "post_pcr_cleanup", "cleanup_pcr", "blue_pippin", "size_select", "selected_library"), SOURCE)
    return graph.require_publishable()

from nascent_protocols import make_spec, role_library

from nascent_protocols import Construct, Definition, InputRole, MaterialKind, NascentEndCapture, ProtocolGraph, StepKind, finish_cdna, ligate_rna_adapter, molecule, role
DEFINITION = Definition("Start-seq — start-site-associated RNA sequencing", "10.1126/science.1181421", "Nechaev et al., Science (2010)", "Short capped nuclear RNAs are sequenced from both ends so their 5′ ends report initiation and their 3′ ends report promoter-proximal pausing.", "Separate compatible libraries read the 5′ and 3′ ends of 25–65-nt capped Start-RNAs.", "01_start_seq.html")

FINAL_LIBRARY = role_library(DEFINITION.title + " final library", (role("assay-derived cDNA"),),
                             indexed=True, platform="Illumina-compatible")
globals().update(vars(make_spec(DEFINITION, FINAL_LIBRARY,
    unavailable='The authoritative paper establishes the library roles and read direction, but does not expose one complete immutable sequencing-primer geometry; adapters are therefore shown by role.')))

def _start_graph(d: Definition, final: Construct) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    nuclear_rna = graph.input(
        "nuclear_rna", "nuclear RNA containing short capped Start-RNAs",
        MaterialKind.BIOLOGICAL,
        state=molecule("nuclear RNA", (role("short capped Start-RNA"),
                                       role("other nuclear RNA"))), primary=True)
    size_rule = graph.input("size_selection", "25–65-nt RNA size-selection rule",
                            MaterialKind.SELECTION)
    short_rna = graph.product("short_rna", "25–65-nt short nuclear RNA",
                              MaterialKind.MOLECULAR, state=nuclear_rna.state)
    long_rna = graph.product("rejected_sizes", "nuclear RNA outside 25–65 nt",
                             MaterialKind.WASTE)
    graph.react("select_size", "Select 25–65-nt short nuclear RNAs",
                StepKind.SELECTION,
                ((nuclear_rna, InputRole.TEMPLATE), (size_rule, InputRole.SELECTION)),
                (short_rna, long_rna))
    graph.discard(long_rna, "falls outside the Start-RNA size interval")
    cap_reagent = graph.input("cap_selection", "capped-RNA retention chemistry",
                              MaterialKind.SELECTION)
    capped = graph.product(
        "capped_start_rna", "capped Start-RNA population", MaterialKind.MOLECULAR,
        state=NascentEndCapture("short capped RNA from isolated nuclei", "both", True).RNA())
    uncapped = graph.product("uncapped_short_rna", "uncapped short nuclear RNA",
                             MaterialKind.WASTE)
    graph.react("select_caps", "Retain capped Start-RNAs", StepKind.SELECTION,
                ((short_rna, InputRole.TEMPLATE), (cap_reagent, InputRole.SELECTION)),
                (capped, uncapped))
    graph.discard(uncapped, "lacks the protected capped end required by Start-seq")
    five_input = graph.product("five_end_input", "Start-RNA assigned to 5-prime library",
                               MaterialKind.MOLECULAR, state=capped.state)
    three_input = graph.product("three_end_input", "Start-RNA assigned to 3-prime library",
                                MaterialKind.MOLECULAR, state=capped.state)
    graph.react("split_end_libraries", "Split Start-RNA into 5-prime and 3-prime libraries",
                StepKind.SPLIT, ((capped, InputRole.TEMPLATE),),
                (five_input, three_input), changes_molecule=False)
    reviewed = [
        "nuclear_rna", "size_selection", "select_size", "cap_selection", "select_caps",
        "split_end_libraries",
    ]
    five3 = ligate_rna_adapter(graph, "five_library_three", five_input,
                                "5-prime-library 3-prime", reviewed)
    five5 = ligate_rna_adapter(graph, "five_library_five", five3,
                                "5-prime-library 5-prime", reviewed)
    finish_cdna(graph, "five_library", five5, final,
                 "5-prime Start-seq library: " + d.endpoint, reviewed)
    three3 = ligate_rna_adapter(graph, "three_library_three", three_input,
                                 "3-prime-library 3-prime", reviewed)
    three5 = ligate_rna_adapter(graph, "three_library_five", three3,
                                 "3-prime-library 5-prime", reviewed)
    finish_cdna(graph, "three_library", three5, final,
                 "3-prime Start-seq library: " + d.endpoint, reviewed)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    key = 'start-seq'
    return _start_graph(d, final)

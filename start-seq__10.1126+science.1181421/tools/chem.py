from nascent_protocols import make_spec, role_library

from nascent_protocols import Construct, Definition, InputRole, MaterialKind, NascentEndCapture, ProtocolGraph, StepKind, finish_cdna, ligate_rna_adapter, molecule, role
DEFINITION = Definition("Start-seq — start-site-associated RNA sequencing", "10.1126/science.1181421", "Nechaev et al., Science (2010)", "Short capped nuclear RNAs are sequenced from both ends so their 5′ ends report initiation and their 3′ ends report promoter-proximal pausing.", "Separate compatible libraries read the 5′ and 3′ ends of 25–100-nt size-selected capped Start-RNAs.", "01_start_seq.html")

FINAL_LIBRARY = role_library(DEFINITION.title + " final library", (role("assay-derived cDNA"),),
                             indexed=False, platform="Illumina-compatible")
globals().update(vars(make_spec(DEFINITION, FINAL_LIBRARY,
    unavailable='The authoritative paper establishes the library roles and read direction, but does not expose one complete immutable sequencing-primer geometry; adapters are therefore shown by role.')))
CAVEAT = ('The original <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM463298">5′-read</a> '
          'and <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM463301">3′-read</a> '
          'GEO methods establish 25–100-nt gel selection and distinct standard versus custom '
          'small-RNA linkers. The original Science supplement with the custom oligo details is '
          'currently inaccessible; exact bases and junctions in the 3′-read branch remain unresolved.')

def _start_graph(d: Definition, final: Construct) -> ProtocolGraph:
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    nuclear_rna = graph.input(
        "nuclear_rna", "nuclear RNA containing short capped Start-RNAs",
        MaterialKind.BIOLOGICAL,
        state=molecule("nuclear RNA", (role("short capped Start-RNA"),
                                       role("other nuclear RNA"))), primary=True)
    size_rule = graph.input("size_selection", "gel size selection: 25–100 nt",
                            MaterialKind.SELECTION)
    short_rna = graph.product("short_rna", "25–100-nt nuclear RNA",
                              MaterialKind.MOLECULAR, state=nuclear_rna.state)
    long_rna = graph.product("rejected_sizes", "longer nuclear RNA",
                             MaterialKind.WASTE)
    graph.react("select_size", "Gel-select 25–100-nt nuclear RNAs",
                StepKind.SELECTION,
                ((nuclear_rna, InputRole.TEMPLATE), (size_rule, InputRole.SELECTION)),
                (short_rna, long_rna))
    graph.discard(long_rna, "falls outside the short-RNA source selection")
    cap_reagent = graph.input("cap_selection", "5′ polyphosphatase then Terminator exonuclease",
                              MaterialKind.SELECTION)
    capped = graph.product(
        "capped_start_rna", "capped Start-RNA population", MaterialKind.MOLECULAR,
        state=NascentEndCapture("short capped RNA from isolated nuclei", "both", True).RNA())
    uncapped = graph.product("uncapped_short_rna", "uncapped short nuclear RNA",
                             MaterialKind.WASTE)
    graph.react("select_caps", "Remove uncapped RNAs after 5′ polyphosphatase and Terminator treatment", StepKind.SELECTION,
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
    alkaline_phosphatase = graph.input(
        "alkaline_phosphatase", "alkaline phosphatase before 5′ adapter ligation",
        MaterialKind.CATALYST)
    tobacco_pyrophosphatase = graph.input(
        "tobacco_acid_pyrophosphatase", "tobacco acid pyrophosphatase for cap removal",
        MaterialKind.CATALYST)
    reviewed.extend(("alkaline_phosphatase", "tobacco_acid_pyrophosphatase"))

    def prepare_five_prime(prefix, three_prime_ligated):
        processed = graph.product(
            prefix + "_decapped", "3′-adapter-ligated RNA ready for 5′ linker",
            MaterialKind.MOLECULAR, state=three_prime_ligated.state)
        graph.react(prefix + "_prepare_five_prime",
                    "Treat with alkaline phosphatase then tobacco acid pyrophosphatase",
                    StepKind.CHEMICAL_CONVERSION,
                    ((three_prime_ligated, InputRole.TEMPLATE),
                     (alkaline_phosphatase, InputRole.CATALYST),
                     (tobacco_pyrophosphatase, InputRole.CATALYST)), (processed,))
        reviewed.append(prefix + "_prepare_five_prime")
        return processed

    five3 = ligate_rna_adapter(graph, "five_library_three", five_input,
                                "5′-read Illumina small-RNA 3′", reviewed)
    five5 = ligate_rna_adapter(graph, "five_library_five", prepare_five_prime("five_library", five3),
                                "5′-read Illumina small-RNA 5′", reviewed)
    finish_cdna(graph, "five_library", five5, final,
                 "5-prime Start-seq library: " + d.endpoint, reviewed,
                 rt_label="Illumina small-RNA kit RT primer", indexed=False)
    three3 = ligate_rna_adapter(graph, "three_library_three", three_input,
                                 "3′-read custom 3′", reviewed)
    three5 = ligate_rna_adapter(graph, "three_library_five", prepare_five_prime("three_library", three3),
                                 "3′-read custom 5′", reviewed)
    finish_cdna(graph, "three_library", three5, final,
                 "3-prime Start-seq library: " + d.endpoint, reviewed,
                 rt_label="custom 3′-read RT primer", indexed=False)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    key = 'start-seq'
    return _start_graph(d, final)

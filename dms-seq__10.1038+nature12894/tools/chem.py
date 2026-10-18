"""Protocol-owned source, oligos, molecular workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, Workflow
from circular import CircularizationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role

DMS_LINKER = Segment("microRNA cloning linker-1", "CTGTAGGCACCATCAAT", "r2")
DMS_BLOCK = Segment("terminal dideoxyC", "C", note="/3ddC/: terminal nucleotide lacks 3′ OH")


TITLE = 'DMS-seq — in-vivo RNA structure probing'
NOTES = '01_dms_seq.html'
SOURCE = 'Rouskin et al., Nature (2014); <a href="https://doi.org/10.1038/nature12894">doi:10.1038/nature12894</a>.'
SUMMARY = 'DMS marks solvent-accessible A and C in living cells; reverse-transcriptase termination sites report nucleotide accessibility after adapter ligation and cDNA circularization.'
CAVEAT = ("The graph follows the source yeast in-vivo 30 °C branch: DMS probing precedes "
          "quenching, RNA extraction and poly(A) selection. Refolded/denatured controls and "
          "mammalian preparations are separate source conditions. Exact RT/PCR primer members, "
          "circle junction and run layout remain unassigned; the cited modified predecessor "
          "recipe is not assumed to establish all DMS-seq oligos.")
FINAL_LIBRARY = Construct((role('P5'), role('DMS-stop cDNA'), role('indexing status unassigned'), role('P7')), name='DMS-seq library')
FINAL_CAPTION = 'Role-level PCR product; linker-1 is source-published; selected PCR and sequencing-primer geometry remains unverified.'
SEQ_PRIMERS = ()
READ_LENGTHS = None
SEQUENCING_INTRO = ''
SEQUENCING_ENDING = 'The defining experiment sequenced RT-termination libraries on an Illumina HiSeq 2000.'
SEQUENCING_UNAVAILABLE = 'PCR and sequencing-primer geometry is not assigned in this role-level model.'

def oligos():
    return ('microRNA cloning linker-1 — 5′-/5rApp/CTGTAGGCACCATCAAT/3ddC/-3′',)

def workflow() -> Workflow:
    return chain(molecule("cells for in-vivo probing", (role("native cellular RNA structure"),)), (
        ("Probe accessible A and C with DMS",
         molecule("DMS-modified RNA", (role("RNA"), role("methylated exposed A/C"))),
         "Parallel in-vivo, refolded in-vitro and denatured samples separate cellular structure from intrinsic reactivity."),
        ("Quench and extract poly(A)-selected DMS-modified RNA",
         molecule("DMS-modified mRNA", (role("poly(A)-selected modified RNA"),)),
         "In-vivo probing occurs before RNA extraction; in-vitro controls are separate preparations."),
        ("Fragment with Zn²⁺, select 60–80 nt, and repair 3′ phosphates",
         molecule("end-repaired DMS RNA fragment", (role("60–80-nt modified RNA"), role("3-prime hydroxyl"))),
         "T4 PNK resolves the 3′ phosphates left by chemical fragmentation."),
        ("Ligate preadenylated cloning linker-1",
         molecule("linker-ligated DMS RNA", (role("modified RNA fragment"), DMS_LINKER, DMS_BLOCK),
                  junctions=(("modified RNA fragment", "microRNA cloning linker-1", "Rnl2tr K227Q ligation"),)),
         "The linker is 5′ preadenylated and 3′ ddC-blocked, so ligation proceeds without ATP."),
        ("Reverse-transcribe and isolate prematurely terminated cDNA",
         molecule("DMS-stop cDNA", (role("adapter-primed cDNA ending before modified A/C"),)),
         "A gel band 25–50 nt above RT-primer size enriches SuperScript III termination products."),
        ("Circularize cDNA and add Illumina adapters by PCR",
         molecule("DMS-seq library", tuple(FINAL_LIBRARY)),
         "CircLigase precedes eight to ten PCR cycles and HiSeq 2000 sequencing."),
    ))



DMS_METHODS = "Rouskin2014 PMC3966492 Methods DMS Modification and Library Generation, lines50–58"


def _gel(graph, identifier, template, label, retained_state, excluded_label):
    gel = graph.input(identifier + "_gel", label, MaterialKind.SELECTION, source=DMS_METHODS)
    kept = graph.product(identifier + "_product", label + " retained fraction", MaterialKind.MOLECULAR,
                         state=retained_state)
    excluded = graph.product(identifier + "_excluded", excluded_label, MaterialKind.WASTE)
    graph.react(identifier, label, StepKind.SELECTION,
                ((template, InputRole.TEMPLATE), (gel, InputRole.SELECTION)), (kept, excluded),
                changes_molecule=False, source=DMS_METHODS)
    graph.discard(excluded, "physically excluded when the selected gel slice is recovered")
    return kept


def linker_reagent():
    sc = Scene()
    sc.strand("adapter", [DMS_LINKER, DMS_BLOCK], label="microRNA cloning linker-1",
              mod5="rApp", mod3="ddC (blocked)")
    return sc


def dms_ligation_geometry():
    sc = linker_reagent()
    sc.strand("RNA", [role("PNK-treated RNA fragment")], label="RNA acceptor", mod3="OH")
    sc.ligation_end("adapter", DMS_LINKER.name, kind="single-strand", polarity="5-prime",
                    owner="adapter", text="5′ preadenylated donor; RNA 3′ OH acceptor")
    return AdapterLigationGeometry.known("single-strand", sc,
        statement="Rouskin2014 Library Generation: 3′ dephosphorylation then ligation to /5rApp/ linker-1 "
                  "with Rnl2tr K227Q; linker terminates in /3ddC/.")


def protocol_graph() -> ProtocolGraph:
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    cells = graph.input("input", "exponentially growing yeast cells (in-vivo 30 °C branch)",
                         MaterialKind.BIOLOGICAL, primary=True,
                         state=molecule("cells before probing", (role("cellular RNA structures"),)))
    dms = graph.input("probe_condition", "300–600 µL DMS in 15 mL yeast culture; 2–4 min at 30 °C",
                       MaterialKind.CONDITION, source=DMS_METHODS)
    probed = graph.product("probe_product", "DMS-probed cells", MaterialKind.BIOLOGICAL,
                           state=molecule("probed cells", (role("RNA with DMS-modified exposed A/C"),)))
    graph.react("probe", "Probe exposed A/C inside living cells", StepKind.CHEMICAL_CONVERSION,
                ((cells, InputRole.TEMPLATE), (dms, InputRole.CONDITION)), (probed,))
    recovery = graph.input("rna_recovery", "BME/isoamyl-alcohol quench, cold washes, hot acid phenol and poly(A)+ Dynal beads",
                            MaterialKind.CONDITION, source=DMS_METHODS)
    mrna = graph.product("mrna", "poly(A)-selected DMS-modified mRNA", MaterialKind.MOLECULAR,
                         state=molecule("modified mRNA", (role("DMS-modified poly(A) RNA"),)))
    removed = graph.product("rna_recovery_waste", "cell debris and fractions excluded during extraction/poly(A) selection",
                            MaterialKind.WASTE)
    graph.react("recover_rna", "Quench, extract RNA and select poly(A) mRNA", StepKind.EXTRACTION,
                ((probed, InputRole.TEMPLATE), (recovery, InputRole.CONDITION)), (mrna, removed))
    graph.discard(removed, "removed by source extraction and poly(A)-selection procedures")
    zinc = graph.input("fragment_condition", "95 °C denature 2 min; Ambion Zn²⁺ fragmentation buffer, 95 °C 2 min; stop and ice",
                        MaterialKind.CONDITION, source=DMS_METHODS)
    fragments = graph.product("fragment_product", "randomly fragmented modified RNA", MaterialKind.MOLECULAR,
                              state=molecule("RNA fragments", (role("modified RNA fragment"), role("3′ phosphate end"))))
    graph.react("fragment", "Fragment the selected DMS-modified mRNA", StepKind.FRAGMENTATION,
                ((mrna, InputRole.TEMPLATE), (zinc, InputRole.CONDITION)), (fragments,))
    selected = _gel(graph, "rna_size", fragments, "10% TBU PAGE: recover 60–80-nt RNA",
                    fragments.state, "RNA outside the 60–80-nt gel fraction")
    pnk = graph.input("repair_enzyme", "T4 PNK in 1× PNK buffer, 37 °C 1 h, no ATP addition specified",
                       MaterialKind.CATALYST, source=DMS_METHODS)
    repaired = graph.product("repair_product", "3′-dephosphorylated modified RNA", MaterialKind.MOLECULAR,
                             state=molecule("RNA acceptor", (role("modified RNA fragment"), role("3′ OH"))))
    graph.react("repair", "Remove fragmentation-derived RNA 3′ phosphates", StepKind.DEPHOSPHORYLATION,
                ((selected, InputRole.TEMPLATE), (pnk, InputRole.CATALYST)), (repaired,))
    adapter = graph.input("ligate_adapter", "1 µg /5rApp/ linker-1 /3ddC/", MaterialKind.ADAPTER,
                           state=linker_reagent(), source=DMS_METHODS)
    ligase = graph.input("ligate_ligase", "T4 RNA ligase 2 truncated K227Q; DTT, PEG, ligase2 buffer; RT 1.5 h",
                          MaterialKind.CATALYST, source=DMS_METHODS)
    tagged_state = molecule("linker-tagged RNA", (role("modified RNA fragment"), DMS_LINKER, DMS_BLOCK),
                            junctions=(("modified RNA fragment", DMS_LINKER.name, "RNA-linker ligation"),))
    tagged = graph.product("ligate_product", "linker-ligated DMS RNA", MaterialKind.MOLECULAR, state=tagged_state)
    graph.react("ligate", "Ligate preadenylated linker to the RNA 3′ OH", StepKind.LIGATION,
                ((repaired, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)),
                (tagged,), geometry=dms_ligation_geometry())
    purified = _gel(graph, "ligation_gel", tagged, "10% TBU PAGE: recover linker-ligated RNA",
                     tagged_state, "unligated excess linker and RNA outside the recovered band")
    primer = graph.input("rt_primer", "DMS-seq RT primer (selected member/modified predecessor applicability unassigned)",
                          MaterialKind.OLIGO, source=DMS_METHODS)
    enzyme = graph.input("rt_enzyme", "SuperScript III, 20 µL reaction at 52 °C", MaterialKind.CATALYST,
                          source=DMS_METHODS)
    cdna_state = molecule("stop-defined single-strand cDNA", (role("RT-primer-derived end (unassigned)"),
                          role("antisense cDNA terminating at a DMS-dependent stop")))
    cdna = graph.product("rt_product", "DMS-stop cDNA population", MaterialKind.MOLECULAR, state=cdna_state)
    graph.react("rt", "Reverse-transcribe with SuperScript III at 52 °C", StepKind.REVERSE_TRANSCRIPTION,
                ((purified, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL), (enzyme, InputRole.CATALYST)), (cdna,))
    short = _gel(graph, "cdna_size", cdna, "Gel select truncated cDNA 25–50 nt above RT-primer size",
                  cdna_state, "cDNA outside the selected stop-product band")
    circ = graph.input("circle_enzyme", "CircLigase (Epicentre)", MaterialKind.CATALYST, source=DMS_METHODS)
    circle = graph.product("circle_product", "circularized DMS-stop cDNA (junction unassembled)",
                            MaterialKind.MOLECULAR,
                            state=molecule("circular cDNA boundary", (role("circularized stop cDNA; exact closing ends unassigned"),)))
    graph.react("circle", "Circularize the selected stop-defined cDNA", StepKind.CIRCULARIZATION,
                ((short, InputRole.TEMPLATE), (circ, InputRole.CATALYST)), (circle,),
                geometry=CircularizationGeometry.unavailable("DMS-seq uses a modified ribosome-profiling protocol; "
                    "the selected RT-primer/circle geometry has not been assigned to this study."))
    forward = graph.input("library_forward", "DMS-seq library forward PCR primer (member unassigned)", MaterialKind.OLIGO)
    reverse = graph.input("library_reverse", "DMS-seq library reverse PCR primer (indexing unassigned)", MaterialKind.OLIGO)
    polymerase = graph.input("library_polymerase", "DMS-seq library PCR mix (identity unassigned)", MaterialKind.CATALYST)
    library = graph.product("library_library", "DMS-seq Illumina sequencing library", MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    graph.react("library_pcr", "Introduce Illumina sequencing adapters by 8–10 PCR cycles", StepKind.PCR,
                ((circle, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL), (reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (library,), geometry=PCRGeometry.unavailable(circle.state,
                    reason="The study-specific RT and PCR primer members/circle are unassigned; predecessor identity "
                           "alone does not establish the actual DMS-seq PCR landing geometry."))
    graph.sequence(library, SequencingPlan(acquisition="unspecified", limitation=SEQUENCING_ENDING +
                    " Read length, read type and index acquisition remain unassigned in this declaration."))
    graph.declare_reviewed_scope(("input", "probe_condition", "probe", "rna_recovery", "recover_rna", "fragment_condition", "fragment",
        "rna_size_gel", "rna_size", "repair_enzyme", "repair", "ligate_adapter", "ligate_ligase", "ligate", "ligation_gel_gel", "ligation_gel",
        "rt_primer", "rt_enzyme", "rt", "cdna_size_gel", "cdna_size", "circle_enzyme", "circle", "library_forward", "library_reverse",
        "library_polymerase", "library_pcr", "library_library"), DMS_METHODS)
    return graph.require_publishable()

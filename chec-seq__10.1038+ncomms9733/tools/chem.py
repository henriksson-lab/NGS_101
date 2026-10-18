"""Protocol-owned source, oligos, workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, feature, panel
import illumina as il
import seqprimers as sp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind)
from protocol_model import molecule, role
from endprep import EndPreparedDuplex, dA_tailed_scene
from protocol_graph import Disposition, Terminal

CHEC_UNIVERSAL_ADAPTER = il.P5[:-4] + il.TRUSEQ_READ1

CHEC_INDEXED_ADAPTER = il.CANONICAL_ARM_READ2 + "N" * 6 + il.P7_RC

CHEC_P5_PCR = il.P5[:21]

CHEC_P7_PCR = il.P7[:22]

SAMPLE_INDEX = feature("orsi_illustrative_i7", "sample_index", "unknown")


def indexed_adapter_segments():
    return [Segment("Orsi Read 2 adapter arm", CHEC_INDEXED_ADAPTER[:33], "r2"),
            Segment("six-position sample index; oligo orientation", "N" * 6, "i7",
                    placeholder=True, feature=SAMPLE_INDEX),
            Segment("Orsi P7 complement", CHEC_INDEXED_ADAPTER[39:], "p7")]


def orsi_adapter_scene():
    # The universal terminal T is a one-base 3′ overhang, not part of the paired core.
    scene = Scene()
    scene.strand("universal", [Segment("Orsi P5/Read1 Y arm", CHEC_UNIVERSAL_ADAPTER[:-13]),
                                Segment("universal paired core", CHEC_UNIVERSAL_ADAPTER[-13:-1]),
                                Segment("adapter terminal dT", CHEC_UNIVERSAL_ADAPTER[-1])],
                 label="Orsi Universal Adapter", mod3="terminal phosphorothioate")
    bottom = [Segment("indexed paired core", CHEC_INDEXED_ADAPTER[:12]),
              Segment("indexed Read2 Y arm", CHEC_INDEXED_ADAPTER[12:33]),
              *indexed_adapter_segments()[1:]]
    scene.anneal("indexed", bottom, to="universal",
                 pair=("indexed paired core", "universal paired core"),
                 unpaired=("indexed Read2 Y arm", "six-position sample index; oligo orientation", "Orsi P7 complement"),
                 label="Orsi Indexed Adapter", mod5="/5Phos/", mod3="terminal phosphorothioate")
    scene.overhang("universal", "adapter terminal dT", polarity="3-prime", owner="adapter",
                   text="unpaired source terminal dT")
    return scene


ORSI_DA_LIGATION = AdapterLigationGeometry.known("cohesive", orsi_adapter_scene(),
    dA_tailed_scene(EndPreparedDuplex(Construct([role("reference ChEC insert")]), five_prime_phosphate=True, three_prime_overhang="A")),
    statement="Separate Orsi Basic Protocol3 supplies both protected Y-adapter oligos and the dA-tailed insert workflow. The universal adapter's terminal T is unpaired and complements the fragment dA; these exact oligos are an illustrated variant, not a demonstrated 2015 sample-index assignment.")


def chec_library():
    library = Construct([Segment("complete Orsi Universal Adapter; P5/Read1 overlap", CHEC_UNIVERSAL_ADAPTER, "r1"),
                         role("DNA released by two nearby protein-tethered MNase cuts"),
                         Segment("fragment-derived dA junction", "A"), *indexed_adapter_segments()],
                        name="ChEC readout illustrated with the separately documented Orsi adapter variant")
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
    errors = sp.verify(library, primers, required_roles=tuple(p.role for p in primers))
    if errors:
        raise ValueError("Orsi-derived ChEC readout: " + "; ".join(errors))
    return library, primers


def oligos():
    return (
        f"TruSeq Universal Adapter — 5′-{CHEC_UNIVERSAL_ADAPTER[:-1]}*{CHEC_UNIVERSAL_ADAPTER[-1]}-3′",
        f"TruSeq Indexed Adapter — 5′-/5Phos/{CHEC_INDEXED_ADAPTER[:-1]}*{CHEC_INDEXED_ADAPTER[-1]}-3′",
        f"P5 library-enrichment primer — 5′-{CHEC_P5_PCR}-3′",
        f"P7 library-enrichment primer — 5′-{CHEC_P7_PCR}-3′",
    )


ORSI_REFERENCE, ORSI_REFERENCE_PRIMERS = chec_library()
FINAL_LIBRARY = None
SEQ_PRIMERS = ()
TOPOLOGY_PANELS = (("Separate later Orsi adapter reference", panel(orsi_adapter_scene().rows(), caption="2015 Basic Protocol3 reagent reference; not a study-specific adapter assignment.")),)
TITLE = 'ChEC-seq — calcium-triggered chromatin endogenous cleavage sequencing'
NOTES = '01_chec_seq.html'
SOURCE = '<a href="https://doi.org/10.1038/ncomms9733">Zentner et al., Nature Communications (2015)</a>; the exact displayed library oligos come from the Henikoff laboratory’s separate <a href="https://doi.org/10.1002/0471142727.mb2131s110">Orsi et al. stepwise protocol (2015)</a>.'
SUMMARY = 'A protein of interest is genetically fused to MNase. Digitonin admits calcium into living yeast, activating cleavage beside occupied binding sites; short two-cut fragments are selected and converted into a paired-end Illumina library.'
CAVEAT = 'ChEC-seq reports calcium-triggered cleavage by a mobile protein–MNase fusion, not an antibody-captured footprint. A free nuclear MNase control and early kinetic samples are essential because slower cleavage contains accessibility, enzyme-bias and transient-collision signal. The defining Methods specify2.5:1 negative AMPure selection: keep the small-DNA supernatant, then re-extract and precipitate it to remove RNase before library preparation; this second recovery is summarized in the graph. The exact oligo strings shown here come from a separate same-laboratory Orsi procedure; the defining ChEC article does not itself print a sample-index assignment.'
FINAL_CAPTION = 'The selected insert is bounded by two MNase-created ends near a protein-binding event. End repair and dA-tailing erase the original overhang geometry but retain both genomic boundary coordinates; the illustrated Orsi variant has one six-base i7 index and no i5 index. This does not identify the exact adapters/index member in the defining 2015 ChEC experiment.'
READ_LENGTHS = {'Read 1': 25, 'Read 2': 25}
SEQUENCING_INTRO = 'The defining data set used 25-base paired-end reads. Standard TruSeq Read 1 and Read 2 sites are computed here on the separately illustrated Orsi-derived construct. Its six-position index is shown with six Index1 cycles; this is not a recovered 2015 instrument mask or sample-index assignment.'
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''

def dna_population(label, *, rna=False, protein=False, carrier=False):
    sc=Scene();sc.strand("dna",[role(label)],label="native DNA population; strand pairing/ends unassigned")
    if rna:sc.strand("rna",[role("independent native RNA/hydrolysates; stage-dependent")],label="separate RNA-derived population")
    if protein:sc.strand("protein",[role("associated protein/proteolysis products; stage-dependent")],label="separate non-DNA material")
    if carrier:sc.strand("carrier",[role("GlycoBlue carrier; component-specific recovery unassigned")],label="independent recovery carrier")
    return sc

def native_library_state(label, *, carrier=False):
    sc=Scene();parts=[role("native library end A; strand/boundary unassigned",inferred=True),role(label,inferred=True),
         role("native indexed end(s); scheme/width unassigned",inferred=True,molecular_feature=feature("native_sample_index","sample_index","unknown"))]
    sc.strand("native",parts,label="INFERRED boundaries — native library molecule")
    sc.strand("opposite",[role("native opposite strand class; copied/paired extent unassigned",inferred=True)],label="native opposite DNA strand")
    if carrier:sc.strand("carrier",[role("possible co-carried GlycoBlue; recovery unassigned")],label="independent carrier")
    if carrier:sc.strand("rna_residual",[role("possible co-carried RNAhydrolysates; class-specific recovery unassigned")],label="independent RNA-derived residuals")
    return sc

def protocol_graph() -> ProtocolGraph:
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    source = ("Zentner et al., Nature Communications (2015); doi:10.1038/ncomms9733; "
              "separate illustrated library variant: Orsi et al., Current Protocols (2015)")
    graph = ProtocolGraph("ChEC-seq", authority=source)
    cells = graph.input(
        "fusion_cells", "living yeast expressing the chromatin-protein–MNase fusion",
        MaterialKind.BIOLOGICAL,
        state=dna_population("native chromatin DNA",rna=True,protein=True),
        primary=True)
    digitonin = graph.input("digitonin", "0.1% digitonin permeabilization condition",
                            MaterialKind.CONDITION)
    chromatin = graph.product(
        "permeabilized_chromatin", "accessible intracellular chromatin",
        MaterialKind.MOLECULAR,
        state=dna_population("accessible intracellular DNA",rna=True,protein=True))
    tethered_mnase = graph.product("tethered_mnase", "inactive protein-tethered MNase",
                                   MaterialKind.CATALYST)
    graph.react("permeabilize", "Permeabilize cells with digitonin",
                StepKind.PERMEABILIZATION,
                ((cells, InputRole.TEMPLATE), (digitonin, InputRole.CONDITION)),
                (chromatin, tethered_mnase))
    calcium = graph.input("calcium", "2 mM calcium activation condition",
                          MaterialKind.CONDITION)
    cut = graph.product(
        "cut_chromatin", "chromatin cut near occupied protein sites",
        MaterialKind.MOLECULAR,
        state=dna_population("MNase-cut DNA fragment population",rna=True,protein=True))
    graph.react("timed_digest", "Activate tethered MNase and remove timed aliquots",
                StepKind.DIGESTION,
                ((chromatin, InputRole.TEMPLATE), (tethered_mnase, InputRole.CATALYST),
                 (calcium, InputRole.CONDITION)), (cut,))
    quench = graph.input("quench_reagents", "EDTA/EGTA, SDS and high-salt quench",
                         MaterialKind.CONDITION)
    stopped = graph.product(
        "stopped_digest", "stopped ChEC digest", MaterialKind.MOLECULAR,
        state=dna_population("quenched MNase-cut DNA",rna=True,protein=True))
    graph.react("quench", "Quench each timed digestion aliquot", StepKind.QUENCHING,
                ((cut, InputRole.TEMPLATE), (quench, InputRole.CONDITION)), (stopped,))
    proteinase = graph.input("proteinase_k", "proteinase K, 55°C for30min", MaterialKind.CATALYST)
    deproteinized = graph.product("deproteinized", "deproteinized ChEC digest", MaterialKind.MOLECULAR, state=dna_population("DNA after proteolysis; digestionproducts remain until extraction",rna=True,protein=True))
    graph.react("deproteinize", "Digest proteins before nucleic-acid extraction", StepKind.DEPROTEINIZATION,
                ((stopped, InputRole.TEMPLATE), (proteinase, InputRole.CATALYST)), (deproteinized,))
    extraction=graph.input("extraction_condition","phenol/chloroform/isoamyl extraction",MaterialKind.CONDITION)
    extracted=graph.product("extracted_nucleic_acids","aqueous extracted ChEC nucleic acid",MaterialKind.MOLECULAR,state=dna_population("native cut DNA after organic extraction",rna=True))
    waste=graph.product("extraction_waste","removed organic/interface fraction",MaterialKind.WASTE)
    graph.react("extract","Extract nucleic acids",StepKind.EXTRACTION,((deproteinized,InputRole.TEMPLATE),(extraction,InputRole.CONDITION)),(extracted,waste))
    graph.discard(waste,"Physical phenol/chloroform extraction fraction; proteolysis itself is not discard.")
    etoh=graph.input("ethanol","ethanol precipitation and source ethanolwash",MaterialKind.SELECTION)
    carrier=graph.input("glycoblue","GlycoBlue recovery carrier",MaterialKind.MOLECULAR)
    precip=graph.product("precipitated","recovered nucleic acids and carrier",MaterialKind.MOLECULAR,state=dna_population("native recovered DNA",rna=True,carrier=True))
    sup=graph.product("precipitation_waste","removed precipitation supernatant/wash fraction",MaterialKind.WASTE)
    graph.react("precipitate","Precipitate and wash extracted nucleic acids",StepKind.CLEANUP,((extracted,InputRole.TEMPLATE),(carrier,InputRole.SPIKE_IN),(etoh,InputRole.SELECTION)),(precip,sup))
    graph.discard(sup,"Physical source ethanol recovery fraction.")
    rnase=graph.input("rnase","RNaseA;37°C20min",MaterialKind.CATALYST)
    extracted=graph.product("extracted_dna","DNA/RNAhydrolysates afterRNase",MaterialKind.MOLECULAR,state=dna_population("native ChEC DNA",rna=True,carrier=True))
    graph.react("digest_rna","Digest RNA after nucleic-acid recovery",StepKind.DIGESTION,((precip,InputRole.TEMPLATE),(rnase,InputRole.CATALYST)),(extracted,))
    beads = graph.input("ampure_selection", "2.5:1 AMPure size-selection rule",
                        MaterialKind.SELECTION)
    selected = graph.product(
        "selected_fragments", "short ChEC DNA fragments", MaterialKind.MOLECULAR,
        state=dna_population("selected short MNase fragments",rna=True,protein=True,carrier=True))
    large = graph.product("large_genomic_dna", "large bead-bound genomic DNA",
                          MaterialKind.WASTE)
    graph.react("size_select", "Retain the small-DNA supernatant", StepKind.SELECTION,
                ((extracted, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (selected, large))
    graph.discard(large, "binds the AMPure beads while the short ChEC fragments remain in supernatant")
    reviewed = [
        "fusion_cells", "digitonin", "permeabilize", "calcium", "timed_digest",
        "proteinase_k", "deproteinize", "quench_reagents", "quench", "extraction_condition", "extract","ethanol","glycoblue","precipitate","rnase","digest_rna",
        "ampure_selection", "size_select",
    ]
    # Defining article explicitly re-extracts and precipitates the negative-selection supernatant.
    recovery = graph.input("postselect_recovery", "phenol/chloroform extraction then ethanol/GlycoBlue precipitation", MaterialKind.CONDITION)
    recovered = graph.product("recovered_chec", "recovered ChEC DNA after RNase-removal extraction", MaterialKind.MOLECULAR,
                              state=dna_population("recovered selected fragments", rna=True, carrier=True))
    removed = graph.product("postselect_removed", "physically removed extraction/precipitation fraction", MaterialKind.WASTE)
    graph.react("recover_selected", "Re-extract small-DNA supernatant and precipitate before library preparation", StepKind.CLEANUP,
                ((selected, InputRole.TEMPLATE), (recovery, InputRole.CONDITION)), (recovered, removed))
    graph.discard(removed, "Physical source recovery removes RNase/reagents; complete species-specific depletion is not asserted.")
    preparation = graph.input("native_method", "native library preparation from cited2011/2013 methods; exact implementation unassigned", MaterialKind.CONDITION)
    prepared = graph.product("native_prepared", "actual native library-PCR precursor; ends unassigned", MaterialKind.MOLECULAR,
                             state=native_library_state("native pre-PCR molecules", carrier=True))
    graph.react("native_prepare", "Native cited-method library preparation boundary", StepKind.MAGIC,
                ((recovered,InputRole.TEMPLATE),(preparation,InputRole.CONDITION)), (prepared,),
                reason="Defining Methods cites Henikoff2011 and Zentner2013 and specifies KAPA amplification. Those detailed predecessor implementations are outside this capped pass; the later Orsi reference does not establish native ends, chemistry or primer identity.")
    forward=graph.input("native_forward", "unassigned native forward library-primer role", MaterialKind.OLIGO, inferred=True)
    reverse=graph.input("native_reverse", "unassigned native opposite library-primer role", MaterialKind.OLIGO, inferred=True)
    polymerase=graph.input("native_kapa", "KAPA polymerase named by defining ChEC Methods", MaterialKind.CATALYST)
    library=graph.product("native_library", "INFERRED — native ChEC amplified library", MaterialKind.LIBRARY,
                          state=native_library_state("native amplified DNA; end/index scheme unassigned", carrier=True), inferred=True)
    residual=graph.product("native_residual", "original pre-PCR molecules/carrier after copying", MaterialKind.MOLECULAR, state=prepared.state)
    graph.react("native_pcr", "INFERRED boundaries — KAPA library amplification", StepKind.PCR,
                ((prepared,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),(reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (library,residual), geometry=PCRGeometry.unavailable(prepared.state,
                    reason="Actual native prepared molecules are the input. Opposing primer roles are unassigned, not two verified source oligos; first-cycle chemistry and ends remain unresolved."))
    graph.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,"Copying does not erase original templates or carrier; source-specific final recovery is unassigned."))
    graph.sequence(library,SequencingPlan(acquired_roles=("Read 1","Read 2"),read_lengths=READ_LENGTHS,limitation=SEQUENCING_ENDING))
    reviewed.extend(("postselect_recovery","recover_selected","native_method","native_prepare","native_forward","native_reverse","native_kapa","native_pcr","native_library"))
    graph.declare_reviewed_scope(reviewed, source)
    graph.annotate(SourceDiscrepancy(
        "fusion_cells", "Published correction: original Reb1 result table",
        (SourceClaim("Original Supplementary Data 1, Reb1 sheet",
                     "The initially published workbook sorted columns E, F, G, H and K independently of the other columns."),
         SourceClaim("2015 corrigendum, doi:10.1038/ncomms10264, sole correction paragraph",
                     "The authors replaced that workbook and state that the incorrect version was not used in the paper's analyses.")),
        "This page describes molecular preparation and uses no peak values from the original result table.",
        "The defining article links the original corrigendum, which identifies the affected sheet and columns.",
        "A historical downloaded result workbook needs its own version check; the correction does not identify an adapter or reaction change."))
    graph.annotate(SourceDiscrepancy(
        "timed_digest", "Published correction: DNA-shape interpretation",
        (SourceClaim("Defining article Fig. 7 and DNA-shape Methods",
                     "High- and low-scoring sites were motif-centered, whereas random controls were not."),
         SourceClaim("2017 corrigendum, doi:10.1038/ncomms15723, correction paragraphs",
                     "With motif-centered random controls, the reported analysis cannot establish shape-dependent, sequence-independent binding at low-scoring sites.")),
        "The schematic retains calcium-triggered MNase cleavage and makes no sequence-independent DNA-shape binding claim.",
        "The authors' correction withdraws that interpretation and states that it does not invalidate ChEC-seq as a method to identify binding sites.",
        "Cleavage kinetics alone do not establish the molecular basis of a particular binding event; the linked correspondence and reply are separate analysis sources."))
    return graph.require_publishable()

# The separate 2015 ordering reference above is not the defining experiment's endpoint.
SHOW_READ_CONTENT=False
READ_SPECS=(sp.read_spec("Read 1",25,"Defining source paired-end genomic boundary read","include"),
            sp.read_spec("Read 2",25,"Defining source paired-end opposite boundary read","include"),
            sp.read_spec("Index 1 (i7)",None,"Native historical index acquisition/width unassigned","undetermined"),
            sp.read_spec("Index 2 (i5)",None,"Native historical index scheme/acquisition unassigned","undetermined"))
SEQUENCING_ENDING="Defining source25-cycle paired-end sequencing. Native adapter/index identity, acquired index cycles and run-primer landing remain unassigned in this capped pass."
SEQUENCING_INTRO=""
CAVEAT += " INFERRED — native library ends and opposing primer roles are unassigned. The exact Orsi oligos are a separate2015 reference and do not assign historical native index geometry or sequencing cycles."

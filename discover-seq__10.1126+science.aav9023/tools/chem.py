"""Protocol-owned source, oligos, molecular workflow and reaction graph."""
import illumina as il
import seqprimers as sp
from adapter_ligation import reported_adapter_end
from chemdraw import Construct, Scene, Segment, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (Disposition, Terminal, Highlight, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind)
from protocol_model import chain, molecule, role

DISCOVER_AR001 = "ATCACG"


def discover_library() -> tuple[Construct, tuple]:
    """Historical single-index TruSeq ChIP product with a six-cycle index read."""
    # The defining procedure delegates to TruSeq ChIP; the kit guide specifies
    # single six-base indexes. Canonical arms are a conditional endpoint model,
    # not disclosed ordering sequences of the physical PCR Primer Cocktail.
    lib=Construct([Segment("P5 overlapping Read1",il.P5[:-4],"p5",inferred=True),
        Segment("Read 1 arm",il.TRUSEQ_READ1,"r1",inferred=True),
        role("MRE11-bound genomic DNA from a Cas9 cut to a random sonication end"),
        Segment("dA junction","A"),Segment("Index 1 / Read 2 arm",il.INDEX1_PRIMER,"r2",inferred=True),
        Segment("i7 reverse complement","I"*6,"cbc",placeholder=True,
                feature=feature("sample_index_i7","sample_index","unknown")),
        Segment("P7 reverse complement",il.P7_RC,"p7",inferred=True)],name="DISCOVER-Seq TruSeq ChIP conditional endpoint")
    primers=(sp.TRUSEQ["R1"],sp.TRUSEQ["I1"],sp.TRUSEQ["R2"])
    errors=sp.verify(lib,primers,required_roles=tuple(p.role for p in primers))
    if errors: raise ValueError("DISCOVER-Seq library: "+"; ".join(errors))
    return lib,primers


def protein_dna_complex_scene(protein_label, dna_parts, label):
    """Display associated protein and DNA as separate molecular rows."""
    scene = Scene()
    scene.strand("dna", list(dna_parts), label=label)
    scene.strand("protein", [role(protein_label)], label=protein_label)
    scene.labels("dna")
    scene.labels("protein")
    return scene


def workflow():
    return chain(molecule("living cells or tissue receiving a CRISPR nuclease", (role("chromosomal on-target and off-target sites"),)), (
        ("Let MRE11 bind newly made double-strand breaks",
         protein_dna_complex_scene("MRE11 at the DNA end", (role("Cas9 cleavage end"), role("surrounding chromatin")), "edited chromatin with repair complexes"),
         "The sampling time is optimized for the delivery system; the detailed protocol suggests approximately 8–16 hours after RNP delivery."),
        ("Crosslink cells in 1% formaldehyde, lyse nuclei and sonicate chromatin",
         molecule("crosslinked chromatin fragments of a few hundred base pairs", (role("random sonication end"), role("MRE11-bound DNA reaching a nuclease cut"), role("Cas9 cleavage end"))),
         "1% formaldehyde for 15 min at room temperature is quenched with 125 mM glycine; a Covaris S2 shears the nuclei. The same detailed protocol states the target chromatin size twice with different windows (see the published-discrepancy notice)."),
        ("Immunoprecipitate MRE11 on antibody-loaded Protein A magnetic beads",
         protein_dna_complex_scene("antibody–MRE11 complex", (role("DNA adjacent to a cellular DSB"), role("Cas9 cleavage boundary")), "MRE11-enriched chromatin"),
         "The detailed protocol pre-binds a cross-reactive anti-MRE11 antibody (Novus NB100-142, or Abcam ab208020) to Protein A Dynabeads. An unedited or non-targeting sample measures background; ChIP enrichment is checked by qPCR before library construction."),
        ("Reverse crosslinks and purify the ChIP DNA",
         molecule("purified break-adjacent DNA", (role("random sonication end"), role("MRE11 footprint region"), role("native nuclease-cleavage end"))),
         "Crosslinks are reversed at 65 °C for at least 6 h in SDS elution buffer; RNase A (37 °C) then proteinase K (55 °C) remove RNA and protein, and the DNA is recovered on a Qiagen MinElute column, retaining the endpoints used by BLENDER."),
        ("End-repair and add one 3-prime dA",
         molecule("A-tailed ChIP fragments", (role("end-repaired random end"), role("break-adjacent insert"), role("end-repaired cleavage end"), role("3-prime dA"))),
         "The historical TruSeq ChIP kit blunts overhangs and then adenylates the repaired fragments."),
        ("Ligate a six-base indexed TruSeq ChIP adapter",
         molecule("adapter-ligated ChIP library", (role("TruSeq Read 1 side"), role("break-adjacent genomic insert"), role("TruSeq Read 2 side with six-base i7 index"))),
         "The discontinued kit supplies single-index adapters; Table 4/5 instructs recording six index bases in the sample sheet."),
        ("PCR-enrich, remove adapter dimers and pool",
         molecule("single-index DISCOVER-Seq library", (role("P5 and Read 1 arm"), role("MRE11-enriched insert"), role("Read 2 arm, i7 and P7"))),
         "Amplification completes the flow-cell ends; individual libraries are pooled in equal amounts and cleaned with 0.8× AMPure XP beads, and the pooled library is expected to average 300–500 bp."),
    ))


def oligos():
    return (f"TruSeq ChIP AR001 sample index example — 5′-{DISCOVER_AR001}-3′",
            "PCR Primer Cocktail PN15026786: physical oligo sequences not established in this model.")


TITLE = 'DISCOVER-Seq — MRE11 ChIP-seq at CRISPR cleavage sites'
NOTES = '01_discover_seq.html'
SOURCE = '<a href="https://doi.org/10.1126/science.aav9023">Wienert et al., Science (2019)</a>; stepwise implementation from the authors’ <a href="https://doi.org/10.1038/s41596-020-0309-5">Nature Protocols procedure</a> and the named TruSeq ChIP kit guide.'
SUMMARY = 'MRE11 recruited to newly cut chromosomal DNA is fixed in place and immunoprecipitated; ordinary paired-end ChIP-seq then preserves the sharp genomic endpoints that BLENDER recognizes as nuclease cuts.'
CAVEAT = 'INFERRED — canonical TruSeq library arms illustrate kit-compatible primer sites; they do not identify the supplied PCR Primer Cocktail sequences. DISCOVER-Seq does not use a custom capture oligo or UMI. Its defining implementation delegates library construction to the now-discontinued single-index TruSeq ChIP kit; index identities are public, while the exact supplied oligos remain unassigned in this model.'
FINAL_CAPTION = 'INFERRED — canonical TruSeq arm model. Either genomic orientation can be captured. One library strand is shown from a random sonication boundary through MRE11-enriched DNA to the cellular nuclease-cleavage boundary; the six-base i7 sample index comes from the historical TruSeq ChIP adapter.'
READ_LENGTHS = {'Read 1': 100, 'Index 1 (i7)': 6, 'Read 2': 100}
SEQUENCING_INTRO = 'The detailed protocol specifies 100-bp paired-end sequencing on a HiSeq 4000, NovaSeq or similar instrument at roughly 30–40 M read pairs per ChIP sample. Standard TruSeq Read 1, Index 1 and Read 2 primers bind the final library; paired fragment endpoints, not a barcode read, provide the cut-site signal.'
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''
FINAL_LIBRARY, SEQ_PRIMERS = discover_library()

def protocol_graph():
    source = SOURCE
    graph = ProtocolGraph(TITLE, authority=source)
    chromatin = graph.input("chromatin", "edited cells or tissue with nuclease-created breaks",
                            MaterialKind.BIOLOGICAL,
                            state=molecule("edited chromatin", (role("chromosomal nuclease cuts"),)),
                            primary=True)
    antibody = graph.input("mre11_antibody", "anti-MRE11 antibody (Novus NB100-142 or Abcam ab208020) on Protein A Dynabeads",
                        MaterialKind.SELECTION)
    mre11 = graph.input("mre11_protein", "endogenous MRE11 repair protein",MaterialKind.MOLECULAR)
    bound = graph.product("mre11_bound", "MRE11-bound break chromatin",
                          MaterialKind.MOLECULAR, state=chromatin.state)
    unbound = graph.product("mre11_unbound", "chromatin without captured MRE11-bound ends",
                            MaterialKind.MOLECULAR)
    graph.react("bind_mre11", "Recognize MRE11 at newly made DNA breaks", StepKind.BINDING,
                ((chromatin, InputRole.TEMPLATE), (mre11, InputRole.STRUCTURAL)),
                (bound, unbound), changes_molecule=False)
    graph.mark_terminal(Terminal(unbound.id,Disposition.UNRESOLVED,
        "MRE11-negative chromatin is not physically separated before fixation. Its later ChIP recovery is not resolved by this representative MRE11-bound branch."))
    formaldehyde = graph.input("formaldehyde", "1% formaldehyde, 15 min at room temperature, quenched with 125 mM glycine",
                               MaterialKind.CONDITION)
    fixed = graph.product("fixed", "cross-linked MRE11–break chromatin",
                          MaterialKind.MOLECULAR, state=bound.state)
    graph.react("crosslink", "Cross-link MRE11–DNA complexes", StepKind.CROSS_LINKING,
                ((bound, InputRole.TEMPLATE), (formaldehyde, InputRole.CONDITION)), (fixed,))
    shear = graph.input("sonication", "Covaris S2 sonication of lysed nuclei", MaterialKind.CONDITION)
    fragments = graph.product("fragments", "cross-linked chromatin fragments of a few hundred base pairs",
                              MaterialKind.MOLECULAR, state=fixed.state)
    graph.react("sonicate", "Sonicate fixed chromatin", StepKind.FRAGMENTATION,
                ((fixed, InputRole.TEMPLATE), (shear, InputRole.CONDITION)), (fragments,))
    beads = graph.input("protein_a_beads", "Protein A magnetic beads", MaterialKind.SELECTION)
    captured = graph.product("captured", "MRE11-enriched chromatin fragments",
                             MaterialKind.MOLECULAR, state=fragments.state)
    rejected = graph.product("chip_flowthrough", "MRE11-depleted ChIP flow-through",
                             MaterialKind.WASTE)
    graph.react("chip", "Immunoprecipitate MRE11-bound chromatin", StepKind.AFFINITY_CAPTURE,
                ((fragments, InputRole.TEMPLATE), (beads, InputRole.SELECTION), (antibody,InputRole.SELECTION)),
                (captured, rejected))
    graph.discard(rejected, "is not retained by MRE11 immunoprecipitation")
    reversal = graph.input("crosslink_reversal", "65 °C cross-link reversal (≥6 h), RNase A, proteinase K and MinElute recovery",
                           MaterialKind.CONDITION)
    dna = graph.product("chip_dna", "purified break-adjacent ChIP DNA",
                        MaterialKind.MOLECULAR, state=captured.state)
    graph.react("purify_dna", "Reverse cross-links and remove protein",
                StepKind.DEPROTEINIZATION,
                ((captured, InputRole.TEMPLATE), (reversal, InputRole.CONDITION)), (dna,))
    repair = graph.input("repair_enzyme", "TruSeq ChIP end-repair enzyme mix",
                         MaterialKind.CATALYST)
    repaired = graph.product("repaired", "end-repaired ChIP fragments",
                             MaterialKind.MOLECULAR, state=dna.state)
    graph.react("repair", "End-repair ChIP DNA", StepKind.END_REPAIR,
                ((dna, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    tailing = graph.input("tailing_enzyme", "single-dA tailing enzyme", MaterialKind.CATALYST)
    tailed = graph.product("tailed", "dA-tailed ChIP fragments",
                           MaterialKind.MOLECULAR, state=repaired.state)
    graph.react("tail", "Add one 3-prime dA", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    adapter = graph.input("truseq_adapter", "six-base indexed TruSeq ChIP adapter",
                          MaterialKind.ADAPTER)
    ligase = graph.input("adapter_ligase", "TruSeq adapter ligase", MaterialKind.CATALYST)
    ligated = graph.product("ligated", "adapter-ligated ChIP fragments",
                            MaterialKind.MOLECULAR, state=tailed.state)
    graph.react("ligate", "Ligate the indexed TruSeq ChIP adapter", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,),
                geometry=reported_adapter_end(
                    "TruSeq ChIP adapter", kind="cohesive", polarity="3-prime",
                    terminal_sequence="T",
                    statement="The named TruSeq ChIP workflow ligates its indexed adapter to dA-tailed fragments."))
    reviewed = ["chromatin", "mre11_antibody", "bind_mre11", "formaldehyde", "crosslink",
                "sonication", "sonicate", "protein_a_beads", "chip", "crosslink_reversal",
                "purify_dna", "repair_enzyme", "repair", "tailing_enzyme", "tail",
                "truseq_adapter", "adapter_ligase", "ligate"]
    forward=graph.input("library_forward","TruSeq ChIP PCR Primer Cocktail forward component",MaterialKind.OLIGO)
    reverse=graph.input("library_reverse","TruSeq ChIP PCR Primer Cocktail reverse component",MaterialKind.OLIGO)
    polymerase=graph.input("library_polymerase","TruSeq ChIP PCR Master Mix PMM PN15026785",MaterialKind.CATALYST)
    library=graph.product("library_library","INFERRED — canonical-arm DISCOVER-Seq library",MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    graph.react("library_pcr","Enrich adapter-ligated ChIP DNA with PPC PN15026786",StepKind.PCR,
        ((ligated,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
         (reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
        geometry=PCRGeometry.unavailable(ligated.state,
            reason="The inspected TruSeq ChIP guide names PCR Primer Cocktail PN15026786; exact physical primer sequences and ligated-template coordinates remain unassigned. Canonical final-library sites do not supply those missing input declarations."))
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS),read_lengths=READ_LENGTHS))
    reviewed.extend(("mre11_protein","library_forward","library_reverse","library_polymerase","library_pcr","library_library"))
    graph.annotate(SourceDiscrepancy(
        "sonicate", "The detailed protocol states two chromatin-size windows",
        (SourceClaim("Wienert et al., Nature Protocols (2020), PMC7305837 'Sonication' overview, line113",
                     "chromatin sheared to a fragment size of approximately 150–500 bp"),
         SourceClaim("Same article, 'Anticipated results', line994 (and the authors' protocols.io v2, step note, line292)",
                     "sheared chromatin should have a median size of 200–500 bp (gel smear 200–500 bp)")),
        "The model leaves the exact shearing window unassigned and draws the fragments as 'a few hundred base pairs'.",
        "Both passages are printed in the same original author protocol; neither is a text-extraction artefact, and the authors state that sonication must be optimized per cell type and sonicator.",
        "Whether 150–500 bp is an informal restatement or a different acceptance window is unresolved; no single numeric range is attributed to the sequenced libraries."))
    graph.annotate(Highlight(
        "library_pcr", "Pooled-library cleanup and size expectation",
        "The detailed protocol pools equal amounts of the individual TruSeq ChIP libraries, removes adapter dimers with 0.8× AMPure XP beads (repeating or using PippinPrep/gel if dimers persist) and expects a 300–500 bp mean library size. Adapter-dimer species are therefore an acknowledged by-product, and only their bead-based removal is represented as a physical separation."))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

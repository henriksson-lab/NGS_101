"""Mumbach et al. 2016 HiChIP molecular model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Segment, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind, Terminal)
from protocol_model import molecule, role

TITLE = "HiChIP — protein-directed chromatin contacts"
NOTES = "01_hichip.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/nmeth.3999">Mumbach et al., <i>Nature Methods</i> (2016)</a>.'
SUMMARY = "Build biotin-marked in situ Hi-C contacts, immunoprecipitate contacts associated with a chosen protein, capture biotin-bearing purified DNA on streptavidin beads, and tagment it into an Illumina library."
CAVEAT = "The contact panel derives the original MboI fill and GATCGATC ligation scar. The opposite-strand biotin is present in the duplex but only the top-strand biotin is separately coloured; the final example uses published Ad1_noMX and Ad2.1_TAAGGCGA, a single-i7 design without i5. Exact tagmentation population geometry remains unresolved."
# Published Supplementary Table 4: representative first indexed primer member.
AD1 = il.P5 + nx.READ1_PRIMER[:21]
AD2_INDEX_WRITTEN = "TCGCCTTA"
AD2 = il.P7 + AD2_INDEX_WRITTEN + nx.READ2_PRIMER[:21]
FINAL_LIBRARY = Construct([
    Segment("P5", il.P5, "p5"), Segment("Read 1 site", nx.READ1_PRIMER, "r1"),
    Segment("protein-enriched contact fragment", "X" * 36, placeholder=True),
    Segment("Read 2 site reverse complement", revcomp(nx.READ2_PRIMER), "r2"),
    Segment("representative i7", revcomp(AD2_INDEX_WRITTEN), "cbc",
            feature=feature("sample_i7", "sample_index", "whitelist", whitelist="HiChIP Supplementary Table 4")),
    Segment("P7 reverse complement", il.P7_RC, "p7"),
], name="HiChIP single-i7 library")
SEQ_PRIMERS = (sp.NEXTERA["R1"], sp.NEXTERA["I1"], sp.NEXTERA["R2"])
if problems := sp.verify(FINAL_LIBRARY, SEQ_PRIMERS, required_roles=tuple(p.role for p in SEQ_PRIMERS)):
    raise ValueError("invalid HiChIP library: " + "; ".join(problems))


def oligos():
    return [f"Ad1_noMX: 5′-{AD1}-3′", f"Ad2.1_TAAGGCGA (representative): 5′-{AD2}-3′"]
FINAL_CAPTION = "Tn5-built paired-end library from immunoprecipitated, biotin-bearing proximity-ligation products."
SEQUENCING_INTRO = "Paired reads recover the two genomic ends of each protein-enriched contact molecule; an internal ligation scar identifies informative contacts."
SHOW_READ_CONTENT = True

def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    cells = graph.input("cells", "input cells", MaterialKind.BIOLOGICAL,
                        state=molecule("cells", (role("native chromatin contacts"),)),
                        primary=True)
    formaldehyde = graph.input("formaldehyde", "formaldehyde cross-linking condition",
                               MaterialKind.CONDITION)
    fixed = graph.product("fixed_chromatin", "cross-linked chromatin",
                          MaterialKind.MOLECULAR, state=cells.state)
    graph.react("crosslink", "Cross-link chromatin contacts", StepKind.CROSS_LINKING,
                ((cells, InputRole.TEMPLATE), (formaldehyde, InputRole.CONDITION)),
                (fixed,))
    mboi = graph.input("mboi", "MboI restriction endonuclease", MaterialKind.CATALYST)
    digested = graph.product("digested_chromatin", "MboI-digested chromatin",
                             MaterialKind.MOLECULAR, state=fixed.state)
    graph.react("digest", "Digest fixed chromatin with MboI", StepKind.DIGESTION,
                ((fixed, InputRole.TEMPLATE), (mboi, InputRole.CATALYST)), (digested,))
    biotin_datp = graph.input("biotin_datp", "biotin-dATP fill-in nucleotide",
                              MaterialKind.MOLECULAR)
    fill_polymerase = graph.input("fill_polymerase", "DNA Polymerase I Large (Klenow) Fragment",
                                  MaterialKind.CATALYST)
    filled = graph.product("filled_ends", "biotin-filled MboI restriction ends",
                           MaterialKind.MOLECULAR, state=digested.state)
    graph.react("fill_ends", "Fill restriction ends with biotin-dATP",
                StepKind.EXTENSION,
                ((digested, InputRole.TEMPLATE), (biotin_datp, InputRole.STRUCTURAL),
                 (fill_polymerase, InputRole.CATALYST)), (filled,))
    ligase = graph.input("proximity_ligase", "T4 DNA ligase",
                         MaterialKind.CATALYST)
    contacts = graph.product(
        "contacts", "biotin-bearing in-situ proximity-ligation products",
        MaterialKind.MOLECULAR,
        state=molecule("HiChIP contact chromatin", (
            role("protein-associated chromatin"), role("biotin-marked GATCGATC junction"))))
    graph.react("proximity_ligation", "Ligate spatially adjacent filled restriction ends",
                StepKind.PROXIMITY_LIGATION,
                ((filled, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (contacts,))
    sonication = graph.input("sonication", "chromatin-sonication condition",
                             MaterialKind.CONDITION)
    sheared = graph.product("sheared_contacts", "sheared proximity-ligated chromatin",
                            MaterialKind.MOLECULAR, state=contacts.state)
    graph.react("shear", "Shear proximity-ligated chromatin", StepKind.FRAGMENTATION,
                ((contacts, InputRole.TEMPLATE), (sonication, InputRole.CONDITION)),
                (sheared,))
    antibody = graph.input("target_antibody", "target-specific antibody and capture beads",
                           MaterialKind.SELECTION)
    enriched = graph.product("protein_contacts", "protein-associated contact chromatin",
                             MaterialKind.MOLECULAR, state=sheared.state)
    flowthrough = graph.product("ip_flowthrough", "unbound chromatin fragments",
                                MaterialKind.WASTE)
    graph.react("immunoprecipitate", "Immunoprecipitate protein-associated contacts",
                StepKind.AFFINITY_CAPTURE,
                ((sheared, InputRole.TEMPLATE), (antibody, InputRole.SELECTION)),
                (enriched, flowthrough))
    graph.discard(flowthrough, "is not retained by the target-directed immunoprecipitation")
    proteinase = graph.input("proteinase_k", "Proteinase K and reverse-crosslinking conditions", MaterialKind.CATALYST)
    purified = graph.product("purified_dna", "eluted, reverse-crosslinked and column-purified post-ChIP DNA", MaterialKind.MOLECULAR)
    graph.react("release_dna", "Elute from IP beads, reverse crosslinks and column-purify DNA",
                StepKind.MAGIC, ((enriched, InputRole.TEMPLATE), (proteinase, InputRole.CATALYST)), (purified,), reason="Elution, proteolysis, reverse crosslinking and column recovery are source-stated; their molecular transitions and removed populations are not reconstructed.")
    streptavidin = graph.input("streptavidin", "streptavidin capture beads",
                               MaterialKind.SELECTION)
    captured = graph.product("captured_contacts", "biotin-bearing purified contact DNA",
                             MaterialKind.MOLECULAR, state=purified.state)
    nonjunction = graph.product("nonjunction_fragments", "post-ChIP DNA not retained by biotin capture",
                                MaterialKind.WASTE)
    graph.react("capture_biotin", "Capture biotin-bearing contact products",
                StepKind.AFFINITY_CAPTURE,
                ((purified, InputRole.TEMPLATE), (streptavidin, InputRole.SELECTION)),
                (captured, nonjunction))
    graph.discard(nonjunction, "does not retain a biotin-marked proximity junction")
    loaded_adapter = graph.input("loaded_nextera_adapter", "Nextera adapter-loaded transposome",
                                 MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagmented_contacts", "on-bead tagmented contact fragments",
                           MaterialKind.MOLECULAR, state=captured.state)
    tn5_waste = graph.product("tn5_waste", "candidate library-ineligible tagmentation products",
                              MaterialKind.MOLECULAR, inferred=True)
    graph.react("tagment", "Tagment streptavidin-bound purified contact DNA",
                StepKind.TAGMENTATION,
                ((captured, InputRole.TEMPLATE), (loaded_adapter, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, tn5_waste))
    graph.mark_terminal(Terminal(tn5_waste.id, Disposition.UNRESOLVED,
        "Candidate end-combination products, not an experimentally established population; "
        "PCR ineligibility does not establish physical removal. Exact end-pool geometry remains unmodelled."))
    forward = graph.input("pcr_forward", "Ad1_noMX universal PCR primer",
                          MaterialKind.OLIGO, state=Construct([Segment("Ad1_noMX", AD1)]))
    reverse = graph.input("pcr_reverse", "Ad2.1 indexed PCR primer (representative)",
                          MaterialKind.OLIGO, state=Construct([Segment("Ad2.1", AD2)]))
    polymerase = graph.input("pcr_polymerase", "Phusion HF polymerase",
                             MaterialKind.CATALYST)
    library = graph.product("amplified_library", FINAL_LIBRARY.name,
                            MaterialKind.MOLECULAR, state=FINAL_LIBRARY)
    graph.react("library_pcr", "Complete and index the tagmented contact library",
                StepKind.PCR,
                ((tagged, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,
                    reason="The pre-PCR on-bead contact population is not represented as one exact construct."))
    gel = graph.input("size_gel", "6% PAGE, selected 300–700 bp interval", MaterialKind.SELECTION)
    selected = graph.product("sequencing_library", FINAL_LIBRARY.name, MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    offsize = graph.product("offsize", "products outside the recovered PAGE interval", MaterialKind.WASTE)
    graph.react("size_selection", "Recover 300–700 bp libraries by the published PAGE branch",
                StepKind.SELECTION, ((library, InputRole.TEMPLATE), (gel, InputRole.SELECTION)),
                (selected, offsize))
    graph.discard(offsize, "outside the recovered gel interval; alternative two-sided AMPure branch not selected")
    graph.sequence(selected, SequencingPlan(
        library=FINAL_LIBRARY, primers=tuple(SEQ_PRIMERS),
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
        read_lengths={"Read 1": 75, "Read 2": 75},
        limitation="Published primer table supports a single eight-base i7; acquired index cycle count is not assigned here."))
    graph.declare_reviewed_scope((
        "cells", "formaldehyde", "crosslink", "mboi", "digest", "biotin_datp",
        "fill_polymerase", "fill_ends", "proximity_ligase", "proximity_ligation",
        "sonication", "shear", "target_antibody", "immunoprecipitate",
        "proteinase_k", "release_dna", "loaded_nextera_adapter", "tn5", "tagment", "streptavidin",
        "capture_biotin", "pcr_forward", "pcr_reverse", "pcr_polymerase",
        "library_pcr", "size_gel", "size_selection", "sequencing_library"), SOURCE)
    return graph.require_publishable()

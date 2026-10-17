
from conformation_replication_specs import _three_c
from protocol_graph import SourceClaim, SourceDiscrepancy
from conformation_replication_specs import AdapterLigationGeometry, Construct, Definition, InputRole, MaterialKind, PCRGeometry, ProtocolGraph, SequencingPlan, StepKind, da_tailed, molecule, role, sp
DEFINITION = Definition(
        "NG Capture-C — next-generation Capture-C", "10.1038/nmeth.3664",
        "Davies et al., Nature Methods (2016)",
        "A DpnII 3C library is sonicated, converted to an indexed TruSeq library and enriched by two sequential rounds of hybridization to 5′-biotinylated 120-nt oligos, streptavidin pull-down and PCR from the beads.",
        (role("captured viewpoint"), role("proximity-ligation junction"), role("reporter fragment")),
        "Paired-end Illumina reads (mostly 150-nt MiSeq; one 100-nt HiSeq run) recover the viewpoint-containing sonication fragment and its ligated reporter. PCR duplicates are removed computationally: after in silico DpnII digestion, reads with identical start and end coordinates of each sub-fragment are collapsed.", True)


SOURCE_METHODS = ("Online Methods of the accepted author manuscript, "
                  "PMC4724891 supplementary file NIHMS65861-supplement-01")
RESTRICTION_ENZYME = "DpnII"


def _ng_capture_c_graph(definition, final: Construct,
                        primers: tuple[sp.SeqPrimer, ...]) -> ProtocolGraph:
    """NG Capture-C graph from the Online Methods (3C, library prep, two captures)."""
    source = f"{definition.citation}; doi:{definition.doi}; {SOURCE_METHODS}"
    graph = ProtocolGraph(definition.title, authority=source)
    fixed = graph.input("fixed_chromatin", "2% formaldehyde-fixed cells, lysed to nuclei",
                        MaterialKind.BIOLOGICAL,
                        state=molecule("formaldehyde-fixed chromatin",
                                       ("viewpoint", "proximal chromatin")), primary=True)
    restriction = graph.input("primary_restriction_enzyme",
                              f"{RESTRICTION_ENZYME} (3 x 500 U, 37 °C, 16–24 h, after SDS/Triton X-100)",
                              MaterialKind.CATALYST)
    digested = graph.product(
        "restriction_fragments", "DpnII-digested fixed chromatin",
        MaterialKind.MOLECULAR,
        state=molecule("restriction fragments", ("viewpoint fragment",
                                                   "proximal reporter fragments")))
    graph.react("primary_digest", f"Digest fixed chromatin with {RESTRICTION_ENZYME}",
                StepKind.DIGESTION,
                ((fixed, InputRole.TEMPLATE),
                 (restriction, InputRole.CATALYST)), (digested,))
    ligase = graph.input("proximity_ligase",
                         "high-concentration T4 DNA ligase (Thermo Scientific), 16 °C overnight",
                         MaterialKind.CATALYST)
    contact_model = _three_c(RESTRICTION_ENZYME,
                             ("captured viewpoint fragment", "reporter fragment"))
    contact_state = contact_model.scene()
    contact = graph.product("contact_product", "proximity-ligated 3C contact",
                            MaterialKind.MOLECULAR, state=contact_state)
    graph.react("proximity_ligation",
                "Heat-inactivate DpnII (65 °C) and proximity-ligate restriction fragments",
                StepKind.PROXIMITY_LIGATION,
                ((digested, InputRole.TEMPLATE),
                 (ligase, InputRole.CATALYST)), (contact,))
    proteinase = graph.input("crosslink_reversal",
                             "Proteinase K at 65 °C overnight, RNase, phenol-chloroform and ethanol precipitation",
                             MaterialKind.CONDITION)
    deproteinized = graph.product("deproteinized_3c_dna", "deproteinized 3C library DNA",
                                  MaterialKind.MOLECULAR, state=contact_state)
    graph.react("reverse_crosslinks", "Reverse cross-links and recover 3C DNA",
                StepKind.DEPROTEINIZATION,
                ((contact, InputRole.TEMPLATE),
                 (proteinase, InputRole.CONDITION)), (deproteinized,))
    sonication = graph.input("sonication",
                             "Covaris S220 sonication of 5 µg 3C DNA to about 200 bp",
                             MaterialKind.CONDITION)
    sonicated = graph.product("sonicated_3c_dna", "sonicated 3C DNA",
                              MaterialKind.MOLECULAR, state=deproteinized.state)
    graph.react("sonicate", "Sonicate recovered 3C DNA", StepKind.FRAGMENTATION,
                ((deproteinized, InputRole.TEMPLATE),
                 (sonication, InputRole.CONDITION)), (sonicated,))
    repair = graph.input("end_repair_enzymes", "NEBNext end-repair reagents",
                         MaterialKind.CATALYST)
    repaired = graph.product("end_repaired_dna", "end-repaired sonicated 3C DNA",
                             MaterialKind.MOLECULAR, state=sonicated.state)
    graph.react("end_repair", "End-repair sonicated 3C DNA", StepKind.END_REPAIR,
                ((sonicated, InputRole.TEMPLATE),
                 (repair, InputRole.CATALYST)), (repaired,))
    tailing = graph.input("a_tailing_enzyme", "NEBNext dA-tailing reagents",
                          MaterialKind.CATALYST)
    tailed_state = da_tailed("dA-tailed 3C fragment", "viewpoint-contact fragment")
    tailed = graph.product("da_tailed_dna", "dA-tailed 3C DNA",
                           MaterialKind.MOLECULAR, state=tailed_state)
    graph.react("da_tail", "Add one 3-prime dA to each end", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE),
                 (tailing, InputRole.CATALYST)), (tailed,))
    adapters = graph.input("illumina_adapters", "Illumina TruSeq sequencing adapters",
                           MaterialKind.ADAPTER, inferred=True)
    adapter_ligase = graph.input("adapter_ligase", "NEBNext adapter-ligation reagents",
                                 MaterialKind.CATALYST)
    adapted = graph.product("adapter_ligated_dna", "adapter-ligated 3C DNA",
                            MaterialKind.MOLECULAR, state=final, inferred=True)
    graph.react("adapter_ligation", "Ligate TruSeq adapters (AMPure XP 1:1.8 cleanups)",
                StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE),
                 (adapters, InputRole.STRUCTURAL),
                 (adapter_ligase, InputRole.CATALYST)), (adapted,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The Online Methods name Illumina TruSeq indexed adapters ligated with NEBnext "
                    "reagents (E6000/E6040/E7335/E7500; E7370 for low cell numbers) but print "
                    "no adapter strands or end geometry."))
    polymerase = graph.input("pcr_polymerase", "Herculase II PCR kit (Agilent), 6–8 cycles",
                             MaterialKind.CATALYST)
    forward = graph.input("index_forward_primer", "library PCR primer 1 (sequence not printed)",
                          MaterialKind.OLIGO)
    reverse = graph.input("index_reverse_primer", "TruSeq index-adding PCR primer 2 (sequence not printed)",
                          MaterialKind.OLIGO)
    indexed = graph.product("indexed_library", "indexed 3C library",
                            MaterialKind.MOLECULAR, state=final)
    graph.react("index_pcr", "Index and amplify the 3C library", StepKind.PCR,
                ((adapted, InputRole.TEMPLATE),
                 (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (indexed,),
                geometry=PCRGeometry.unavailable(
                    final, reason=(
                        "The Online Methods say libraries were indexed with TruSeq indexed adapters "
                        "using NEBNext reagents for adaptor ligation and PCR indexing, and that 6–8 "
                        "Herculase II cycles added the TruSeq indices; whether the index enters on "
                        "the adapter or the PCR primer, and the primer sequences, are not printed.")))
    reviewed = [
        "fixed_chromatin", "primary_restriction_enzyme", "primary_digest",
        "proximity_ligase", "proximity_ligation", "crosslink_reversal",
        "reverse_crosslinks", "sonication", "sonicate", "end_repair_enzymes",
        "end_repair", "a_tailing_enzyme", "da_tail", "illumina_adapters",
        "adapter_ligase", "adapter_ligation", "pcr_polymerase",
        "index_forward_primer", "index_reverse_primer", "index_pcr",
    ]

    beads = graph.input("streptavidin_beads",
                        "M270 streptavidin beads with SeqCap EZ wash buffers (47 °C stringent washes)",
                        MaterialKind.SELECTION)
    lm_polymerase = graph.input(
        "post_capture_pcr_mix",
        "SeqCap EZ Post-Capture LM-PCR master mix (18 cycles) or KAPA mix from SeqCap EZ accessory kit v2 (14 cycles)",
        MaterialKind.CATALYST)
    lm_forward = graph.input("post_capture_primer_1", "SeqCap EZ Post LM-PCR oligo, first of pair (sequence not printed)",
                             MaterialKind.OLIGO)
    lm_reverse = graph.input("post_capture_primer_2", "SeqCap EZ Post LM-PCR oligo, second of pair (sequence not printed)",
                             MaterialKind.OLIGO)
    template = indexed
    for n, hours in ((1, "64–72 h"), (2, "24 h")):
        probes = graph.input(f"capture_oligo_round{n}",
                             f"5′-biotinylated 120-nt viewpoint capture oligo pool, round {n} (13 pmol)",
                             MaterialKind.OLIGO)
        blockers = graph.input(
            f"hybridization_blockers_round{n}",
            "COT DNA, NimbleGen HE universal and index-specific blocking oligos, NimbleGen hybridization buffer",
            MaterialKind.CONDITION)
        hybrid = graph.product(f"capture_hybrid_round{n}",
                               f"round-{n} library hybridized to biotinylated probes",
                               MaterialKind.MOLECULAR, state=final)
        substrate = ("1.5–2 µg indexed library, or pooled indexed libraries"
                     if n == 1 else "75% of round-1 amplified material, up to 2 µg")
        graph.react(f"capture_hybridization_round{n}",
                    f"Denature and hybridize {substrate} to capture oligos at 47 °C for {hours}",
                    StepKind.HYBRIDIZATION,
                    ((template, InputRole.TEMPLATE),
                     (probes, InputRole.STRUCTURAL),
                     (blockers, InputRole.CONDITION)), (hybrid,))
        bound = graph.product(f"capture_round{n}", f"round-{n} streptavidin-bead-bound library",
                              MaterialKind.MOLECULAR, state=final)
        reject = graph.product(f"capture_round{n}_reject", f"uncaptured round-{n} library",
                               MaterialKind.WASTE)
        graph.react(f"capture_round{n}_step",
                    f"Round-{n} streptavidin pull-down and SeqCap EZ washes",
                    StepKind.AFFINITY_CAPTURE,
                    ((hybrid, InputRole.TEMPLATE),
                     (beads, InputRole.SELECTION)), (bound, reject))
        graph.discard(reject, "not bound to a biotinylated probe; removed with the "
                              "bead-binding and wash supernatants")
        last = n == 2
        amplified = graph.product(
            "sequencing_library" if last else f"capture_round{n}_amplified",
            definition.title + " captured library" if last else "round-1 amplified captured library",
            MaterialKind.LIBRARY if last else MaterialKind.MOLECULAR, state=final)
        graph.react(f"post_capture_pcr_round{n}",
                    f"Amplify round-{n} captured material from the beads (two reactions, AMPure XP cleanup)",
                    StepKind.PCR,
                    ((bound, InputRole.TEMPLATE),
                     (lm_forward, InputRole.STRUCTURAL),
                     (lm_reverse, InputRole.STRUCTURAL),
                     (lm_polymerase, InputRole.CATALYST)), (amplified,),
                    geometry=PCRGeometry.unavailable(
                        final, reason=(
                            "Fig. 1a legend: PCR from beads uses the adapter P5 and P7 sequences; "
                            "the SeqCap EZ post-LM-PCR oligo sequences are vendor reagents not printed in the source.")))
        reviewed.extend((
            f"capture_oligo_round{n}", f"hybridization_blockers_round{n}",
            f"capture_hybridization_round{n}", f"capture_round{n}_step",
            f"post_capture_pcr_round{n}",
        ))
        template = amplified
    reviewed.extend(("streptavidin_beads", "post_capture_pcr_mix",
                                "post_capture_primer_1", "post_capture_primer_2",
                                "sequencing_library"))
    graph.sequence(template, SequencingPlan(
        final, tuple(primers), tuple(p.role for p in primers)))
    graph.annotate(SourceDiscrepancy(
        target="capture_round2_step",
        title="Reported double-capture enrichment",
        claims=(SourceClaim("Results paragraph on single versus double capture (PMC4724891 JATS)",
                            "two sequential capture steps resulted in up to 1,000,000-fold enrichment"),
                SourceClaim("Figure 1a and Figure 2a legends (PMC4724891 JATS)",
                            "double capture produced enrichments up to 3,000,000-fold")),
        model_choice="The model encodes no fold-enrichment value; the two capture rounds are drawn without a yield claim.",
        evidence="Both values are printed in the open-access accepted manuscript; the paywalled publisher version of the main text was not inspected.",
        uncertainty="Which maximum belongs to which experiment, and whether the published version harmonizes the numbers, is unresolved."))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _ng_capture_c_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)


from types import SimpleNamespace
from conformation_replication_specs import chain, _three_c
from chromosome_capture import ThreeCLibrary
from replication_labeling import CopyNumberReplication

def workflow():
    contact = _three_c(RESTRICTION_ENZYME, ("captured viewpoint fragment", "reporter fragment"))
    adapted = molecule("indexed 3C sequencing library",
                       ("left Illumina adapter", "viewpoint", "contact junction", "reporter", "right Illumina adapter"),
                       junctions=(("left Illumina adapter", "viewpoint", "adapter ligation"),
                                  ("viewpoint", "contact junction", "proximity ligation"),
                                  ("reporter", "right Illumina adapter", "adapter ligation")))
    return chain(molecule("formaldehyde-fixed chromatin", ("viewpoint", "proximal chromatin")), (
        ("Digest with DpnII and proximity-ligate to make a 3C library", contact.scene(),
         "No biotin fill-in is used to select the contact junction."),
        ("Sonicate to about 200 bp, end-repair, dA-tail, ligate TruSeq adapters and index by 6–8 PCR cycles", adapted,
         "Random sonication ends let identical-coordinate PCR duplicates be removed computationally."),
        ("Two rounds of hybridization to biotinylated oligos, streptavidin capture and PCR from the beads", adapted,
         "Uncaptured library is washed away; each round ends with post-capture PCR."),
    ))

d = DEFINITION
if d.standard_truseq:
    from batch_ngs import truseq_library
    FINAL_LIBRARY, SEQ_PRIMERS = truseq_library(list(d.final_parts), d.title + " final library", dual_index=False)
else:
    FINAL_LIBRARY = Construct([role("Illumina-compatible left end"), *d.final_parts, role("Illumina-compatible right end")], name=d.title + " final library")
    SEQ_PRIMERS = ()
OLIGOS = (
    "<b>Published Hba-1 viewpoint capture example</b> — 5′-/5Biosg/GAT CCA TGG TAG CAC AGG GCA GCT AAG ATG CAA GTC TGA AGG AGG AGT CTG GCG AGC TGC TCC TGC AGT TCC CTG GAC CCA GAA GGA TGA GCT AGC AGA TTC ACT TGA GCC AAA GGA TTC-3′. The 120-base probe and 5′ biotin notation are from the original Supplementary Data workbook, ‘Oligonucleotide sequences’ sheet, row 94; capture probes vary by viewpoint.",
)
SPEC = SimpleNamespace(
    TITLE=d.title, NOTES="01_ng_capture_c.html",
    SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
    SUMMARY=d.summary, CAVEAT=d.caveat, FINAL_LIBRARY=FINAL_LIBRARY,
    SEQ_PRIMERS=SEQ_PRIMERS, FINAL_CAPTION=d.sequencing,
    SEQUENCING_INTRO=d.sequencing, READ_LENGTHS=None,
    SEQUENCING_ENDING=d.sequencing,
    SEQUENCING_UNAVAILABLE=("The defining protocol uses locus-specific or custom primer sequence that cannot be represented by one universal binding sequence."
                            if not SEQ_PRIMERS else ""),
    workflow=workflow, oligos=lambda: OLIGOS,
)
SPEC.protocol_graph = protocol_graph

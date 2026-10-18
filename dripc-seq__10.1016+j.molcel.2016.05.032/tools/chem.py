"""Source-owned 2016 DRIPc-seq pathway; the 2019 procedure is a comparator."""
from types import SimpleNamespace
from chemdraw import Construct, Scene, Segment
from pcr import PCRGeometry
from protocol_graph import (Disposition, Highlight, InputRole, MaterialKind,
                            ProtocolGraph, SequencingPlan, StepKind, Terminal)
from protocol_model import molecule, role

TITLE = "DRIPc-seq — strand-specific sequencing of R-loop RNA"
NOTES = "01_dripc_seq.html"
SOURCE = ('<a href="https://doi.org/10.1016/j.molcel.2016.05.032">'
          'Sanz et al., Molecular Cell (2016)</a>, Experimental Procedures, '
          'DRIP-seq and DRIPc-seq (PMC4955522, paragraph 137).')
SUMMARY = ("After S9.6 immunoprecipitation, DNase treatment and RNA recovery precede "
           "reverse transcription and dUTP second-strand synthesis. UDG treatment before "
           "PCR provides strand specificity without establishing physical disposal of the second strand.")
CAVEAT = ("INFERRED — the 2016 molecular end sequences, adapter implementation and read acquisition "
          "remain unassigned. The graph follows the defining 2016 method; detailed end-repair, "
          "ligation, primer and selection choices from the later 2019 protocol are not substituted. "
          "The 2019 primers below are comparator inventory only. UDG excises uracil bases; "
          "the subsequent physical integrity and recovery of the affected strand are unresolved.")
PCR_PRIMER_P5 = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGA"
PCR_PRIMER_P7 = "CAAGCAGAAGACGGCATACGAGAT"
FINAL_LIBRARY = Construct([
    role("2016 library left end; unassigned", inferred=True),
    role("amplified first-strand cDNA from R-loop RNA"),
    role("2016 library right end; unassigned", inferred=True),
], name="2016 DRIPc-seq role-level sequencing library")
SEQ_PRIMERS = ()
READ_LENGTHS = None
SEQUENCING_INTRO = ("The defining article reports HiSeq 2000 sequencing. Read roles, lengths, "
                    "index acquisition and primer direction relative to the original RNA remain unassigned. "
                    "The later 2019 protocol's SE100 example is not an acquisition declaration for 2016.")
SEQUENCING_ENDING = SEQUENCING_INTRO
SEQUENCING_UNAVAILABLE = "Exact 2016 sequencing-primer sites and acquisition mask await the applicable originals."
FINAL_CAPTION = "INFERRED — unknown-length 2016 library ends; strand selection retains first-strand cDNA."


def pcr_primer_p5():
    """Printed 2019 comparator, not a 2016 molecular input."""
    return Construct([Segment("2019 P5 primer", PCR_PRIMER_P5)])


def pcr_primer_p7():
    """Printed 2019 comparator, not a 2016 molecular input."""
    return Construct([Segment("2019 P7 primer", PCR_PRIMER_P7)])


def oligos():
    return (
        f"<b>2019 comparator PCR primer 1.0 P5</b> — 5′-{PCR_PRIMER_P5}-3′",
        f"<b>2019 comparator PCR primer 2.0 P7</b> — 5′-{PCR_PRIMER_P7}-3′",
    )


def rloop_state():
    scene = Scene()
    for name in ("DNA template in RNA:DNA hybrid", "hybrid RNA", "displaced DNA strand"):
        scene.strand(name, [role(name)], label=name + "; bases/extent unassigned")
    return scene


def first_strand_state():
    scene = Scene()
    scene.strand("RNA template", [role("recovered R-loop RNA")], label="RNA template population")
    scene.strand("first cDNA", [role("first-strand cDNA")], label="reverse-transcribed DNA; paired extent unassigned")
    return scene


def marked_duplex():
    scene = Scene.duplex([role("RNA-derived cDNA insert")], label="first/second cDNA duplex; bases unassigned", unpaired=("RNA-derived cDNA insert'",))
    scene.note("top", "first-strand cDNA")
    scene.note("bottom", "second strand incorporates dUTP instead of dTTP; no base-level U positions assigned")
    return scene


def retained_precursor():
    return molecule("retained first-strand cDNA entering PCR", (
        role("2016 first-strand library end; unassigned", inferred=True),
        role("first-strand cDNA insert"),
        role("opposite library end; unassigned", inferred=True)))


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    sample = g.input("rloop_chromatin", "R-loop-containing cellular nucleic acid", MaterialKind.BIOLOGICAL,
                     state=rloop_state(), primary=True)
    extraction = g.input("extraction_condition", "DRIP extraction per the cited Ginno 2012 method; current detailed conditions unassigned", MaterialKind.CONDITION)
    extracted = g.product("extracted", "extracted R-loop material", MaterialKind.MOLECULAR, state=rloop_state())
    g.react("extract", "Prepare nucleic acid for DRIP", StepKind.EXTRACTION,
            ((sample, InputRole.TEMPLATE), (extraction, InputRole.CONDITION)), (extracted,))
    restriction = g.input("restriction_mix", "referenced DRIP fragmentation preparation; 2016 enzyme assignment unassigned", MaterialKind.CATALYST)
    fragments = g.product("fragments", "R-loop material for immunoprecipitation", MaterialKind.MOLECULAR, state=rloop_state())
    g.react("restriction_digest", "Prepare DRIP fragments under the referenced method", StepKind.DIGESTION,
            ((extracted, InputRole.TEMPLATE), (restriction, InputRole.CATALYST)), (fragments,))
    antibody = g.input("s96", "S9.6 hybrid antibody", MaterialKind.SELECTION)
    support = g.input("capture_support", "DRIP capture support; detailed formulation unassigned", MaterialKind.SELECTION)
    captured = g.product("captured", "captured RNA:DNA hybrids and displaced DNA", MaterialKind.MOLECULAR, state=rloop_state())
    unbound = g.product("unbound", "non-retained immunoprecipitation fraction", MaterialKind.WASTE)
    g.react("immunoprecipitate", "Capture hybrids by the referenced DRIP method", StepKind.AFFINITY_CAPTURE,
            ((fragments, InputRole.TEMPLATE), (antibody, InputRole.SELECTION), (support, InputRole.SELECTION)), (captured, unbound))
    g.discard(unbound, "non-retained material is separated in immunoprecipitation; detailed washes remain unassigned")
    dnase = g.input("dnase", "2016: 4 U DNase I, 35 min at 37 °C", MaterialKind.CATALYST)
    rna = g.product("released_rna", "R-loop RNA after DNA digestion", MaterialKind.MOLECULAR,
                    state=molecule("released RNA", (role("R-loop RNA"),)))
    degraded = g.product("digested_dna", "DNA digestion products in the reaction", MaterialKind.MOLECULAR)
    g.react("digest_dna", "Digest the captured DNA component", StepKind.DIGESTION,
            ((captured, InputRole.TEMPLATE), (dnase, InputRole.CATALYST)), (rna, degraded))
    g.mark_terminal(Terminal(degraded.id, Disposition.UNRESOLVED,
        "Hydrolysis does not physically discard DNA-derived material. Its recovery through the later precipitation is unassigned."))
    precipitation = g.input("rna_recovery", "ethanol precipitation", MaterialKind.SELECTION)
    recovered = g.product("recovered_rna", "ethanol-recovered R-loop RNA", MaterialKind.MOLECULAR, state=rna.state)
    g.react("recover_rna", "Recover the RNA by ethanol precipitation", StepKind.GENERIC,
            ((rna, InputRole.TEMPLATE), (precipitation, InputRole.SELECTION)), (recovered,), changes_molecule=False)
    primers = g.input("rt_primers", "poly(T) and random hexamers in iScript RT supermix", MaterialKind.OLIGO)
    rt = g.input("reverse_transcriptase", "iScript reverse transcription supermix (Bio-Rad)", MaterialKind.CATALYST)
    first = g.product("first_cdna", "RNA/first-cDNA products", MaterialKind.MOLECULAR, state=first_strand_state())
    g.react("reverse_transcription", "Reverse-transcribe recovered R-loop RNA", StepKind.REVERSE_TRANSCRIPTION,
            ((recovered, InputRole.TEMPLATE), (primers, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (first,))
    second_mix = g.input("second_strand_mix", "2016 second-strand synthesis preparation; enzyme identities unassigned", MaterialKind.CATALYST)
    dutp = g.input("dutp", "dUTP replaces dTTP in second-strand synthesis", MaterialKind.CONDITION)
    duplex = g.product("marked_duplex", "first cDNA paired to dUTP-containing second DNA", MaterialKind.MOLECULAR, state=marked_duplex())
    other_rna = g.product("rt_rna_fate", "RNA-derived material after second-strand synthesis", MaterialKind.MOLECULAR)
    g.react("second_strand", "Synthesize the dUTP-marked second strand", StepKind.EXTENSION,
            ((first, InputRole.TEMPLATE), (second_mix, InputRole.CATALYST), (dutp, InputRole.CONDITION)), (duplex, other_rna))
    g.mark_terminal(Terminal(other_rna.id, Disposition.UNRESOLVED,
        "The defining paragraph does not assign the detailed RNA-processing chemistry; 2019 enzyme choices are not transferred."))
    qc_condition = g.input("qc_sampling", "aliquot for R-loop-positive versus R-loop-negative locus qPCR; volume unassigned", MaterialKind.CONDITION)
    bulk = g.product("bulk_cdna", "cDNA retained for library construction", MaterialKind.MOLECULAR, state=duplex.state)
    qc = g.product("qc_cdna", "cDNA enrichment QC aliquot", MaterialKind.QC, state=duplex.state)
    g.react("sample_qc", "Check enrichment on newly synthesized cDNA", StepKind.SAMPLING,
            ((duplex, InputRole.TEMPLATE), (qc_condition, InputRole.CONDITION)), (bulk, qc), changes_molecule=False)
    g.mark_terminal(Terminal(qc.id, Disposition.QC, "2016 locus-enrichment qPCR; molecular assay not reconstructed"))
    sonication = g.input("sonication", "sonication to approximately 200 bp", MaterialKind.CONDITION)
    sheared = g.product("sheared_cdna", "sonicated marked cDNA", MaterialKind.MOLECULAR, state=duplex.state)
    g.react("fragment_cdna", "Sonicate cDNA before library construction", StepKind.FRAGMENTATION,
            ((bulk, InputRole.TEMPLATE), (sonication, InputRole.CONDITION)), (sheared,))
    construction = g.input("construction_reagents", "2016 library-construction reagents; identities and ordering unassigned", MaterialKind.CONDITION)
    tagged = g.product("udg_input", "library-processed cDNA before UDG", MaterialKind.MOLECULAR, inferred=True)
    g.react("library_construction", "INFERRED — unresolved 2016 library-end construction", StepKind.MAGIC,
            ((sheared, InputRole.TEMPLATE), (construction, InputRole.CONDITION)), (tagged,),
            reason="The 2016 method specifies sonication, library construction and size selection but the applicable detailed 2016 supplement was not available in the local bundle for this bounded pass. The 2019 end-repair, dA-tailing, ligation, selected adapter and size-window choices are not assignments for this 2016 graph.")
    udg = g.input("udg", "UDG DNA glycosylase before PCR", MaterialKind.CATALYST)
    retained = g.product("retained_cdna", "retained first-strand library precursor", MaterialKind.MOLECULAR,
                         state=retained_precursor(), inferred=True)
    affected = g.product("affected_second_strand", "UDG-affected dUTP second strand and excised uracil", MaterialKind.MOLECULAR)
    g.react("udg_treatment", "INFERRED ends — retain first cDNA after UDG treatment", StepKind.CHEMICAL_CONVERSION,
            ((tagged, InputRole.TEMPLATE), (udg, InputRole.CATALYST)), (retained, affected))
    g.mark_terminal(Terminal(affected.id, Disposition.UNRESOLVED,
        "Base excision and suppression of second-strand amplification do not establish physical strand disposal; fragmentation and recovery remain unassigned."))
    forward = g.input("pcr_forward", "2016 library forward primer; identity unassigned", MaterialKind.OLIGO)
    reverse = g.input("pcr_reverse", "2016 library opposite-end primer; identity unassigned", MaterialKind.OLIGO)
    pol = g.input("pcr_polymerase", "2016 library amplification preparation; identity unassigned", MaterialKind.CATALYST)
    library = g.product("library", "strand-specific 2016 DRIPc library", MaterialKind.LIBRARY, state=FINAL_LIBRARY, inferred=True)
    g.react("library_pcr", "INFERRED ends — amplify the retained cDNA library strand", StepKind.PCR,
            ((retained, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL), (reverse, InputRole.STRUCTURAL), (pol, InputRole.CATALYST)), (library,),
            geometry=PCRGeometry.unavailable(retained.state, reason="Actual preceding first-strand precursor; 2016 adapter and PCR binding sites remain unassigned. The printed 2019 pair is comparison inventory, not evidence of the historical pair."))
    g.sequence(library, SequencingPlan(library=FINAL_LIBRARY, acquisition="unspecified",
                                     limitation=SEQUENCING_INTRO))
    g.annotate(Highlight("digest_dna", "Separate protocol versions",
        "The 2016 graph uses its source value: 4 U DNase I for 35 min. The later 2019 protocol states 6 U for 45 min; that later procedure and its exact P5/P7 primers are comparators, not a conflict or substitutes for the 2016 implementation."))
    g.declare_reviewed_scope(("rloop_chromatin", "extraction_condition", "extract", "restriction_mix", "restriction_digest",
        "s96", "capture_support", "immunoprecipitate", "dnase", "digest_dna", "rna_recovery", "recover_rna",
        "rt_primers", "reverse_transcriptase", "reverse_transcription", "second_strand_mix", "dutp", "second_strand",
        "qc_sampling", "sample_qc", "sonication", "fragment_cdna", "construction_reagents", "library_construction",
        "udg", "udg_treatment", "pcr_forward", "pcr_reverse", "pcr_polymerase", "library_pcr", "library"), SOURCE)
    return g.require_publishable()


SPEC = SimpleNamespace(TITLE=TITLE, NOTES=NOTES, SOURCE=SOURCE, SUMMARY=SUMMARY, CAVEAT=CAVEAT,
    FINAL_LIBRARY=FINAL_LIBRARY, SEQ_PRIMERS=SEQ_PRIMERS, READ_LENGTHS=READ_LENGTHS,
    FINAL_CAPTION=FINAL_CAPTION, SEQUENCING_INTRO=SEQUENCING_INTRO,
    SEQUENCING_ENDING=SEQUENCING_ENDING, SEQUENCING_UNAVAILABLE=SEQUENCING_UNAVAILABLE,
    oligos=oligos, protocol_graph=protocol_graph)

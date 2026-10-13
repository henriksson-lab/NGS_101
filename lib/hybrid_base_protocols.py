"""Reusable models for R-loop capture and differential cytosine readout protocols."""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

import seqprimers as sp
import illumina as il
from adapter_ligation import AdapterLigationGeometry
from base_conversion import ace_path, oxbs_path, tab_path
from batch_ngs import seg, truseq_library
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow
from endprep import (ThreePrimeOverhangAdapter, cohesive_ligation_geometry,
                     repair_and_dA_tail)
from pcr import PCRGeometry
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)


TRUSEQ_Y_END = ThreePrimeOverhangAdapter(
    "TruSeq Y adapter", il.STEM + "T", il.STEM_COMPLEMENT)


def _part(name: str, symbol: str = "X", n: int = 24, **kw) -> Segment:
    return Segment(name, symbol * n, placeholder=True, **kw)


@dataclass(frozen=True)
class RLoopCapture:
    """A three-stranded R-loop plus a capture route and sequenced component."""
    probe: str
    sequenced_component: str
    in_vivo: bool = False

    def __post_init__(self) -> None:
        if self.sequenced_component not in {"DNA", "RNA"}:
            raise ValueError("an R-loop capture must sequence either its DNA or RNA component")
        if not self.probe.strip():
            raise ValueError("an R-loop capture requires an affinity probe")

    def rloop(self) -> Scene:
        sc = Scene()
        sc.strand("template DNA", [_part("DNA in RNA:DNA hybrid")], label="R-loop")
        sc.anneal("hybrid RNA", [_part("RNA in RNA:DNA hybrid", "x")],
                  to="template DNA", pair=("RNA in RNA:DNA hybrid", "DNA in RNA:DNA hybrid"),
                  label="R-loop")
        sc.strand("displaced DNA", [_part("displaced single-stranded DNA", "Y")],
                  label="displaced DNA")
        sc.note("hybrid RNA", f"captured by {self.probe}")
        return sc

    def captured_component(self) -> Scene:
        sc = Scene()
        sc.strand("selected molecule", [_part(f"R-loop {self.sequenced_component} component")],
                  label=f"selected {self.sequenced_component}")
        sc.note("selected molecule", f"{self.probe}; " +
                ("captured in fixed cells" if self.in_vivo else "captured after gentle extraction"))
        return sc


def _workflow(initial: Scene, steps) -> Workflow:
    w = Workflow(MolecularState("Starting material", tuple(initial.rows())))
    for step in steps:
        action, scene, note, *geometry = step
        if geometry:
            w.ligate_adapter(action, geometry[0], scene.rows(), name=action, note=note)
        else:
            w.react(action, scene.rows(), name=action, note=note)
    return w


def _role_library(name: str, payload: tuple[str, ...]) -> Construct:
    return Construct([_part("Illumina-compatible left adapter", "L", 18),
                      *[_part(x, "X" if i % 2 == 0 else "Y", 26)
                        for i, x in enumerate(payload)],
                      _part("Illumina-compatible right adapter", "R", 18)], name=name)


def _conversion_scene(states: tuple[tuple[str, str], ...], label: str) -> Scene:
    sc = Scene()
    for i, (origin, state) in enumerate(states):
        sc.strand(origin, [Segment(state, state, "me" if "m" in state else None,
                                   placeholder=True, length_bp=1)], label=f"{origin} → {state}", at=0)
    return sc


def _conversion_workflow(key: str) -> Workflow:
    initial = _conversion_scene((("original C", "C"), ("original 5mC", "5mC"),
                                 ("original 5hmC", "5hmC")), "genomic cytosines")
    if key == "oxbs-seq":
        end_prep = repair_and_dA_tail(Construct([_part("unconverted genomic insert")]))
        truseq_ends = cohesive_ligation_geometry(
            end_prep, TRUSEQ_Y_END,
            statement="The protocol specifies TruSeq adapters: their public 3′-T end pairs with the dA-tailed insert.")
        tagged = _role_library_scene("adapter-ligated oxBS precursor", ("unconverted genomic insert",))
        oxidized = _conversion_scene((("original C", "C"), ("original 5mC", "5mC"),
                                      ("original 5hmC", "5fC")), "after oxidation")
        readout = _conversion_scene(tuple((b, oxbs_path(b).states[-1])
                                          for b in ("C", "5mC", "5hmC")), "PCR readout")
        return _workflow(initial, (
            ("Sonicate, end-repair and add one 3-prime dA to each end",
             _role_library_scene("dA-tailed oxBS insert", ("unconverted genomic insert",), ligated=False),
             "The dA-tailed insert is the substrate for the following adapter ligation."),
            ("Ligate TruSeq adapters", tagged,
             "The Nature Protocols workflow installs the library adapters before oxidation and bisulfite.",
             truseq_ends),
            ("Oxidize 5hmC selectively with potassium perruthenate", oxidized,
             "5hmC becomes 5fC; 5mC remains protected."),
            ("Bisulfite-convert, then amplify with uracil-tolerant polymerase", readout,
             "C and oxidized 5hmC read as T; only 5mC reads as C."),
        ))
    if key == "tab-seq":
        protected = _conversion_scene((("original C", "C"), ("original 5mC", "5mC"),
                                       ("original 5hmC", "5gmC")), "after β-GT")
        oxidized = _conversion_scene((("original C", "C"), ("original 5mC", "5caC"),
                                      ("original 5hmC", "5gmC")), "after mTet1")
        readout = _conversion_scene(tuple((b, tab_path(b).states[-1])
                                          for b in ("C", "5mC", "5hmC")), "PCR readout")
        return _workflow(initial, (
            ("Sonicate genomic DNA and add conversion-control molecules", initial,
             "C, 5mC and 5hmC spike-ins measure conversion, oxidation and protection."),
            ("Glucosylate 5hmC with T4 β-glucosyltransferase", protected,
             "5gmC is protected from the following TET oxidation."),
            ("Oxidize 5mC to 5caC with excess mTet1", oxidized,
             "The protected 5hmC branch remains 5gmC."),
            ("Prepare the Illumina WGBS library, bisulfite-convert and amplify", readout,
             "Only original 5hmC reads as C; C and original 5mC read as T.",
             AdapterLigationGeometry.unavailable(
                 "The defining TAB-seq source names its WGBS library preparation but does not print the adapter-end oligos.")),
        ))
    protected = _conversion_scene((("original C", "C"), ("original 5mC", "5mC"),
                                   ("original 5hmC", "5gmC")), "after β-GT")
    readout = _conversion_scene(tuple((b, ace_path(b).states[-1])
                                      for b in ("C", "5mC", "5hmC")), "PCR readout")
    return _workflow(initial, (
        ("Shear genomic DNA and add C/5mC/5hmC controls", initial,
         "The defining whole-genome experiment used approximately 300-bp input fragments."),
        ("Glucosylate 5hmC with T4 β-glucosyltransferase", protected,
         "Glucosylation strengthens the resistance of 5hmC to APOBEC3A."),
        ("Denature, then deaminate C and 5mC with APOBEC3A", readout,
         "C and 5mC read as T; protected 5hmC remains C without bisulfite."),
        ("Install Accel Methyl-NGS adapters and index by PCR",
         _role_library_scene("ACE-seq library", ("ACE-treated insert",), ligated=False),
         "Adapter construction follows deamination; the commercial oligo sequences are not public.",
         AdapterLigationGeometry.unavailable(
             "Accel Methyl-NGS adapter-end geometry is proprietary and not disclosed by the defining source.")),
    ))


def _conversion_graph(key: str, definition, final: Construct,
                      primers: tuple) -> ProtocolGraph:
    """Typed source-provenance graph for the three cytosine-conversion methods."""
    source = f"{definition.source}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    initial_state = _conversion_scene((
        ("original C", "C"), ("original 5mC", "5mC"),
        ("original 5hmC", "5hmC")), "genomic cytosines")
    dna = graph.input("genomic_dna", "genomic DNA", MaterialKind.BIOLOGICAL,
                      state=initial_state, primary=True)
    shear = graph.input("shearing_condition", "mechanical shearing condition",
                        MaterialKind.CONDITION)
    fragmented = graph.product(
        "fragmented_dna", "fragmented genomic DNA", MaterialKind.MOLECULAR,
        state=initial_state)
    graph.react("fragment", "Shear genomic DNA", StepKind.FRAGMENTATION,
                ((dna, InputRole.TEMPLATE), (shear, InputRole.CONDITION)), (fragmented,))

    if key == "oxbs-seq":
        repair = graph.input("end_repair_enzymes", "end-repair enzyme mix",
                             MaterialKind.CATALYST)
        repaired = graph.product(
            "end_repaired_dna", "end-repaired genomic fragments",
            MaterialKind.MOLECULAR,
            state=_role_library_scene("end-repaired oxBS insert",
                                      ("unconverted genomic insert",), ligated=False))
        graph.react("end_repair", "End-repair fragmented DNA", StepKind.END_REPAIR,
                    ((fragmented, InputRole.TEMPLATE), (repair, InputRole.CATALYST)),
                    (repaired,))
        tailing = graph.input("a_tailing_enzyme", "dA-tailing enzyme",
                              MaterialKind.CATALYST)
        da_scene = _role_library_scene(
            "dA-tailed oxBS insert", ("unconverted genomic insert",), ligated=False)
        tailed = graph.product("da_tailed_dna", "dA-tailed oxBS fragments",
                               MaterialKind.MOLECULAR, state=da_scene)
        graph.react("da_tail", "Add one 3-prime dA to each end", StepKind.TAILING,
                    ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)),
                    (tailed,))
        adapters = graph.input("truseq_adapters", "TruSeq Y adapters",
                               MaterialKind.ADAPTER)
        ligase = graph.input("dna_ligase", "DNA ligase", MaterialKind.CATALYST)
        ligated = graph.product(
            "ligated_precursor", "TruSeq-adapter-ligated oxBS precursor",
            MaterialKind.MOLECULAR,
            state=_role_library_scene("adapter-ligated oxBS precursor",
                                      ("unconverted genomic insert",)))
        end_prep = repair_and_dA_tail(Construct([_part("unconverted genomic insert")]))
        geometry = cohesive_ligation_geometry(
            end_prep, TRUSEQ_Y_END,
            statement="The protocol specifies TruSeq adapters: their public 3′-T end pairs with the dA-tailed insert.")
        graph.react("adapter_ligation", "Ligate TruSeq adapters", StepKind.LIGATION,
                    ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                     (ligase, InputRole.CATALYST)),
                    (ligated,), geometry=geometry)
        oxidant = graph.input("perruthenate", "potassium perruthenate",
                              MaterialKind.CONDITION)
        oxidized = graph.product(
            "oxidized_dna", "selectively oxidized adapter-ligated DNA",
            MaterialKind.MOLECULAR,
            state=_conversion_scene((
                ("original C", "C"), ("original 5mC", "5mC"),
                ("original 5hmC", "5fC")), "after oxidation"))
        graph.react("oxidize_5hmc", "Oxidize 5hmC selectively",
                    StepKind.CHEMICAL_CONVERSION,
                    ((ligated, InputRole.TEMPLATE), (oxidant, InputRole.CONDITION)),
                    (oxidized,))
        conversion_input = oxidized
        reviewed = [
            "genomic_dna", "shearing_condition", "fragment", "end_repair_enzymes",
            "end_repair", "a_tailing_enzyme", "da_tail", "truseq_adapters",
            "dna_ligase", "adapter_ligation", "perruthenate", "oxidize_5hmc",
        ]
    elif key == "tab-seq":
        controls = graph.input("conversion_controls", "C/5mC/5hmC control molecules",
                               MaterialKind.MOLECULAR)
        beta_gt = graph.input("beta_gt", "T4 beta-glucosyltransferase",
                              MaterialKind.CATALYST)
        udp_glucose = graph.input("udp_glucose", "UDP-glucose", MaterialKind.CONDITION)
        protected = graph.product(
            "glucosylated_dna", "5hmC-glucosylated genomic DNA",
            MaterialKind.MOLECULAR,
            state=_conversion_scene((
                ("original C", "C"), ("original 5mC", "5mC"),
                ("original 5hmC", "5gmC")), "after beta-GT"))
        graph.react("glucosylate_5hmc", "Glucosylate 5hmC", StepKind.CHEMICAL_CONVERSION,
                    ((fragmented, InputRole.TEMPLATE), (controls, InputRole.SPIKE_IN),
                     (beta_gt, InputRole.CATALYST),
                     (udp_glucose, InputRole.CONDITION)), (protected,))
        mtet = graph.input("mtet1", "mTet1 dioxygenase", MaterialKind.CATALYST)
        tet_conditions = graph.input("tet_cofactors", "mTet1 oxidation cofactors",
                                     MaterialKind.CONDITION)
        oxidized = graph.product(
            "tet_oxidized_dna", "mTet1-oxidized protected DNA",
            MaterialKind.MOLECULAR,
            state=_conversion_scene((
                ("original C", "C"), ("original 5mC", "5caC"),
                ("original 5hmC", "5gmC")), "after mTet1"))
        graph.react("oxidize_5mc", "Oxidize 5mC to 5caC", StepKind.CHEMICAL_CONVERSION,
                    ((protected, InputRole.TEMPLATE), (mtet, InputRole.CATALYST),
                     (tet_conditions, InputRole.CONDITION)), (oxidized,))
        conversion_input = oxidized
        reviewed = [
            "genomic_dna", "shearing_condition", "fragment", "conversion_controls",
            "beta_gt", "udp_glucose", "glucosylate_5hmc", "mtet1",
            "tet_cofactors", "oxidize_5mc",
        ]
    else:
        controls = graph.input("conversion_controls", "C/5mC/5hmC control molecules",
                               MaterialKind.MOLECULAR)
        beta_gt = graph.input("beta_gt", "T4 beta-glucosyltransferase",
                              MaterialKind.CATALYST)
        udp_glucose = graph.input("udp_glucose", "UDP-glucose", MaterialKind.CONDITION)
        protected = graph.product(
            "glucosylated_dna", "5hmC-glucosylated genomic DNA",
            MaterialKind.MOLECULAR,
            state=_conversion_scene((
                ("original C", "C"), ("original 5mC", "5mC"),
                ("original 5hmC", "5gmC")), "after beta-GT"))
        graph.react("glucosylate_5hmc", "Glucosylate 5hmC", StepKind.CHEMICAL_CONVERSION,
                    ((fragmented, InputRole.TEMPLATE), (controls, InputRole.SPIKE_IN),
                     (beta_gt, InputRole.CATALYST),
                     (udp_glucose, InputRole.CONDITION)), (protected,))
        denaturing = graph.input("denaturing_condition", "heat-denaturing condition",
                                 MaterialKind.CONDITION)
        apobec = graph.input("apobec3a", "APOBEC3A", MaterialKind.CATALYST)
        treated = graph.product(
            "ace_treated_dna", "APOBEC3A-deaminated protected DNA",
            MaterialKind.MOLECULAR,
            state=_conversion_scene(tuple(
                (base, ace_path(base).states[-1])
                for base in ("C", "5mC", "5hmC")), "ACE readout"))
        graph.react("apobec_deamination", "Denature and deaminate C and 5mC",
                    StepKind.CHEMICAL_CONVERSION,
                    ((protected, InputRole.TEMPLATE),
                     (denaturing, InputRole.CONDITION),
                     (apobec, InputRole.CATALYST)), (treated,))
        adapters = graph.input("accel_adapters", "Accel Methyl-NGS adapter reagents",
                               MaterialKind.ADAPTER, inferred=True)
        ligase = graph.input("adapter_installation_enzymes",
                             "Accel adapter-installation enzyme mix",
                             MaterialKind.CATALYST, inferred=True)
        ligated = graph.product(
            "adapter_tagged_dna", "Accel-adapter-tagged ACE-treated DNA",
            MaterialKind.MOLECULAR,
            state=_role_library_scene("ACE-seq adapter-tagged template",
                                      ("ACE-treated insert",), ligated=False),
            inferred=True)
        unavailable = AdapterLigationGeometry.unavailable(
            "Accel Methyl-NGS adapter-end geometry is proprietary and not disclosed by the defining source.")
        graph.react("adapter_installation", "Install Accel Methyl-NGS adapters",
                    StepKind.LIGATION,
                    ((treated, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                     (ligase, InputRole.CATALYST)), (ligated,), geometry=unavailable)
        conversion_input = ligated
        reviewed = [
            "genomic_dna", "shearing_condition", "fragment", "conversion_controls",
            "beta_gt", "udp_glucose", "glucosylate_5hmc", "denaturing_condition",
            "apobec3a", "apobec_deamination", "accel_adapters",
            "adapter_installation_enzymes", "adapter_installation",
        ]

    if key in {"oxbs-seq", "tab-seq"}:
        if key == "tab-seq":
            repair = graph.input("end_repair_enzymes", "WGBS end-repair enzyme mix",
                                 MaterialKind.CATALYST)
            repaired = graph.product("end_repaired_dna", "end-repaired TAB-seq DNA",
                                     MaterialKind.MOLECULAR, state=conversion_input.state)
            graph.react("end_repair", "End-repair treated DNA", StepKind.END_REPAIR,
                        ((conversion_input, InputRole.TEMPLATE),
                         (repair, InputRole.CATALYST)), (repaired,))
            tailing = graph.input("a_tailing_enzyme", "WGBS dA-tailing enzyme",
                                  MaterialKind.CATALYST)
            tailed = graph.product("da_tailed_dna", "dA-tailed TAB-seq DNA",
                                   MaterialKind.MOLECULAR, state=repaired.state)
            graph.react("da_tail", "Add one 3-prime dA to each end", StepKind.TAILING,
                        ((repaired, InputRole.TEMPLATE),
                         (tailing, InputRole.CATALYST)), (tailed,))
            adapters = graph.input("wgbs_adapters", "WGBS library adapters",
                                   MaterialKind.ADAPTER, inferred=True)
            ligase = graph.input("dna_ligase", "DNA ligase", MaterialKind.CATALYST)
            ligated = graph.product("ligated_precursor", "adapter-ligated TAB-seq DNA",
                                    MaterialKind.MOLECULAR, state=tailed.state,
                                    inferred=True)
            graph.react("adapter_ligation", "Ligate WGBS adapters", StepKind.LIGATION,
                        ((tailed, InputRole.TEMPLATE),
                         (adapters, InputRole.STRUCTURAL),
                         (ligase, InputRole.CATALYST)), (ligated,),
                        geometry=AdapterLigationGeometry.unavailable(
                            "The defining source names its WGBS preparation but does not print the adapter-end oligos."))
            conversion_input = ligated
            reviewed.extend(("end_repair_enzymes", "end_repair", "a_tailing_enzyme",
                             "da_tail", "wgbs_adapters", "dna_ligase",
                             "adapter_ligation"))
        bisulfite = graph.input("bisulfite_reagent", "bisulfite conversion reagent",
                                MaterialKind.CONDITION)
        converted = graph.product(
            "bisulfite_converted_dna", "bisulfite-converted library template",
            MaterialKind.MOLECULAR,
            state=_conversion_scene(tuple(
                (base, (oxbs_path(base) if key == "oxbs-seq" else tab_path(base)).states[-1])
                for base in ("C", "5mC", "5hmC")), "converted readout"))
        graph.react("bisulfite_conversion", "Bisulfite-convert the treated DNA",
                    StepKind.CHEMICAL_CONVERSION,
                    ((conversion_input, InputRole.TEMPLATE),
                     (bisulfite, InputRole.CONDITION)), (converted,))
        conversion_input = converted
        reviewed.extend(("bisulfite_reagent", "bisulfite_conversion"))

    forward = graph.input("library_forward_primer", "library PCR forward primer",
                          MaterialKind.OLIGO, inferred=not bool(primers))
    reverse = graph.input("library_reverse_primer", "indexed library PCR reverse primer",
                          MaterialKind.OLIGO, inferred=not bool(primers))
    polymerase = graph.input("library_polymerase", "library amplification polymerase",
                             MaterialKind.CATALYST)
    library = graph.product("sequencing_library", definition.title + " sequencing library",
                            MaterialKind.LIBRARY, state=final,
                            inferred=not bool(primers))
    graph.react("library_pcr", "Index and amplify the sequencing library", StepKind.PCR,
                ((conversion_input, InputRole.TEMPLATE),
                 (forward, InputRole.STRUCTURAL), (reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(
                    conversion_input.state if isinstance(conversion_input.state, Construct)
                    else final,
                    reason="The source identifies the library PCR but does not print one immutable selected primer pair."))
    plan = (SequencingPlan(final, tuple(primers), tuple(p.role for p in primers))
            if primers else SequencingPlan(limitation=definition.caveat or definition.sequencing))
    graph.sequence(library, plan)
    reviewed.extend(("library_forward_primer", "library_reverse_primer",
                     "library_polymerase", "library_pcr", "sequencing_library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _hybrid_workflow(key: str) -> Workflow:
    if key == "drip-seq":
        cap = RLoopCapture("S9.6 antibody", "DNA")
        return _workflow(cap.rloop(), (
            ("Gently extract nucleic acid and restriction-digest genomic DNA", cap.rloop(),
             "The five-enzyme cocktail preserves RNA:DNA hybrids while reducing fragment size."),
            ("Immunoprecipitate RNA:DNA hybrids with S9.6", cap.captured_component(),
             "An RNase H-treated aliquot is the hybrid-specificity control."),
            ("Sonicate captured DNA, end-repair and add one 3-prime dA to each end",
             _role_library_scene("dA-tailed captured DNA", ("captured DNA fragment",), ligated=False),
             "The dA-tailed captured fragment is shown before adapter ligation."),
            ("Ligate adapters to the dA-tailed captured DNA",
             _role_library_scene("non-strand-specific DRIP-seq library", ("captured DNA fragment",)),
             "The DNA component is sequenced; the original implementation used a standard dsDNA library.",
             AdapterLigationGeometry.unavailable(
                 "The detailed DRIP-seq source specifies a standard dsDNA library but does not print the adapter-end oligos.")),
        ))
    if key == "dripc-seq":
        cap = RLoopCapture("S9.6 antibody", "RNA")
        return _workflow(cap.rloop(), (
            ("Gently extract, restriction-digest and immunoprecipitate with S9.6", cap.rloop(),
             "The same hybrid capture precedes a different molecular readout."),
            ("Digest DNA and reverse-transcribe the retained R-loop RNA", cap.captured_component(),
             "Poly(dT) and random hexamers copy the RNA component into cDNA."),
            ("Make the second strand with dUTP, sonicate and ligate indexed adapters",
             _role_library_scene("dUTP-marked DRIPc-seq library", ("R-loop RNA-derived cDNA", "dUTP strand")),
             "UNG removes the dUTP-containing strand before PCR, preserving RNA strand identity.",
             AdapterLigationGeometry.unavailable(
                 "The DRIPc-seq source specifies indexed adapters but does not print their ligation-end geometry.")),
        ))
    cap = RLoopCapture("V5-tagged catalytically inactive RNASEH1 D210N", "DNA", in_vivo=True)
    primed = Scene(); primed.strand("copied molecule", [
        _part("N9-primer handle"), _part("N9 random primer", "N", 9),
        _part("R-loop DNA-derived copy")], label="unidirectional R-ChIP copy")
    primed.note("copied molecule", "the 5′ primer overhang blocks adapter ligation at this end")
    return _workflow(cap.rloop(), (
        ("Crosslink cells, shear chromatin and immunoprecipitate V5–RNASEH1", cap.captured_component(),
         "The D210N mutant binds R-loops; the binding-defective WKKD mutant is the control."),
        ("Denature and copy captured DNA with the tailed N9 primer", primed,
         "Random priming converts the DNA component to dsDNA and creates a blocking 5′ overhang."),
        ("Add one 3-prime dA to the unobstructed end", primed,
         "The tailed N9-primer end remains blocked; only the opposite end becomes ligation-ready."),
        ("Ligate one Illumina adapter and PCR-amplify",
         _rchip_library_scene(),
         "Only the unobstructed end accepts the adapter, so read orientation is retained.",
         AdapterLigationGeometry.unavailable(
             "The R-ChIP source prints the tailed N9 primer but identifies the opposite reagent only as an Illumina adapter.")),
    ))


def _hybrid_graph(key: str, definition, final: Construct) -> ProtocolGraph:
    """Typed graph for the DNA:RNA-hybrid capture protocols in this module."""
    source = f"{definition.source}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    if key == "drip-seq":
        cap = RLoopCapture("S9.6 antibody", "DNA")
    elif key == "dripc-seq":
        cap = RLoopCapture("S9.6 antibody", "RNA")
    else:
        cap = RLoopCapture("V5-tagged catalytically inactive RNASEH1 D210N",
                           "DNA", in_vivo=True)
    chromatin = graph.input("rloop_chromatin", "chromatin containing R-loops",
                            MaterialKind.BIOLOGICAL, state=cap.rloop(), primary=True)

    if key in {"drip-seq", "dripc-seq"}:
        extraction = graph.input("gentle_extraction", "gentle nucleic-acid extraction",
                                 MaterialKind.CONDITION)
        extracted = graph.product("extracted_nucleic_acid",
                                  "gently extracted R-loop-containing nucleic acid",
                                  MaterialKind.MOLECULAR, state=cap.rloop())
        graph.react("extract", "Gently extract nucleic acid", StepKind.EXTRACTION,
                    ((chromatin, InputRole.TEMPLATE),
                     (extraction, InputRole.CONDITION)), (extracted,))
        restriction = graph.input("restriction_cocktail",
                                  "R-loop-preserving restriction-enzyme cocktail",
                                  MaterialKind.CATALYST)
        digested = graph.product("restriction_fragments",
                                 "restriction-fragmented R-loop nucleic acid",
                                 MaterialKind.MOLECULAR, state=cap.rloop())
        graph.react("restriction_digest", "Restriction-digest genomic DNA",
                    StepKind.DIGESTION,
                    ((extracted, InputRole.TEMPLATE),
                     (restriction, InputRole.CATALYST)), (digested,))
        antibody = graph.input("s96_antibody", "S9.6 antibody",
                               MaterialKind.SELECTION)
        beads = graph.input("capture_beads", "antibody-capture beads",
                            MaterialKind.SELECTION)
        captured = graph.product("captured_rloops", "S9.6-captured R-loops",
                                 MaterialKind.MOLECULAR,
                                 state=cap.captured_component())
        rejected = graph.product("uncaptured_nucleic_acid", "uncaptured nucleic acid",
                                 MaterialKind.WASTE)
        graph.react("immunoprecipitate", "Immunoprecipitate RNA:DNA hybrids",
                    StepKind.AFFINITY_CAPTURE,
                    ((digested, InputRole.TEMPLATE),
                     (antibody, InputRole.SELECTION),
                     (beads, InputRole.SELECTION)), (captured, rejected))
        graph.discard(rejected, "not retained by S9.6 immunoprecipitation")
        reviewed = [
            "rloop_chromatin", "gentle_extraction", "extract",
            "restriction_cocktail", "restriction_digest", "s96_antibody",
            "capture_beads", "immunoprecipitate",
        ]
    else:
        formaldehyde = graph.input("formaldehyde", "formaldehyde cross-linking reagent",
                                   MaterialKind.CONDITION)
        fixed = graph.product("fixed_chromatin", "cross-linked R-loop chromatin",
                              MaterialKind.MOLECULAR, state=cap.rloop())
        graph.react("crosslink", "Cross-link cells", StepKind.CROSS_LINKING,
                    ((chromatin, InputRole.TEMPLATE),
                     (formaldehyde, InputRole.CONDITION)), (fixed,))
        sonication = graph.input("sonication", "chromatin sonication condition",
                                 MaterialKind.CONDITION)
        sheared = graph.product("sheared_chromatin", "sheared fixed chromatin",
                                MaterialKind.MOLECULAR, state=cap.rloop())
        graph.react("shear", "Shear cross-linked chromatin", StepKind.FRAGMENTATION,
                    ((fixed, InputRole.TEMPLATE),
                     (sonication, InputRole.CONDITION)), (sheared,))
        v5_antibody = graph.input("v5_antibody", "anti-V5 antibody",
                                  MaterialKind.SELECTION)
        inactive_rnaseh = graph.input("inactive_rnaseh1",
                                      "V5–RNASEH1 D210N R-loop binder",
                                      MaterialKind.SELECTION)
        captured = graph.product("captured_rloop_dna", "captured R-loop chromatin",
                                 MaterialKind.MOLECULAR,
                                 state=cap.captured_component())
        rejected = graph.product("uncaptured_chromatin", "uncaptured chromatin",
                                 MaterialKind.WASTE)
        graph.react("immunoprecipitate", "Immunoprecipitate V5–RNASEH1-bound R-loops",
                    StepKind.AFFINITY_CAPTURE,
                    ((sheared, InputRole.TEMPLATE),
                     (inactive_rnaseh, InputRole.SELECTION),
                     (v5_antibody, InputRole.SELECTION)), (captured, rejected))
        graph.discard(rejected, "not retained by anti-V5 immunoprecipitation")
        reviewed = [
            "rloop_chromatin", "formaldehyde", "crosslink", "sonication", "shear",
            "v5_antibody", "inactive_rnaseh1", "immunoprecipitate",
        ]

    if key == "drip-seq":
        sonication = graph.input("sonication", "post-capture sonication condition",
                                 MaterialKind.CONDITION)
        fragmented = graph.product("captured_dna_fragments", "sonicated captured DNA",
                                   MaterialKind.MOLECULAR,
                                   state=cap.captured_component())
        graph.react("fragment_captured_dna", "Sonicate captured DNA",
                    StepKind.FRAGMENTATION,
                    ((captured, InputRole.TEMPLATE),
                     (sonication, InputRole.CONDITION)), (fragmented,))
        repair = graph.input("end_repair_enzymes", "end-repair enzyme mix",
                             MaterialKind.CATALYST)
        repaired = graph.product("end_repaired_dna", "end-repaired captured DNA",
                                 MaterialKind.MOLECULAR, state=fragmented.state)
        graph.react("end_repair", "End-repair captured DNA", StepKind.END_REPAIR,
                    ((fragmented, InputRole.TEMPLATE),
                     (repair, InputRole.CATALYST)), (repaired,))
        tailing = graph.input("a_tailing_enzyme", "dA-tailing enzyme",
                              MaterialKind.CATALYST)
        tailed = graph.product(
            "da_tailed_dna", "dA-tailed captured DNA", MaterialKind.MOLECULAR,
            state=_role_library_scene("dA-tailed captured DNA",
                                      ("captured DNA fragment",), ligated=False))
        graph.react("da_tail", "Add one 3-prime dA to each end", StepKind.TAILING,
                    ((repaired, InputRole.TEMPLATE),
                     (tailing, InputRole.CATALYST)), (tailed,))
        adapters = graph.input("library_adapters", "standard dsDNA library adapters",
                               MaterialKind.ADAPTER, inferred=True)
        ligase = graph.input("dna_ligase", "DNA ligase", MaterialKind.CATALYST)
        tagged = graph.product(
            "adapter_ligated_dna", "adapter-ligated captured DNA",
            MaterialKind.MOLECULAR,
            state=_role_library_scene("non-strand-specific DRIP-seq library",
                                      ("captured DNA fragment",)), inferred=True)
        graph.react("adapter_ligation", "Ligate adapters", StepKind.LIGATION,
                    ((tailed, InputRole.TEMPLATE),
                     (adapters, InputRole.STRUCTURAL),
                     (ligase, InputRole.CATALYST)), (tagged,),
                    geometry=AdapterLigationGeometry.unavailable(
                        "The detailed source specifies a standard dsDNA library but does not print the adapter-end oligos."))
        pcr_input = tagged
        reviewed.extend((
            "sonication", "fragment_captured_dna", "end_repair_enzymes",
            "end_repair", "a_tailing_enzyme", "da_tail", "library_adapters",
            "dna_ligase", "adapter_ligation",
        ))
    elif key == "dripc-seq":
        dnase = graph.input("dnase", "DNase", MaterialKind.CATALYST)
        retained_rna = graph.product("retained_rloop_rna", "retained R-loop RNA",
                                     MaterialKind.MOLECULAR,
                                     state=cap.captured_component())
        dna_waste = graph.product("digested_dna", "digested DNA component",
                                  MaterialKind.WASTE)
        graph.react("digest_dna", "Digest the captured DNA component",
                    StepKind.DIGESTION,
                    ((captured, InputRole.TEMPLATE), (dnase, InputRole.CATALYST)),
                    (retained_rna, dna_waste))
        graph.discard(dna_waste, "DNase removes the DNA component before RNA copying")
        rt_primers = graph.input("rt_primer_mix", "poly(dT) and random-hexamer RT primers",
                                 MaterialKind.OLIGO)
        rtase = graph.input("reverse_transcriptase", "reverse transcriptase",
                            MaterialKind.CATALYST)
        cdna = graph.product("first_strand_cdna", "R-loop RNA-derived cDNA",
                             MaterialKind.MOLECULAR, state=cap.captured_component())
        graph.react("reverse_transcription", "Reverse-transcribe retained R-loop RNA",
                    StepKind.REVERSE_TRANSCRIPTION,
                    ((retained_rna, InputRole.TEMPLATE),
                     (rt_primers, InputRole.STRUCTURAL),
                     (rtase, InputRole.CATALYST)), (cdna,))
        second_primer = graph.input("second_strand_primer", "second-strand primer",
                                    MaterialKind.OLIGO)
        second_polymerase = graph.input("second_strand_polymerase",
                                        "dUTP second-strand polymerase",
                                        MaterialKind.CATALYST)
        dutp = graph.input("dutp", "dUTP-containing nucleotide mix",
                           MaterialKind.CONDITION)
        duplex = graph.product("dutp_cdna_duplex", "dUTP-marked cDNA duplex",
                               MaterialKind.MOLECULAR,
                               state=_role_library_scene("dUTP-marked cDNA",
                                                         ("R-loop RNA-derived cDNA",
                                                          "dUTP strand"), ligated=False))
        graph.react("second_strand", "Synthesize the dUTP-marked second strand",
                    StepKind.EXTENSION,
                    ((cdna, InputRole.TEMPLATE),
                     (second_primer, InputRole.STRUCTURAL),
                     (second_polymerase, InputRole.CATALYST),
                     (dutp, InputRole.CONDITION)), (duplex,))
        sonication = graph.input("sonication", "cDNA sonication condition",
                                 MaterialKind.CONDITION)
        fragmented = graph.product("fragmented_cdna", "fragmented dUTP-marked cDNA",
                                   MaterialKind.MOLECULAR, state=duplex.state)
        graph.react("fragment_cdna", "Sonicate the cDNA duplex", StepKind.FRAGMENTATION,
                    ((duplex, InputRole.TEMPLATE),
                     (sonication, InputRole.CONDITION)), (fragmented,))
        adapters = graph.input("indexed_adapters", "indexed library adapters",
                               MaterialKind.ADAPTER, inferred=True)
        ligase = graph.input("dna_ligase", "DNA ligase", MaterialKind.CATALYST)
        tagged = graph.product("adapter_ligated_cdna", "adapter-ligated cDNA",
                               MaterialKind.MOLECULAR, state=fragmented.state,
                               inferred=True)
        graph.react("adapter_ligation", "Ligate indexed adapters", StepKind.LIGATION,
                    ((fragmented, InputRole.TEMPLATE),
                     (adapters, InputRole.STRUCTURAL),
                     (ligase, InputRole.CATALYST)), (tagged,),
                    geometry=AdapterLigationGeometry.unavailable(
                        "The source specifies indexed adapters but does not print their ligation-end geometry."))
        ung = graph.input("ung", "uracil-DNA glycosylase", MaterialKind.CATALYST)
        retained = graph.product("retained_cdna_strand", "retained cDNA library strand",
                                 MaterialKind.MOLECULAR, state=tagged.state,
                                 inferred=True)
        uracil_waste = graph.product("removed_dutp_strand", "removed dUTP-marked strand",
                                     MaterialKind.WASTE)
        graph.react("remove_dutp_strand", "Remove the dUTP-containing strand",
                    StepKind.DIGESTION,
                    ((tagged, InputRole.TEMPLATE), (ung, InputRole.CATALYST)),
                    (retained, uracil_waste))
        graph.discard(uracil_waste, "UNG removes the dUTP-containing strand before PCR")
        pcr_input = retained
        reviewed.extend((
            "dnase", "digest_dna", "rt_primer_mix", "reverse_transcriptase",
            "reverse_transcription", "second_strand_primer",
            "second_strand_polymerase", "dutp", "second_strand", "sonication",
            "fragment_cdna", "indexed_adapters", "dna_ligase", "adapter_ligation",
            "ung", "remove_dutp_strand",
        ))
    else:
        denaturing = graph.input("denaturing_condition", "captured-DNA denaturation",
                                 MaterialKind.CONDITION)
        single = graph.product("single_stranded_captured_dna",
                               "single-stranded captured R-loop DNA",
                               MaterialKind.MOLECULAR,
                               state=cap.captured_component())
        partner = graph.product("released_complement", "released complementary strand",
                                MaterialKind.WASTE)
        graph.react("denature", "Denature captured R-loop DNA",
                    StepKind.STRAND_SEPARATION,
                    ((captured, InputRole.TEMPLATE),
                     (denaturing, InputRole.CONDITION)), (single, partner))
        graph.discard(partner, "the copied captured strand defines the directional library")
        n9 = graph.input("tailed_n9_primer", "tailed N9 random primer",
                         MaterialKind.OLIGO)
        polymerase = graph.input("copy_polymerase", "DNA polymerase",
                                 MaterialKind.CATALYST)
        primed_scene = Scene(); primed_scene.strand("copied molecule", [
            _part("N9-primer handle"), _part("N9 random primer", "N", 9),
            _part("R-loop DNA-derived copy")], label="unidirectional R-ChIP copy")
        primed = graph.product("n9_copied_dna", "tailed-N9 copied R-loop DNA",
                               MaterialKind.MOLECULAR, state=primed_scene)
        graph.react("random_primer_copy", "Copy captured DNA with the tailed N9 primer",
                    StepKind.EXTENSION,
                    ((single, InputRole.TEMPLATE), (n9, InputRole.STRUCTURAL),
                     (polymerase, InputRole.CATALYST)), (primed,))
        tailing = graph.input("a_tailing_enzyme", "dA-tailing enzyme",
                              MaterialKind.CATALYST)
        tailed = graph.product("one_end_da_tailed", "one-end dA-tailed copied DNA",
                               MaterialKind.MOLECULAR, state=primed_scene)
        graph.react("da_tail", "Add one 3-prime dA to the unobstructed end",
                    StepKind.TAILING,
                    ((primed, InputRole.TEMPLATE),
                     (tailing, InputRole.CATALYST)), (tailed,))
        adapter = graph.input("illumina_adapter", "one Illumina adapter",
                              MaterialKind.ADAPTER, inferred=True)
        ligase = graph.input("dna_ligase", "DNA ligase", MaterialKind.CATALYST)
        tagged = graph.product("one_end_ligated_dna", "one-end-adapter-ligated DNA",
                               MaterialKind.MOLECULAR, state=_rchip_library_scene(),
                               inferred=True)
        graph.react("adapter_ligation", "Ligate one Illumina adapter",
                    StepKind.LIGATION,
                    ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                     (ligase, InputRole.CATALYST)), (tagged,),
                    geometry=AdapterLigationGeometry.unavailable(
                        "The source prints the tailed N9 primer but identifies the opposite reagent only as an Illumina adapter."))
        pcr_input = tagged
        reviewed.extend((
            "denaturing_condition", "denature", "tailed_n9_primer", "copy_polymerase",
            "random_primer_copy", "a_tailing_enzyme", "da_tail",
            "illumina_adapter", "dna_ligase", "adapter_ligation",
        ))

    forward = graph.input("library_forward_primer", "library PCR forward primer",
                          MaterialKind.OLIGO, inferred=True)
    reverse = graph.input("library_reverse_primer", "library PCR reverse primer",
                          MaterialKind.OLIGO, inferred=True)
    pcr_polymerase = graph.input("pcr_polymerase", "library PCR polymerase",
                                 MaterialKind.CATALYST)
    library = graph.product("sequencing_library", definition.title + " sequencing library",
                            MaterialKind.LIBRARY, state=final, inferred=True)
    graph.react("library_pcr", "PCR-amplify the sequencing library", StepKind.PCR,
                ((pcr_input, InputRole.TEMPLATE),
                 (forward, InputRole.STRUCTURAL), (reverse, InputRole.STRUCTURAL),
                 (pcr_polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(
                    pcr_input.state if isinstance(pcr_input.state, Construct) else final,
                    reason="The authoritative source does not expose one exact selected library-primer pair."))
    graph.sequence(library, SequencingPlan(
        limitation=definition.caveat or definition.sequencing))
    reviewed.extend(("library_forward_primer", "library_reverse_primer",
                     "pcr_polymerase", "library_pcr", "sequencing_library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _role_library_scene(name: str, payload: tuple[str, ...], *, ligated: bool = True) -> Scene:
    lib = _role_library(name, payload)
    sc = Scene.duplex(list(lib), label=name)
    if ligated:
        sc.junction("top", "Illumina-compatible left adapter", payload[0], "adapter ligation")
        sc.junction("top", payload[-1], "Illumina-compatible right adapter", "adapter ligation")
    return sc


def _rchip_library() -> Construct:
    return Construct([
        Segment("tailed N9 handle", "CAAGCAGAAGACGGCATACGAG", "p7"),
        _part("N9 random-primer remnant", "N", 9),
        _part("R-loop DNA-derived insert"),
        _part("one ligated Illumina adapter", "R", 18),
    ], name="strand-specific R-ChIP library")


def _rchip_library_scene() -> Scene:
    sc = Scene.duplex(list(_rchip_library()), label="strand-specific R-ChIP library")
    sc.junction("top", "R-loop DNA-derived insert", "one ligated Illumina adapter",
                "adapter ligation")
    return sc


@dataclass(frozen=True)
class Definition:
    title: str
    doi: str
    source: str
    summary: str
    payload: tuple[str, ...]
    sequencing: str
    notes: str
    truseq: bool = False
    caveat: str = ""


DEFS = {
    "drip-seq": Definition("DRIP-seq — DNA:RNA immunoprecipitation sequencing",
        "10.1016/j.molcel.2012.01.017", "Ginno et al., Molecular Cell (2012)",
        "Gentle extraction and restriction digestion preserve R-loops for S9.6 immunoprecipitation; a conventional dsDNA library sequences the captured DNA component.",
        ("captured R-loop DNA",), "Single-end Illumina reads report captured restriction fragments without strand identity.", "01_drip_seq.html",
        caveat="The defining paper names the library reagents but does not expose complete adapter or sequencing-primer sequences; adapter regions are shown by role."),
    "dripc-seq": Definition("DRIPc-seq — strand-specific sequencing of R-loop RNA",
        "10.1016/j.molcel.2016.05.032", "Sanz et al., Molecular Cell (2016)",
        "After S9.6 capture, DNase releases the R-loop RNA for reverse transcription; dUTP second-strand synthesis and UNG make the final library strand-specific.",
        ("R-loop RNA-derived cDNA",), "Single-end reads derive from the retained cDNA strand and therefore identify the RNA strand in the original hybrid.", "01_dripc_seq.html",
        caveat="The paper and detailed protocol specify indexed Illumina construction but do not print complete adapter sequences; adapter regions are shown by role."),
    "r-chip": Definition("R-ChIP — in-vivo R-loop capture with inactive RNASEH1",
        "10.1016/j.molcel.2017.10.008", "Chen et al., Molecular Cell (2017)",
        "V5-tagged RNASEH1 D210N binds R-loops in fixed cells; a tailed N9 primer and one-ended adapter ligation produce a strand-specific DNA-component library.",
        ("tailed N9 end", "R-loop DNA-derived insert"), "Single-end Illumina reads preserve the orientation of the DNA strand captured from each R-loop.", "01_r_chip.html",
        caveat="The defining paper prints the tailed N9 primer, but identifies the opposite adapter only as an Illumina adapter; that adapter is shown by role."),
    "oxbs-seq": Definition("oxBS-seq — oxidative bisulfite sequencing",
        "10.1126/science.1220671", "Booth et al., Science (2012)",
        "Potassium perruthenate converts 5hmC to 5fC before bisulfite, so 5hmC and unmodified C read as T while 5mC alone reads as C.",
        ("oxidative-bisulfite-converted insert",), "TruSeq reads interrogate the converted insert; an i7 read identifies the library.", "01_oxbs_seq.html", truseq=True),
    "tab-seq": Definition("TAB-seq — TET-assisted bisulfite sequencing",
        "10.1016/j.cell.2012.04.027", "Yu et al., Cell (2012)",
        "β-GT protects 5hmC as 5gmC, then mTet1 oxidizes 5mC to 5caC; after bisulfite, only original 5hmC reads as C.",
        ("TAB-treated insert",), "TruSeq reads interrogate the converted insert; an i7 read identifies the library.", "01_tab_seq.html", truseq=True),
    "ace-seq": Definition("ACE-seq — APOBEC-coupled epigenetic sequencing",
        "10.1038/nbt.4204", "Schutsky et al., Nature Biotechnology (2018)",
        "β-GT protects 5hmC before APOBEC3A deaminates C and 5mC, giving a bisulfite-free readout in which only original 5hmC remains C.",
        ("ACE-treated insert",), "The defining experiment used single-end NextSeq sequencing of Accel Methyl-NGS libraries.", "01_ace_seq.html",
        caveat="Accel Methyl-NGS adapter oligos are proprietary and are shown by role; no exact sequencing-primer sequence is asserted."),
}


def get(key: str):
    d = DEFS[key]
    if d.truseq:
        final, primers = truseq_library([seg(d.payload[0], "X" * 34, placeholder=True)],
                                        d.title + " final library", dual_index=False)
    else:
        final, primers = ((_rchip_library() if key == "r-chip" else
                           _role_library(d.title + " final library", d.payload)), ())
    spec = SimpleNamespace(
        TITLE=d.title, NOTES=d.notes,
        SOURCE=(f'Defining source: <a href="https://doi.org/{d.doi}">{d.source}</a>.' +
                (' Detailed chemistry was checked against the authors’ later protocol.'
                 if key in {"drip-seq", "dripc-seq", "oxbs-seq", "tab-seq", "r-chip"} else '')),
        SUMMARY=d.summary, CAVEAT=d.caveat, FINAL_LIBRARY=final, SEQ_PRIMERS=primers,
        FINAL_CAPTION=d.sequencing, SEQUENCING_INTRO=d.sequencing, READ_LENGTHS=None,
        SEQUENCING_ENDING=d.sequencing,
        SEQUENCING_UNAVAILABLE=("The authoritative source does not expose the exact sequencing-primer sequence for this library implementation."
                                if not primers else ""),
        workflow=lambda: (_conversion_workflow(key) if key in {"oxbs-seq", "tab-seq", "ace-seq"}
                          else _hybrid_workflow(key)),
        oligos=(lambda: ('5′-/invddt/CAAGCAGAAGACGGCATACGAGNNNNNNNNN-3′ — tailed N9 random primer',)
                if key == "r-chip" else ()),
    )
    if key in {"oxbs-seq", "tab-seq", "ace-seq"}:
        spec.protocol_graph = lambda: _conversion_graph(key, d, final, primers)
    elif key in {"drip-seq", "dripc-seq", "r-chip"}:
        spec.protocol_graph = lambda: _hybrid_graph(key, d, final)
    return spec

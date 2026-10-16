"""Source-owned CUT&RUN2017 soluble-yeast branch; unresolved library ends stay unassembled."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct
from pcr import PCRGeometry
from protocol_graph import (Disposition, Highlight, InputRole, MaterialKind,
                            ProtocolGraph, SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind, Terminal)
from protocol_model import role

TITLE="CUT&RUN — antibody-tethered MNase chromatin profiling"
NOTES="01_cut-run.html"
SOURCE='Defining source: <a href="https://doi.org/10.7554/eLife.21856">Skene &amp; Henikoff, <i>eLife</i> (2017)</a>, Methods and Appendix 2; library-preparation paragraph and its cited predecessors.'
SUMMARY=("An antibody tethers Protein A–MNase to native chromatin. This view selects the "
         "yeast soluble-fragment branch: stop digestion, release fragments, physically separate "
         "the supernatant, then extract DNA and prepare a sequencing library. The pellet has an optional QC route.")
CAVEAT=("The Kasinathan 2014 predecessor's Supplementary Fig. 1c names Illumina TruSeq adapters; "
        "that is predecessor context, not a selected 2017 adapter, index width/member or ordered primer. "
        "The CUT&RUN source names a KAPA library kit, 8 or more PCR cycles and 25-base paired reads on HiSeq 2500. "
        "The actual 2017 adapter/cap/index assembly and sequencing-primer sites remain unassigned "
        "in this declaration. Bracketed roles do not represent inferred exact bases. "
        "The mammalian, total-DNA and high-salt variants are distinguished in source annotations.")
PAIRED_READ_CYCLES=25
LIBRARY_MIN_PCR_CYCLES=8
READ_LENGTHS={"Read 1":PAIRED_READ_CYCLES,"Read 2":PAIRED_READ_CYCLES}
ACQUIRED_ROLES=tuple(READ_LENGTHS)
SEQ_PRIMERS=()
SEQUENCING_ENDING="The source specifies 25 cycles paired-end on HiSeq 2500. Exact molecular primer sites and index-read acquisition are unassigned; no single-index eight-base shell is supplied by a shared recipe."
SEQUENCING_UNAVAILABLE="The selected 2017 adapter, physical oligos and complete ordered library are not assigned by this model. This is not accepted publication absence."
READ_SPECS=(sp.read_spec("Read 1",25,"released DNA from one end; exact primer boundary unassigned","include"),
            sp.read_spec("Read 2",25,"released DNA from the other end; exact primer boundary unassigned","include"),
            sp.read_spec("Index 1 (i7)",None,"acquisition, width, location and member unassigned","undetermined"),
            sp.read_spec("Index 2 (i5)",None,"acquisition, width, location and member unassigned","undetermined"))


def state(name: str, *parts: str) -> Construct:
    return Construct([role(part) for part in parts],name=name)


RELEASED_DNA=state("released CUT&RUN DNA","target-associated fragments, soluble background and optional spike-in DNA")
LIGATED_LIBRARY=state("adapter-ligated CUT&RUN DNA","unassembled ligation adapter region",
                      "released CUT&RUN fragment","unassembled opposite ligation adapter region")
FINAL_LIBRARY=state("CUT&RUN library — ordered ends unassembled","unassigned adapter/cap/index region",
                    "released CUT&RUN fragment","unassigned opposite adapter/cap/index region")
FINAL_CAPTION="Source-declared library roles only. No exact TruSeq shell, fixed insert length, eight-base index or physical primer sequence is asserted."


def protocol_graph() -> ProtocolGraph:
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    nuclei=graph.input("nuclei","native yeast nuclei (selected bead-based branch)",MaterialKind.BIOLOGICAL,
                        primary=True,state=state("native yeast nuclei","chromatin in permeabilized nuclei"))
    beads=graph.input("cona_beads","concanavalin-A magnetic beads",MaterialKind.SELECTION)
    antibody=graph.input("primary_antibody","target-specific primary antibody",MaterialKind.SELECTION)
    secondary=graph.input("secondary_antibody","secondary IgG for the selected mouse-antibody example",MaterialKind.SELECTION)
    mnase=graph.input("pa_mnase","Protein A–MNase",MaterialKind.CATALYST)
    tethered=graph.product("tethered_nuclei","washed antibody/pA–MNase nuclei",MaterialKind.MOLECULAR)
    wash=graph.product("unbound_targeting_reagents","unbound targeting reagents removed by washes",MaterialKind.WASTE)
    graph.react("target_chromatin","Immobilize nuclei, bind antibodies and pA–MNase, then wash",StepKind.AFFINITY_CAPTURE,
                ((nuclei,InputRole.TEMPLATE),(beads,InputRole.SELECTION),(antibody,InputRole.SELECTION),
                 (secondary,InputRole.SELECTION),(mnase,InputRole.CATALYST)),(tethered,wash))
    graph.discard(wash,"The selected targeted branch washes before activation; the free-MNase control is a separate source branch.")
    calcium=graph.input("calcium","calcium activation on ice",MaterialKind.CONDITION)
    cut=graph.product("cleaved_nuclei","digested nuclei before stopping/release",MaterialKind.MOLECULAR)
    graph.react("targeted_cleavage","Activate antibody-tethered MNase",StepKind.DIGESTION,
                ((tethered,InputRole.TEMPLATE),(mnase,InputRole.CATALYST),(calcium,InputRole.CONDITION)),(cut,))
    stop=graph.input("stop_buffer","EDTA/EGTA STOP buffer",MaterialKind.CONDITION)
    rnase=graph.input("rnase_a","RNase A in the selected yeast STOP workflow",MaterialKind.CATALYST)
    spike=graph.input("spike_in_dna","optional heterologous DNA spike-in (shown)",MaterialKind.MOLECULAR,
                       state=state("optional spike-in","heterologous DNA fragments"))
    stopped=graph.product("stopped_nuclei","quenched digest with optional DNA spike-in",MaterialKind.MOLECULAR)
    graph.react("stop_digest","Stop MNase with chelators; optional heterologous DNA is added here",StepKind.QUENCHING,
                ((cut,InputRole.TEMPLATE),(stop,InputRole.CONDITION),(rnase,InputRole.CATALYST),
                 (spike,InputRole.SPIKE_IN)),(stopped,))
    release=graph.input("release_condition","37°C release in the yeast soluble branch",MaterialKind.CONDITION)
    slurry=graph.product("released_slurry","released fragments plus residual nuclei",MaterialKind.MOLECULAR)
    graph.react("release","Allow excised fragments to enter the soluble phase",StepKind.THERMAL_TREATMENT,
                ((stopped,InputRole.TEMPLATE),(release,InputRole.CONDITION)),(slurry,),changes_molecule=False)
    separator=graph.input("fractionation","centrifugation and magnetic separation",MaterialKind.SELECTION)
    soluble=graph.product("soluble_fraction","recovered supernatant",MaterialKind.MOLECULAR,state=RELEASED_DNA)
    pellet=graph.product("residual_pellet","residual bead/nuclear pellet",MaterialKind.MOLECULAR)
    graph.react("fractionate","Recover supernatant; retain the separate pellet for optional QC",StepKind.SELECTION,
                ((slurry,InputRole.TEMPLATE),(separator,InputRole.SELECTION)),(soluble,pellet),changes_molecule=False)
    pk=graph.input("proteinase_k","Proteinase K",MaterialKind.CATALYST)
    extraction=graph.input("extraction_conditions","SDS, phenol/chloroform and ethanol recovery",MaterialKind.CONDITION)
    recovered=graph.product("recovered_dna","deproteinized supernatant DNA",MaterialKind.MOLECULAR,state=RELEASED_DNA)
    graph.react("recover_dna","Deproteinize and extract supernatant DNA for library preparation",StepKind.DEPROTEINIZATION,
                ((soluble,InputRole.TEMPLATE),(pk,InputRole.CATALYST),(extraction,InputRole.CONDITION)),(recovered,))
    qc=graph.product("pellet_qc","optional extracted pellet DNA for gel analysis",MaterialKind.QC)
    graph.react("pellet_qc_extraction","Optional — extract pellet DNA for gel QC",StepKind.DEPROTEINIZATION,
                ((pellet,InputRole.TEMPLATE),(pk,InputRole.CATALYST),(extraction,InputRole.CONDITION)),(qc,))
    graph.mark_terminal(Terminal(qc.id,Disposition.QC,"Appendix 2 explicitly offers pellet-fraction gel analysis; this optional branch is not mandatory physical discard."))
    repair=graph.input("end_repair_enzymes","source-referenced library end-repair enzymes",MaterialKind.CATALYST)
    repaired=graph.product("end_repaired_dna","end-repaired DNA; exact end chemistry unassigned",MaterialKind.MOLECULAR)
    graph.react("end_repair","Perform source-referenced library end repair",StepKind.END_REPAIR,
                ((recovered,InputRole.TEMPLATE),(repair,InputRole.CATALYST)),(repaired,))
    tailing=graph.input("a_tailing_enzyme","source-referenced library A-addition enzyme",MaterialKind.CATALYST)
    tailed=graph.product("da_tailed_dna","A-addition product; ligation substrate unassembled",MaterialKind.MOLECULAR)
    graph.react("a_tail","Perform source-referenced A-addition",StepKind.TAILING,
                ((repaired,InputRole.TEMPLATE),(tailing,InputRole.CATALYST)),(tailed,))
    adapters=graph.input("library_adapters","2017 library adapters — ordered strands unassigned",MaterialKind.ADAPTER)
    ligase=graph.input("dna_ligase","library DNA ligase",MaterialKind.CATALYST)
    ligated=graph.product("adapter_ligated_dna",LIGATED_LIBRARY.name,MaterialKind.MOLECULAR,state=LIGATED_LIBRARY)
    graph.react("adapter_ligation","Ligate library adapters; exact ends remain unassembled",StepKind.LIGATION,
                ((tailed,InputRole.TEMPLATE),(adapters,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(ligated,),
                geometry=AdapterLigationGeometry.unavailable("The declaration has not assigned the actual 2017 ordered adapters or ligation-end geometry. Predecessor TruSeq labeling is not proof of the selected 2017 adapter; no publication absence is accepted."))
    forward=graph.input("library_forward_primer","forward library-PCR reagent — ordered strand unassigned",MaterialKind.OLIGO)
    reverse=graph.input("library_reverse_primer","reverse library-PCR reagent — ordered strand unassigned",MaterialKind.OLIGO)
    polymerase=graph.input("library_polymerase","KAPA DNA polymerase library preparation kit",MaterialKind.CATALYST)
    library=graph.product("cut_run_library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    graph.react("library_pcr",f"Amplify the soluble-branch library ({LIBRARY_MIN_PCR_CYCLES} or more cycles; no size-selection)",StepKind.PCR,
                ((ligated,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),(reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.unavailable(LIGATED_LIBRARY,
                    reason="The current declaration does not assign the 2017 ordered library-PCR strands or their binding geometry. This is a model/source-assignment gap, not a claim of unpublished oligos."))
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=(),acquired_roles=ACQUIRED_ROLES,
                   read_lengths=READ_LENGTHS,limitation=SEQUENCING_UNAVAILABLE))
    graph.annotate(Highlight("fractionate","Other extraction branches are distinct",
        "Appendix 3 uses a separate K562 protocol with 4°C leaching or whole-reaction DNA extraction. "
        "Its Option B and the yeast total-DNA route use 0.5× AMPure, retain the SUPERNATANT and exclude bead-bound large DNA. "
        "Those physical-selection variants are not applied to this selected soluble yeast route."))
    graph.annotate(Highlight("library_pcr","Amplification preference is not physical discard",
        "Methods line 376 names 14 PCR cycles for total-DNA samples to favor short-fragment amplification. "
        "Appendix 2 also notes some soluble high-MW DNA is not appreciably amplified. "
        "Neither statement proves physical removal of every long or background fragment; no analysis length bin is modeled as a discarded molecule."))
    graph.annotate(SourceDiscrepancy("fractionate", "Deposited total-DNA protocol label and body differ",
        (SourceClaim("eLife v3 Appendix 2, total-DNA paragraph", "Extract the total digestion reaction before negative AMPure selection."),
         SourceClaim("GSE84474_Cut-and-Run_total_DNA_protocol.rtf, complete body", "The file labeled total DNA instead instructs centrifugation and recovery of the supernatant before extraction and negative AMPure selection.")),
        "This page selects the article's soluble yeast branch; the deposited total-DNA file is not used to define whole-reaction extraction.",
        "Original article Appendix 2 and deposited GEO protocol, independently read.",
        "The filename/body difference remains unresolved; no typographical error is assumed."))
    graph.declare_reviewed_scope((
        "nuclei","cona_beads","primary_antibody","secondary_antibody","pa_mnase",
        "target_chromatin","calcium","targeted_cleavage","stop_buffer","rnase_a",
        "spike_in_dna","stop_digest","release_condition","release","fractionation",
        "fractionate","proteinase_k","extraction_conditions","recover_dna","pellet_qc_extraction",
        "end_repair_enzymes","end_repair","a_tailing_enzyme","a_tail","library_adapters",
        "dna_ligase","adapter_ligation","library_forward_primer","library_reverse_primer",
        "library_polymerase","library_pcr","cut_run_library"),SOURCE)
    return graph.require_publishable()

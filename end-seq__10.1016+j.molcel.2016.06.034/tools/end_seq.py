"""Canela et al. END-seq model."""
from copy import deepcopy
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, feature
from end_capture import UserHairpin
from protocol_model import molecule, role
import illumina as il
from pcr import PCRGeometry
from protocol_graph import (Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from topology import HairpinTopology
import seqprimers as sp

TITLE = "END-seq — direct capture of DNA double-strand-break ends"
NOTES = "01_end-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1016/j.molcel.2016.06.034">Canela et al., <i>Molecular Cell</i> (2016)</a>.'
SUMMARY = "Native DNA ends are blunted and trapped in agarose with a biotinylated P5 hairpin; after shearing and enrichment, a P7 hairpin captures the opposite library end and USER opens both loops for indexed PCR."
PCR_CYCLES = 16
SHOW_READ_CONTENT = False
READ_SPECS = (sp.read_spec("Read 1", "50 (HiSeq2500) or75 (NextSeq500)", "processed break-derived genomic end; source alternative runs", "include"),)
CAVEAT = "Source-retained USER arms are modelled directly. Actual outer PCR ends and index members/widths are unassigned; the source specifies single-end sequencing, not acquired dual-index or paired-end reads. Blunting deliberately erases the original overhang chemistry. Read 1 starts at the first base retained at the processed break end."
A1 = UserHairpin("ENDseq adaptor 1", "GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGUUTUTUUACACTCTTTCCCTACACGACGCTCTTCCGATCT", (33,35))
A2 = UserHairpin("ENDseq adaptor 2", "GATCGGAAGAGCACACGTCUUUUUUUUAGACGTGTGCTCTTCCGATCT")
A1_STEM_NT, A1_LOOP_NT = 31, 8
A2_STEM_NT, A2_LOOP_NT = 19, 9
A1_TOPOLOGY = HairpinTopology(A1.name, A1.sequence, 0,
                              A1_STEM_NT + A1_LOOP_NT, A1_STEM_NT,
                              colour="var(--c-r1)")
A2_TOPOLOGY = HairpinTopology(A2.name, A2.sequence, 0,
                              A2_STEM_NT + A2_LOOP_NT, A2_STEM_NT,
                              colour="var(--c-r2)")
# The retained A1 arm is exactly the source-compatible Read1 priming site.
# Other PCR-added ends/index architecture remain unassigned, not a generic shell.
A1_RETAINED = A1.sequence.rsplit("U", 1)[1]
A2_RETAINED = A2.sequence.split("U", 1)[0]
if A1_RETAINED != il.TRUSEQ_READ1:
    raise ValueError("source A1 retained arm does not match the declared Read1 reference")
FINAL_LIBRARY = Construct([
    role("PCR-added P5-side outer end; physical primer unassigned"),
    Segment("source USER-retained adaptor1 Read1 arm", A1_RETAINED, "r1"),
    role("processed-break-to-sonication genomic insert; variable length"),
    Segment("distal dA junction", "A"),
    Segment("source USER-retained adaptor2 arm", A2_RETAINED, "r2"),
    role("index-bearing PCR end; architecture/member/width unassigned", tag="cbc",
         molecular_feature=feature("sample_index", "sample_index", "unknown")),
], name="END-seq source-retained-arm endpoint; outer PCR ends unassigned")
SEQ_PRIMERS = (sp.TRUSEQ["R1"],)
if sp.verify(FINAL_LIBRARY, SEQ_PRIMERS, required_roles=("Read 1",)):
    raise ValueError("source-retained END-seq Read1 site does not land")
FINAL_CAPTION = "Source-retained USER arms flank the processed DSB-to-sonication insert. PCR-added outer ends and index architecture remain unassigned; no canonical dual-index shell is supplied."
SEQUENCING_INTRO = "The source reports50-base HiSeq2500 or75-base NextSeq500 single-end runs. A compatible Read1 reference lands in the source-retained adaptor1 arm; exact physical sequencing reagent and acquired index cycles remain unassigned."
SEQUENCING_ENDING = ""
SEQUENCING_UNAVAILABLE = "Index architecture/width and acquired index cycles remain unassigned; alternative biological read lengths are source-reported, not merged into one run."


def retained_state(label, *names):
    return molecule(label, tuple(role(n) for n in names))


def hairpin_state(adapter):
    scene = Scene()
    scene.strand("hairpin", [Segment(adapter.name, adapter.sequence, "r1",
                 note="U is deoxyuridine; biotin-dT and terminal phosphorothioate are source-written modifications")],
                 label=adapter.name + " (one folded oligo)", mod5="P")
    scene.note("hairpin", "Source terminal C*T phosphorothioate; " +
               ("two internal biotin-dT residues at source loop positions" if adapter.biotin_dt else "no biotin-dT assigned"))
    return scene


def opened_template():
    scene = molecule("actual USER-opened END-seq precursor", (
        Segment("source retained adaptor1 arm", A1_RETAINED, "r1"),
        role("processed-break-to-sonication genomic DNA"),
        Segment("distal dA junction", "A"),
        Segment("source retained adaptor2 arm", A2_RETAINED, "r2")))
    scene.note("molecule", "Representative retained strand path; complementary ends and physical PCR primer identities unassigned.")
    return scene


def oligos():
    """Source-written phosphate, biotin-dT, U and terminal phosphorothioate retained."""
    rows=[]
    for adapter in (A1,A2):
        written="".join("[Biotin-dT]" if i in adapter.biotin_dt else b
                        for i,b in enumerate(adapter.sequence))
        written=written[:-1]+"*"+written[-1]
        rows.append(adapter.name+" — 5′-Phos-"+written+"-3′; U denotes USER-sensitive deoxyuridine")
    return tuple(rows)

def prepared_end_state(label, region):
    scene = molecule(label, (role(region), Segment("source DNA3prime dA", "A")))
    scene.overhang("molecule", "source DNA3prime dA", polarity="3-prime", owner="substrate", text="source-described3′A")
    scene.note("molecule", "Representative prepared end of a source-described duplex; opposite strand/base geometry unassigned.")
    return scene


def hairpin_ligation_geometry(adapter, substrate_state):
    scene = deepcopy(substrate_state)
    scene.strand("hairpin", [Segment("USER hairpin body", adapter.sequence[:-1], "r1"),
                             Segment("hairpin3primeT", adapter.sequence[-1], "r1")], label=adapter.name, mod5="P")
    scene.overhang("hairpin", "hairpin3primeT", polarity="3-prime", owner="adapter", text="source3′T")
    return AdapterLigationGeometry.known("cohesive", scene,
        statement="The source-printed adapter3′T pairs with the actual source-described dA-prepared duplex end; full capped topology remains unassigned.")


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    dna=graph.input("genomic_dna","agarose-embedded genomic DNA with native DSBs",
                    MaterialKind.MOLECULAR,primary=True)
    deproteinize=graph.input("deproteinization_condition","ProteinaseK in agarose,50°C1h then37°C7h;wash,RNaseA,wash",
                             MaterialKind.CONDITION)
    clean=graph.product("clean","deproteinized genomic DNA",MaterialKind.MOLECULAR)
    graph.react("deproteinize","Deproteinize genomic DNA in agarose",
                StepKind.DEPROTEINIZATION,
                ((dna,InputRole.TEMPLATE),(deproteinize,InputRole.CONDITION)),(clean,))
    repair=graph.input("repair_enzymes","source in-plug DSB-end blunting enzymes; identities deferred to Figure1A",MaterialKind.CATALYST)
    blunt=graph.product("blunted","blunted native DSB ends",MaterialKind.MOLECULAR)
    graph.react("blunt","Blunt native DSB ends",StepKind.END_REPAIR,
                ((clean,InputRole.TEMPLATE),(repair,InputRole.CATALYST)),(blunt,))
    tailing=graph.input("tailing_polymerase","in-plug dA-tailing polymerase; exact Figure1A assignment deferred",MaterialKind.CATALYST)
    tailed1=graph.product("tailed1","dA-tailed captured DSB ends",MaterialKind.MOLECULAR, state=prepared_end_state("dA-tailed native break", "blunted native genomic break"))
    graph.react("tail_break","Add 3-prime dA at processed break ends",StepKind.TAILING,
                ((blunt,InputRole.TEMPLATE),(tailing,InputRole.CATALYST)),(tailed1,))
    a1=graph.input("a1","biotinylated P5 USER hairpin adaptor 1",MaterialKind.ADAPTER, state=hairpin_state(A1), chemical_state=A1_TOPOLOGY)
    ligase1=graph.input("ligase1","P5 hairpin ligase",MaterialKind.CATALYST)
    tagged=graph.product("tagged","P5-hairpin-tagged DSB ends",MaterialKind.MOLECULAR,
                         state=retained_state("adaptor1-capped processed genomic break", "source adaptor1 hairpin cap", "break-adjacent genomic DNA"))
    graph.react("ligate_a1","Ligate biotinylated P5 hairpin adaptor 1",
                StepKind.LIGATION,
                ((tailed1,InputRole.TEMPLATE),(a1,InputRole.STRUCTURAL),
                 (ligase1,InputRole.CATALYST)),(tagged,),geometry=hairpin_ligation_geometry(A1,tailed1.state))
    gelase=graph.input("gelase","GELase for agarose digestion",MaterialKind.CATALYST)
    melting=graph.input("agarose_melting","agarose melting condition",MaterialKind.CONDITION)
    recovered=graph.product("recovered","DNA recovered from the agarose plug",MaterialKind.MOLECULAR)
    graph.react("recover_plug","Melt and digest agarose to recover tagged DNA",StepKind.EXTRACTION,
                ((tagged,InputRole.TEMPLATE),(gelase,InputRole.CATALYST),(melting,InputRole.CONDITION)),(recovered,))
    shear=graph.input("shear_condition","Covaris sonication to150–200bp",MaterialKind.CONDITION)
    fragments=graph.product("fragments","sheared P5-hairpin-tagged fragments",
                            MaterialKind.MOLECULAR)
    graph.react("shear","Shear distal genomic DNA",StepKind.FRAGMENTATION,
                ((recovered,InputRole.TEMPLATE),(shear,InputRole.CONDITION)),(fragments,))
    streptavidin=graph.input("streptavidin","MyOneC1 streptavidin beads (Invitrogen)",
                            MaterialKind.SELECTION)
    captured=graph.product("captured","captured P5-hairpin DSB fragments",
                           MaterialKind.MOLECULAR)
    untagged=graph.product("untagged","fragments lacking the biotinylated break end",
                           MaterialKind.WASTE)
    graph.react("capture","Enrich biotinylated break-end fragments",
                StepKind.AFFINITY_CAPTURE,
                ((fragments,InputRole.TEMPLATE),(streptavidin,InputRole.SELECTION)),
                (captured,untagged))
    graph.discard(untagged,"fragments lacking adaptor 1 are washed away")
    distal_repair=graph.input("distal_repair_enzymes","T4DNApolymerase15U,Klenow5U,T4PNK15U",MaterialKind.CATALYST)
    repaired2=graph.product("repaired2","distal-end-repaired captured fragments",
                            MaterialKind.MOLECULAR)
    graph.react("repair_distal","Repair the sonicated distal ends",StepKind.END_REPAIR,
                ((captured,InputRole.TEMPLATE),(distal_repair,InputRole.CATALYST)),(repaired2,))
    distal_tail=graph.input("distal_tailing_enzyme","Klenow exo-minus fragment15U",MaterialKind.CATALYST)
    tailed2=graph.product("tailed2","distal-dA-tailed captured fragments",
                          MaterialKind.MOLECULAR, state=prepared_end_state("distal dA-tailed captured DNA", "adaptor1-capped break-to-sonication DNA"))
    graph.react("tail_distal","Add 3-prime dA at distal fragment ends",StepKind.TAILING,
                ((repaired2,InputRole.TEMPLATE),(distal_tail,InputRole.CATALYST)),(tailed2,))
    a2=graph.input("a2","P7 USER hairpin adaptor 2",MaterialKind.ADAPTER, state=hairpin_state(A2), chemical_state=A2_TOPOLOGY)
    ligase2=graph.input("ligase2","NEB Quick Ligation kit for adaptor2",MaterialKind.CATALYST)
    double_hairpin=graph.product("double_hairpin","two-hairpin captured fragment",
                                 MaterialKind.MOLECULAR, state=retained_state("two un-opened hairpin caps", "adaptor1 hairpin", "processed-break-to-sonication DNA", "adaptor2 hairpin"))
    graph.react("ligate_a2","Ligate P7 hairpin adaptor 2",StepKind.LIGATION,
                ((tailed2,InputRole.TEMPLATE),(a2,InputRole.STRUCTURAL),
                 (ligase2,InputRole.CATALYST)),(double_hairpin,),geometry=hairpin_ligation_geometry(A2,tailed2.state))
    user=graph.input("user","USER enzyme",MaterialKind.CATALYST)
    opened=graph.product("opened","USER-opened PCR-addressable library precursor",
                         MaterialKind.MOLECULAR, state=opened_template())
    loops=graph.product("loops","USER hairpin-cleavage products; recovery unassigned",MaterialKind.MOLECULAR)
    graph.react("open_hairpins","Open both uracil-containing hairpins",
                StepKind.DIGESTION,
                ((double_hairpin,InputRole.TEMPLATE),(user,InputRole.CATALYST)),
                (opened,loops))
    graph.mark_terminal(Terminal(loops.id,Disposition.UNRESOLVED,"USER opens the hairpins. Products are outside the retained adapter arms, but their subsequent physical removal is not established by cleavage itself."))
    forward=graph.input("pcr_forward","TruSeq PCR forward primer",MaterialKind.OLIGO)
    reverse=graph.input("pcr_reverse","indexed TruSeq PCR reverse primer",MaterialKind.OLIGO)
    polymerase=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    graph.react("index_pcr",f"Amplify the opened END-seq library for {PCR_CYCLES} source cycles",StepKind.PCR,
                ((opened,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
                 (reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (library,),geometry=PCRGeometry.unavailable(
                    opened.state,reason="Actual USER-opened retained-arm DNA enters PCR. The source names TruSeq index adapters but exact physical primer members/widths and full complementary precursor ends remain unassigned."))
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,
                                          primers=tuple(SEQ_PRIMERS), acquired_roles=("Read 1",), limitation="The original reports alternative50-base HiSeq2500 and75-base NextSeq500 single-end runs; exact index-cycle recipe is unassigned. Reference primer geometry does not assign actual loaded oligos."))
    graph.declare_reviewed_scope((
        "genomic_dna","deproteinization_condition","deproteinize","repair_enzymes",
        "blunt","tailing_polymerase","tail_break","a1","ligase1","ligate_a1",
        "gelase","agarose_melting","recover_plug","shear_condition","shear","streptavidin","capture","repair_distal",
        "distal_repair_enzymes","distal_tailing_enzyme","tail_distal","a2","ligase2","ligate_a2","user","open_hairpins",
        "pcr_forward","pcr_reverse","pcr_polymerase","index_pcr","library"),SOURCE)
    return graph.require_publishable()

# The bundle scraper splits this bracketed base modification into separate hits.
# selftest.py checks the complete source-written modified oligo and its biotin sites.
SOURCE_OLIGO_EXCEPTIONS = {
    "a1": "Parser limitation: published PMC6299834 Experimental Procedures line159 embeds two [Biotin-dT] bases in one hairpin; scraper splits the whole oligo at these brackets. tools/selftest.py directly parses the complete original modified sequence and checks all bases plus biotin-dT positions."
}

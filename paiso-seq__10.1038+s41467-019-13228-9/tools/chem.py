"""PAIso-seq end extension, USER release and PacBio SMRTbell conversion."""
from batch_ngs import joined_scene, seg
from chemdraw import Construct, Row, feature
from dumbbell import smrtbell_circuit, smrtbell_hairpin_geometry, smrtbell_sequencing_plan
from endprep import dA_tailed_scene, repair_and_dA_tail
from longread_footprinting import conventional_smrtbell, pacbio_entry, smrtbell_rows
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_extensions import capped_full_length_cdna
TITLE="PAIso-seq — poly(A)-inclusive full-length isoform sequencing"
NOTES="01_paiso-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-019-13228-9">Liu et al., <i>Nature Communications</i> (2019)</a>; detailed protocol <a href="https://doi.org/10.1038/s41596-022-00704-8">Nat. Protocols (2022)</a>.'
SUMMARY="A dU-interrupted oligo templates end extension across the complete RNA poly(A) tail, USER releases the copied product, template switching captures the 5-prime end, and PacBio CCS reads the entire isoform plus native tail."
INSERT=Construct([seg("5-prime cDNA handle","AAGCAGTGGTATCAACGCAGAGT","tso"),seg("full-length cDNA","X"*38,placeholder=True),seg("copied native poly(A) tail","T"*18),seg("sample barcode","B"*16,"cbc",placeholder=True,feature=feature("sample_barcode","sample_index","whitelist",whitelist="published PAIso-seq sample-barcode set"))],name="PAIso-seq full-length cDNA")
SMRTBELL=conventional_smrtbell(INSERT)
FINAL_LIBRARY=None; SEQ_PRIMERS=(); SEQUENCING_ENDING=pacbio_entry(SMRTBELL)
def oligos():
    return ["<b>TSO</b> — 5′-AAGCAGTGGTATCAACGCAGAGTACATrGrG+G-3′","<b>RT primer</b> — 5′-AAGCAGTGGTATCAACGCAGAGTAC-3′","<b>PCR primer</b> — 5′-AAGCAGTGGTATCAACGCAGAGT-3′"]
def sections():
    return [("Template an end extension across the native poly(A) tail",[Row(chunks=[("RNA body—poly(A) ← dU-interrupted barcoded oligo", "r1", False)])],"The designed oligo lets polymerase copy the complete tail rather than prime within it."),("Open the dU sites with USER and reverse-transcribe",capped_full_length_cdna(retain_poly_a=True).rows(),"USER removes the templating oligo; RT and template switching place a known handle at the transcript 5-prime end."),("Amplify full-length cDNA",[Row(chunks=[("TSO handle—full-length cDNA—copied native tail—barcode", "tso", False)])],"The tail remains inside the amplicon rather than being replaced by an oligo(dT) tract."),("Ligate PacBio hairpins",smrtbell_rows(SMRTBELL),"Damage repair, end repair and hairpin ligation produce a closed CCS template."),]


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    rna=graph.input("rna","poly(A) RNA with its native tail",MaterialKind.BIOLOGICAL,
                    primary=True)
    extension_oligo=graph.input(
        "extension_oligo","barcoded dU-interrupted poly(A)-tail extension oligo",
        MaterialKind.OLIGO,inferred=True,
        source=SOURCE+" (the protocol defines a barcode family rather than one universal member)")
    extension_polymerase=graph.input("extension_polymerase","tail end-extension polymerase",
                                     MaterialKind.CATALYST)
    extended=graph.product("extended","end-extended RNA–oligo complex",
                           MaterialKind.MOLECULAR)
    graph.react("extend_tail","Template end extension across the complete native poly(A) tail",
                StepKind.EXTENSION,
                ((rna,InputRole.TEMPLATE),(extension_oligo,InputRole.STRUCTURAL),
                 (extension_polymerase,InputRole.CATALYST)),(extended,))
    user=graph.input("user","USER enzyme",MaterialKind.CATALYST)
    released=graph.product("released","USER-released tail-extended RNA template",
                           MaterialKind.MOLECULAR)
    oligo_waste=graph.product("oligo_waste","cleaved dU-templating oligo",MaterialKind.WASTE)
    graph.react("open_du","Open dU sites and release the copied tail end",StepKind.DIGESTION,
                ((extended,InputRole.TEMPLATE),(user,InputRole.CATALYST)),
                (released,oligo_waste))
    graph.discard(oligo_waste,"USER cleavage removes the tail-templating oligo")
    rt_primer=graph.input("rt_primer","published PAIso-seq RT primer",MaterialKind.OLIGO)
    rt=graph.input("reverse_transcriptase","template-switching reverse transcriptase",
                   MaterialKind.CATALYST)
    first=graph.product("first_strand","tail-inclusive first-strand cDNA",
                        MaterialKind.MOLECULAR)
    graph.react("reverse_transcribe","Reverse-transcribe through the copied native tail",
                StepKind.REVERSE_TRANSCRIPTION,
                ((released,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(first,))
    tso=graph.input("tso","published PAIso-seq template-switch oligo",MaterialKind.OLIGO)
    full=graph.product("full_length","full-length tail-inclusive cDNA with both handles",
                       MaterialKind.MOLECULAR,state=INSERT)
    graph.react("template_switch","Capture the transcript 5-prime end by template switching",
                StepKind.TEMPLATE_SWITCHING,
                ((first,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(full,))
    pcr_primer=graph.input("pcr_primer","published PAIso-seq PCR primer",MaterialKind.OLIGO)
    pcr_reverse=graph.input(
        "barcoded_reverse","barcode-family reverse PCR primer",MaterialKind.OLIGO,
        inferred=True,source=SOURCE+" (barcode family varies by sample)")
    pcr_polymerase=graph.input("pcr_polymerase","full-length cDNA PCR polymerase",
                               MaterialKind.CATALYST)
    amplified=graph.product("amplified",INSERT.name,MaterialKind.MOLECULAR,state=INSERT)
    graph.react("pcr","Amplify the tail-inclusive full-length cDNA",StepKind.PCR,
                ((full,InputRole.TEMPLATE),(pcr_primer,InputRole.STRUCTURAL),
                 (pcr_reverse,InputRole.STRUCTURAL),(pcr_polymerase,InputRole.CATALYST)),
                (amplified,),geometry=PCRGeometry.unavailable(
                    INSERT,reason="The source provides a barcode family; one sample-specific reverse-primer member is not universal."))
    end_mix=graph.input("end_mix","PacBio damage-repair, end-repair and dA-tailing mix",
                        MaterialKind.CATALYST)
    end_product=repair_and_dA_tail(INSERT,inferred=True)
    tailed=graph.product("tailed","repaired dA-tailed full-length PAIso-seq cDNA",
                         MaterialKind.MOLECULAR,state=dA_tailed_scene(end_product))
    graph.react("end_prep","Repair and dA-tail full-length cDNA",StepKind.END_REPAIR,
                ((amplified,InputRole.TEMPLATE),(end_mix,InputRole.CATALYST)),(tailed,))
    hairpins=graph.input("hairpins","T-overhang PacBio SMRTbell hairpins",
                         MaterialKind.ADAPTER,inferred=True,
                         source=SOURCE+" (hairpin-body bases proprietary)")
    ligase=graph.input("smrt_ligase","SMRTbell ligase",MaterialKind.CATALYST)
    circuit=smrtbell_circuit(SMRTBELL)
    library=graph.product("library",SMRTBELL.name,MaterialKind.LIBRARY,
                          state=tuple(smrtbell_rows(SMRTBELL)),inferred=True,
                          source=SOURCE+" (hairpin-body bases proprietary)")
    graph.react("ligate_hairpins","Ligate hairpins to both full-length cDNA ends",
                StepKind.LIGATION,
                ((tailed,InputRole.TEMPLATE),(hairpins,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(library,),
                geometry=smrtbell_hairpin_geometry(
                    SMRTBELL,statement="PacBio SMRTbell preparation ligates T-overhang hairpins to dA-tailed cDNA."))
    graph.sequence(library,smrtbell_sequencing_plan(
        SMRTBELL,circuit,
        limitation="The PacBio sequencing-primer site lies in a proprietary hairpin adapter."))
    graph.declare_reviewed_scope((
        "rna","extension_oligo","extension_polymerase","extend_tail","user","open_du",
        "rt_primer","reverse_transcriptase","reverse_transcribe","tso",
        "template_switch","pcr_primer","barcoded_reverse","pcr_polymerase","pcr",
        "end_mix","end_prep","hairpins","smrt_ligase","ligate_hairpins","library"),SOURCE)
    return graph.require_publishable()

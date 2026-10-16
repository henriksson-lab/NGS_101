"""PAIso-seq end extension, USER release and PacBio SMRTbell conversion."""
from batch_ngs import joined_scene, seg
from chemdraw import Construct, Row, feature, revcomp
from dumbbell import Dumbbell, smrtbell_circuit, smrtbell_hairpin_geometry, smrtbell_sequencing_plan
from longread_footprinting import pacbio_entry, smrtbell_rows
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_extensions import capped_full_length_cdna
import rt
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "smart-seq2__10.1038+nmeth.2639" / "tools"))
import smartseq2 as s2
TITLE="PAIso-seq — poly(A)-inclusive full-length isoform sequencing"
NOTES="01_paiso-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-019-13228-9">Liu et al., <i>Nature Communications</i> (2019)</a>; detailed protocol <a href="https://doi.org/10.1038/s41596-022-00704-8">Nat. Protocols (2022)</a>.'
SUMMARY="A dU-interrupted oligo templates end extension across the complete RNA poly(A) tail, USER releases the copied product, template switching captures the 5-prime end, and PacBio CCS reads the entire isoform plus native tail."
CAVEAT=('The 2019 paper (Methods, “PAIso-seq library construction”) names a PacBio '
        'SMRTbell Template Prep Kit without its version. The later <a href="https://doi.org/10.1038/s41596-022-00704-8">'
        '2022 author protocol</a> specifies kit 1.0-SPv3 (Materials and step 29), '
        'but that does not establish the 2019 experiment\'s kit. The <a href="https://www.pacb.com/wp-content/uploads/2015/09/Insert-SMRTbell-Template-Prep-Kit-1.0.pdf">'
        'PacBio 1.0 insert, p. 1</a> lists a Blunt Adapter, while the <a href="https://www.pacb.com/wp-content/uploads/SMRTbell-Express-TPK-2.0-Overview-Large-Insert-gDNA-Library-Prep-Customer-Training.pdf">'
        'Express 2.0 overview, pp. 6 and 8</a> specifies a T-overhang adapter and A/T ligation. '
        'The model preserves the closed SMRTbell topology but leaves this experiment\'s '
        'end-preparation and insert–adapter junction geometry unresolved. Drawn N hairpins '
        'have schematic length; their bases and length are not established. '
        'The later protocol also adds 5′ iCiGiC and 3′ amine blocking to its TSO '
        '(Table 2); this model uses the 2019 Supplementary Table 1 TSO instead.')
# PAIso-seq reuses the SMART-seq2 handle system, so the shared model supplies the bases
# instead of a second transcription of them.  Supplementary Table 1 of the defining paper
# ("Oligos for PAIso-seq library construction") prints TSO, RT primer and PCR primer in
# full, and its TSO carries the four-base ACAT linker of smart-seq2__10.1038+nmeth.2639/tools/smartseq2.TSO_LINKER — the
# Nature Protocols form of SMART-seq2 — so no base here is inherited by assumption.  The
# self-test pins the three ordering strings against that table.
PCR_HANDLE = rt.SMART_HANDLE
GUIDE_LINKER = s2.DT_LINKER
TSO_LINKER = s2.TSO_LINKER
TSO_ORDERING = s2.TSO
RT_PRIMER = PCR_HANDLE + GUIDE_LINKER
PCR_PRIMER = PCR_HANDLE
BARCODES = (
    "TCAGACGATGCGTCAT", "CTATACATGACTCTGC", "TACTAGAGTAGCACTC",
    "TGTGTATCAGTACATG", "GATCTCTACTATATGC", "ACAGTCTATACTGCTG",
    "ATGATGTGCTACATCT", "CTGCGTGCTCTACGAC", "GCGCGATACGATGACT",
    "CGCGCTCAGCTGATCG", "GCGCACGCACTACAGA", "ACACTGACGTCGCGAC",
    "CGTCTATATACGTATA", "ATAGAGACTCAGAGCT", "TAGATGCGAGAGTAGA",
)
DU_TAIL = "dUTTTTTTTTdUTTTTTTTT"


def extension_oligo(sample: int) -> str:
    """Ordering string for the exact Supplementary Table 1 2dU-BC family."""
    if not 1 <= sample <= len(BARCODES):
        raise ValueError("PAIso-seq barcode number must be 1..15")
    return RT_PRIMER + BARCODES[sample - 1] + DU_TAIL


# BC3 is the bulk-library guide named in the defining paper's Methods.  The top
# strand is drawn from the template-switch end toward the RT-primer end, hence
# the copied BC3/linker and terminal handle appear as reverse complements.
REPRESENTATIVE_BARCODE = BARCODES[2]
INSERT=Construct([
    seg("5-prime PCR handle", PCR_HANDLE, "tso"),
    seg("TSO linker", TSO_LINKER),
    seg("template-switch G contribution", "GGG", "tso"),
    seg("full-length transcript", "X"*38, placeholder=True),
    seg("native poly(A) tail", "A"*18),
    seg("BC3 reverse complement", revcomp(REPRESENTATIVE_BARCODE), "cbc",
        feature=feature("sample_barcode", "sample_index", "whitelist",
                        whitelist="published PAIso-seq 2dU-BC1..BC15 set")),
    seg("guide linker reverse complement", revcomp(GUIDE_LINKER)),
    seg("3-prime PCR handle reverse complement", revcomp(PCR_HANDLE), "tso"),
], name="PAIso-seq BC3 full-length cDNA")
SMRTBELL=Dumbbell(INSERT,"N"*18,"N"*18,
                   name=INSERT.name+" SMRTbell",insert_overhang="",
                   adapter_overhang="",junction_unavailable=True)
FINAL_LIBRARY=None; SEQ_PRIMERS=(); SEQUENCING_ENDING=pacbio_entry(SMRTBELL)
def oligos():
    fixed = [
        f"<b>TSO</b> — 5′-{TSO_ORDERING}-3′",
        f"<b>RT primer</b> — 5′-{RT_PRIMER}-3′",
        f"<b>PCR primer</b> — 5′-{PCR_PRIMER}-3′",
    ]
    return fixed + [
        f"<b>2dU-BC{i}</b> — 5′-{extension_oligo(i)}-3′"
        for i in range(1, len(BARCODES) + 1)
    ]
def sections():
    return [("Template an end extension across the native poly(A) tail",[Row(chunks=[("RNA body—poly(A) ← dU-interrupted barcoded oligo", "r1", False)])],"The designed oligo lets polymerase copy the complete tail rather than prime within it."),("Open the dU sites with USER and reverse-transcribe",capped_full_length_cdna(retain_poly_a=True).rows(),"USER removes the templating oligo; RT and template switching place a known handle at the transcript 5-prime end."),("Amplify full-length cDNA",[Row(chunks=[("TSO handle—full-length cDNA—copied native tail—barcode", "tso", False)])],"The tail remains inside the amplicon rather than being replaced by an oligo(dT) tract."),("Convert cDNA to a PacBio SMRTbell",smrtbell_rows(SMRTBELL),"The paper establishes a closed CCS template but does not identify the kit version or its insert–adapter junction chemistry."),]


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    rna=graph.input("rna","poly(A) RNA with its native tail",MaterialKind.BIOLOGICAL,
                    primary=True)
    extension_oligo=graph.input(
        "extension_oligo","published 2dU-BC1..BC15 poly(A)-tail extension oligo family",
        MaterialKind.OLIGO,
        source=SOURCE+" (Supplementary Table 1; BC3 is the bulk-library member)")
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
    pcr_polymerase=graph.input("pcr_polymerase","full-length cDNA PCR polymerase",
                               MaterialKind.CATALYST)
    amplified=graph.product("amplified",INSERT.name,MaterialKind.MOLECULAR,state=INSERT)
    graph.react("pcr","Amplify the tail-inclusive full-length cDNA",StepKind.PCR,
                ((full,InputRole.TEMPLATE),(pcr_primer,InputRole.STRUCTURAL),
                 (pcr_polymerase,InputRole.CATALYST)),
                (amplified,),geometry=PCRGeometry.symmetric(
                    INSERT, PCRPrimer("PCR primer", PCR_PRIMER),
                    statement=("The defining paper states that amplified cDNA uses the "
                               "single published PCR primer; the common handle occurs at "
                               "both inward-facing ends.")))
    end_mix=graph.input("end_mix","PacBio SMRTbell kit end-preparation reagents (version unspecified)",
                        MaterialKind.CATALYST)
    prepared=graph.product("prepared","kit-prepared full-length PAIso-seq cDNA; end geometry unresolved",
                           MaterialKind.MOLECULAR,state=INSERT)
    graph.react("end_prep","Prepare cDNA ends with the version-unspecified SMRTbell kit",StepKind.END_REPAIR,
                ((amplified,InputRole.TEMPLATE),(end_mix,InputRole.CATALYST)),(prepared,))
    hairpins=graph.input("hairpins","PacBio SMRTbell hairpin adapters; ligating ends unresolved",
                         MaterialKind.ADAPTER,inferred=True,
                         source=SOURCE+" (kit version and hairpin ends unspecified)")
    ligase=graph.input("smrt_ligase","SMRTbell ligase",MaterialKind.CATALYST)
    circuit=smrtbell_circuit(SMRTBELL)
    library=graph.product("library",SMRTBELL.name,MaterialKind.LIBRARY,
                          state=tuple(smrtbell_rows(SMRTBELL)),inferred=True,
                          source=SOURCE+" (hairpin-body bases and junction geometry unresolved)")
    graph.react("ligate_hairpins","Ligate hairpins to both full-length cDNA ends",
                StepKind.LIGATION,
                ((prepared,InputRole.TEMPLATE),(hairpins,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(library,),
                geometry=smrtbell_hairpin_geometry(
                    SMRTBELL,statement="PAIso-seq names an unversioned SMRTbell Template Prep Kit; contemporaneous PacBio 1.0 and Express 2.0 documents describe different adapter ends, so the actual ligating ends are unresolved."))
    graph.sequence(library,smrtbell_sequencing_plan(
        SMRTBELL,circuit,
        limitation="The exact hairpin adapter and sequencing-primer binding site are not assigned to the unversioned 2019 kit in this model."))
    graph.declare_reviewed_scope((
        "rna","extension_oligo","extension_polymerase","extend_tail","user","open_du",
        "rt_primer","reverse_transcriptase","reverse_transcribe","tso",
        "template_switch","pcr_primer","pcr_polymerase","pcr",
        "end_mix","end_prep","hairpins","smrt_ligase","ligate_hairpins","library"),SOURCE)
    return graph.require_publishable()

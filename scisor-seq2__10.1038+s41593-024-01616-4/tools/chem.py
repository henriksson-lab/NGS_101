"""ScISOr-seq2 barcoded 10x cDNA enrichment for PacBio and ONT."""
from batch_ngs import seg
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, MolecularState, Row, Scene, Workflow, feature, revcomp
from endprep import dA_tailed_scene, repair_and_dA_tail
from longread_footprinting import conventional_smrtbell, pacbio_entry, smrtbell_rows
TITLE="ScISOr-seq2 — enriched single-cell isoforms by long reads"
NOTES="01_scisor-seq2.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41593-024-01616-4">Joglekar et al., <i>Nature Neuroscience</i> (2024)</a>.'
SUMMARY="10x 3-prime cDNA is strand-selectively re-amplified to retain cell barcodes, exon-capture probes enrich spliced molecules, and the same enriched cDNA can enter PacBio HiFi or ONT ligation sequencing."
INSERT=Construct([seg("Partial Read 1 handle","CTACACGACGCTCTTCCGATCT","r1"),seg("cell barcode","B"*16,"cbc",placeholder=True,feature=feature("cell_barcode","cell_barcode","whitelist",whitelist="10x Chromium 3-prime v3 whitelist")),seg("UMI","U"*12,"umi",placeholder=True,feature=feature("umi","umi","random")),seg("full-length spliced cDNA","X"*42,placeholder=True),seg("Partial TSO handle","AAGCAGTGGTATCAACGCAGAGTACAT","tso")],name="LAP-CAP enriched barcoded cDNA")
SMRTBELL=conventional_smrtbell("ScISOr-seq2 enriched cDNA")
FINAL_LIBRARY=None; SEQ_PRIMERS=(); SEQUENCING_ENDING=pacbio_entry(SMRTBELL)+" The same enriched cDNA may instead receive an ONT motor adapter, in which case no synthesis sequencing primer is used."
def oligos():
    return ["<b>Partial Read 1</b> — 5′-CTACACGACGCTCTTCCGATCT-3′","<b>Partial TSO</b> — 5′-AAGCAGTGGTATCAACGCAGAGTACAT-3′"]
def workflow():
    w=Workflow(MolecularState("10x 3′ cDNA",(Row(chunks=[("Read1 handle—cell barcode—UMI—transcript cDNA—TSO","cbc",False)]),)))
    w.react("Perform linear/asymmetric then exponential amplification",[Row(chunks=[("Partial Read1-only cycles → strand bias; + Partial TSO → full amplicon","r1",False)])],name="LAP amplicons",note="LAP removes non-barcoded products while retaining both terminal handles.")
    w.react("Capture spliced cDNA with exon probes",[Row(chunks=[("blocked cDNA + biotin exon probes ● streptavidin → CAP enrichment","w1",False)])],name="CAP-enriched cDNA",note="Hybrid capture enriches multi-exonic transcripts.")
    tailed=repair_and_dA_tail(Construct([seg("ScISOr-seq2 enriched cDNA","X"*42,placeholder=True)]))
    w.react("Repair and add 3′ dA",dA_tailed_scene(tailed).rows(),name="dA-tailed enriched cDNA",note="PacBio end prep creates the two hairpin-ligation substrates.")
    adapter=Scene(); adapter.strand("hairpin end",[seg("proprietary PacBio hairpin core","[proprietary hairpin]",placeholder=True),seg("3′ T","T")],label="PacBio hairpin ligating end"); adapter.overhang("hairpin end","3′ T",polarity="3-prime",owner="adapter")
    geometry=AdapterLigationGeometry.known("cohesive",dA_tailed_scene(tailed),adapter,statement="The PacBio SMRTbell workflow uses proprietary hairpin sequences with a disclosed 3′-T ligating overhang complementary to each insert 3′-dA.")
    w.ligate_adapter("Ligate PacBio hairpins",geometry,smrtbell_rows(SMRTBELL),name="closed SMRTbell",note="The ONT alternative instead uses a proprietary motor-bearing ligation adapter whose end bases are not disclosed here.")
    return w

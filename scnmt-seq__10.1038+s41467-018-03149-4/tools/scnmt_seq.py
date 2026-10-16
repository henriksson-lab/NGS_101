"""Clark et al. 2018 scNMT-seq chemistry and its two library branches."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/"lib"), str(ROOT/"scbs-seq__10.1038+nmeth.3035"/"tools")]
sys.path.insert(0,str(ROOT/"smart-seq2__10.1038+nmeth.2639"/"tools"))

from base_conversion import bisulfite_path
from batch_ngs import nextera_library, seg
from chemdraw import Construct, Row
from dataclasses import replace
import seqprimers as sp
import rt
import smartseq2 as s2
from multimodal_spatial import finish_pcr_branch
from pcr import PCRGeometry
from protocol_graph import Disposition, Terminal, SequencingPlan, InputRole, MaterialKind, ProtocolGraph, StepKind
import scbs

TITLE = "scNMT-seq — RNA, DNA methylation and accessibility"
NOTES = "01_scnmt-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/s41467-018-03149-4">Clark et al., <i>Nature Communications</i> (2018)</a>.'
SUMMARY = "M.CviPI marks accessible GpC sites, then bead capture separates RNA for Smart-seq2/Nextera XT from DNA for scBS-seq."
CAVEAT = ("<b>INFERRED DNA scaffold.</b> Clark et al., Methods, Library preparation names 6NF for both synthesis stages. The exact 6NF assignment is deferred; the inherited N9 scBS-seq scaffold below is provisional, not established for this experiment. The graph summarizes recovery stages. G–C–G and C–C–G positions are excluded from methylation interpretation (Results and DNA processing) because GCG is ambiguous and CCG can reflect off-target marking. "
          "The scBS-seq reaction and primer scaffold follow the cited original; "
          "the experiment-specific iTAG index member is unresolved. The defining paper explicitly uses "
          "Smart-seq2 followed by Nextera XT for RNA; its cited Smart-seq2 protocol and "
          "the canonical Nextera chemistry supply the reusable primer geometry. It does "
          "not print the template-switch oligo, though, and Smart-seq2 has two published "
          "TSO linkers: <code>AC</code> for the chosen <code>rGrG+G</code> oligo in "
          "Picelli et al. 2013 Supplementary Table 2, <code>ACAT</code> in the "
          "<i>Nature Protocols</i> form. Clark et al. cite both papers at once, so the "
          "longer modelled form is an inherited choice and the whole TSO linker is "
          "drawn dotted.")


def oligos():
    return [
        f"Smart-seq2 oligo-dT30VN: 5′-{s2.OLIGO_DT}-3′",
        f"Smart-seq2 TSO (linker inherited, not printed by this paper): "
        f"5′-{rt.SMART_HANDLE}<inf>{s2.TSO_LINKER}</inf>{rt.TSO_G_TAIL_LNA}-3′",
        f"Smart-seq2 ISPCR primer: 5′-{s2.ISPCR}-3′",
        "Nextera XT indexing primers: P5–i5–s5 and P7–i7–s7; index bases vary by sample",
    ]


def rna_library():
    return nextera_library([seg("Smart-seq2 cDNA insert", "X"*36, placeholder=True)],
                           "scNMT-seq RNA library")


def dna_library():
    ref = scbs.final_library()
    return Construct([replace(s, inferred=True) for s in ref], name="INFERRED scBS-seq reference; scNMT 6NF assignment pending"), tuple(scbs.SEQ_PRIMERS)


RNA_LIBRARY, RNA_PRIMERS = rna_library()
RNA_PREAMP = s2.preamp_template()
DNA_LIBRARY, DNA_PRIMERS = dna_library()
FINAL_LIBRARIES = (
    ("RNA branch", RNA_LIBRARY, RNA_PRIMERS,
     "Nextera XT library made from Smart-seq2 cDNA.",
     "Nextera sequencing reads the RNA-derived insert."),
    ("DNA branch", DNA_LIBRARY, DNA_PRIMERS,
     "INFERRED — inherited scBS-seq reference; study-specific 6NF assignment pending.",
     "The dedicated historical iPCRTag primer reads i7; there is no i5."),
)


def sections():
    return [
        ("Mark accessible chromatin", [
            Row(chunks=[("nucleosome-free GpC + M.CviPI + SAM → Gp5mC", "w1", False)]),
            Row(chunks=[("nucleosome-protected GpC → unmarked", None, False)]),
        ], "GpC methylation records accessibility before the cell is lysed."),
        ("Separate RNA and DNA", [
            Row(chunks=[("poly(A) RNA → oligo-dT magnetic beads → Smart-seq2 cDNA", "r1", False)]),
            Row(chunks=[("supernatant genomic DNA → bisulfite conversion → scBS-seq", "r2", False)]),
        ], "Repeated bead washes are transferred to the DNA plate to maximize recovery."),
        ("Interpret the DNA branch", [
            Row(chunks=[(bisulfite_path(protected=False).text() + "  (closed unmarked GpC or unmethylated CpG)", None, False)]),
            Row(chunks=[(bisulfite_path(protected=True).text() + "  (accessible M.CviPI-marked GpC or endogenous 5mCpG)", "w1", False)]),
        ], "Sequence context separates exogenous GpC accessibility marks from endogenous CpG methylation."),
    ]


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cell=g.input("cell","one permeabilized cell containing chromatin and poly(A) RNA",MaterialKind.BIOLOGICAL,state=sections()[0][1],primary=True)
    methylase=g.input("mcvi_pi","M.CviPI GpC methyltransferase",MaterialKind.CATALYST)
    sam=g.input("sam","S-adenosyl methionine methyl donor",MaterialKind.CONDITION)
    marked=g.product("marked","cell with accessible GpC sites marked as Gp5mC",MaterialKind.MOLECULAR,state=sections()[0][1])
    g.react("mark_accessibility","Methylate accessible GpC sites with M.CviPI",StepKind.CHEMICAL_CONVERSION,((cell,InputRole.TEMPLATE),(methylase,InputRole.CATALYST),(sam,InputRole.CONDITION)),(marked,))
    beads=g.input("oligodt_beads","oligo(dT) magnetic beads",MaterialKind.SELECTION)
    rna=g.product("rna_fraction","bead-bound poly(A) RNA",MaterialKind.MOLECULAR,state=sections()[1][1][:1])
    dna=g.product("dna_fraction","GpC-marked genomic-DNA supernatant plus washes",MaterialKind.MOLECULAR,state=sections()[1][1][1:])
    g.react("separate","Separate poly(A) RNA from marked genomic DNA",StepKind.AFFINITY_CAPTURE,((marked,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(rna,dna),changes_molecule=False)

    rt_primer=g.input("rt_primer","Smart-seq2 anchored oligo(dT) primer",MaterialKind.OLIGO,
                      state=(seg("SMART handle",rt.SMART_HANDLE,"tso"),
                             seg("oligo-dT linker",s2.DT_LINKER),
                             seg("dT30VN","T"*30+"VN",placeholder=True)))
    rt_enzyme=g.input("rt","SuperScript II reverse transcriptase",MaterialKind.CATALYST)
    copied=g.product("copied","first-strand Smart-seq2 cDNA",MaterialKind.MOLECULAR,state=RNA_PREAMP)
    g.react("reverse_transcribe","Reverse-transcribe bead-captured RNA",StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt_enzyme,InputRole.CATALYST)),(copied,))
    tso=g.input("tso","Smart-seq2 template-switch oligo",MaterialKind.OLIGO,
                state=(seg("SMART handle",rt.SMART_HANDLE,"tso"),
                       seg("TSO linker",s2.TSO_LINKER,inferred=True),
                       seg("rGrG+G","GGG","tso")))
    switched=g.product("switched","template-switched full-length cDNA",MaterialKind.MOLECULAR,state=RNA_PREAMP)
    g.react("template_switch","INFERRED — Template-switch the nascent cDNA (the TSO linker is inherited from Smart-seq2; this paper prints no TSO)",StepKind.TEMPLATE_SWITCHING,((copied,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),(rt_enzyme,InputRole.CATALYST)),(switched,))
    rna_f=g.input("rna_preamp_forward","Smart-seq2 ISPCR primer at the 5′ handle",MaterialKind.OLIGO,state=(seg("ISPCR",s2.ISPCR,"tso"),))
    rna_r=g.input("rna_preamp_reverse","same Smart-seq2 ISPCR species at the 3′ handle",MaterialKind.OLIGO,state=(seg("ISPCR",s2.ISPCR,"tso"),))
    rna_pol=g.input("rna_preamp_polymerase","Smart-seq2 cDNA PCR polymerase",MaterialKind.CATALYST)
    preamp=g.product("rna_preamp","amplified full-length cDNA",MaterialKind.MOLECULAR,state=RNA_PREAMP)
    g.react("preamp_rna","Amplify full-length Smart-seq2 cDNA: 18 cycles",StepKind.PCR,((switched,InputRole.TEMPLATE),(rna_f,InputRole.STRUCTURAL),(rna_r,InputRole.STRUCTURAL),(rna_pol,InputRole.CATALYST)),(preamp,),geometry=s2.ispcr_geometry(RNA_PREAMP))
    transposome=g.input("rna_transposome","Nextera XT transposome",MaterialKind.ADAPTER)
    tn5=g.input("rna_tn5","Tn5 transposase",MaterialKind.CATALYST)
    rna_tagged=g.product("rna_tagged","Nextera-tagged Smart-seq2 cDNA",MaterialKind.MOLECULAR,state=RNA_LIBRARY)
    rwaste=g.product("rna_tn5_waste","candidate non-library cDNA tagmentation products",MaterialKind.MOLECULAR)
    g.react("tagment_rna","Tagment Smart-seq2 cDNA with Nextera XT",StepKind.TAGMENTATION,((preamp,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(rna_tagged,rwaste));g.mark_terminal(Terminal(rwaste.id,Disposition.UNRESOLVED,"Candidate insertion products and their physical recovery remain unresolved."))

    bisulfite=g.input("bisulfite","bisulfite-conversion condition",MaterialKind.CONDITION)
    converted=g.product("converted","bisulfite-converted GpC-marked genomic DNA",MaterialKind.MOLECULAR,state=DNA_LIBRARY)
    g.react("bisulfite_convert","Bisulfite-convert the marked genomic-DNA fraction",StepKind.CHEMICAL_CONVERSION,((dna,InputRole.TEMPLATE),(bisulfite,InputRole.CONDITION)),(converted,))
    oligo1=g.input("oligo1","6NF primer, first synthesis; sequence assignment pending",MaterialKind.OLIGO,
                   state=None,
                   source="Clark et al. 2018 Methods name 6NF; exact assignment deferred")
    strand_pol=g.input("strand_polymerase","scBS-seq strand-synthesis polymerase",MaterialKind.CATALYST)
    first=g.product("first_strand","oligo1-primed converted-DNA copies",MaterialKind.MOLECULAR,state=DNA_LIBRARY)
    g.react("prime_oligo1","Random-prime converted DNA: five synthesis rounds",StepKind.EXTENSION,((converted,InputRole.TEMPLATE),(oligo1,InputRole.STRUCTURAL),(strand_pol,InputRole.CATALYST)),(first,))
    oligo2=g.input("oligo2","6NF primer as named for second synthesis; sequence assignment pending",MaterialKind.OLIGO,
                   state=None,
                   source="Clark et al. 2018 Methods name 6NF; exact assignment deferred")
    second=g.product("second_strand","oligo1/oligo2-completed scBS templates",MaterialKind.MOLECULAR,state=DNA_LIBRARY)
    g.react("prime_oligo2","Random-prime the complementary strand with oligo2",StepKind.EXTENSION,((first,InputRole.TEMPLATE),(oligo2,InputRole.STRUCTURAL),(strand_pol,InputRole.CATALYST)),(second,))

    reviewed=["cell","mcvi_pi","sam","mark_accessibility","oligodt_beads","separate","rt_primer","rt","reverse_transcribe","tso","template_switch","rna_preamp_forward","rna_preamp_reverse","rna_preamp_polymerase","preamp_rna","rna_transposome","rna_tn5","tagment_rna","bisulfite","bisulfite_convert","oligo1","strand_polymerase","prime_oligo1","oligo2","prime_oligo2"]
    rna_geometry=s2.nextera_index_geometry(RNA_LIBRARY)
    finish_pcr_branch(g,"rna",rna_tagged,RNA_LIBRARY,RNA_PRIMERS,reviewed,label="scNMT-seq RNA",geometry=rna_geometry)
    dna_geometry=PCRGeometry.unavailable(DNA_LIBRARY, reason="Study-specific 6NF and PCR primer assignment remains deferred; inherited scBS-seq endpoint is only a reference.")
    finish_pcr_branch(g,"dna",second,DNA_LIBRARY,DNA_PRIMERS,reviewed,label="scNMT-seq DNA",geometry=dna_geometry)
    g.steps["dna_pcr"] = replace(g.steps["dna_pcr"], label="Amplify scNMT-seq DNA library: 14 PCR cycles")
    for key, lib, primers in (("rna_library", RNA_LIBRARY, RNA_PRIMERS), ("dna_library", DNA_LIBRARY, DNA_PRIMERS)):
        g.sequencing[key] = SequencingPlan(library=lib, primers=primers,
            acquired_roles=("Read 1", "Read 2"),
            limitation="EL16 serum cohort uses paired 125-base reads for both branches. E14 EB uses DNA paired 75-base and RNA single-end 75-base reads. Index acquisition remains unassigned; DNA primer geometry is a provisional inherited reference.")
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()

READ_SPECS = {
    "RNA branch": tuple(sp.read_spec(r, "125 (EL16); 75 (E14 EB, Read 1 only)", "RNA-derived insert; paired EL16, single-end E14 EB", "include") for r in ("Read 1", "Read 2")),
    "DNA branch": tuple(sp.read_spec(r, "125 (EL16); 75 (E14 EB)", "bisulfite DNA; EL16 QC subset also 50-base single-end", "include") for r in ("Read 1", "Read 2")),
}

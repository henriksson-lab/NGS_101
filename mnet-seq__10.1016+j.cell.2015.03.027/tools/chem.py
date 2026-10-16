"""mNET-seq2015 owns its capture route and unassigned library states."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import seqprimers as sp
import illumina as il
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, feature
from pcr import PCRGeometry
from protocol_graph import Highlight, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind
from protocol_model import role

TITLE="mNET-seq — mammalian native elongating transcript sequencing"
NOTES="01_mnet_seq.html"
SOURCE=('Defining source: <a href="https://doi.org/10.1016/j.cell.2015.03.027">Nojima et al., Cell (2015)</a>, '
        'Methods and original Extended Experimental Procedures; '
        '<a href="https://doi.org/10.1038/nprot.2016.012">2016 successor protocol</a> distinguished below.')
SUMMARY=("MNase solubilizes native chromatin. Pol II immunoprecipitation recovers protected RNA; on-bead phosphorylation "
         "prepares its 5′ end before extraction and gel selection. TruSeq Small RNA library preparation supports paired reads whose "
         "Read 2 first base reports the RNA 3′ end, including cotranscriptional processing intermediates.")
CAVEAT=("This declaration selects the defining 2015 mNET-seq route. Ordered library adapters, physical sequencing primers, "
        "selected index member and acquisition cycle counts remain unassigned model data. Dated 2011 and current manufacturer "
        "documents provide positive kit-family information; they do not establish the exact revision used in this experiment. "
        "The 2016 successor specifies reagents and gives additional transfer/gel details that are not automatically assigned to 2015. "
        "No publication absence is accepted from these model gaps.")
SEQ_PRIMERS=()
ACQUIRED_ROLES=("Read 1","Read 2")
READ_LENGTHS={}
READ_SPECS=(sp.read_spec("Read 1",None,"RNA strand direction; exact molecular primer site and cycles unassigned","include"),
            sp.read_spec("Read 2",None,"first base gives the RNA 3′ end; cycles unassigned","include"),
            sp.read_spec("Index 1 (i7)",None,"selected acquisition and index member unassigned","undetermined"),
            sp.read_spec("Index 2 (i5)",None,"acquisition unassigned","undetermined"))
SEQUENCING_ENDING="The defining source reports paired sequencing on HiSeq 2000/2500. Read 2 starts at the informative RNA 3′ end; Read 1 supplies strand direction."
SEQUENCING_UNAVAILABLE="The selected 2015 ordered library and sequencing-primer sites remain unassembled; source-linked kit-family sequences do not establish a run-specific assignment."


def state(name,*parts):
    return Construct([role(p) for p in parts],name=name)


RNA=state("selected RNA","variable protected RNA with informative 3-prime hydroxyl")
LIGATED=state("two-adapter RNA","unassembled 5-prime adapter","variable protected RNA","unassembled 3-prime adapter")
CDNA=state("first-strand cDNA","unassembled adapter-derived cDNA tract","variable RNA-derived cDNA","unassembled opposite adapter-derived tract")
FINAL_LIBRARY=Construct([role("unassigned outer adapter/primer region"),
                         role("variable protected-RNA-derived interval"),
                         role("unassigned sample-index region",molecular_feature=feature("sample_index","sample_index","unknown")),
                         role("unassigned opposite adapter/primer region")],name="mNET-seq2015 library — ordered ends unassembled")
FINAL_CAPTION="Source-known library roles only. No generic eight-base index, fixed biological insert or exact read-primer oligo is assigned."


def oligos():
    # These are literal positive vendor-family references, not assigned 2015 inputs.
    return ["Manufacturer family comparator only — not selected 2015 oligos: RA5 (RNA, 5′→3′) "
            +il.TRUSEQ_SMALL_RNA_RA5.replace('T','U'),
            "Manufacturer family comparator only: RA3 (DNA, 5′→3′) "+il.TRUSEQ_SMALL_RNA_RA3]


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    chromatin=g.input("native_chromatin","native chromatin containing Pol II/RNA complexes",MaterialKind.BIOLOGICAL,
                      primary=True,state=state("native chromatin","chromatin","Pol II-protected RNA"))
    mnase=g.input("mnase","micrococcal nuclease",MaterialKind.CATALYST)
    digest=g.product("digest","digested chromatin mixture",MaterialKind.MOLECULAR)
    g.react("mnase_digest","Digest native chromatin with MNase",StepKind.DIGESTION,
            ((chromatin,InputRole.TEMPLATE),(mnase,InputRole.CATALYST)),(digest,))
    egta=g.input("egta","EGTA to stop MNase",MaterialKind.CONDITION)
    stopped=g.product("stopped_digest","MNase-inactivated mixture",MaterialKind.MOLECULAR)
    g.react("stop_mnase","Stop digestion with EGTA",StepKind.QUENCHING,
            ((digest,InputRole.TEMPLATE),(egta,InputRole.CONDITION)),(stopped,))
    spin=g.input("fractionation","centrifugation",MaterialKind.SELECTION)
    soluble=g.product("soluble_chromatin","soluble chromatin supernatant",MaterialKind.MOLECULAR)
    pellet=g.product("insoluble_pellet","insoluble chromatin pellet",MaterialKind.WASTE)
    g.react("recover_soluble","Collect the soluble supernatant",StepKind.SELECTION,
            ((stopped,InputRole.TEMPLATE),(spin,InputRole.SELECTION)),(soluble,pellet),changes_molecule=False)
    g.discard(pellet,"The selected mNET-seq route proceeds with the soluble chromatin supernatant.")
    antibodies=g.input("polii_beads","Pol II antibody-conjugated beads",MaterialKind.SELECTION)
    selected=g.product("polii_ip","washed Pol II/RNA immunoprecipitate on beads",MaterialKind.MOLECULAR)
    unbound=g.product("unbound","unbound and washed material",MaterialKind.WASTE)
    g.react("immunoprecipitate","Capture the selected Pol II CTD state and wash",StepKind.AFFINITY_CAPTURE,
            ((soluble,InputRole.TEMPLATE),(antibodies,InputRole.SELECTION)),(selected,unbound))
    g.discard(unbound,"Material not retained in the selected Pol II immunoprecipitate is removed by washes.")
    kinase=g.input("rna_kinase","T4 PNK, 3-prime phosphatase minus",MaterialKind.CATALYST)
    atp=g.input("atp","ATP for the on-bead kinase reaction",MaterialKind.CONDITION)
    phosphorylated=g.product("phosphorylated_ip","on-bead RNA with phosphorylated 5-prime ends",MaterialKind.MOLECULAR)
    g.react("phosphorylate","Phosphorylate RNA 5′ ends on the IP beads",StepKind.PHOSPHORYLATION,
            ((selected,InputRole.TEMPLATE),(kinase,InputRole.CATALYST),(atp,InputRole.CONDITION)),(phosphorylated,))
    extraction=g.input("rna_extraction","post-kinase wash and Trizol RNA extraction",MaterialKind.CONDITION)
    rna=g.product("extracted_rna","extracted immunoprecipitated RNA",MaterialKind.MOLECULAR,state=RNA)
    g.react("extract_rna","Wash the beads, then extract RNA",StepKind.EXTRACTION,
            ((phosphorylated,InputRole.TEMPLATE),(extraction,InputRole.CONDITION)),(rna,))
    size=g.input("rna_size_rule","35–100-nt denaturing RNA gel-selection interval",MaterialKind.SELECTION)
    sized=g.product("sized_rna","35–100-nt protected RNA",MaterialKind.MOLECULAR,state=RNA)
    excluded=g.product("excluded_rna","gel-excluded RNA sizes",MaterialKind.WASTE)
    g.react("select_rna","Gel-select the source 35–100-nt RNA interval",StepKind.SELECTION,
            ((rna,InputRole.TEMPLATE),(size,InputRole.SELECTION)),(sized,excluded))
    g.discard(excluded,"RNA outside the selected gel interval is excluded from the excised band.")
    adapter3=g.input("adapter3","TruSeq Small RNA 3-prime adapter; actual ordered reagent unassigned",MaterialKind.ADAPTER)
    ligase3=g.input("ligase3","kit-family 3-prime RNA-adapter ligase",MaterialKind.CATALYST)
    linked3=g.product("linked3","3-prime adapter-bearing RNA",MaterialKind.MOLECULAR)
    g.react("ligate3","Ligate the kit-family 3′ RNA adapter",StepKind.LIGATION,
            ((sized,InputRole.TEMPLATE),(adapter3,InputRole.STRUCTURAL),(ligase3,InputRole.CATALYST)),(linked3,),
            geometry=AdapterLigationGeometry.unavailable("The model has not assigned the exact 2015 adapter strands/modifications. Manufacturer-family bases are displayed separately as a comparator, not proof of this selected reagent."))
    adapter5=g.input("adapter5","TruSeq Small RNA 5-prime adapter; actual ordered reagent unassigned",MaterialKind.ADAPTER)
    ligase5=g.input("ligase5","kit-family 5-prime RNA-adapter ligase with ATP",MaterialKind.CATALYST)
    linked5=g.product("linked5",LIGATED.name,MaterialKind.MOLECULAR,state=LIGATED)
    g.react("ligate5","Ligate the kit-family 5′ RNA adapter",StepKind.LIGATION,
            ((linked3,InputRole.TEMPLATE),(adapter5,InputRole.STRUCTURAL),(ligase5,InputRole.CATALYST)),(linked5,),
            geometry=AdapterLigationGeometry.unavailable("Actual 2015 ordered 5-prime adapter and end chemistry are unassigned; a shared small-RNA recipe supplies no replacement."))
    rtprimer=g.input("rt_primer","kit RNA RT primer; selected ordered strand unassigned",MaterialKind.OLIGO)
    rt=g.input("rt_enzyme","kit reverse transcriptase; actual 2015 identity unassigned",MaterialKind.CATALYST)
    cdna=g.product("cdna",CDNA.name,MaterialKind.MOLECULAR,state=CDNA)
    g.react("reverse_transcribe","Reverse-transcribe the two-adapter RNA",StepKind.REVERSE_TRANSCRIPTION,
            ((linked5,InputRole.TEMPLATE),(rtprimer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,))
    rp1=g.input("pcr_forward","kit RP1-family PCR primer; selected strand unassigned",MaterialKind.OLIGO)
    rpi=g.input("pcr_index","kit RPI-family PCR primer; selected member unassigned",MaterialKind.OLIGO)
    polymerase=g.input("pcr_enzyme","kit PCR mixture; actual 2015 formulation unassigned",MaterialKind.CATALYST)
    library=g.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    g.react("library_pcr","Amplify and index the cDNA using the source-named kit",StepKind.PCR,
            ((cdna,InputRole.TEMPLATE),(rp1,InputRole.STRUCTURAL),(rpi,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
            geometry=PCRGeometry.unavailable(CDNA,reason="Actual 2015 PCR strands and binding geometry are unassembled; named manufacturer-family primers are positive references, not a selected sample-member assignment."))
    g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=(),acquired_roles=ACQUIRED_ROLES,
                                     limitation=SEQUENCING_UNAVAILABLE))
    g.annotate(Highlight("phosphorylate","On-bead kinase preserves the informative end",
        "The 2015 source explicitly uses T4 PNK 3′ phosphatase minus plus ATP before RNA extraction. The Pol II-protected RNA 3′ hydroxyl is retained. "
        "No generic RNA end-repair or 3′ dephosphorylation step is substituted."))
    g.annotate(Highlight("library_pcr","Vendor family and later protocol remain distinct",
        "The March 2011 manual and current sequence table describe a six-base RPI family, not the previous model's eight-base default. "
        "The 2011 guide supplies SuperScript II; the 2016 successor instead names SuperScript III and gives a 150–230-bp final gel interval with 120–125-bp background exclusion. "
        "Those later transfer amounts, enzyme selection and gel cutoffs are not assigned to the defining 2015 experiment."))
    g.annotate(Highlight("library","Paired read semantics and computational filtering",
        "Read 2's first base gives the RNA 3′ end; Read 1 provides directionality. Processing trims adapters, drops short/unpaired reads and keeps the informative nucleotide computationally. "
        "These operations do not physically remove RNA or library molecules. Processing-intermediate signals are not all interpreted as polymerase positions. "
        "The exact 2015 cycle counts and physical primer sites remain unassigned."))
    g.annotate(SourceDiscrepancy("phosphorylate", "PNKT wash-buffer concentration differs between deposited sources",
        (SourceClaim("2015 article JATS line 201 and original Document S1 text line 2502", "The PNKT bead wash contains 0.1% Triton X-100."),
         SourceClaim("GSE60358 original extraction metadata, repeated mNET-seq paragraph", "The PNKT bead wash contains 0.05% Triton X-100.")),
        "The graph preserves the wash and on-bead phosphorylation order without selecting a Triton concentration.",
        "Defining article, original supplemental methods and deposited GEO extraction metadata, independently read.",
        "The actual experimental wash concentration is unresolved; no typographical error or dose correction is assumed."))
    g.declare_reviewed_scope(("native_chromatin","mnase","mnase_digest","egta","stop_mnase","fractionation","recover_soluble",
        "polii_beads","immunoprecipitate","rna_kinase","atp","phosphorylate","rna_extraction","extract_rna","rna_size_rule","select_rna",
        "adapter3","ligase3","ligate3","adapter5","ligase5","ligate5","rt_primer","rt_enzyme","reverse_transcribe",
        "pcr_forward","pcr_index","pcr_enzyme","library_pcr","library"),SOURCE)
    return g.require_publishable()

"""Molecular model for LAST-seq (Lyu and Chen, 2023)."""
from __future__ import annotations

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import (Construct, MolecularState, Row, Scene, Segment, Workflow,
                      annotation_rows, feature, oligo, complement_segments)
from pcr import PCRGeometry, PCRPrimer
from chemdraw import revcomp
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role
from topology import HairpinTopology


def seg(name, bases, tag=None, **kw):
    return Segment(name, bases, tag, **kw)


T7_PROMOTER = "TAATACGACTCACTATAGGG"
RANDOM_RT_HANDLE = "TACACGACGCTCTTCCGATCT"
RANDOM_RT = RANDOM_RT_HANDLE + "NNNNNN"
DTRU_LINKER = "CCACCTTTCATTCACCCT"
DTRU_DNA = "TTTTTT"
DTRU_RNA = "U" * 19
DTRU_TAIL = DTRU_LINKER + DTRU_DNA + DTRU_RNA
HAIRPIN_LINKER_ARM = "GGGTGAATGAAAGGTGG"
HAIRPIN_CORE = "CTCTAATACGACTCACTATAGGGAGACGTGTGCTCTTCCGATCT"
HAIRPIN_CLOSING_ARM = "TCGTATTAGAG"
CELL_BARCODE_LEN = 6
UMI_LEN = 8
BC2 = "GTATCC"
BC2_HAIRPIN = HAIRPIN_LINKER_ARM + HAIRPIN_CORE + "N" * UMI_LEN + BC2 + HAIRPIN_CLOSING_ARM
BC2_HAIRPIN_TOPOLOGY = HairpinTopology(
    "LAST-seq BC2 hairpin module", BC2_HAIRPIN,
    first_stem_start=len(HAIRPIN_LINKER_ARM),
    second_stem_start=len(BC2_HAIRPIN) - len(HAIRPIN_CLOSING_ARM),
    stem_length=len(HAIRPIN_CLOSING_ARM), colour="var(--c-cbc)")

# Additional file 2, Table S3 library primers (NEBNext1 shown; index window CGTGAT).
NEBNEXT_UNIVERSAL = il.NEBNEXT_UNIVERSAL_PRIMER          # printed with a 3' phosphorothioate
NEBNEXT_INDEX_5 = "CAAGCAGAAGACGGCATACGAGAT"
NEBNEXT_INDEX_3 = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"  # printed ...GATC*T
NEBNEXT_INDEX_LEN = 6
LIB_P5 = "AATGATACGGCGACCACCGAGATCTACAC"
LIB_P7 = "CAAGCAGAAGACGGCATACGAGAT"


def _overlap_3prime_site(primer: str, site: str) -> int:
    """Longest primer 3' end found verbatim in ``site``."""
    for k in range(len(primer), 0, -1):
        if primer[-k:] in site:
            return k
    return 0


# Computed annealing of the Table S3 primers on the RT-product ends: the NEBNext
# primers have 21-nt 3' sites on the RandomRT handle / hairpin PCR handle, whereas the
# short P5/P7 primers have only 5 and 3 nt there.  Methods name "P5/P7 primers" for
# the first 8-cycle PCR; see CAVEAT.
NEB_SITE_NT = (_overlap_3prime_site(NEBNEXT_UNIVERSAL, RANDOM_RT_HANDLE),
               _overlap_3prime_site(NEBNEXT_INDEX_3, HAIRPIN_CORE))
SHORT_SITE_NT = (_overlap_3prime_site(LIB_P5, RANDOM_RT_HANDLE),
                 _overlap_3prime_site(LIB_P7, HAIRPIN_CORE))
if min(NEB_SITE_NT) < 18 or max(SHORT_SITE_NT) >= 10:
    raise ValueError("LAST-seq PCR primer site computation changed: %r %r"
                     % (NEB_SITE_NT, SHORT_SITE_NT))

CELL_FEATURE = feature("last_cell_barcode", "cell_barcode", "whitelist",
                       whitelist="Supplementary Table S3")
UMI_FEATURE = feature("last_umi", "umi", "random")


def hairpin_parts():
    before=HAIRPIN_LINKER_ARM+HAIRPIN_CORE
    return [seg("hairpin fixedDNA",before),seg("UMI", "N"*8,"umi",placeholder=True,feature=UMI_FEATURE),seg("cell barcode",BC2,"cbc",feature=CELL_FEATURE),seg("closingarm",HAIRPIN_CLOSING_ARM)]

def tail_parts():
    return [seg("sourceDNA linker",DTRU_LINKER),seg("dT6",DTRU_DNA),seg("rU19",DTRU_RNA,note="19 physicalRNA uridines"),seg("physicalterminalddC","C",note="source/3ddC/ addsphysicalC; no3′OH")]

def capture_primer():return [*hairpin_parts(),*tail_parts()]

def capture_scene():
    sc=Scene();sc.strand("primer",capture_primer(),mod3="ddC-blocked",label="actual covalent assembledcaptureprimer; nativehairpinfold validatedseparately")
    sc.strand("RNA",[role("independentRNAbody; physicalU"),role("originalRNApolyAtail; annealing/nickcoordinateunassigned")],label="originalRNA; no RT beforeIVT")
    sc.strand("genomicDNA",[role("possible independent genomicDNA/othercellularpolymers; puritynotasserted")])
    return sc

def primer_ligation_geometry():
    sc=Scene();sc.strand("hairpin",[seg("hairpin linker arm",HAIRPIN_LINKER_ARM),role("remaining foldedhairpin; native3′OH nickend")])
    sc.anneal("rU-dTmodule",[seg("paired linker",DTRU_LINKER[:-1]),seg("unpairedlinkerT",DTRU_LINKER[-1]),seg("dT6",DTRU_DNA),seg("rU19",DTRU_RNA),seg("terminalddC","C")],to="hairpin",pair=("paired linker","hairpin linker arm"),unpaired=("unpairedlinkerT","dT6","rU19","terminalddC"),mod5="p",mod3="ddC")
    sc.ligation_end("rU-dTmodule","paired linker",kind="splinted",polarity="5-prime",owner="adapter",text="sourceEcoliDNAligase seals hairpin3′OH→module5′P nick")
    return AdapterLigationGeometry.known("splinted",sc,statement="SourceTableS3 gives17bp overlap plusoneunpairedlinkerdT; originalFigS3 showsEcoli nickligation. AdditionalphysicalterminalddC blocks3′extension, not an inertannotation.")

def promoter_scene():
    sc=Scene();p=seg("T7promoter",T7_PROMOTER,"t7");sc.strand("nontemplate",[p],label="sourcepromoterDNA inhairpin")
    sc.anneal("template",[seg("promotercomplement",revcomp(T7_PROMOTER))],to="nontemplate",pair=("promotercomplement","T7promoter"),label="new DNAextension copieshairpinpromoter")
    sc.strand("templatecontext",[role("originalRNA5′body/polyA remnant covalentlyjoined tonew3′DNAextension; nick/extensionextentunassigned"),seg("copiedUMI","N"*8,"umi",placeholder=True,feature=UMI_FEATURE),seg("copiedBC2",revcomp(BC2),"cbc",feature=CELL_FEATURE)],label="context ofsameRNA–DNAtemplate; diagramnot a secondequivalenttemplate or fullfold")
    sc.strand("originalprimer",capture_primer(),mod3="ddC",label="independent originalhairpin/rU-dT captureprimer")
    sc.strand("residuals",[role("RNaseH-cleaved originalRNA3′fragments, carrier/spikeRNA, residualprimerRNA and independentgenomicDNA; hydrolysisnotphysicalremoval")])
    return sc

def arna_scene():
    sc=Scene();sc.strand("aRNA",[role("T7initiatedRNA5′footprint fromsourcehairpin; nativeendpointunassigned"),seg("RNAcopiedUMI","N"*8,"umi",placeholder=True,feature=UMI_FEATURE),seg("RNAcopiedBC2",BC2,"cbc",feature=CELL_FEATURE),seg("RNAcopiedclosingarm",HAIRPIN_CLOSING_ARM.replace("T","U")),role("RNAcopyoflinker/polyAcontext and antisensetranscript; physicalU, junctionextentunassigned")],label="antisenseRNA copies; linear amplification")
    sc.strand("originals",[role("originalRNA–DNAtemplate/promoter and physicalddCcaptureprimer persist; not consumedbytranscription")])
    sc.strand("residuals",[role("carrier/spike/otherRNA, independentgenomicDNA/possibleofftargetcopies and primer/RNaseHfragments; unknownclassrecovery")])
    return sc

def random_rt_scene():
    sc=Scene();sc.strand("cDNA",[seg("RandomRT handle",RANDOM_RT_HANDLE),seg("randomN6","N"*6,placeholder=True,note="randomprimingpositions, not a newUMI"),role("actualssDNAcopyofantisenseRNA; nativeendpointunassigned"),role("copiedsourceUMI; orientationunassigned",length_bp=8,molecular_feature=UMI_FEATURE),role("copiedsourcecellbarcode; orientationunassigned",length_bp=6,molecular_feature=CELL_FEATURE)],label="actualincoming firstssDNA, not anannealedunextendedRTprimer")
    sc.strand("aRNA",[role("independent originalaRNAtemplates; physicalU; noRNAremovaloperation assigned")])
    sc.strand("template_residuals",[role("possible originalssRNA–DNAtemplate/captureprimer/carrier/spike/cleavedRNA/oligos and independentgenomicDNA/possibleofftargetcopies; physicalrecoveryunassigned")])
    return sc


FINAL_LIBRARY = Construct([
    seg("P5", LIB_P5, "p5"),
    seg("Read 1 arm", NEBNEXT_UNIVERSAL[len(LIB_P5):], "r1"),
    seg("random-primer N6", "N" * 6, placeholder=True, length_bp=6),
    seg("aRNA-derived cDNA", "X" * 38, placeholder=True),
    role("T-Acontext; sourceFigureS8 labels26bp, exactnickposition notassigned"),
    seg("rU-dT linker complement", revcomp(DTRU_LINKER[:-1]), "tso"),
    seg("hairpin closing-arm complement", revcomp(HAIRPIN_CLOSING_ARM)),
    seg("cell barcode", revcomp(BC2), "cbc",
        feature=CELL_FEATURE),
    seg("UMI", "N" * UMI_LEN, "umi", placeholder=True, feature=UMI_FEATURE),
    seg("Read 2 arm (NEBNext index primer)", revcomp(NEBNEXT_INDEX_3), "r2"),
    seg("i7 reverse complement", "I" * NEBNEXT_INDEX_LEN, "cbc", placeholder=True,
        feature=feature("last_sample_index", "sample_index", "whitelist",
                        whitelist="NEBNext1-9 and NEBNext12 index primers, Table S3")),
    seg("P7 reverse complement", revcomp(LIB_P7), "p7"),
], name="LAST-seq library")

SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
_errors = sp.verify(FINAL_LIBRARY, SEQ_PRIMERS,
                    required_roles=tuple(p.role for p in SEQ_PRIMERS))
if _errors:
    raise ValueError("LAST-seq final library: " + "; ".join(_errors))


def final_rows():
    sc = Scene.duplex(list(FINAL_LIBRARY), label="LAST-seq",unpaired=tuple(p.name+"'" for p in FINAL_LIBRARY if p.is_role_token()))
    sc.labels("top")
    return [*sc.rows(), *annotation_rows(FINAL_LIBRARY)]




TITLE = "LAST-seq — direct linear amplification of single-stranded RNA"
NOTES = "01_last-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1186/s13059-023-03025-5">'
          'Lyu and Chen, <i>Genome Biology</i> (2023)</a>.')
SUMMARY = ("A ligated hairpin/rU-dT primer installs a short duplex T7 promoter at an mRNA’s "
           "3′ end. T7 polymerase then amplifies the original RNA directly, before any reverse transcription.")
CAVEAT = ("Published PCR-stage ambiguity: Methods names shortP5/P7 primers for8-cycle firstPCR, then furtherPCR, while sourceFigureS8 displays NEBNextUniversal/Index endpoints. TableS3 suppliesboth families. Source shortP5/P7 have only%d/%d nt3′matches onknownRT/hairpinends versus%d/%d forNEBNext. This graph preservesnamedP5/P7 reagentinput and actualssDNAprecursor while leaving firstcycle/productassembly unavailable; it doesnot silentlysubstituteNEBNext asactualfirstPCR. The unassignedsecondPCR is not assigned to a particularfamily/cyclecount. Endpointreference followsoriginalFigS8 without claimingstageidentity. SourceBC5 andBC7 bothprintCAGGAA despite80selectedbarcodes; BC2(GTATCC) is the selectedexample, no correctioninvented. " % (SHORT_SITE_NT+NEB_SITE_NT))
FINAL_CAPTION = ("Representative PCR-completed LAST-seq molecule. Read 2 reads the hairpin-derived "
                 "UMI, then the cell barcode, then the 11-nt hairpin closing arm; Read 1 enters "
                 "the randomly primed aRNA-derived insert. Native completed endpointreference fromsourceFigS8; exactfirstPCRproductisnotasserted. The P7 side is drawn from the "
                 "NEBNext index primer (6-nt index) printed in Table S3; the poly(A) length is variable.")
SEQUENCING_INTRO = ("NextSeq 550, High Output v2.5 kit: Read 1, 50 cycles; Read 2, 25 cycles. "
                    "The paper reports 6-nt cellular barcodes and 8-nt UMIs; index-read cycles "
                    "are not stated.")
READ_LENGTHS = {"Read 1": 50, "Read 2": 25}
SHOW_READ_CONTENT = False
READ_SPECS=(sp.read_spec("Read 1",50,"aRNA-derived transcriptDNA; sourceFigS8compatible endpointreference","include"),sp.read_spec("Read 2",25,"sourcehairpin UMI8 + cellbarcode6 + closingarm11","extract"),sp.read_spec("Index 1 (i7)","unassigned","native6ntNEBNextindexfamily; acquiredcyclesunassigned","undetermined"),sp.read_spec("Index 2 (i5)","unassigned","nativei5presence/acquisitionunassigned","undetermined"))


def oligos():
    yield oligo("rU-dT module", [seg("linker", DTRU_LINKER, "tso"),
                                  seg("dT6", DTRU_DNA), seg("rU19", DTRU_RNA), seg("physicalterminalddC","C")],
                mods="5′ phosphate; 3′ ddC")
    yield oligo("NEBNext Universal (Table S3)",
                [seg("P5", LIB_P5, "p5"),
                 seg("Read 1 arm", NEBNEXT_UNIVERSAL[len(LIB_P5):], "r1")],
                mods="3′-terminal phosphorothioate (printed C*T)")
    yield oligo("NEBNext index primer (Table S3)",
                [seg("P7", NEBNEXT_INDEX_5, "p7"),
                 seg("i7 index", "I" * NEBNEXT_INDEX_LEN, "cbc", placeholder=True),
                 seg("Read 2 arm", NEBNEXT_INDEX_3, "r2")],
                mods="3′-terminal phosphorothioate (printed C*T)")
    yield oligo("P5 (Table S3)", [seg("P5", LIB_P5, "p5")])
    yield oligo("P7 (Table S3)", [seg("P7", LIB_P7, "p7")])
    yield oligo("RandomRT primer6", [seg("Read 1 handle", RANDOM_RT_HANDLE, "r1"),
                                      seg("N6", "N" * 6, placeholder=True)])


def native_pcr_state(stage):
    sc=Scene();sc.strand("PCRproducts",[role("source"+stage+"PCR products; actualfirstcycle/primerfootprints unassigned becauseMethods/TableS3/FigS8 stageidentityisambiguous"),role("retainedsourceUMI; strandextentunassigned",length_bp=8,molecular_feature=UMI_FEATURE),role("retainedsourcewellbarcode; strandextentunassigned",length_bp=6,molecular_feature=CELL_FEATURE)],label="actualnative PCRproductpopulation")
    sc.strand("originals",[role("possible originalssDNA/aRNA, RNA–DNAtemplate/promoter, physicalcaptureprimer/ddC, carrier/spikeRNA, independentgenomicDNA/possibleofftargetcopies; classrecoveryunassigned")])
    return sc


def protocol_graph():
    from protocol_graph import Terminal,Disposition
    g=ProtocolGraph(TITLE,authority=SOURCE)
    M,I,K=MaterialKind,InputRole,StepKind
    def inp(i,l,k=M.CONDITION,state=None):return g.input(i,l,k,state=state)
    def op(i,l,k,x,extra,state,kind=M.MOLECULAR):
        y=g.product(i+"_out",l,kind,state=state);g.react(i,l,k,((x,I.TEMPLATE),*extra),(y,));return y
    def clean(i,l,x):
        sel=inp(i+"_selection",l,M.SELECTION);y=g.product(i+"_out",l+"retainedmixture",x.kind,state=x.state);w=g.product(i+"_waste","physicalnonretainedfraction",M.WASTE)
        g.react(i,l,K.CLEANUP,((x,I.TEMPLATE),(sel,I.SELECTION)),(y,w));g.discard(w,"Actualphysicalnonretainedfraction; no completeclassclearanceasserted.");return y
    h=inp("hairpin","originalBC2hairpinDNA; N8UMI/6ntwellbarcode",M.OLIGO,Construct(hairpin_parts()))
    tailstate=Scene();tailstate.strand("rU-dT",tail_parts(),mod5="p",mod3="ddC")
    tail=inp("ru_dt","source5′P mixedDNA/rU19/physicalterminalddCmodule",M.OLIGO,tailstate)
    lig=inp("ecoli_ligase","EcoliDNAligase16°C30min",M.CATALYST)
    assembled=Scene();assembled.strand("primer",capture_primer(),mod3="ddC-blocked",label="assembledsourcecaptureprimer; no cellularRNA introduced yet")
    capture=g.product("captureprimer","assembledcaptureprimer",M.OLIGO,state=assembled)
    g.react("assemble","Anneal andnickligate originalmodules",K.LIGATION,((h,I.TEMPLATE),(tail,I.STRUCTURAL),(lig,I.CATALYST)),(capture,),geometry=primer_ligation_geometry())
    capture=clean("primer_page","10%ureaPAGE;300µMNaOAc sourceunit;ethanolprecipitation",capture)
    sample=Scene();sample.strand("RNA",[role("singlecellRNA; physicalU andnativepolyA")]);sample.strand("other",[role("independent cellularDNA/otherpolymers")])
    cell=g.input("cell","onecell in2µllysisbuffer; selected10wellpoolroute",M.BIOLOGICAL,state=sample,primary=True)
    carrier=inp("carrier","60pg physicalRNAcarrier",M.MOLECULAR,Construct([seg("sourcecarrierRNA","AUCCGUCAGGUUUUCCUGAC")]))
    spike=inp("spike","independentERCCRNA source1:1000000; A60is separatealternative",M.MOLECULAR)
    ly=inp("lysis","0.1%Triton/RNaseinhibitor/dNTP;65°C5min/ice")
    lysed=g.product("lysed","nativeRNA plusindependentcarrier/spike andcaptureprimer",M.MOLECULAR,state=capture_scene())
    lysed.state.strand("carrier",[role("independent60pgsourcecarrierRNA; no covalentjoining")]);lysed.state.strand("spike",[role("independentERCCRNA; samecaptureoperations applytoeligibleRNA")])
    g.react("lyse","Lyse independentlybarcodedcellwell",K.THERMAL_TREATMENT,((cell,I.TEMPLATE),(capture,I.STRUCTURAL),(carrier,I.SPIKE_IN),(spike,I.SPIKE_IN),(ly,I.CONDITION)),(lysed,))
    dry=inp("evaporation","Vacufuge30°C25min");x=op("dry","Evaporate eachsourcewelllysate",K.THERMAL_TREATMENT,lysed,((dry,I.CONDITION),),lysed.state)
    rh=inp("rnase_h","RNaseH0.02U/µl; nicksoriginalRNA withinRNA/dT hybrid",M.CATALYST);kl=inp("klenow","Klenowexo−0.3U/µl plusdNTP;37°C30min",M.CATALYST)
    x=op("install_promoter","RNaseHnicking plusDNAextension oforiginalRNA3′end; noRT",K.EXTENSION,x,((rh,I.CATALYST),(kl,I.CATALYST)),promoter_scene())
    peers=inp("other_wells","nineother independent sourcebarcodedcelllysates/promotertemplates",M.MOLECULAR,promoter_scene());pooled=g.product("pooled","tenindependentwellmixtures pooled; not concatenated",M.MOLECULAR,state=x.state)
    g.react("pool","Physicallypool10barcodedlysates;80wellroute is separate",K.POOL,((x,I.TEMPLATE),(peers,I.TEMPLATE)),(pooled,),changes_molecule=False)
    x=clean("pool_cleanup","0.5×Alinebeads+1.3×buffer;rehydratebeadsin25µlIVTmix",pooled)
    t7=inp("t7","T7RNApolymerase2.5U/µl",M.CATALYST);ntp=inp("ivt_mix","sourceNTP1.25mM/Mg/DTT;37°Covernight≤16h")
    x=op("ivt","DirectlinearIVT from originalRNA–DNApromotertemplate",K.TRANSCRIPTION,x,((t7,I.CATALYST),(ntp,I.CONDITION)),arna_scene())
    x=clean("arna_cleanup","0.8×beads+1×buffer;7µlwater",x)
    rp=inp("random_rt","RandomRTprimer6",M.OLIGO,Construct([seg("RandomRT",RANDOM_RT)]));rt=inp("ssiv","SSIV;23°C10min50°C15min80°C10min",M.CATALYST)
    x=op("rt","Random-prime actualaRNAcopies; RTafterIVT",K.REVERSE_TRANSCRIPTION,x,((rp,I.STRUCTURAL),(rt,I.CATALYST)),random_rt_scene())
    p5=inp("methods_p5","MethodsfirstPCRP5;TableS3shortP5 stageapplicabilityunresolved",M.OLIGO,Construct([seg("shortP5",LIB_P5)]));p7=inp("methods_p7","MethodsfirstPCRP7;TableS3shortP7 stageapplicabilityunresolved",M.OLIGO,Construct([seg("shortP7",LIB_P7)]));pol=inp("q5","Q5HighFidelityMasterMix",M.CATALYST)
    amp=g.product("first_pcr","source8cyclePCRoutput; firstcycleassemblyunassigned",M.MOLECULAR,state=native_pcr_state("first"))
    g.react("first_pcr_step","Methods8cyclePCR withnamedP5/P7; sourcephaseconflictbounded",K.PCR,((x,I.TEMPLATE),(p5,I.STRUCTURAL),(p7,I.STRUCTURAL),(pol,I.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(x.state,reason="ActualrandomRTssDNA input; TableS3shortP5/P7 lack fullnativefirstcyclesites. Sourceknown8cyclePCRpreserved without silentlyusingNEBNext endpointsastemplate or assertedreactionreagents."))
    x=clean("pcr1_cleanup","0.9×Alinebeads",amp)
    f=inp("further_forward","sourcefurtherPCRforwardprimer; familyassignmentunresolved",M.OLIGO);r=inp("further_reverse","sourcefurtherPCRreverseprimer; familyassignmentunresolved",M.OLIGO)
    lib=g.product("amplified_library","nativefurtherPCR librarymixture; endpointreference separate",M.LIBRARY,state=native_pcr_state("further"))
    g.react("further_pcr","FurtherPCR; cycles/primerstageidentity unassigned",K.PCR,((x,I.TEMPLATE),(f,I.STRUCTURAL),(r,I.STRUCTURAL),(pol,I.CATALYST)),(lib,),geometry=PCRGeometry.unavailable(x.state,reason="ActualrecoveredfirstPCRmixture enterssourcefurtherPCR. Nativeassemblyphase remainsbounded; sourceFigS8finishedNEBNextcomponentreference isnot actualinputgeometry."))
    x=clean("final_cleanup","Sourcefurtherbeadpurification;11µl0.1×TE",lib)
    g.sequence(x,SequencingPlan(acquired_roles=("Read 1","Read 2"),read_lengths={"Read 1":50,"Read 2":25},read_uses={"Read 1":"include","Read 2":"extract"},limitation="SourceNextSeq55050/25mask; R2sourcehairpin8ntUMI+6ntwellbarcode+11ntclosingarm; native6ntNEBindexfamily separatefromunassignedindexacquisition. Physicalcarrier/spike/originalRNA do not disappear by nonamplification or computationalQC."))
    g.declare_reviewed_scope(REVIEWED,SOURCE)
    return g.require_publishable()
REVIEWED=("hairpin","ru_dt","ecoli_ligase","assemble","primer_page_selection","primer_page","cell","carrier","spike","lysis","lyse","evaporation","dry","rnase_h","klenow","install_promoter","other_wells","pool","pool_cleanup_selection","pool_cleanup","t7","ivt_mix","ivt","arna_cleanup_selection","arna_cleanup","random_rt","ssiv","rt","methods_p5","methods_p7","q5","first_pcr_step","pcr1_cleanup_selection","pcr1_cleanup","further_forward","further_reverse","further_pcr","final_cleanup_selection","final_cleanup","final_cleanup_out")

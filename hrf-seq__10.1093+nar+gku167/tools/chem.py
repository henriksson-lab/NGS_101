"""Protocol-owned source, oligos, molecular workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, feature, revcomp
from circular import CircularizationGeometry
import illumina as il
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, Disposition, Terminal
from protocol_model import molecule, role
import seqprimers as sp


HRF_RT_RANDOM = "AGACGTGTGCTCTTCCGATCT" + "N" * 8 + "S"

HRF_RT_CASSETTE = "AGACGTGTGCTCTTCCGATCTGAACCGGACCGAAGCCCG"

HRF_LIGATION_ADAPTER = "N" * 7 + "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

HRF_PCR_FORWARD = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCT"

HRF_INDEXES = ("GGAACT", "GGACGG", "CGTACG", "GCTACC")

def _hrf_reverse(index: str) -> str:
    return il.P7 + index + il.TRUSEQ_READ2


BARCODE_FEATURE=feature("hrf_molecular_barcode","umi","random")

def rt_product_scene():
    sc=Scene();sc.strand("cDNA",[
        Segment("RT handle",HRF_RT_RANDOM[:-9]),
        Segment("random priming coordinate","N"*8+"S",placeholder=True,length_bp=9),
        Segment("break-terminated cDNA","X"*42,placeholder=True)],
        label="ribosome route: tracked single-strand cDNA target component",mod3="OH")
    sc.ligation_end("cDNA","break-terminated cDNA",kind="single-strand",polarity="3-prime",owner="substrate",text="cDNA3′OH")
    return sc

def ligation_adapter_scene():
    sc=Scene();sc.strand("adapter",[
        Segment("7-nt molecular barcode","N"*7,placeholder=True,length_bp=7,feature=BARCODE_FEATURE),
        Segment("ligation adapter arm",HRF_LIGATION_ADAPTER[7:])],label="LIGATION_ADAPTER_RB",mod5="P",mod3="3NHC3")
    sc.ligation_end("adapter","7-nt molecular barcode",kind="single-strand",polarity="5-prime",owner="adapter",text="source5′phosphate")
    sc.note("adapter","3′aminomodifier is a nonnucleotide terminal block, not a copied cDNA base")
    return sc

def ligated_cdna_scene():
    sc=Scene();sc.strand("ligated cDNA",[*rt_product_scene().strands["cDNA"].segs,
        *ligation_adapter_scene().strands["adapter"].segs],label="actual single-stranded cDNA–adapter precursor",mod3="3NHC3")
    sc.junction("ligated cDNA","break-terminated cDNA","7-nt molecular barcode","cDNA3′OH–adapter5′P ligation")
    return sc

HRF_PCR_TEMPLATE=ligated_cdna_scene()
HRF_PCR_GEOMETRY=PCRGeometry.unavailable(HRF_PCR_TEMPLATE,
    reason="Tracked PCR target is the blocked3′ single-stranded RT/ligation product; other reaction material, including residual RNA, has separately unresolved recovery. Source primers are known, but the forward-first copying and subsequent reverse-primer phase are not implemented as native first-cycle geometry.")

HRF_FINAL_LIBRARY = Construct((
    Segment("P5 / partial Read 1 PCR primer", HRF_PCR_FORWARD, "p5"),
    Segment("remaining Read 1 site", il.TRUSEQ_READ1[-10:], "r1"),
    Segment("7-nt molecular barcode", "N" * 7, "umi", placeholder=True,
            feature=BARCODE_FEATURE),
    Segment("accessibility-defined cDNA", "X" * 42, placeholder=True),
    Segment("random-primer coordinate reverse complement", "S" + "N" * 8, placeholder=True,length_bp=9),
    Segment("Index 1 primer site", revcomp(il.TRUSEQ_READ2), "r2"),
    Segment("sample index", revcomp(HRF_INDEXES[0]), "cbc",
            feature=feature("hrf_sample_index", "sample_index", "whitelist",
                            whitelist="Kielpinski and Vinther Supplementary Table 1")),
    Segment("P7 reverse complement", il.P7_RC, "p7"),
), name="HRF-Seq representative indexed library")

HRF_SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])


TITLE = 'HRF-Seq — hydroxyl-radical RNA accessibility'
NOTES = '01_hrf_seq.html'
SOURCE = 'Kielpinski and Vinther, Nucleic Acids Research (2014), <a href="https://doi.org/10.1093/nar/gku167">doi:10.1093/nar/gku167</a>.'
SUMMARY = 'Hydroxyl radicals cleave solvent-exposed ribose backbones; random-primed RT, a seven-base molecular barcode and paired-end sequencing recover both break and priming coordinates.'
CAVEAT = 'The selected ribosome/random-RT route probes solvent-exposed ribose. RNaseP/cassette-specific RT is a distinct source variant. Actual ligated cDNA is single stranded and3′aminoblocked; first-cycle PCR geometry remains unavailable rather than a manufactured duplex.'
FINAL_LIBRARY = HRF_FINAL_LIBRARY
PCR_GEOMETRY = HRF_PCR_GEOMETRY
FINAL_CAPTION = 'Representative indexed member of the published four-primer family. Read 1 begins with the seven-base molecular barcode; Read 2 reports the random-priming end.'
SEQ_PRIMERS = HRF_SEQ_PRIMERS
ACQUIRED_ROLES = ('Read 1', 'Read 2')
READ_LENGTHS = {'Read 1': 100, 'Read 2': 100}
READ_USES = {'Read 1': 'extract_then_include', 'Index 1 (i7)': 'demultiplex', 'Read 2': 'include'}
SHOW_READ_CONTENT = False
READ_SPECS = (
    sp.read_spec('Read 1', 100, 'Seven-base molecular barcode followed by accessibility-defined cDNA', 'extract_then_include'),
    sp.read_spec('Index 1 (i7)', 'unassigned', 'Physical six-base sample index; acquisition and cycle count unassigned', 'undetermined'),
    sp.read_spec('Read 2', 100, 'Random-primer coordinate followed by cDNA', 'include'),
)
READ_SPECS_NOTE = 'Methods establish paired100. A compatible index-primer site and a physical sample index do not establish whether an index read was acquired.'
SEQUENCING_INTRO = 'The source uses standard paired-end Illumina primers. Their sites, the first-read molecular barcode, sample index and second-read random-primer coordinate are derived from the published oligos.'
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''

def oligos():
    return ('RT_random_primer — 5′-AGACGTGTGCTCTTCCGATCTNNNNNNNNS-3′', 'RT_structure_cassette — 5′-AGACGTGTGCTCTTCCGATCTGAACCGGACCGAAGCCCG-3′', 'LIGATION_ADAPTER_RB — 5′-PHO-NNNNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT-3NHC3-3′', 'PCR_forward — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCT-3′', *tuple((f'PCR reverse {index} — 5′-{_hrf_reverse(index)}-3′' for index in HRF_INDEXES)))

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    start=g.input("input","folded E.coli ribosomes; random-RT route selected",MaterialKind.BIOLOGICAL,state=molecule("ribosomal RNA",(role("folded ribosomal RNA"),)),primary=True)
    treated=g.product("treated","hydroxyl-radical probing aliquot",MaterialKind.MOLECULAR,state=start.state)
    control=g.product("control","parallel water-control aliquot",MaterialKind.MOLECULAR,state=start.state)
    g.react("split","Prepare probing and untreated water-control samples",StepKind.SPLIT,((start,InputRole.TEMPLATE),),(treated,control),changes_molecule=False)
    g.mark_terminal(Terminal(control.id,Disposition.UNRESOLVED,"Parallel6µLwater control follows the same library route; repeated physical library graph is not drawn, not a discard."))
    radicals=g.input("fenton","peroxidative Fenton reagents; three2µLdroplets:5mM ferrous ammonium sulfate–EDTA,50mM ascorbate,1.5%H2O2",MaterialKind.CONDITION)
    cut=g.product("cleaved","backbone-cleaved RNA population",MaterialKind.MOLECULAR,state=molecule("hydroxyl-radical RNA products",(role("cleaved RNA; chemical end classes unresolved"),)))
    g.react("cleave","Mix Fenton reagents for60s to probe solvent-exposed ribose",StepKind.FRAGMENTATION,((treated,InputRole.TEMPLATE),(radicals,InputRole.CONDITION)),(cut,))
    precip=g.input("ethanol","318µLice-cold ethanol+10µg glycogen;−80°C30min,centrifuge,resuspend12.5µLwater",MaterialKind.SELECTION)
    probed=g.product("probed","ethanol-recovered probed RNA",MaterialKind.MOLECULAR,state=cut.state)
    supernatant=g.product("probing_supernatant","ethanol precipitation supernatant",MaterialKind.WASTE)
    g.react("recover_probed","Stop probing and recover RNA by ethanol precipitation",StepKind.CLEANUP,((cut,InputRole.TEMPLATE),(precip,InputRole.SELECTION)),(probed,supernatant),changes_molecule=False)
    g.discard(supernatant,"physically excluded after centrifugation; exact small-product recovery unassigned")
    primer=g.input("rt_primer","RT_random_primer; ribosome branch",MaterialKind.OLIGO,state=Construct(list(rt_product_scene().strands["cDNA"].segs)[:2]),source="SupplementaryTableS1")
    rt=g.input("rt_enzyme","PrimeScript;30s25°C,30min42°C,10min each50/56/60°C",MaterialKind.CATALYST)
    cdna=g.product("cdna","break-terminated first-strand cDNA",MaterialKind.MOLECULAR,state=rt_product_scene())
    residual_rna=g.product("rt_template_rna","RNA template population remaining after RT; downstream recovery unassigned",MaterialKind.MOLECULAR,state=probed.state)
    g.react("rt","Random-prime and copy probed ribosomal RNA",StepKind.REVERSE_TRANSCRIPTION,((probed,InputRole.TEMPLATE),(primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,residual_rna))
    g.mark_terminal(Terminal(residual_rna.id,Disposition.UNRESOLVED,"RT copies RNA without physically removing it. The graph follows the cDNA target component through RNAClean; RNA recovery/degradation through that delegated purification and later reactions is unassigned, not an asserted RNA-free tube."))
    beads=g.input("rt_cleanup","RNACleanXP recovery;25µL5mM Na-citrate pH6;200×dilution inwater",MaterialKind.SELECTION)
    recovered=g.product("rt_clean","diluted recovered single-strand cDNA",MaterialKind.MOLECULAR,state=cdna.state)
    rt_waste=g.product("rt_cleanup_waste","RNAClean supernatant/washes",MaterialKind.WASTE)
    g.react("rt_recovery","Recover and dilute ribosome-route cDNA before ligation",StepKind.CLEANUP,((cdna,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(recovered,rt_waste),changes_molecule=False)
    g.discard(rt_waste,"excluded by source bead purification; exact molecular partition not reconstructed")
    adapter=g.input("adapter","5′PHO-LIGATION_ADAPTER_RB-3′NHC3",MaterialKind.ADAPTER,state=ligation_adapter_scene(),source="SupplementaryTableS1")
    ligase=g.input("ligase","CircLigase withATP,MnCl2,PEG6000,betaine;2h60°C,1h68°C,10min80°C",MaterialKind.CATALYST)
    ligated=g.product("ligated","3′adapter-ligated cDNA",MaterialKind.MOLECULAR,state=HRF_PCR_TEMPLATE)
    g.react("ligate","Ligate the5′phosphorylated barcode adapter to cDNA3′OH",StepKind.LIGATION,((recovered,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(ligated,),geometry=AdapterLigationGeometry.known("single-strand",rt_product_scene(),ligation_adapter_scene(),statement="Methods CircLigase and TableS1 identify5′PHO adapter/3′NHC3 block; native cDNA3′OH joins adapter5′P. The nonnucleotide3′block remains terminal."))
    ampure=g.input("ligation_beads","AMPureXP cleanup;16µLwater elution",MaterialKind.SELECTION)
    clean=g.product("pcr_template","purified ligated single-strand cDNA",MaterialKind.MOLECULAR,state=ligated.state)
    lig_waste=g.product("ligation_waste","AMPure supernatant/washes",MaterialKind.WASTE)
    g.react("ligation_cleanup","Purify ligated cDNA",StepKind.CLEANUP,((ligated,InputRole.TEMPLATE),(ampure,InputRole.SELECTION)),(clean,lig_waste),changes_molecule=False)
    g.discard(lig_waste,"physically excluded in bead purification")
    f=g.input("library_forward","PCR_forward",MaterialKind.OLIGO,state=HRF_PCR_FORWARD,source="SupplementaryTableS1")
    r=g.input("library_reverse","PCR_REVERSE_INDEX.14_AGTTCC representative",MaterialKind.OLIGO,state=_hrf_reverse(HRF_INDEXES[0]),source="SupplementaryTableS1")
    pol=g.input("library_polymerase","Phusion;4cycles98/64/72°C then20cycles98/72°C;5min72°C",MaterialKind.CATALYST)
    amp=g.product("amplified","indexed PCR product",MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
    g.react("library_pcr","Amplify and index actual ligated cDNA (4+20cycles)",StepKind.PCR,((clean,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(amp,),geometry=HRF_PCR_GEOMETRY)
    size=g.input("library_size_selection","AMPure,pool;2%E-gel200–600bp;QiagenPCR column;finalAMPure",MaterialKind.SELECTION)
    library=g.product("library_library","purified HRF sequencing library",MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    waste=g.product("library_cleanup_waste","out-of-window gel material and purification flowthrough/supernatants",MaterialKind.WASTE)
    g.react("library_cleanup","Purify,pool and size-select200–600bp library",StepKind.CLEANUP,((amp,InputRole.TEMPLATE),(size,InputRole.SELECTION)),(library,waste),changes_molecule=False)
    g.discard(waste,"physically excluded by source gel/column/bead transfers")
    g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=SEQ_PRIMERS,acquired_roles=ACQUIRED_ROLES,read_lengths=READ_LENGTHS,read_uses=READ_USES,limitation="Methods establish paired100. Six-base sample-index sequence comes from TableS1; exact index-read acquisition is not separately specified."))
    g.declare_reviewed_scope(("input","split","fenton","cleave","ethanol","recover_probed","rt_primer","rt_enzyme","rt","rt_cleanup","rt_recovery","adapter","ligase","ligate","ligation_beads","ligation_cleanup","library_forward","library_reverse","library_polymerase","library_pcr","library_size_selection","library_cleanup","library_library"),SOURCE)
    return g.require_publishable()

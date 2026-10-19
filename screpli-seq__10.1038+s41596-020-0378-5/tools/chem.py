"""Miura2020 source-owned scRepli-seq: independent cells and native library states."""
from types import SimpleNamespace
from chemdraw import Construct, Scene, Segment, feature
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (Disposition, InputRole as I, MaterialKind as M,
                            ProtocolGraph, SequencingPlan, StepKind as K, Terminal)
import seqprimers as sp

TITLE = "scRepli-seq — single-cell DNA replication sequencing"
AUTHORITY = "Miura et al., Nature Protocols2020, doi:10.1038/s41596-020-0378-5; published Procedure17–79, Reagents, WGA overview and Library preparation/sequencing"
SOURCE = '<a href="https://doi.org/10.1038/s41596-020-0378-5">Miura et al., Nature Protocols (2020)</a>; published Procedure17–79 and Reagents.'
SUMMARY = "Individual S-phase and G1 control cells undergo SeqPlex WGA, Eco57I adaptor digestion, KAPA LTP preparation, NEXTflex index-adaptor ligation and TPC1/TPC2 library PCR. Genomic copy number reports replication timing."
CAVEAT = "Native SeqPlex and NEXTflex strands, cut coordinates, residual handles and first-cycle PCR geometry remain unassigned in this bounded scope. The predicted SEQXE and single-index trimming patterns are analysis observations, not synthesized oligos or a measured index width. No canonical TruSeq shell is substituted. Independent original DNA and modified-primer-derived PCR daughters remain distinct; component-specific cleanup recovery is unresolved."
TPC1 = "AATGATACGGCGACCACCGAG"
TPC2 = "CAAGCAGAAGACGGCATACGAG"
def sample_index(label):
    key = {"S-phase cell": "s_phase_index", "G1 control": "g1_index"}[label]
    return feature(key, "sample_index", "unknown", group=key,
                   note="Native ligated NEXTflex index identifies this individual cell library; placement, width and member unassigned")


def role(name, *, molecular_feature=None):
    return Segment(name, "[" + name + "]", placeholder=True, feature=molecular_feature)


def tpc(name):
    seq = {"TPC1": TPC1, "TPC2": TPC2}[name]
    return Construct([Segment(name, seq, note="Source Reagents: phosphorothioate bond immediately before terminal G; HPLC grade")],name=name)


def state(label, stage):
    sc = Scene()
    if stage in {"cell", "fragmented"}:
        sc.strand("genomic", [role(label+" genomic DNA" if stage=="cell" else label+" heat-fragmented genomic DNA; native ends unassigned")],label="genomic population within this one cell lineage")
    elif stage in {"preamp", "wga"}:
        # Source places universal SEQXE adaptor sequence at the ends of WGA products
        # (overview); pre-amplification ends and any end asymmetry are not reported.
        end = ("pre-amplification end; SEQXE presence/structure unassigned" if stage=="preamp"
               else "universal SEQXE adaptor end; native strands unassigned")
        sc.strand("copies", [role("left "+end),role(label+" DNA copies"),role("right "+end)],label="pre-amplification DNA" if stage=="preamp" else "SeqPlex WGA daughter population")
        sc.strand("originals", [role("possible carried original genomic DNA and primers; native recovery unassigned")],label="independent possible carried populations")
    elif stage in {"digested", "repaired", "tailed"}:
        sc.strand("dna", [role("possible residual SEQXE; left end unassigned"),role(label+" WGA genomic fragment"),role("possible residual SEQXE; right end unassigned")],label=stage+" native DNA population")
        sc.strand("carried", [role("possible carried genomic originals and residual primers; class recovery unassigned")],label="independent possible carried originals")
        if stage=="tailed":
            sc.note("dna", "Source on-bead A-tailing is declared; individual strand lengths, overhang bases and phosphorylation geometry are not reconstructed from kit name.")
    elif stage in {"ligated", "pcr"}:
        if stage=="ligated":
            sc.strand("native", [role(label+" native NEXTflex-indexed WGA DNA; full strand ends, junctions, residual SEQXE and index position/width unassigned", molecular_feature=sample_index(label))],label="native ligated DNA; feature belongs to this polymer population, coordinates unassigned")
        else:
            for name in ("TPC1","TPC2"):
                sc.strand(name,list(tpc(name))+[role(label+" extended ordinary DNA copy; native insert/distal end and index position unassigned",molecular_feature=sample_index(label))],label=name+"-initiated physical daughter; terminal primer bond remains phosphorothioate")
            sc.note("TPC1", "Physical supplied primers retain their phosphorothioate bond upon extension. Copying the opposite template makes ordinary DNA and does not copy backbone modifications; these rows do not assert annealed strand geometry.")
        sc.strand("carried", [role("possible carried original DNA, unligated products and residual primers/adaptors; class recovery unassigned")],label="independent possible co-carried molecular populations")
    else:
        raise ValueError(stage)
    return sc


def protocol_graph():
    g=ProtocolGraph("scRepli-seq",authority=AUTHORITY)
    # Explicit source mapping stays local and is not derived from whatever nodes exist.
    reviewed=["sorted_cells","separate_cells"]
    inputs=Scene()
    inputs.strand("s",[role("one sorted S-phase cell genome")],label="independent S-phase cell")
    inputs.strand("g1",[role("one separately sorted G1 control genome")],label="independent G1 control cell")
    cells=g.input("sorted_cells","separately sorted S-phase and G1 cells (Steps17–18)",M.BIOLOGICAL,state=inputs,primary=True)
    s=g.product("s_cell","individual S-phase cell",M.MOLECULAR,state=state("S-phase","cell"))
    g1=g.product("g1_cell","individual G1 control cell",M.MOLECULAR,state=state("G1","cell"))
    g.react("separate_cells","Keep independently sorted cells in separate wells throughout preparation",K.SPLIT,((cells,I.TEMPLATE),),(s,g1),changes_molecule=False)

    def branch(prefix,cell,label):
        def inp(key,text,kind,**kw):
            key=prefix+"_"+key; reviewed.append(key)
            return g.input(key,label+": "+text,kind,**kw)
        def prod(key,text,stage=None,kind=M.MOLECULAR):
            return g.product(prefix+"_"+key,label+": "+text,kind,state=state(label,stage) if stage else None)
        def react(key,text,kind,ins,outs,**kw):
            key=prefix+"_"+key; reviewed.append(key)
            return g.react(key,text+" ("+label+")",kind,ins,outs,**kw)
        def unresolved(x,why):
            g.mark_terminal(Terminal(x.id,Disposition.UNRESOLVED,why))
        def cleanup(key,x,description,stage):
            beads=inp(key+"_reagent",description,M.SELECTION)
            y=prod(key+"_retained","physical retained cleanup fraction",stage)
            waste=prod(key+"_waste","removed supernatant/washes; molecular composition not fully assigned",kind=M.WASTE)
            react(key,"Recover source bead-bound DNA fraction",K.CLEANUP,((x,I.TEMPLATE),(beads,I.SELECTION)),(y,waste),changes_molecule=False)
            g.discard(waste,"Physical supernatant/washes removed. Individual original/daughter/primer class recovery is not inferred.")
            return y
        def aliquot(key,x,amount,stage):
            y=prod(key,"source aliquot "+amount,stage)
            rest=prod(key+"_remainder","remaining source material; quantity and subsequent use unassigned",stage)
            react(key+"_split","Allocate "+amount+" for the next source operation",K.SPLIT,((x,I.TEMPLATE),),(y,rest),changes_molecule=False)
            unresolved(rest,"Remainder is physically separate from selected aliquot; subsequent storage/use unassigned.")
            return y
        lysis=inp("lysis","6µl Single Cell Lysis & Fragmentation Buffer plus proteinase K;55°C1h,99°C4min",M.CONDITION)
        x=prod("fragmented","heat-fragmented genomic DNA","fragmented")
        react("lyse","Lyse and heat-fragment single-cell genomic DNA",K.FRAGMENTATION,((cell,I.TEMPLATE),(lysis,I.CONDITION)),(x,))
        prep=inp("seqplex_prep","SeqPlex Library Preparation Buffer/Enzyme;95°C2min then16/24/37°C20min each,75°C5min",M.CATALYST)
        pre=prod("preamp","native pre-amplification population","preamp")
        react("preamp_step","Perform source SeqPlex pre-amplification",K.MAGIC,((x,I.TEMPLATE),(prep,I.CATALYST)),(pre,),reason="Steps20–22 specify kit preparation but native ordered oligos, priming geometry and intermediates remain unassigned; no final library is substituted.")
        f=inp("wga_forward","native SeqPlex (DOP-PCR-based) priming role A within kit reagents; sequence and number of distinct species unassigned (not asserted distinct from role B)",M.OLIGO)
        r=inp("wga_reverse","native SeqPlex (DOP-PCR-based) priming role B within kit reagents; sequence and number of distinct species unassigned (may be the same species as role A)",M.OLIGO)
        enzyme=inp("wga_polymerase","DNA Polymerase for SeqPlex and5×Amplification Mix",M.CATALYST)
        wga=prod("wga","SeqPlex WGA population","wga")
        react("wga_step","Whole-genome amplify 94°C2min;24cycles94°C15s/70°C5min;70°C30min (Steps23–25)",K.PCR,((pre,I.TEMPLATE),(f,I.STRUCTURAL),(r,I.STRUCTURAL),(enzyme,I.CATALYST)),(wga,),geometry=PCRGeometry.unavailable(pre.state,reason="Actual native pre-amplification input. SeqPlex primer species and first-cycle sites unassigned; end roles do not establish distinct ordered primer species."))
        x=cleanup("wga_cleanup",wga,"1.7×AMPure XP;Steps26–32","wga")
        x=aliquot("digestion_input",x,"500–800ng (Step34)","wga")
        eco=inp("eco57i","Primer Removal Enzyme Eco57I / Buffer / Solution;37°C1h,65°C20min",M.CATALYST)
        dna=prod("digested","Eco57I-digested WGA genomic population","digested")
        cuts=prod("cut_adaptors","released~40bp SEQXE adaptor polymers; exact strands/cut coordinates unassigned")
        react("digest","Digest universal SEQXE adaptors",K.DIGESTION,((x,I.TEMPLATE),(eco,I.CATALYST)),(dna,cuts))
        beads=inp("adaptor_removal_reagent","1.8×AMPure XP;repeat if~40bp fragments remain;Steps36–43",M.SELECTION)
        x=prod("adaptor_removed","retained100–500bp WGA population; residual handles possible","digested")
        waste=prod("adaptor_removal_waste","physical supernatant/washes including digested SEQXE adaptors",kind=M.WASTE)
        react("adaptor_cleanup","Physically deplete short digested adaptors",K.CLEANUP,((dna,I.TEMPLATE),(cuts,I.TEMPLATE),(beads,I.SELECTION)),(x,waste),changes_molecule=False)
        g.discard(waste,"Source short-adaptor depletion; incomplete recovery is allowed (Step43 repeat). Computational SEQXE trimming later does not remove molecules.")
        x=aliquot("end_repair_input",x,"20–40ng (Step44)","digested")
        repair=inp("repair_mix","KAPA LTP End Repair Enzyme Mix;20°C30min",M.CATALYST)
        y=prod("repaired","native end-repaired WGA DNA","repaired")
        react("end_repair","End repair source WGA aliquot",K.END_REPAIR,((x,I.TEMPLATE),(repair,I.CATALYST)),(y,))
        x=cleanup("repair_cleanup",y,"1.7×AMPure XP;retain beads for on-bead A-tailing;Steps46–50","repaired")
        tail=inp("tail_mix","KAPA LTP A-Tailing Enzyme/Buffer on beads;30°C30min",M.CATALYST)
        y=prod("tailed","native A-tailed population; per-strand geometry unassigned","tailed")
        react("a_tail","A-tail DNA on the retained beads",K.TAILING,((x,I.TEMPLATE),(tail,I.CATALYST)),(y,))
        x=cleanup("tail_cleanup",y,"1.8×PEG/NaCl,retain beads;Steps53–57","tailed")
        adapter=inp("adaptor","NEXTflex DNA Barcodes–48 index adaptor,1µl of5µM stock; native strands/index width/member unassigned",M.ADAPTER)
        ligase=inp("ligase","KAPA LTP DNA Ligase;20°C15min",M.CATALYST)
        y=prod("ligated","native NEXTflex-ligated reaction","ligated")
        react("ligation","Ligate native source single-index NEXTflex adaptor",K.LIGATION,((x,I.TEMPLATE),(adapter,I.STRUCTURAL),(ligase,I.CATALYST)),(y,),geometry=AdapterLigationGeometry.unavailable("Step58 source NEXTflex DNA Barcodes–48; native adapter strands/ends and precise junction chemistry remain unassigned. Allowed dual/unique-dual alternatives are separate variants, not defaults."))
        x=cleanup("ligation_cleanup",y,"1.0×PEG/NaCl;elute11µlBufferEB;Steps60–65","ligated")
        q=prod("qpcr_input","1.5µl ligated eluate for real-time cycle determination","ligated")
        bulk=prod("bulk_input","7µl ligated eluate for library PCR","ligated")
        reserve=prod("eluate_remainder","remaining ligated eluate; subsequent use unassigned","ligated")
        react("pcr_allocation","Allocate distinct1.5µl diagnostic and7µl bulk-PCR aliquots from11µl eluate",K.SPLIT,((x,I.TEMPLATE),),(q,bulk,reserve),changes_molecule=False)
        unresolved(reserve,"Unused eluate after the two specified aliquots; no disposal is asserted.")
        pf=inp("tpc1","TPC1:AATGATACGGCGACCACCGA*G;terminal phosphorothioate",M.OLIGO,state=tpc("TPC1"))
        pr=inp("tpc2","TPC2:CAAGCAGAAGACGGCATACGA*G;terminal phosphorothioate",M.OLIGO,state=tpc("TPC2"))
        qmix=inp("qpcr_mix","KK2702 HiFi HotStart real-time mix;TPC primers;standards1–4 measured in separate wells",M.CATALYST)
        qc=prod("qpcr_product","real-time PCR product for cycle determination","pcr",kind=M.QC)
        react("qpcr","Determine cycle count from least cycle exceedingStd2;5–7cycles for bulkPCR",K.PCR,((q,I.TEMPLATE),(pf,I.STRUCTURAL),(pr,I.STRUCTURAL),(qmix,I.CATALYST)),(qc,),geometry=PCRGeometry.unavailable(q.state,reason="Actual1.5µl native ligated aliquot; TPC strands are known but their native NEXTflex sites are unassigned."))
        g.mark_terminal(Terminal(qc.id,Disposition.QC,"Physically separate diagnostic PCR aliquot; not returned to bulk library."))
        mix=inp("pcr_mix","KAPA HiFi HotStart ReadyMix KK2602 7.5µl, 0.5µl 10µM TPC primer mix;98°C2min, n×(98°C15s/60°C30s/72°C30s), 72°C1min; n=5–7, source datasets5or6 (Steps70–71)",M.CATALYST)
        amplified=prod("pcr_product","TPC-initiated library daughters and possible carried originals","pcr")
        react("library_pcr","Amplify separate bulk-library aliquot with TPC1/TPC2",K.PCR,((bulk,I.TEMPLATE),(pf,I.STRUCTURAL),(pr,I.STRUCTURAL),(mix,I.CATALYST)),(amplified,),geometry=PCRGeometry.unavailable(bulk.state,reason="Actual7µl native ligated aliquot; source TPC sequences alone do not establish complete native NEXTflex sites or first-cycle product ends."))
        x=cleanup("post_pcr_cleanup",amplified,"0.9×AMPure XP;repeat if adaptor dimers remain;Steps72–78","pcr")
        library=prod("library","individual index-adaptor-bearing library for pooled sequencing","pcr",kind=M.LIBRARY)
        react("library_recovery","Carry recovered library to sequencing;analytical QC volumes unassigned",K.GENERIC,((x,I.TEMPLATE),),(library,),changes_molecule=False)
        reviewed.append(library.id)
        g.sequence(library,SequencingPlan(acquired_roles=("Read 1","Index 1 (i7)"),read_lengths={"Read 1":80},limitation="Published HiSeq1500 80bp single-read example;Step79 recommends>80bp. This graph selects source NEXTflex–48 single-index variant; native primer sites, index width/member and index cycle count remain unassigned. Separate cell libraries are pooled(~40cells/lane); full pool membership/volumes are unexpanded. Read trimming and G1 normalization are computational, not molecular disposal."))
    branch("s_phase",s,"S-phase cell")
    branch("g1",g1,"G1 control")
    g.declare_reviewed_scope(reviewed,AUTHORITY)
    return g.require_publishable()


def oligos():
    return ("TPC1:5′-AATGATACGGCGACCACCGA*G-3′ (*terminal phosphorothioate;Reagents)",
            "TPC2:5′-CAAGCAGAAGACGGCATACGA*G-3′ (*terminal phosphorothioate;Reagents)")

SPEC=SimpleNamespace(TITLE=TITLE,NOTES="01_screpli_seq.html",SOURCE=SOURCE,SUMMARY=SUMMARY,
    CAVEAT=CAVEAT,FINAL_LIBRARY=None,SEQ_PRIMERS=(),protocol_graph=protocol_graph,oligos=oligos,
    SEQUENCING_ENDING="Native individual NEXTflex–48 libraries;published80bp single-read HiSeq1500 example; index cycles unassigned.",
    READ_SPECS=(sp.read_spec("Read 1",80,"Genomic WGA DNA;possible residual SEQXE;computational trimming then G1-normalized copy number","include"),sp.read_spec("Index 1 (i7)",None,"Ligated sample index identifies one cell library;native width and cycle count unassigned","demultiplex")))

"""Protocol-owned source, oligos, molecular workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Segment, Scene, feature
import illumina as il
from pcr import PCRGeometry
from protocol_graph import InputRole as I, MaterialKind as M, ProtocolGraph, SequencingPlan, StepKind as K, Terminal, Disposition
from protocol_model import role
import seqprimers as sp


# SI p.1 prints the 58-nt A_adapter_t; P5 and the Read 1 site share ACAC, present once.
# (Model correction 2026-10-17: the earlier il.P5 + il.TRUSEQ_READ1 duplicated it, 62 nt.)
SHAPE_A_T = il.P5 + il.TRUSEQ_READ1[4:]

SHAPE_A_B = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG"
SHAPE_PE_READ2 = il.PE_READ2_PRIMER

# SI p.1: RT primers carry a 5' handle RRRY for (+) and YYYR for (-) reactions; an inline
# degenerate condition tag read first in Read 1, not a designed whitelist barcode.
SHAPE_HANDLE_FEATURE = feature("shape_treatment_handle", "inline_barcode", "degenerate",
                               note="RRRY = 1M7-treated (+); YYYR = solvent control (-)")

SHAPE_RNA_FEATURE = feature("shape_rna_identity", "feature_barcode", "degenerate",
                            note="Four-nucleotide degenerate RNA barcode generated via DNA template primers; selected seven-RNA pool is barcoded")

SHAPE_FINAL_LIBRARY = Construct((
    Segment("P5", il.P5, "p5"),
    Segment("A_adapter_t / Read 1 site remainder", SHAPE_A_T[len(il.P5):], "r1"),
    Segment("plus/minus handle", "RRRY", "cbc", placeholder=True,
            feature=SHAPE_HANDLE_FEATURE),
    # SI Fig. S1 / Table S1: the 3' structure cassette ends 5'-linker hairpin-[4-nt bar
    # code]-RT primer binding site-3'.  The RT primer therefore extends straight into
    # the bar code, so the bar code sits next to the RT primer (Read 1 end), and the
    # SHAPE stop lies at the opposite, A_adapter_b-ligated end (main text: the two
    # pieces of information "are on opposite ends of the cDNA molecules").
    role("target-specific RT primer (cassette RT-primer binding site)"),
    role("cDNA of RNA identity barcode", tag="cbc", length_bp=4,
         molecular_feature=SHAPE_RNA_FEATURE),
    role("cDNA of 3-prime linker hairpin"),
    role("cDNA ending before the SHAPE-modified nucleotide"),
    Segment("A_adapter_b sequencing arm", SHAPE_A_B[:-len(il.P7_RC)], "r2"),
    Segment("P7 reverse complement", il.P7_RC, "p7"),
), name="SHAPE-Seq v1 library")


TITLE = 'SHAPE-Seq — multiplexed RNA flexibility probing'
NOTES = '01_shape_seq.html'
SOURCE = 'Lucks et al., PNAS (2011), <a href="https://doi.org/10.1073/pnas.1106501108">doi:10.1073/pnas.1106501108</a>.'
SUMMARY = 'SHAPE chemistry acylates flexible riboses, reverse transcriptase stops at the adducts, and paired-end sequencing links each stop to a pre-attached RNA identity barcode.'
CAVEAT = 'Selected 1M7-treated channel of the seven-RNA in-vitro pool. The separate larger RNaseP library includes an unbarcoded exception and is outside this seven-RNA selection. The physical solvent-control allocation is shown, with its downstream processing and eventual pooling timing outside this selected branch. RT enzyme and target-specific primer, adapter 5′ phosphorylation and PCR-primer identities remain unassigned. Native original RNA/adducts, RNA hydrolysates, blocked adapter and DNA daughters remain distinct populations.'
FINAL_LIBRARY = SHAPE_FINAL_LIBRARY
FINAL_CAPTION = 'Source-exact v1 adapter arms flank a target-specific RT product. Read 1 begins with the four-base treatment handle, then the RT primer and the cDNA copy of the 3′-cassette RNA bar code; Read 2 starts at the A_adapter_b ligation junction and locates the SHAPE-dependent RT stop. Model correction (2026-10-17): an earlier version drew the RNA bar code beside A_adapter_b; SI Fig. S1/Table S1 place it between the 3′ linker and the RT-primer binding site. Target-specific RT-primer bases remain unassigned pending the cited Mortimer procedure. The final duplex is a component/read-site reference; first-cycle PCR assembly remains unresolved. Original adapter C3 is not copied into daughter DNA. The cited Kozarewa method prints its own adapter_b beginning GATC and a phosphorothioate-marked adapter_t; the defining SHAPE SI instead prints AGATC and leaves that RT-tail modification unmarked. This page uses SHAPE-specific oligos and does not transfer predecessor modifications.'
SEQ_PRIMERS = (
    sp.TRUSEQ['R1'],
    sp.custom('Read 2', 'Illumina PE Read 2', SHAPE_PE_READ2,
              'Illumina 11 January 2011 adapter and sequencing-primer letter, p. 6'),
)
READ_LENGTHS = {'Read 1': 50, 'Read 2': 50}
SEQUENCING_INTRO = ''
SEQUENCING_ENDING = 'The defining method amplified the library for 9 or 12 PCR cycles (SI; the main text reports nine) with primers binding the Illumina adapter sequences (bases not printed), followed by AMPure cleanup without size selection, and used 2×50-cycle paired-end sequencing on a Genome Analyzer IIx. Handle binning, handle/adapter trimming, unique zero-mismatch alignment and stop counting are computational operations, not physical separation.'
SEQUENCING_UNAVAILABLE = ''

def oligos():
    return (f'A_adapter_t / PE_forward — 5′-{SHAPE_A_T}-3′', f'A_adapter_b — 5′-{SHAPE_A_B}/3SpC3/-3′ (5′ phosphorylation is inferred for ligation, not printed in the SHAPE-Seq SI)', 'treated RT primer pattern — 5′-[A_adapter_t]-RRRY-[target-specific RT primer]-3′', 'control RT primer pattern — 5′-[A_adapter_t]-YYYR-[target-specific RT primer]-3′')

def rt_parts():
    return [Segment('A_adapter_t tail',SHAPE_A_T),
            Segment('treated handle','RRRY',placeholder=True,feature=SHAPE_HANDLE_FEATURE),
            role('target-specific RT-primer DNA; sequence/length unassigned')]


def cdna_parts():
    return [*rt_parts(), role('DNA copy of RNA identity barcode',length_bp=4,molecular_feature=SHAPE_RNA_FEATURE),
            role('DNA copy of 3-prime linker hairpin'),
            role('DNA extension ending at RT stop; exact coordinate unassigned')]


def adapter_state():
    sc=Scene();sc.strand('adapter',[Segment('A_adapter_b',SHAPE_A_B)],mod3='C3 spacer',label='native ssDNA adapter; 5′ phosphate unverified, terminal C3 is nonnucleotide')
    return sc


def native(stage):
    sc=Scene()
    if stage in ('rna','modified','rt'):
        sc.strand('original_RNA',[role('RNA target; physical ribonucleotides'),role('RNA 3-prime linker hairpin'),role('transcribed RNA identity barcode',length_bp=4,molecular_feature=SHAPE_RNA_FEATURE),role('RNA RT-primer binding site')],label='original RNA; 1M7 2′-O-adducts are substituents, not extra nucleotides' if stage!='rna' else 'folded seven-RNA pool')
    else:
        sc.strand('RNA_hydrolysates',[role('original RNA alkaline hydrolysis products; adduct chemical fate and physical recovery unresolved')],label='RNA chemical degradation does not itself discard material')
    if stage not in ('rna','modified'):
        parts=cdna_parts()
        if stage in ('ligated','library'):parts += [Segment('ligated A_adapter_b',SHAPE_A_B)]
        sc.strand('native_cDNA',parts,mod3='C3 spacer on original adapter' if stage in ('ligated','library') else '',label='actual native single-stranded DNA; RNA adducts not copied')
        if stage in ('ligated','library'):
            sc.junction('native_cDNA','DNA extension ending at RT stop; exact coordinate unassigned','ligated A_adapter_b','proposed 3′ DNA–5′ adapter junction; adapter 5′P unverified')
            sc.strand('unligated_DNA',cdna_parts(),label='possible original unligated stop cDNA; recovery unresolved')
            sc.strand('residual_adapter',[Segment('possible free native A_adapter_b',SHAPE_A_B)],mod3='C3 spacer',label='possible free blocked adapter and off-products; cleanup not class-pure')
        sc.strand('residual_RT_primer',rt_parts(),label='possible original unextended primer; recovery unresolved')
        if stage=='library':
            sc.strand('PCR_daughters',[role('new amplified DNA with source adapter-derived ends; exact primer-bounded assembly unassigned')],label='DNA daughters; native adapter C3 is not copied into DNA')
            sc.strand('PCR_residuals',[role('possible original PCR primers and off-target DNA; class-specific recovery unresolved')])
    return sc


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    def inp(i,l,k=M.CONDITION,state=None):return g.input(i,l,k,state=state)
    def clean(i,l,x):
        reagent=inp(i+'_selection',l,M.SELECTION)
        y=g.product(i+'_out',l+' retained mixture',x.kind,state=x.state)
        w=g.product(i+'_waste',l+' nonretained fraction',M.WASTE)
        g.react(i,l,K.CLEANUP,((x,I.TEMPLATE),(reagent,I.SELECTION)),(y,w))
        g.discard(w,'Physical separation; exact hydrolysate/primer/adapter/off-product recovery is unresolved.')
        return y
    pool=g.input('rna_pool','folded pool of seven designed bar-coded RNAs',M.BIOLOGICAL,state=native('rna'),primary=True)
    plus=g.product('plus_aliquot','physical RNA aliquot for 1M7 treatment',M.MOLECULAR,state=native('rna'))
    minus=g.product('minus_aliquot','physical RNA aliquot for independent solvent control',M.MOLECULAR,state=native('rna'))
    g.react('split_channels','Split folded RNA into treated and control aliquots',K.SPLIT,((pool,I.TEMPLATE),),(plus,minus))
    g.mark_terminal(Terminal(minus.id,Disposition.UNRESOLVED,'Solvent channel uses YYYR and independent RT/ligation; downstream processing and eventual pooling timing outside selected treated-channel graph. This is not physical discard.'))
    condition=inp('one_m7','source 1M7 modification condition')
    modified=g.product('modified_rna','1M7-modified RNA mixture',M.MOLECULAR,state=native('modified'))
    g.react('modify','Acylate flexible original RNA 2′ hydroxyls',K.CHEMICAL_CONVERSION,((plus,I.TEMPLATE),(condition,I.CONDITION)),(modified,))
    recovery=inp('modified_recovery','source recovered modified RNA; detailed recovery chemistry not reconstructed')
    recovered=g.product('recovered_rna','recovered modified RNA',M.MOLECULAR,state=modified.state)
    g.react('recover_modified','Recover modified RNA before RT',K.MAGIC,((modified,I.TEMPLATE),(recovery,I.CONDITION)),(recovered,),reason='SI says recovered modified RNA; cited modification protocol details remain outside this source cap. No purification reagent or purity is invented.')
    prim=Scene();prim.strand('primer',rt_parts(),label='native treated RT primer; exact target-specific 3′ sequence unassigned')
    rp=inp('rt_primer','A_adapter_t–RRRY–target-specific DNA RT primer',M.OLIGO,prim)
    rt=inp('rt_enzyme','RT catalyst from cited Mortimer procedure; identity unassigned in inspected SHAPE SI',M.CATALYST)
    cond=inp('rt_condition','source heat step and primer-extension condition; cited procedure remains deferred')
    cdna=g.product('stop_cdna','RNA plus native primer-extended DNA mixture',M.MOLECULAR,state=native('rt'))
    g.react('reverse_transcribe','Extend actual handle-tailed primer to RNA stops',K.REVERSE_TRANSCRIPTION,((recovered,I.TEMPLATE),(rp,I.STRUCTURAL),(rt,I.CATALYST),(cond,I.CONDITION)),(cdna,))
    naoh=inp('naoh','source NaOH/heat RNA hydrolysis')
    hydro=g.product('hydrolysed','native ssDNA plus original RNA hydrolysis products',M.MOLECULAR,state=native('hydrolysed'))
    g.react('hydrolyse','Chemically hydrolyse RNA without asserting physical removal',K.CHEMICAL_CONVERSION,((cdna,I.TEMPLATE),(naoh,I.CONDITION)),(hydro,))
    x=clean('cdna_precipitation','Ethanol precipitation after RNA hydrolysis',hydro)
    ad=inp('adapter','source A_adapter_b with 3′ C3 spacer; 5′P unverified',M.ADAPTER,adapter_state())
    lig=inp('circligase','CircLigase ssDNA ligase; source ligation condition',M.CATALYST)
    adapted=g.product('adapted','native adapter-ligated ssDNA mixture',M.MOLECULAR,state=native('ligated'))
    g.react('ligate','Intermolecular ssDNA adapter ligation, not cDNA circularization',K.LIGATION,((x,I.TEMPLATE),(ad,I.STRUCTURAL),(lig,I.CATALYST)),(adapted,),geometry=AdapterLigationGeometry.unavailable('Source prints adapter DNA and terminal C3, but adapter 5′ phosphorylation remains unverified. Native ssDNA input retained; no cohesive/blunt adapter default.'))
    heat=inp('ligation_stop_condition','source heat stop')
    stopped=g.product('ligation_stopped','heat-stopped ligation mixture',M.MOLECULAR,state=adapted.state)
    g.react('stop_ligation','Stop ligation by heat',K.THERMAL_TREATMENT,((adapted,I.TEMPLATE),(heat,I.CONDITION)),(stopped,),changes_molecule=False)
    x=clean('ligation_precipitation','Ethanol precipitation after ligation',stopped)
    x=clean('adapter_cleanup','AMPure XP removal of excess adapter',x)
    f=inp('pcr_forward','adapter-binding PCR forward primer; bases unassigned',M.OLIGO)
    r=inp('pcr_reverse','adapter-binding PCR reverse primer; bases unassigned',M.OLIGO)
    pol=inp('pcr_catalyst','PCR catalyst; identity deferred to cited procedure',M.CATALYST)
    library=g.product('library','amplified DNA and native residual mixture',M.LIBRARY,state=native('library'))
    g.react('library_pcr','PCR of actual native ligated ssDNA mixture',K.PCR,((x,I.TEMPLATE),(f,I.STRUCTURAL),(r,I.STRUCTURAL),(pol,I.CATALYST)),(library,),geometry=PCRGeometry.unavailable(x.state,reason='Actual native ssDNA, terminal C3 and residual classes enter PCR; primer identities/sites and first-cycle assembly remain unassigned. Final read-site reference is not the incoming template.'))
    x=clean('library_cleanup','AMPure XP post-PCR cleanup; no size selection',library)
    g.sequence(x,SequencingPlan(library=FINAL_LIBRARY,primers=SEQ_PRIMERS,read_lengths=READ_LENGTHS,limitation='Component/read-site reference from published adapter arms; exact PCR-primer-bounded assembly is unresolved. '+SEQUENCING_ENDING+' Bioanalyzer QC sampling and eventual treated/control pooling allocations remain unassigned.'))
    g.declare_reviewed_scope(REVIEWED,SOURCE)
    return g.require_publishable()

REVIEWED=('rna_pool','split_channels','one_m7','modify','modified_recovery','recover_modified','rt_primer','rt_enzyme','rt_condition','reverse_transcribe','naoh','hydrolyse','cdna_precipitation_selection','cdna_precipitation','adapter','circligase','ligate','ligation_stop_condition','stop_ligation','ligation_precipitation_selection','ligation_precipitation','adapter_cleanup_selection','adapter_cleanup','pcr_forward','pcr_reverse','pcr_catalyst','library_pcr','library_cleanup_selection','library_cleanup','library_cleanup_out')

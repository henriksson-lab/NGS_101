"""Wang2019/Cheng2018 selected single-transposon3′branch stLFR native graph."""
from types import SimpleNamespace
from chemdraw import Construct, Scene, Segment, complement_segments, revcomp, panel
from protocol_model import role
from pcr import PCRGeometry
from protocol_graph import Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind
import seqprimers as sp
import stlfr_chem as S
if not S.BEAD_COMMON_T.endswith(S.PCR1):raise ValueError('workbook PCR1 doesnot land at common-primer boundary')
if revcomp(S.BRANCH_T)!=S.PCR2:raise ValueError('workbook PCR2/BranchT orientation mismatch')
if revcomp(S.BRANCH_T[:12])!=S.BRANCH_B+'G':raise ValueError('physical BranchB3′ddG doesnot pair at blunt5′BranchT end')
TITLE='stLFR — single-tube long-fragment reads'
NOTES='01_stlfr.html'
SOURCE='Wang2019 GenomeResearch doi:10.1101/gr.245126.118 Methods; Cheng2018 original detailedprotocol doi:10.1038/protex.2018.116 steps61–81 and original orderedworkbook.'
SUMMARY='Selected singleTransposon1T/3′branch route forstLFR1–3: capture on preassembled3setbarcoded beads, RNA-splinted barcode ligation, source digestions/washes, recessed3′branch ligation, PCR and reportedBGI circle/DNB readout.'
CAVEAT='ThreeN10barcode sets and source modifications remain local; multiple long molecules can share a bead barcode. Native genomic ends/gaps and physical BGIcircle assembly remain unassigned. Original branchPCR and common-adapter variants are disclosed below. No two-transposon/I llumina alternative is substituted. Hydrolysis/collision/nonamplification are not physical discard classes.'
READOUT_REFERENCE=Construct([role('BGI native readout flank; physical circle geometry unassigned'),*S.barcode_read().segments,role('genomic subfragment'),role('native circle/DNB assembly unassigned')],name='stLFR source readout reference; physical circle not constructed')
FINAL_CAPTION='Source42position barcode architecture and genomic readout; this is not a drawn physical DNAcircle.'
FINAL_LIBRARY=None
SEQ_PRIMERS=()
SEQUENCING_INTRO='Original2018 detailed protocol describes PE100 genomic reads and42position barcode/common-adapter string appended toRead1. Currentauthor code42/54positionRead2 conventions are distinct historical/source scopes.'
SOURCE_LIMITATION='Workbook lists sequencing oligos, but source-matched physical circle/native primer landing is not constructed; study circle and read/index geometry remain bounded.'
SEQUENCING_ENDING=SOURCE_LIMITATION
SEQUENCING_UNAVAILABLE=''  # Graph draws native circle/DNB scope; do not invent a duplex endpoint.
_readout_view=Scene();_readout_view.strand('reported string',READOUT_REFERENCE.segments,label='source barcode/readout reference; not native DNB topology')
TOPOLOGY_PANELS=(('Reported barcode/readout architecture',panel(_readout_view.rows(),caption='Source-reported single string reference only. It is not the physical circle or DNB repeat topology.')),)
SHOW_READ_CONTENT=False
READ_SPECS=(sp.read_spec('Read 1',100,'genomic read plus separately reported42positionbarcode/common-adapter string','include'),sp.read_spec('Read 2',100,'opposite genomic read','include'))
oligos=S.oligos

def native_readout_state(kind):
 sc=Scene()
 parts=[role('reported '+kind+' physical polymer/closure context; strand andrepeat arrangement unassigned')]
 for source_segment in S.barcode_read().segments:
  if source_segment.feature is not None:
   parts.append(role(source_segment.name+' logical identity; native placement/multiplicity unassigned',molecular_feature=source_segment.feature))
 parts.append(role('native genomic fragment and opposite context; physical arrangement unassigned'))
 sc.strand('native scope',parts,label=kind+' scope; logical barcode lineages, not base-resolved topology')
 return sc


def transposome_scene():
 sc=Scene()
 sc.strand('top',[Segment('capture donor prefix',S.TRANSPOSON_1_T[:-19]),Segment('donor ME',S.TRANSPOSON_1_T[-19:])],label='Transposon1T',mod5='p')
 sc.anneal('bottom',[Segment('uracil-bearing DNA strand',S.TRANSPOSON_B,note='Source U residues are DNA deoxyuridines; UDG substrate, not RNA')],'top',pair=('uracil-bearing DNA strand','donor ME'),label='TransposonB DNA/dU',mod5='p')
 return sc

def native_state(stage):
 sc=Scene()
 if stage in ('tagged','captured'):
  sc.strand('genomic',[Segment('transferred Transposon1T',S.TRANSPOSON_1_T),role('long genomic DNA; native insertion/break geometry unassigned')],label='representative tagged genomic lineage; multiple long molecules may share barcode',mod5='p')
  sc.strand('bottom',[Segment('TransposonB dU-bearing DNA',S.TRANSPOSON_B,note='DNA deoxyuridines; post-digestion survival unassigned')],label='independent nontransferred donor-strand population')
  if stage=='captured':
   sc.strand('bead',S.assembled_bead_oligo().segments,label='independent captured bead adapter before ligation',mod5='dual-biotin')
   sc.strand('bridge',[Segment('RNA bridge',S.BRIDGE_RNA,note='This bridge is RNA')],label='independent RNA splint; paired span unassigned')
  return sc
 parts=[*S.assembled_bead_oligo().segments,Segment('ligated Transposon1T',S.TRANSPOSON_1_T),role('genomic subfragment; native extent/ends unassigned')]
 if stage=='branch':parts += [Segment('BranchT paired region',S.BRANCH_T[:12]),Segment('BranchT3′tail',S.BRANCH_T[12:])]
 sc.strand('barcoded',parts,label='source-ordered barcoded lineage reference; native opposite/gap extent unassigned',mod5='dual-biotin')
 if stage=='branch':
  sc.anneal('BranchB',[Segment('BranchB DNA',S.BRANCH_B),Segment('terminal3′ddG','G',note='Source /3′dd/G is a physical terminal dideoxyG; not copied as an added internalG')],'barcoded',pair=('terminal3′ddG','BranchT paired region'),label='source blocked opposite adapter strand',mod3='ddG')
 else:sc.strand('other',[role('native opposite donor/DNA/RNA-splint populations; post-digestion/wash survival unassigned')],label='independent molecular classes')
 return sc

def branch_scene():
 sc=Scene();sc.strand('BranchT',[Segment('BranchT paired region',S.BRANCH_T[:12]),Segment('BranchT3′tail',S.BRANCH_T[12:])],label='source BranchT',mod5='p')
 sc.anneal('BranchB',[Segment('BranchB DNA',S.BRANCH_B),Segment('terminal3′ddG','G',note='Physical dideoxyG terminal residue')],'BranchT',pair=('terminal3′ddG','BranchT paired region'),label='source BranchB',mod3='ddG')
 return sc

def wash(g,key,template,label):
 reagent=g.input(key+'_condition',label,MaterialKind.SELECTION)
 kept=g.product(key+'_retained','source washed bead-associated DNA',MaterialKind.MOLECULAR,state=template.state)
 waste=g.product(key+'_waste','physically removed supernatant/washes',MaterialKind.WASTE)
 g.react(key,label,StepKind.CLEANUP,((template,InputRole.TEMPLATE),(reagent,InputRole.SELECTION)),(kept,waste),changes_molecule=False)
 g.discard(waste,'Source magnetic/bead wash fractions; no barcode-collision selection or complete molecular-class depletion inferred.')
 return kept

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    dna=g.input('dna','long genomic DNA',MaterialKind.BIOLOGICAL,state=Construct([role('long genomic DNA; native extent unassigned')],name='long DNA'),primary=True)
    transposome=g.input('transposome','singleTransposon1T +5′P DNA/dUTransposonB',MaterialKind.ADAPTER,state=transposome_scene())
    tn5=g.input('tn5','Tn5 insertion, source55°C10min',MaterialKind.CATALYST)
    tagged=g.product('tagged','capture-tagged long genomic DNA',MaterialKind.MOLECULAR,state=native_state('tagged'))
    g.react('tagment','Insert selected singlecapture transposon',StepKind.TAGMENTATION,((dna,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,))
    bead_sc=Scene();bead_sc.strand('bead',S.assembled_bead_oligo().segments,label='orderedworkbook bead variant',mod5='dual-biotin')
    bead=g.input('barcoded_bead','preassembled3setN10barcoded beads; statistical occupancy',MaterialKind.SELECTION,state=bead_sc)
    bridge=g.input('bridge','source RNA capture bridge',MaterialKind.OLIGO,state=Construct([Segment('RNA bridge',S.BRIDGE_RNA,note='RNA')],name='RNA bridge'))
    captured=g.product('captured','bead-associated tagged long DNA; multiple occupancy allowed',MaterialKind.MOLECULAR,state=native_state('captured'))
    g.react('capture','Hybridize tagged DNA to bead barcodes through RNAbridge',StepKind.HYBRIDIZATION,((tagged,InputRole.TEMPLATE),(bead,InputRole.SELECTION),(bridge,InputRole.STRUCTURAL)),(captured,))
    ligase=g.input('ligase','sourceT4DNA ligase',MaterialKind.CATALYST)
    barcoded=g.product('barcoded','bead-barcode ligated DNA native reference',MaterialKind.MOLECULAR,state=native_state('barcoded'))
    g.react('ligate_barcode','Ligate bead barcode at source RNA-splinted nick',StepKind.LIGATION,((captured,InputRole.TEMPLATE),(captured,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(barcoded,),geometry=S.barcode_transfer_geometry())
    current=wash(g,'wash_capture',barcoded,'Afterbarcode ligation:magnet/supernatant removal/low-saltwash')
    exo1=g.input('exo1','90UExonucleaseI',MaterialKind.CATALYST)
    exo3=g.input('exo3','100UExonucleaseIII;37°C10min',MaterialKind.CATALYST)
    digested=g.product('adapter_digested','adapter-digested native bead-associated DNA',MaterialKind.MOLECULAR,state=native_state('processed'))
    debris=g.product('adapter_hydrolysate','enzyme-hydrolysed excess adapter products',MaterialKind.MOLECULAR)
    g.react('digest_adapters','Digest adapter oligos while Tn5 remains bound',StepKind.DIGESTION,((current,InputRole.TEMPLATE),(exo1,InputRole.CATALYST),(exo3,InputRole.CATALYST)),(digested,debris))
    g.mark_terminal(Terminal(debris.id,Disposition.UNRESOLVED,'Hydrolysis is chemical, not physical removal; survival/recovery through source washes unassigned.'))
    sds=g.input('sds','11µl1%SDS;RT10min',MaterialKind.CONDITION)
    deprot=g.product('deproteinized','native bead DNA after Tn5 dissociation',MaterialKind.MOLECULAR,state=digested.state)
    g.react('remove_tn5','Stopdigestion/remove Tn5 withSDS',StepKind.DEPROTEINIZATION,((digested,InputRole.TEMPLATE),(sds,InputRole.CONDITION)),(deprot,))
    current=wash(g,'wash_sds',deprot,'Magnet;low-saltwash thenNEB2wash')
    udg=g.input('udg','10UUDG',MaterialKind.CATALYST)
    ape1=g.input('ape1','30UAPE1;37°C30min',MaterialKind.CATALYST)
    processed=g.product('capture_digested','capture-oligo-digested native lineage',MaterialKind.MOLECULAR,state=native_state('processed'))
    hydro=g.product('capture_hydrolysate','dU/base-excision cleavage products; survival unassigned',MaterialKind.MOLECULAR)
    g.react('digest_capture','Digest excesscapture oligos usingUDG/APE1',StepKind.DIGESTION,((current,InputRole.TEMPLATE),(udg,InputRole.CATALYST),(ape1,InputRole.CATALYST)),(processed,hydro))
    g.mark_terminal(Terminal(hydro.id,Disposition.UNRESOLVED,'Base excision/cleavage is not physical discard; selected route doesnot add dual-transposonExoI at this stage.'))
    current=wash(g,'wash_udg',processed,'Magnet;high-saltwash thenlow-saltwash')
    ssb=g.input('ssb','ETSSB;37°C30min',MaterialKind.MOLECULAR)
    bound=g.product('ssb_bound','SSB-associated native DNA; protein not a DNA strand segment',MaterialKind.MOLECULAR,state=current.state)
    g.react('bind_ssb','BindETSSB before3′branch ligation',StepKind.BINDING,((current,InputRole.TEMPLATE),(ssb,InputRole.STRUCTURAL)),(bound,))
    branch=g.input('branch_adapter','source5′PBranchT / BranchB withphysical3′ddG',MaterialKind.ADAPTER,state=branch_scene())
    joined=g.product('branch_ligated','native3′branch ligation precursor',MaterialKind.MOLECULAR,state=native_state('branch'))
    g.react('branch_ligation','LigateBranchT5′phosphate to recessedgenomic3′OH;RT2h',StepKind.LIGATION,((bound,InputRole.TEMPLATE),(branch,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(joined,),geometry=S.branch_adapter_geometry())
    current=wash(g,'wash_branch',joined,'Magnet;high-saltwash thenPCRbufferwash')
    p1=g.input('pcr1','sourceworkbook PCR1',MaterialKind.OLIGO,state=Construct([Segment('PCR1',S.PCR1)],name='PCR1'))
    p2=g.input('pcr2','sourcebranch PCR2',MaterialKind.OLIGO,state=Construct([Segment('PCR2',S.PCR2)],name='PCR2'))
    pol=g.input('polymerase','sourcePfuCx; programme variants disclosed separately',MaterialKind.CATALYST)
    amplicon=Construct([Segment('PCR1 copied boundary',S.PCR1),*S.assembled_bead_oligo().segments[1:],Segment('sourcecapture donor',S.TRANSPOSON_1_T),role('native genomic subfragment'),Segment('ligatedBranchT',S.BRANCH_T)],name='stLFR amplified lineage reference')
    amplified=g.product('amplified','source co-barcoded PCR products',MaterialKind.MOLECULAR,state=amplicon)
    originals=g.product('originals','bead-tethered originals/RNAbridge/blockedBranchB and other unassigned polymers',MaterialKind.MOLECULAR,state=current.state)
    g.react('pcr','Amplify actual native3′branch-ligated lineage withPCR1/PCR2',StepKind.PCR,((current,InputRole.TEMPLATE),(p1,InputRole.STRUCTURAL),(p2,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(amplified,originals),geometry=PCRGeometry.unavailable(current.state,reason='Actual source branch-ligated native state; complete genomicbreak/gap andfirst-cycle duplex proof unassigned. Physical BranchB3′ddG is not an insertedinternal base.'))
    g.mark_terminal(Terminal(originals.id,Disposition.UNRESOLVED,'PCRcopies do not establish disposal ofbiotin-tethered originals orblocked/ribose strands; physical recovery unassigned.'))
    amplified=wash(g,'pcr_cleanup',amplified,'Originalstep81→59–60:AgencourtXPpurification; source no barcode-ambiguity rejection')
    circle=g.product('native_circle','reportedBGIcircle; native strand/closure geometry unassigned',MaterialKind.MOLECULAR,state=native_readout_state('BGIcircle'))
    g.react('circle','Reported circularization—nativeassembly boundary',StepKind.MAGIC,((amplified,InputRole.TEMPLATE),),(circle,),reason='Source reports circularization/DNB and citesDrmanac2010. Native circle precursor strand selection/closure sequence is outside current3-original source scope; no completed readout is used as a circle template.')
    priming=g.input('dnb_primers','nativeDNB priming reagents; selectedbranch-compatible species unassigned',MaterialKind.OLIGO)
    dnb_enzyme=g.input('dnb_enzyme','reportedDNB amplification enzyme; identity unassigned',MaterialKind.CATALYST)
    library=g.product('library','reportedDNB sequencing material; physicalassembly unassigned',MaterialKind.LIBRARY,state=native_readout_state('DNA nanoball/repeated-copy population'))
    g.react('dnb','Amplify reportedcircle intoDNB material',StepKind.EXTENSION,((circle,InputRole.TEMPLATE),(priming,InputRole.STRUCTURAL),(dnb_enzyme,InputRole.CATALYST)),(library,))
    g.sequence(library,SequencingPlan(acquisition='unspecified',limitation=SOURCE_LIMITATION))
    g.declare_reviewed_scope(REVIEWED,SOURCE)
    g.annotate(SourceDiscrepancy("barcoded_bead", "Published common-adapter variants",
        (SourceClaim("Original Protocol Exchange Supplementary Figure 1",
                     "/52-Bio/AAAAAAAACCCGTAGCCATGTCGTTCTG; 8A linker and a different common-primer region."),
         SourceClaim("Original Supplement_0.xlsx, Common adapter sequences, BeadCommonT/PCR1",
                     "/52-Bio/AAAAAAAAAATGTGAGCCAAGGAGTTG; 10A linker, with matching PCR1 TGTGAGCCAAGGAGTTG.")),
        "The model uses the ordered workbook variant and its matching PCR1.",
        "The workbook supplies explicit oligo order sequences and a self-consistent primer pair; the protocol also allows changing the common-primer sequence for different sequencing technologies.",
        "The figure and workbook do not explain which physical preparations used each variant. This is an unresolved source variant, not a demonstrated author error."))
    g.annotate(SourceDiscrepancy("pcr", "Published branch-route PCR programme",
        (SourceClaim("Wang2019 Methods, 'stLFR with 3′ branched adapter ligation' (JATS line 212) referring to the two-transposon PCR (line 208)",
                     "PCR mix and conditions 'the same as the two-transposon process': PfuCx, 400 pmol each primer, 10 min at 72°C then 12 cycles of 95°C 10 s / 58°C 30 s / 72°C 2 min."),
         SourceClaim("Cheng2018 original manuscript step 80 (p. 9), single-transposon 3′ branch route",
                     "95°C 3 min, then 95°C 10 s / 58°C 30 s / 72°C 2 min repeated 10–12 times; no initial 72°C 10 min step (that step appears in the two-transposon step 58).")),
        "The model draws one PCR with source PCR1/PCR2 and does not select either programme or a cycle number.",
        "Both are authoritative originals for the selected route; the primers, enzyme scale and annealing/extension cycle agree, while the initial step and cycle count differ.",
        "The sources do not state which programme was used for libraries stLFR-1–3 or whether a 72°C extension before denaturation is needed after 3′ branch ligation."))
    return g.require_publishable()

REVIEWED=('dna','transposome','tn5','tagment','barcoded_bead','bridge','capture','ligase','ligate_barcode','wash_capture_condition','wash_capture','exo1','exo3','digest_adapters','sds','remove_tn5','wash_sds_condition','wash_sds','udg','ape1','digest_capture','wash_udg_condition','wash_udg','ssb','bind_ssb','branch_adapter','branch_ligation','wash_branch_condition','wash_branch','pcr1','pcr2','polymerase','pcr','pcr_cleanup_condition','pcr_cleanup','circle','dnb_primers','dnb_enzyme','dnb','library')
SPEC=SimpleNamespace(**{k:v for k,v in globals().copy().items() if k.isupper()},protocol_graph=protocol_graph,oligos=oligos)

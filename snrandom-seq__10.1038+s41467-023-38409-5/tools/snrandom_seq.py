"""Xu et al. 2023; source-owned snRandom reaction and identifier declarations."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]; sys.path[:0]=[str(ROOT/'lib')]
from chemdraw import Construct, Scene, Segment, feature, oligo, revcomp
from protocol_model import role
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
import seqprimers as sp
TITLE='snRandom-seq — FFPE single-nucleus total RNA'
NOTES='01_snrandom-seq.html'
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-023-38409-5">Xu et al., Nature Communications (2023)</a>; published Methods and SI Table1.'
SUMMARY='Nuclear DNA blocking, random/oligo(dT) RT and variable dA tailing precede droplet bead extension, PFO recovery, cDNA PCR and ND607-01 library preparation.'
CAVEAT='Native Vazyme adapter strands, PCR ends and index allocation remain unassigned in this bounded pass. No TruSeq shell or fixed index lengths are supplied. The molecular release site at the bead deoxyuridine is not assigned from the inspected Methods; no cleavage reaction is invented.'
BLOCK_PRIMER="GAGAATGTGAGTGAAGATGTATGGTGANNNNNNN"
RT_HANDLE="GGAGTTGGAGTGAGTGGATGAGTGATGGAAGGAAT"
BEAD_LINKER_12="TGGT"
BEAD_LINKER_23="GAGA"
BEAD_ANCHOR_WRITTEN="/5Acryd/ATTATATATATUGTGAGTGATGGTTGAGGATGTGTGGAGATA"


PCR1='GGAGTTGGAGTGAGTGGATGAGTGATG'
PCR2='GTGAGTGATGGTTGAGGATGTGTGGAGATA'
cbc1=Segment('cell barcode 1','B'*10,'cbc',placeholder=True,feature=feature('cell_barcode_1','cell_barcode','combinatorial',group='cell_barcode',part='round 1'))
cbc2=Segment('cell barcode 2','C'*10,'cbc',placeholder=True,feature=feature('cell_barcode_2','cell_barcode','combinatorial',group='cell_barcode',part='round 2'))
cbc3=Segment('cell barcode 3','D'*10,'cbc',placeholder=True,feature=feature('cell_barcode_3','cell_barcode','combinatorial',group='cell_barcode',part='round 3'))
umi=Segment('UMI','U'*8,'umi',placeholder=True,feature=feature('umi','umi','random'))
def barcode_parts():
 return [cbc1,Segment('TGGT linker',BEAD_LINKER_12),cbc2,Segment('GAGA linker',BEAD_LINKER_23),cbc3,umi,Segment('poly(dT)','T'*20)]
FINAL_LIBRARY=Construct([role('Vazyme native end A; strand geometry unassigned'),*barcode_parts(),role('RNA-derived insert; fragment boundary unassigned'),role('Vazyme native end B; strand geometry unassigned')],name='snRandom-seq source-limited native endpoint')
SOURCE_OLIGO_EXCEPTIONS={'barcoded_bead':'Published SI Table1 prints /5Acryd/ATTATATATAT U GTG AGT GAT GGT TGA GGA TGT GTG GAGATA followed by three [10 bases barcode] slots, TGGT/GAGA linkers, N8 and T20. Verified both lines directly; the guard cannot match the assembled class-expanded string against those named placeholders. Molecular dU and 5′ acrydite are retained.'}
SEQ_PRIMERS=()
FINAL_CAPTION='Published barcode and UMI structure, with native vendor-end roles; no canonical Illumina outer sequence assignment.'
SEQUENCING_INTRO='NovaSeq6000 S4 PE150. The paper extracts 30 barcode bases and eight UMI bases from Read1; Read2 supplies expression sequence. Positions relative to unassigned native primer landing are source-declared, not computed.'
SEQUENCING_ENDING=SEQUENCING_INTRO
SHOW_READ_CONTENT=False
READ_SPECS=(sp.read_spec('Read 1',150,'Three10-base barcode parts, intervening linkers and8-base UMI; exact cycle offsets unassigned','extract'),sp.read_spec('Read 2',150,'Transcript sequence; primers and dA-tailing sequence computationally trimmed','include'),sp.read_spec('Index 1 (i7)',None,'Native index acquisition unassigned','undetermined'),sp.read_spec('Index 2 (i5)',None,'Native index acquisition unassigned','undetermined'))
def nuclear_scene(stage):
 sc=Scene(); sc.strand('genome',[role('genomic DNA; '+('block-primer extension present' if stage!='initial' else 'accessible ssDNA'))],label='genomic DNA; retained nuclear population'); sc.strand('RNA',[role('nuclear RNA template')],label='independent RNA population')
 if stage in ('rt','tail'):
  parts=[Segment('representative variant1 RT handle',RT_HANDLE),role('RNA-derived first-strand cDNA')]
  if stage=='tail': parts.append(role('variable-length 3′ dA tail'))
  sc.strand('cDNA',parts,label='cDNA; other four RT handle variants are distinct molecules')
 return sc

def barcoded_scene():
 sc=Scene(); sc.strand('barcoded',[role('upstream acrydite–ATTATATATAT–dU boundary; release and retention unassigned'),Segment('bead-side PCR2 handle',PCR2),*barcode_parts(),role('RNA-derived copied insert'),Segment('copy of representative RT handle',revcomp(RT_HANDLE))],label='barcoded DNA product; release-site boundary unassigned'); sc.strand('genome',[role('genomic DNA remaining in droplet mixture')],label='independent genomic DNA population'); sc.strand('RNA',[role('RNA remaining in droplet mixture')],label='independent RNA population'); return sc

def oligos():
 yield oligo('DNA blocking primer',[Segment('block',BLOCK_PRIMER)])
 yield oligo('Representative random RT primer1',[Segment('RT handle',RT_HANDLE),Segment('N7','N'*7,placeholder=True)])
 yield oligo('Bead anchor as written',[Segment('acrydite/deoxyuridine anchor',BEAD_ANCHOR_WRITTEN,placeholder=True)])
 yield oligo('PCR primer1',[Segment('primer1',PCR1)])
 yield oligo('PCR primer2',[Segment('primer2',PCR2)])

def native_duplex(parts,*,label):
 sc=Scene(); sc.strand('DNA',parts,label=label+'; declared strand roles'); sc.strand('opposite',[role('opposite DNA strand; exact paired geometry unassigned')],label='separate complementary strand role'); return sc

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 nucleus=g.input('nucleus','qualified fixed permeabilized FFPE nuclei',MaterialKind.BIOLOGICAL,state=nuclear_scene('initial'),primary=True)
 block=g.input('block_primer','N7 genomic DNA block primer',MaterialKind.OLIGO,state=Construct([Segment('block primer',BLOCK_PRIMER)]))
 pol=g.input('block_polymerase','M20 DNA Polymerase kit',MaterialKind.CATALYST)
 blocked=g.product('blocked','DNA-blocked nuclei',MaterialKind.BIOLOGICAL,state=nuclear_scene('blocked'))
 g.react('block_dna','Block accessible ssDNA: 37°C 30 min',StepKind.EXTENSION,((nucleus,InputRole.TEMPLATE),(block,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(blocked,))
 wash=g.input('wash_buffer','PBS+0.05% Tween20',MaterialKind.CONDITION)
 washed=g.product('washed_block','washed DNA-blocked nuclei',MaterialKind.BIOLOGICAL,state=blocked.state)
 waste=g.product('wash_block_waste','PBST washes with residual blocking primers/dimers',MaterialKind.WASTE)
 g.react('wash_block','Wash nuclei three times; remove residual blocking primers/dimers',StepKind.CLEANUP,((blocked,InputRole.TEMPLATE),(wash,InputRole.CONDITION)),(washed,waste))
 g.mark_terminal(Terminal(waste.id,Disposition.DISCARD,reason='Non-retained physical fraction at this source-described purification; no species-specific depletion asserted.'))
 rtprim=g.input('rt_primers','five randomN7 plus five oligo(dT) primer variants, SI Table1',MaterialKind.OLIGO)
 rt=g.input('reverse_transcriptase','M20 reverse transcriptase kit',MaterialKind.CATALYST)
 cdna=g.product('cdna','in-nucleus RNA/cDNA mixture',MaterialKind.MOLECULAR,state=nuclear_scene('rt'))
 g.react('reverse_transcribe','12 annealing ramps8–42°C then42°C30min',StepKind.REVERSE_TRANSCRIPTION,((washed,InputRole.TEMPLATE),(rtprim,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,))
 rtwashed=g.product('washed_rt','washed in-nucleus cDNA',MaterialKind.MOLECULAR,state=cdna.state)
 waste=g.product('wash_rt_waste','PBST washes with residual RT primers/dimers',MaterialKind.WASTE)
 g.react('wash_rt','Three PBST washes remove residual RT primers/dimers',StepKind.CLEANUP,((cdna,InputRole.TEMPLATE),(wash,InputRole.CONDITION)),(rtwashed,waste))
 g.mark_terminal(Terminal(waste.id,Disposition.DISCARD,reason='Non-retained physical fraction at this source-described purification; no species-specific depletion asserted.'))
 tdt=g.input('tdt','TdT+dATP+CoCl2',MaterialKind.CATALYST)
 tailed=g.product('tailed','variable-length dA-tailed cDNA in nuclei',MaterialKind.MOLECULAR,state=nuclear_scene('tail'))
 g.react('da_tail','TdT dA tailing:37°C30min',StepKind.TAILING,((rtwashed,InputRole.TEMPLATE),(tdt,InputRole.CATALYST)),(tailed,))
 tailwashed=g.product('washed_tail','washed dA-tailed nuclei',MaterialKind.MOLECULAR,state=tailed.state)
 waste=g.product('wash_tail_waste','tailing wash fractions; composition unassigned',MaterialKind.WASTE)
 g.react('wash_tail','Three PBS/Tween washes',StepKind.CLEANUP,((tailed,InputRole.TEMPLATE),(wash,InputRole.CONDITION)),(tailwashed,waste))
 g.mark_terminal(Terminal(waste.id,Disposition.DISCARD,reason='Non-retained physical fraction at this source-described purification; no species-specific depletion asserted.'))
 bead_sc=Scene(); bead_sc.strand('bead',[Segment('anchor before dU','ATTATATATAT'),Segment('deoxyuridine','U',note='source deoxyuridine, not ribouridine'),Segment('bead PCR handle',PCR2),*barcode_parts()],label='bead-linked primer; split-pool assembly documented separately',mod5='acrydite')
 bead=g.input('barcoded_bead','bead primers assembled by three split-pool ligations',MaterialKind.OLIGO,state=bead_sc)
 extension=g.input('extension_mix','M20 2X DNA extension reaction mix',MaterialKind.CATALYST)
 barcoded=g.product('barcoded','droplet barcoding mixture',MaterialKind.MOLECULAR,state=barcoded_scene())
 g.react('bead_extend','Encapsulate nuclei/bead;37°C1h,50°C30min,60°C30min,75°C20min',StepKind.EXTENSION,((tailwashed,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(extension,InputRole.CATALYST)),(barcoded,))
 pfo=g.input('pfo','PFO buffer',MaterialKind.CONDITION)
 aqueous=g.product('aqueous','recovered aqueous phase',MaterialKind.MOLECULAR,state=barcoded.state)
 g.react('break_emulsion','Break droplets with PFO; recover aqueous phase',StepKind.EXTRACTION,((barcoded,InputRole.TEMPLATE),(pfo,InputRole.CONDITION)),(aqueous,))
 beads=g.input('ampure','AMPure XP',MaterialKind.SELECTION)
 purified=g.product('purified','bead-purified barcoded DNA precursor',MaterialKind.MOLECULAR,state=barcoded.state)
 waste=g.product('purify_preamp_waste','pre-PCR AMPure non-retained fraction; composition unassigned',MaterialKind.WASTE)
 g.react('purify_preamp','AMPure XP before cDNA PCR',StepKind.CLEANUP,((aqueous,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(purified,waste))
 g.mark_terminal(Terminal(waste.id,Disposition.DISCARD,reason='Non-retained physical fraction at this source-described purification; no species-specific depletion asserted.'))
 primers=g.input('preamp_primers','SI Table1 PCR Primer1',MaterialKind.OLIGO,state=Construct([Segment('PCR1',PCR1)]))
 reverse=g.input('preamp_reverse','SI Table1 PCR Primer2',MaterialKind.OLIGO,state=Construct([Segment('PCR2',PCR2)]))
 amp_sc=native_duplex([Segment('PCR2-side handle',PCR2),*barcode_parts(),role('copied transcript insert'),Segment('PCR1 opposite end',revcomp(PCR1))],label='barcoded cDNA amplicon')
 amp=g.product('amplified','barcoded cDNA PCR product',MaterialKind.MOLECULAR,state=amp_sc)
 other=g.product('nonamplified','nonamplified RNA/genomic DNA; partition unassigned',MaterialKind.MOLECULAR)
 g.react('preamp','Amplify purified product with Primer1/2; cycle count unassigned',StepKind.PCR,((purified,InputRole.TEMPLATE),(primers,InputRole.STRUCTURAL),(reverse,InputRole.STRUCTURAL),(extension,InputRole.CATALYST)),(amp,other),geometry=PCRGeometry.unavailable(purified.state,reason='Actual recovered barcoding mixture; complete first-copy strand and release geometry unassigned.'))
 g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,reason='Non-amplification or cleanup does not establish a genomic DNA/RNA discard fraction.'))
 clean=g.product('clean_amp','AMPure-purified cDNA amplicons',MaterialKind.MOLECULAR,state=amp_sc)
 waste=g.product('clean_amplicon_waste','post-PCR AMPure non-retained fraction; composition unassigned',MaterialKind.WASTE)
 g.react('clean_amplicon','AMPure purification after cDNA amplification',StepKind.CLEANUP,((amp,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(clean,waste))
 g.mark_terminal(Terminal(waste.id,Disposition.DISCARD,reason='Non-retained physical fraction at this source-described purification; no species-specific depletion asserted.'))
 frag_sc=native_duplex([*barcode_parts(),role('barcoded amplified cDNA fragment; exact cut boundary unassigned')],label='fragmented DNA')
 fragmented=g.product('fragmented','50ng fragmented DNA for named kit',MaterialKind.MOLECULAR,state=frag_sc)
 g.react('fragment','Fragmentation boundary before kit end preparation',StepKind.MAGIC,((clean,InputRole.TEMPLATE),),(fragmented,),reason='Published library-preparation paragraph names fragmented DNA but does not assign the fragmentation mechanism in the inspected scope; no generic nuclease is substituted.')
 endmix=g.input('endprep_reagents','ND607-01 end-repair/adenylation reagents',MaterialKind.CATALYST)
 repaired=g.product('endprepped','end-repaired and adenylated barcoded DNA fragments',MaterialKind.MOLECULAR,state=frag_sc)
 g.react('endprep','End repair+adenylation:30°C30min,65°C30min',StepKind.END_REPAIR,((fragmented,InputRole.TEMPLATE),(endmix,InputRole.CATALYST)),(repaired,))
 adapter=g.input('adapter','working Vazyme adapter; native strands unassigned',MaterialKind.ADAPTER)
 ligase=g.input('ligase','kit ligation enzymes',MaterialKind.CATALYST)
 adapted=g.product('adapted','adapter-ligated DNA',MaterialKind.MOLECULAR,state=native_duplex([role('native adapter end A'),*barcode_parts(),role('barcoded cDNA fragment'),role('native adapter end B')],label='ligation product; exact junction unassigned'))
 g.react('ligate_adapter','20°C15min adapter ligation',StepKind.LIGATION,((repaired,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(adapted,),geometry=AdapterLigationGeometry.unavailable('Native vendor strands and end geometry are unassigned in inspected source scope.'))
 selected=g.product('selected','AMPure size-selected adapted DNA',MaterialKind.MOLECULAR,state=adapted.state)
 waste=g.product('select_library_waste','AMPure size-selection rejected fraction; size boundaries unassigned',MaterialKind.WASTE)
 g.react('select_library','Purify and size-select ligated DNA with AMPure XP',StepKind.SELECTION,((adapted,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(selected,waste))
 g.mark_terminal(Terminal(waste.id,Disposition.DISCARD,reason='Non-retained physical fraction at this source-described purification; no species-specific depletion asserted.'))
 kitpcr=g.input('kit_pcr','native kit PCR reagents/primers; identities unassigned',MaterialKind.CATALYST)
 libraryamp=g.product('library_amp','amplified native library',MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 kitf=g.input('kit_forward','native kit PCR primer A',MaterialKind.OLIGO)
 kitr=g.input('kit_reverse','native kit PCR primer B',MaterialKind.OLIGO)
 g.react('library_pcr','Kit library amplification',StepKind.PCR,((selected,InputRole.TEMPLATE),(kitf,InputRole.STRUCTURAL),(kitr,InputRole.STRUCTURAL),(kitpcr,InputRole.CATALYST)),(libraryamp,),geometry=PCRGeometry.unavailable(selected.state,reason='Actual size-selected native adapted precursor; no generic final duplex substituted.'))
 library=g.product('library',FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 waste=g.product('final_cleanup_waste','final AMPure non-retained fraction; composition unassigned',MaterialKind.WASTE)
 g.react('final_cleanup','Final AMPure XP purification; Qubit/Qsep QC',StepKind.CLEANUP,((libraryamp,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(library,waste))
 g.mark_terminal(Terminal(waste.id,Disposition.DISCARD,reason='Non-retained physical fraction at this source-described purification; no species-specific depletion asserted.'))
 g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=(),read_lengths={'Read 1':150,'Read 2':150},limitation='PE150 and identifier extraction source-declared; exact native sequencing primer landing and index allocation unassigned.'))
 g.declare_reviewed_scope(('preamp_reverse','kit_forward','kit_reverse','nucleus','block_primer','block_polymerase','block_dna','blocked','wash_buffer','wash_block','washed_block','rt_primers','reverse_transcriptase','reverse_transcribe','cdna','wash_rt','washed_rt','tdt','da_tail','tailed','wash_tail','washed_tail','barcoded_bead','extension_mix','bead_extend','barcoded','pfo','break_emulsion','aqueous','ampure','purify_preamp','purified','preamp_primers','preamp','amplified','clean_amplicon','clean_amp','fragment','fragmented','endprep','endprep_reagents','endprepped','adapter','ligase','ligate_adapter','adapted','select_library','selected','kit_pcr','library_pcr','library_amp','final_cleanup','library'),SOURCE)
 return g.require_publishable()

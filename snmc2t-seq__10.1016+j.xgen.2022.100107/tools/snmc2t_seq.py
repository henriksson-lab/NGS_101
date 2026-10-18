"""Luo2022 snmC2T nuclear RNA/methylome branch; selected human-brain RT conditions."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Scene, Segment, feature, oligo
from protocol_model import role
from pcr import PCRGeometry
from protocol_graph import ProtocolGraph, MaterialKind, InputRole, StepKind, SequencingPlan, Terminal, Disposition
import seqprimers as sp
TITLE='snmC2T-seq — transcriptome and methylome in one bisulfite library'
NOTES='01_snmc2t-seq.html'
SOURCE='<a href="https://doi.org/10.1016/j.xgen.2022.100107">Luo et al., Cell Genomics2022</a>, published PMC9004682 Results215/Methods540–554; nuclear finish explicitly inherits Luo2018 snmC-seq2.'
SUMMARY='Selected nuclear RNA/methylome subworkflow: methyl-dCTP RT/cDNA amplification, brain UDG-only treatment, mixed-DNA bisulfite conversion and native snmC-seq2 library finishing. RNA/DNA read partitioning is computational.'
CAVEAT='The broader snmCAT workflow adds NOMe marking; that accessibility front end is outside this snmC2T subworkflow. Whole cultured cells use the earlier snmC finish and12cDNA-PCR cycles; brain nuclei use14cycles and UDG without EndoVIII. Native 2018-family finishing order is inherited only for nuclear libraries. Study-specific random-primer/index members and Adaptase end geometry remain unassigned; no finished canonical duplex is substituted for an actual PCR input.'
DT30VN = "AAGCAGUGGUAUCAACGCAGAGUACUTTTTTUTTTTTUTTTTTUTTTTTUTTTTTVN"
N6 = "AAGCAGUGGUAUCAACGCAGAGUACNNNNNN"
TSO = "AAGCAGUGGUAUCAACGCAGAGUGAAUrGrG+G"
ISPCR23 = "AAGCAGUGGUAUCAACGCAGAGU"



RP_HANDLE='TTCCCTACACGACGCTCTTCCGATCT'
CELL_FEATURE=feature('snmc2t_inline','cell_barcode','unknown',group='cell_identity',note='2018 nuclear-finish inline family; experiment-specific member unassigned.')
I5_FEATURE=feature('snmc2t_i5','sample_index','unknown',group='cell_identity',note='Native dual-index family; selected member/cycle allocation unassigned.')
I7_FEATURE=feature('snmc2t_i7','sample_index','unknown',group='cell_identity',note='Native dual-index family; selected member/cycle allocation unassigned.')
def random_parts():
 return [Segment('2018 RP-H handle',RP_HANDLE),role('inline cell barcode; native member unassigned',molecular_feature=CELL_FEATURE),Segment('2018 RP-H random tract','H'*9,placeholder=True)]
FINAL_LIBRARY=Construct([role('native P5/indexed end; member/geometry unassigned',molecular_feature=I5_FEATURE),*random_parts()[1:],role('genomic-DNA-derived OR RNA-derived insert; methylation encoded in base sequence'),role('Adaptase-dependent span and native P7/index end; member/geometry unassigned',molecular_feature=I7_FEATURE)],name='snmC2T-seq mixed DNA/cDNA library')
SEQ_PRIMERS=()
FINAL_CAPTION='One pool of independent DNA and RNA-derived molecules; RNA/DNA identity is a read-level computation, not a covalent junction or extra modality barcode.'
SEQUENCING_INTRO='Selected brain nuclear finish: NovaSeq6000 S4 PE150. Index cycle allocation remains unassigned.'
SEQUENCING_ENDING=SEQUENCING_INTRO
SHOW_READ_CONTENT=False
READ_SPECS=(sp.read_spec('Read 1',150,'Mixed genomic-DNA/RNA-derived sequence with inline identifier/priming region; exact native cycle offsets unassigned','extract_then_include'),sp.read_spec('Read 2',150,'Mixed genomic-DNA/RNA-derived sequence; modalities partitioned computationally by mCH','include'),sp.read_spec('Index 1 (i7)','unassigned','Native pool index acquisition unassigned','undetermined'),sp.read_spec('Index 2 (i5)','unassigned','Native pool index acquisition unassigned','undetermined'))
def oligos():
 yield '<b>Brain dT30VN_4</b> — 5′-/5SpC3/'+DT30VN+'-3′'
 yield '<b>N6_2</b> — 5′-/5SpC3/'+N6+'-3′'
 yield '<b>TSO_3</b> — 5′-/5SpC3/'+TSO+'-3′'
 yield '<b>ISPCR23_2</b> — 5′-/5SpC3/'+ISPCR23+'-3′'

def source_oligo(sequence,label):
 sc=Scene()
 if sequence.endswith('rGrG+G'):
  parts=[Segment(label+' DNA/dU tract',sequence[:-6],placeholder=True,note='U is deoxyuridine in this DNA tract'),Segment('two terminal RNA guanosines','GG',note='two riboguanosine residues, written rGrG'),Segment('terminal LNA guanosine','G',note='source terminal +G is LNA')]
 else:
  parts=[Segment(label,sequence,placeholder=True,note='U is deoxyuridine in the DNA primer')]
 sc.strand('oligo',parts,label=label,mod5='SpC3')
 return sc

def mixture(stage):
 sc=Scene();sc.strand('genome',[role('genomic DNA; endogenous methylation')],label='independent genomic DNA population')
 sc.strand('RNA',[role('nuclear RNA')],label='independent RNA templates; later physical fate unassigned')
 if stage!='raw':
  ends='native primer-derived ends with dU/C3; terminal retention unassigned' if stage!='treated' else 'primer-derived end chemistry after UDG; abasic/cleavage boundary unassigned'
  sc.strand('cDNA',[role(ends),role('RNA-derived DNA body; newly synthesized cytosines are5mC'),role('TSO-derived end; '+ends if stage in ('switched','amplified','treated') else 'nascent3′ cDNA end')],label=('amplified double-strand cDNA population' if stage in ('amplified','treated') else 'single-strand cDNA'))
 return sc

def converted_state():
 sc=Scene()
 sc.strand('genome',[role('converted original genomic DNA: unmethylated C→physical U; methyl-C retained')],label='genomic DNA population')
 sc.strand('cDNA',[role('RNA-derived original cDNA body:5mC retained through bisulfite; primer-end abasic cleavage unassigned')],label='separate cDNA population; no covalent genomic/RNA join')
 return sc

def copied_state(*,adapted=False):
 sc=converted_state()
 for label in ('genomic DNA','RNA-derived DNA'):
  parts=[*random_parts(),role('new copy of '+label+'; base sequence encodes converted-parent methylation, new copies not assumed chemically methylated')]
  if adapted:parts.append(role('native Adaptase3′ end; sequence/extent unassigned'))
  sc.strand(label+' copy',parts,label='expected library-copy population; original/copy occupancy unassigned' if not adapted else 'single-strand native Adaptase/PCR precursor',mod5='SpC3')
 return sc

def post_pcr_state():
 sc=Scene.duplex(list(FINAL_LIBRARY),label='expected amplified library population',unpaired=tuple(p.name+"'" for p in FINAL_LIBRARY if p.is_role_token()))
 sc.strand('residual0',[role('possible co-carried original genomic DNA; U/C chemistry/integrity/partition unassigned')],label='co-carried molecular lineage; no physical separation asserted')
 sc.strand('residual1',[role('possible co-carried original RNA-derived methylated cDNA; chemistry/integrity/partition unassigned')],label='co-carried molecular lineage; no physical separation asserted')
 sc.strand('residual2',[role('possible co-carried5′C3 RP-H/Adaptase genomic-copy templates; integrity/partition unassigned')],label='co-carried molecular lineage; no physical separation asserted')
 sc.strand('residual3',[role('possible co-carried5′C3 RP-H/Adaptase RNA-derived-copy templates; integrity/partition unassigned')],label='co-carried molecular lineage; no physical separation asserted')
 return sc

POST_PCR_STATE=post_pcr_state()

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 nucleus=g.input('nucleus','brain nucleus; RNA/methylome subworkflow before optional NOMe front end',MaterialKind.BIOLOGICAL,state=mixture('raw'),primary=True)
 primer_sc=Scene();primer_sc.strand('dT',[Segment('dT30VN_4',DT30VN,placeholder=True,note='DNA deoxyuridines as printed')],label='brain dT primer',mod5='SpC3');primer_sc.strand('N6',[Segment('N6_2',N6,placeholder=True,note='DNA deoxyuridines and randomN6 as printed')],label='nuclear random RT primer',mod5='SpC3')
 primers=g.input('rt_primers','source dT30VN_4 plus nuclear N6_2',MaterialKind.OLIGO,state=primer_sc)
 rt=g.input('rt','10U SuperScriptII',MaterialKind.CATALYST)
 mdctp=g.input('methyl_dctp','5-methyl-dCTP replaces dCTP',MaterialKind.CONDITION)
 first=g.product('first','methylated first-strand cDNA alongside RNA and genome',MaterialKind.MOLECULAR,state=mixture('first'))
 g.react('reverse_transcribe','RT25°C5min,42°C90min,70°C15min',StepKind.REVERSE_TRANSCRIPTION,((nucleus,InputRole.TEMPLATE),(primers,InputRole.STRUCTURAL),(rt,InputRole.CATALYST),(mdctp,InputRole.CONDITION)),(first,))
 tso=g.input('tso','source TSO_3 with5′C3 and rGrG+G tail',MaterialKind.OLIGO,state=source_oligo(TSO,'TSO_3 as written'))
 switched=g.product('switched','template-switched methylated cDNA reaction',MaterialKind.MOLECULAR,state=mixture('switched'))
 g.react('template_switch','Template switch nascent cDNA during RT',StepKind.TEMPLATE_SWITCHING,((first,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(switched,))
 ispcr=g.input('ispcr','source ISPCR23_2',MaterialKind.OLIGO,state=source_oligo(ISPCR23,'ISPCR23_2 as written'))
 kapa2g=g.input('kapa2g','KAPA2G Robust HotStart',MaterialKind.CATALYST)
 amplified=g.product('amplified','amplified methylated cDNA beside genome/RNA',MaterialKind.MOLECULAR,state=mixture('amplified'))
 g.react('cdna_pcr','Brain nuclear cDNA amplification14cycles; cultured nuclei15 is a separate condition',StepKind.PCR,((switched,InputRole.TEMPLATE),(ispcr,InputRole.STRUCTURAL),(kapa2g,InputRole.CATALYST),(mdctp,InputRole.CONDITION)),(amplified,),geometry=PCRGeometry.unavailable(switched.state,reason='Actual switched single-strand cDNA mixture; dU-containing primer-end topology not fully reconstructed.',symmetric_species=True))
 udg=g.input('udg','brain UDG-only primer treatment; no added EndoVIII',MaterialKind.CATALYST)
 treated=g.product('treated','UDG-treated mixed genome/cDNA reaction',MaterialKind.MOLECULAR,state=mixture('treated'))
 abasic=g.product('abasic_oligos','unused-oligo abasic products; later alkaline cleavage',MaterialKind.MOLECULAR)
 g.react('udg_treat','Treat unused dU-containing DNA oligos with UDG',StepKind.DIGESTION,((amplified,InputRole.TEMPLATE),(udg,InputRole.CATALYST)),(treated,abasic))
 g.mark_terminal(Terminal(abasic.id,Disposition.UNRESOLVED,reason='UDG creates abasic sites; later alkaline desulfonation cleaves them. Chemical cleavage does not assign physical removal.'))
 bs=g.input('bisulfite','25µl CT conversion reagent; nuclear modified2018 method',MaterialKind.CONDITION)
 converted=g.product('converted','mixed converted genomic DNA and protected methylated cDNA',MaterialKind.MOLECULAR,state=converted_state())
 rna_residual=g.product('rna_residual','RNA templates/hydrolysis products after conversion; chemical and physical fate unassigned',MaterialKind.MOLECULAR)
 g.react('convert','Bisulfite convert mixed DNA; original genome carries physical U',StepKind.CHEMICAL_CONVERSION,((treated,InputRole.TEMPLATE),(bs,InputRole.CONDITION)),(converted,rna_residual))
 g.mark_terminal(Terminal(rna_residual.id,Disposition.UNRESOLVED,reason='RNA template/hydrolysis fate and partition during conversion/recovery are unassigned; no physical discard is asserted.'))
 binding=g.input('binding','90µl M-binding buffer; native2018 wash/desulfonation recovery',MaterialKind.SELECTION)
 recovered=g.product('recovered','desulfonated/recovered mixed original DNA',MaterialKind.MOLECULAR,state=converted.state)
 bs_waste=g.product('bs_waste','conversion-recovery non-retained fraction',MaterialKind.WASTE)
 g.react('recover_bs','Bind, desulfonate/alkaline-cleave abasic oligos, wash and elute',StepKind.CLEANUP,((converted,InputRole.TEMPLATE),(binding,InputRole.SELECTION)),(recovered,bs_waste));g.discard(bs_waste,'Source-described physical conversion-recovery fraction; exact species partition unassigned.')
 rp=g.input('random_primer','native2018 RP-H inline-indexed random primer; selected member unassigned',MaterialKind.OLIGO)
 heat=g.input('denature','inherited2018 denaturation95°C3min/ice2min',MaterialKind.CONDITION)
 denatured=g.product('denatured','denatured mixed original DNA with RP-H',MaterialKind.MOLECULAR,state=recovered.state)
 g.react('denature_originals','Denature recovered DNA before native random extension',StepKind.THERMAL_TREATMENT,((recovered,InputRole.TEMPLATE),(rp,InputRole.STRUCTURAL),(heat,InputRole.CONDITION)),(denatured,))
 klenow=g.input('klenow','native2018 Klenow exo-minus extension mix',MaterialKind.CATALYST)
 primed=g.product('primed','native RP-H copies of both molecule classes',MaterialKind.MOLECULAR,state=copied_state())
 g.react('random_prime','Inherited2018 random extension4°C5min,25°C5min,37°C60min',StepKind.EXTENSION,((denatured,InputRole.TEMPLATE),(klenow,InputRole.CATALYST)),(primed,))
 exo=g.input('exo','ExonucleaseI',MaterialKind.CATALYST);rsap=g.input('rsap','rSAP',MaterialKind.CATALYST)
 exo_treated=g.product('exo_treated','treated expected mixed DNA copies',MaterialKind.MOLECULAR,state=primed.state)
 products=g.product('exo_products','ExoI/rSAP chemical products; physical recovery unassigned',MaterialKind.MOLECULAR)
 g.react('exo_rsap','Native2018 ExoI/rSAP37°C30min',StepKind.DIGESTION,((primed,InputRole.TEMPLATE),(exo,InputRole.CATALYST),(rsap,InputRole.CATALYST)),(exo_treated,products));g.mark_terminal(Terminal(products.id,Disposition.UNRESOLVED,reason='Chemical inactivation/hydrolysis is distinct from SPRI removal.'))
 others=g.input('other_reactions','other native inline-indexed reactions for eight-way2018 pooling',MaterialKind.MOLECULAR)
 pooled=g.product('pooled','pooled mixed libraries; selected lineage shown',MaterialKind.MOLECULAR,state=exo_treated.state)
 g.react('pool_eight','Inherited2018 eight-way inline-indexed reaction pooling',StepKind.POOL,((exo_treated,InputRole.TEMPLATE),(others,InputRole.TEMPLATE)),(pooled,))
 spri=g.input('spri','SPRI beads',MaterialKind.SELECTION)
 purified=g.product('purified','native copies after0.8×SPRI recovery',MaterialKind.MOLECULAR,state=pooled.state)
 waste=g.product('priming_waste','SPRI non-retained fraction',MaterialKind.WASTE)
 g.react('purify_primed','Inherited0.8×SPRI before Adaptase',StepKind.CLEANUP,((pooled,InputRole.TEMPLATE),(spri,InputRole.SELECTION)),(purified,waste));g.discard(waste,'Physical source-method cleanup fraction.')
 adapt_heat=g.input('adapt_heat','95°C3min then ice2min',MaterialKind.CONDITION)
 ss=g.product('ss_copies','denatured copies before Adaptase',MaterialKind.MOLECULAR,state=purified.state)
 g.react('denature_copies','Denature recovered copies before Adaptase',StepKind.THERMAL_TREATMENT,((purified,InputRole.TEMPLATE),(adapt_heat,InputRole.CONDITION)),(ss,))
 adaptase=g.input('adaptase','native Adaptase module',MaterialKind.CATALYST);adapter=g.input('adapter','native Adaptase3′ completion end; bases unassigned',MaterialKind.ADAPTER)
 adapted=g.product('adapted','native Adaptase-completed ss copies',MaterialKind.MOLECULAR,state=copied_state(adapted=True))
 g.react('adaptase_complete','Inherited Adaptase37°C30min,95°C2min',StepKind.ADAPTER_ATTACHMENT,((ss,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(adaptase,InputRole.CATALYST)),(adapted,))
 f=g.input('pcr_forward','native2018-family P5 indexed PCR primer; member unassigned',MaterialKind.OLIGO);r=g.input('pcr_reverse','native2018-family P7 indexed PCR primer; member unassigned',MaterialKind.OLIGO);kapa=g.input('kapa','KAPA HiFi library PCR',MaterialKind.CATALYST)
 amp=g.product('library_amp','native indexed mixed-library pool',MaterialKind.LIBRARY,state=POST_PCR_STATE)
 g.react('library_pcr','Inherited2018 library amplification15cycles',StepKind.PCR,((adapted,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(kapa,InputRole.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(adapted.state,reason='Actual native ss Adaptase products of two molecule classes; selected primer members and adapter-site geometry unassigned.'))
 clean1=g.product('clean1','first recovered library pool',MaterialKind.LIBRARY,state=POST_PCR_STATE);w1=g.product('waste1','first cleanup non-retained fraction',MaterialKind.WASTE)
 g.react('cleanup1','Inherited final0.8×SPRI cleanup1',StepKind.CLEANUP,((amp,InputRole.TEMPLATE),(spri,InputRole.SELECTION)),(clean1,w1));g.discard(w1,'Physical cleanup fraction.')
 clean2=g.product('clean2','second recovered library pool',MaterialKind.LIBRARY,state=POST_PCR_STATE);w2=g.product('waste2','second cleanup non-retained fraction',MaterialKind.WASTE)
 g.react('cleanup2','Inherited final0.8×SPRI cleanup2',StepKind.CLEANUP,((clean1,InputRole.TEMPLATE),(spri,InputRole.SELECTION)),(clean2,w2));g.discard(w2,'Physical cleanup fraction.')
 library=g.product('library',FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=POST_PCR_STATE);w3=g.product('waste3','third cleanup non-retained fraction',MaterialKind.WASTE)
 g.react('cleanup3','Inherited final0.8×SPRI cleanup3',StepKind.CLEANUP,((clean2,InputRole.TEMPLATE),(spri,InputRole.SELECTION)),(library,w3));g.discard(w3,'Physical cleanup fraction.')
 g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=(),read_lengths={'Read 1':150,'Read 2':150},limitation='Brain nuclear NovaSeq6000 S4 PE150; native sequencing primer/index allocation unassigned. RNA/DNA partitioning by mCH is computational.'))
 g.declare_reviewed_scope(('nucleus','rt_primers','rt','methyl_dctp','reverse_transcribe','first','tso','template_switch','switched','ispcr','kapa2g','cdna_pcr','amplified','udg','udg_treat','treated','bisulfite','convert','converted','binding','recover_bs','recovered','random_primer','denature','denature_originals','denatured','klenow','random_prime','primed','exo','rsap','exo_rsap','exo_treated','other_reactions','pool_eight','pooled','spri','purify_primed','purified','adapt_heat','denature_copies','ss_copies','adaptase','adapter','adaptase_complete','adapted','pcr_forward','pcr_reverse','kapa','library_pcr','library_amp','cleanup1','clean1','cleanup2','clean2','cleanup3','library'),SOURCE)
 return g.require_publishable()

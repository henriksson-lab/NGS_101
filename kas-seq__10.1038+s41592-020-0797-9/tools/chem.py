"""Source-qualified bulk and low-input KAS-seq; no imported kit recipe."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
from chemdraw import Construct, Scene, Segment, feature
from pcr import PCRGeometry
from protocol_model import molecule, role
from protocol_graph import InputRole as R, MaterialKind as K, ProtocolGraph, SequencingPlan, StepKind as S
import seqprimers as sp
PAGE_FILE='kas-seq.html'
TITLE='KAS-seq — reversible guanine tagging and DNA enrichment'
NOTES='01_kas-seq.html'
SOURCE='<a href="https://doi.org/10.1038/s41592-020-0797-9">Wu et al. (2020)</a>, publisher Supplementary Protocol sections1–3; associated accepted manuscript PMC7205578; <a href="https://doi.org/10.1038/s41592-020-0881-1">Author Correction</a>.'
SUMMARY='N3-kethoxal tags accessible guanines in single-stranded DNA. Azide–DBCO coupling adds biotin for streptavidin enrichment. Bulk sonication/Swift library preparation and low-input tagmentation/on-bead PCR remain distinct branches.'
CAVEAT=('INFERRED — the displayed guanine represents a possible accessible site, not an assigned genomic coordinate or complete labeling. '
'Only fragment enrichment is mapped; the sequencing read does not directly identify the adduct position. '
'The graph follows the published Supplementary Protocol: bulk click input2µg and RNaseA15min; low-input D4074 and Tn5 15027865. '
'The associated accepted-manuscript Methods instead describes1µg/5min, D4068 and10527865. Those source/version differences are not silently merged. '
'The published correction changes the synthetic-oligo characterization concentration from100mM to100µM; live-cell5mM N3-kethoxal is unchanged. '
'Exact historical kit internals and adapters are unassembled; the bulk kit remains an explicit unvalidated boundary.')

@dataclass(frozen=True)
class TaggedDNA:
    """Local chemical label state, separate from unassigned genomic positions.

    These declarations carry the source-described functional group, not an atom-level
    reaction proof or a claim that this position is known in a biological molecule.
    """
    guanine_state: str
    genomic_positions: str='unassigned; representative accessible G only'
    end_architecture: str='unassigned genomic ends'
    rna_component: str='not represented; no absence claim'
    def __post_init__(self):
        if self.guanine_state not in {'accessible G','azido-kethoxal G adduct','biotin-clicked G adduct','heat-deprotected DNA'}:
            raise ValueError('unknown KAS functional-group state')
    def rows(self):
        sc=Scene()
        segments=[role('genomic flank'),Segment('representative accessible guanine','G',inferred=True),role('genomic flank after selected site')]
        if self.end_architecture=='Tn5-tagged; polymers unresolved':
            segments=[role('unassembled loaded-adapter-derived end'),*segments,role('unassembled opposite loaded-adapter-derived end')]
        sc.strand('dna',segments,label='INFERRED — representative genomic site')
        sc.note('dna','End architecture: '+self.end_architecture)
        if self.rna_component!='not represented; no absence claim':
            sc.strand('rna',[role(self.rna_component)],label='possible co-isolated RNA material; amount unassigned')
        sc.note('dna',self.guanine_state+'; '+self.genomic_positions)
        if self.guanine_state=='azido-kethoxal G adduct':sc.note('dna','Source: kethoxal adduct at exposed guanine N1/N2, with azide click handle. Watson–Crick-protected guanines are not asserted labeled.')
        if self.guanine_state=='biotin-clicked G adduct':sc.note('dna','Source: DBCO-PEG4-biotin couples to azide; reversible kethoxal linkage remains on guanine until source heating.')
        if self.guanine_state=='heat-deprotected DNA':sc.note('dna','95 C treatment removes the reversible label operationally; released tag chemistry and conversion yield are not assigned.')
        return sc.rows()

BULK_LIBRARY=Construct((role('unassembled Swift30024 library end'),role('genomic DNA from selected fraction'),role('unassembled opposite Swift library end')),name='bulk KAS library; kit architecture unresolved')
LOW_LIBRARY=Construct((role('unassembled Tn5/i5-primer-derived arm'),Segment('i5 index; length and bases unassigned','[i5 index]',placeholder=True,feature=feature('kas_low_i5','sample_index','unknown')),role('genomic DNA from selected fraction'),Segment('i7 index; length and bases unassigned','[i7 index]',placeholder=True,feature=feature('kas_low_i7','sample_index','unknown')),role('unassembled i7-primer-derived arm')),name='low-input KAS library; adapter/index members unresolved')
FINAL_LIBRARY=Construct((role('branch-specific unassembled library end'),role('genomic DNA from input or enriched fraction'),role('branch-specific unassembled opposite end')),name='role-only KAS library overview')
FINAL_CAPTION='Role-only overview; bulk Swift and low-input Tn5 endpoints are kept separate in the graph. No shared adapter architecture is asserted.'
SEQ_PRIMERS=()
SEQUENCING_UNAVAILABLE='Exact source-version adapter and sequencing-primer polymers remain unassembled. Bulk is source-associated NextSeq500 single-end80bp; publisher low-input protocol states single-end without assigning a cycle length. No low-input cycle mask is borrowed from the bulk description.'
READ_SPECS=(sp.read_spec('Read 1','bulk80; low-input unassigned','genomic sequence of selected DNA fragments; adduct is not a base-resolved readout','include'),)
SEQUENCING_ENDING='Both input controls and streptavidin-enriched fractions are sequenced; source aims for30million single-end reads per sample. Reads extended to150bp in analysis are computational intervals, not150 sequenced bases.'


def _cleanup(g,prefix,substrate,label,*,output_kind=K.MOLECULAR,retained_state=None):
 condition=g.input(prefix+'_condition',label,K.CONDITION)
 state=substrate.state if retained_state is None else retained_state
 recovered=g.product(prefix+'_dna','recovered DNA',output_kind,state=state)
 waste=g.product(prefix+'_waste','removed recovery fractions',K.WASTE)
 g.react(prefix,label,S.CLEANUP,((substrate,R.TEMPLATE),(condition,R.CONDITION)),(recovered,waste))
 g.discard(waste,'Source purification separates discarded fractions; their exact composition and recovery efficiency are unassigned.')
 return recovered


def _front(g,prefix,cells,n3,dbco,rnase,*,extraction,rnase_minutes,click_amount):
 labelled=g.product(prefix+'_labelled','N3-kethoxal-tagged DNA in cells',K.MOLECULAR,state=TaggedDNA('azido-kethoxal G adduct'))
 g.react(prefix+'_label','Label accessible guanines:5mM N3-kethoxal,10min37C',S.CHEMICAL_CONVERSION,((cells,R.TEMPLATE),(n3,R.CONDITION)),(labelled,))
 kit=g.input(prefix+'_extraction',extraction,K.CONDITION)
 extracted=g.product(prefix+'_extracted','isolated adduct-bearing DNA',K.MOLECULAR,state=labelled.state)
 waste=g.product(prefix+'_extraction_waste','DNA-extraction removed fractions',K.WASTE)
 g.react(prefix+'_extract','Isolate DNA into25mM borate pH7',S.EXTRACTION,((labelled,R.TEMPLATE),(kit,R.CONDITION)),(extracted,waste));g.discard(waste,'Source DNA-isolation kit removes unretained material; label-containing RNA contamination is not assumed absent.')
 clicked=g.product(prefix+'_clicked','biotin-clicked DNA preparation with possible RNA carryover',K.MOLECULAR,state=TaggedDNA('biotin-clicked G adduct',rna_component='possible co-isolated RNA carrying RNaseA-sensitive phosphodiester bonds'))
 g.react(prefix+'_click',click_amount+'; DBCO-PEG4-biotin37C90min',S.CHEMICAL_CONVERSION,((extracted,R.TEMPLATE),(dbco,R.CONDITION)),(clicked,))
 treated=g.product(prefix+'_rnase_treated','RNaseA-treated DNA preparation',K.MOLECULAR,state=TaggedDNA('biotin-clicked G adduct',rna_component='RNA digestion fragments; completeness not quantified'))
 g.react(prefix+'_rnase','RNaseA37C for'+str(rnase_minutes)+'min',S.DIGESTION,((clicked,R.TEMPLATE),(rnase,R.CATALYST)),(treated,))
 return _cleanup(g,prefix+'_click_cleanup',treated,'DNA Clean & Concentrator D4043; borate elution',retained_state=TaggedDNA('biotin-clicked G adduct'))


def _capture(g,prefix,dna,*,reserve,bead_volume):
 aliquot=g.input(prefix+'_split_condition',reserve,K.CONDITION)
 control=g.product(prefix+'_input','reserved input DNA control',K.MOLECULAR,state=dna.state)
 capture=g.product(prefix+'_capture_input','DNA aliquot for streptavidin enrichment',K.MOLECULAR,state=dna.state)
 g.react(prefix+'_split',reserve,S.SPLIT,((dna,R.TEMPLATE),(aliquot,R.CONDITION)),(control,capture))
 beads=g.input(prefix+'_beads',str(bead_volume)+'µL MyOne streptavidinC1;15min room temperature;5washes',K.SELECTION)
 enriched=g.product(prefix+'_enriched','bead-bound enriched DNA',K.MOLECULAR,state=dna.state)
 waste=g.product(prefix+'_unretained','capture supernatant and wash fractions',K.WASTE)
 g.react(prefix+'_capture','Capture biotin-bearing DNA; wash five times',S.AFFINITY_CAPTURE,((capture,R.TEMPLATE),(beads,R.SELECTION)),(enriched,waste));g.discard(waste,'Physical unretained capture/wash fractions; not a claim of perfect exclusion of untagged molecules.')
 return control,enriched


def _bulk_library(g,prefix,dna,cycles):
 kit=g.input(prefix+'_kit','Swift Accel-NGS Methyl-seq30024; historical formulation unassembled',K.CONDITION)
 library=g.product(prefix+'_library','bulk library; source '+cycles+'PCR cycles;200–500bp',K.LIBRARY,state=BULK_LIBRARY)
 g.react(prefix+'_library_prep','UNRESOLVED KIT BOUNDARY — Swift30024 library construction;'+cycles+'PCR cycles',S.MAGIC,((dna,R.TEMPLATE),(kit,R.CONDITION)),(library,),reason='Source names kit30024 and resulting sequencing library; exact historical adapter/end-processing/PCR geometry not reconstructed in this capped pass. No bisulfite or Adaptase recipe imported.')
 g.sequence(library,SequencingPlan(library=BULK_LIBRARY,limitation='Bulk kit boundary unresolved; exact run primers and index mask unassigned. Associated accepted manuscript specifies NextSeq500 SE80; publisher protocol specifies single-end.',read_lengths={'Read 1':80},acquired_roles=('Read 1',)))


def _low_library(g,prefix,dna,polymerase,i5,i7):
 filled=g.product(prefix+'_extended','initially extended and heated low-input template',K.MOLECULAR,state=dna.state)
 g.react(prefix+'_initial_extension','Initial72C5min extension then95C10min',S.EXTENSION,((dna,R.TEMPLATE),(polymerase,R.CATALYST)),(filled,))
 amplified=g.product(prefix+'_amplified','low-input indexed PCR products',K.MOLECULAR,state=LOW_LIBRARY)
 g.react(prefix+'_pcr','i5/i7 PCR:15cycles98C10s/60C30s/72C1min',S.PCR,((filled,R.TEMPLATE),(i5,R.STRUCTURAL),(i7,R.STRUCTURAL),(polymerase,R.CATALYST)),(amplified,),geometry=PCRGeometry.unavailable(molecule('low-input PCR substrate',(role('source-named Tn5-tagged genomic DNA; exact end polymers unresolved'),)),reason='Historical Tn5 adapter and source20027213 primer members remain unassembled; binding sites are not invented.'))
 library=_cleanup(g,prefix+'_cleanup',amplified,'Remove beads if enriched; MinElute PCR purification;200–600bp',output_kind=K.LIBRARY)
 g.sequence(library,SequencingPlan(library=LOW_LIBRARY,limitation='Publisher low-input protocol specifies single-end; exact cycle lengths, Tn5 adapter polymers, i5/i7 members and loaded run primers remain unassigned.',acquired_roles=('Read 1',)))


def protocol_graph():
 g=ProtocolGraph(TITLE,authority='Published Supplementary Protocol, MOESM1 sections1–3; associated manuscript Results/Methods; correction2020')
 bulk=g.input('bulk_cells','bulk1–5million mammalian cells',K.BIOLOGICAL,state=TaggedDNA('accessible G'),primary=True)
 low=g.input('low_cells','alternative low-input1000/5000/10000cells',K.BIOLOGICAL,state=TaggedDNA('accessible G'),primary=True)
 n3=g.input('n3_kethoxal','5mM N3-kethoxal live-cell labeling condition',K.CONDITION)
 dbco=g.input('dbco','DBCO-PEG4-biotin click reagent in borate/PBS',K.CONDITION)
 rnase=g.input('rnase_a','RNaseA',K.CATALYST)
 bulk_dna=_front(g,'bulk',bulk,n3,dbco,rnase,extraction='PureLink K182001;25mM borate elution',rnase_minutes=15,click_amount='2µg DNA in100µL click mixture')
 low_dna=_front(g,'low',low,n3,dbco,rnase,extraction='Quick-DNA Microprep Plus D4074;43µL borate elution',rnase_minutes=5,click_amount='42.5µL DNA in50µL click mixture')
 sonication=g.input('bulk_sonication_condition','Bioruptor Pico30cycles30s-on/off;150–350bp',K.CONDITION)
 fragmented=g.product('bulk_fragmented','sonicated biotinylated genomic fragments',K.MOLECULAR,state=bulk_dna.state)
 g.react('bulk_fragment','Sonicate1µg tagged DNA in100µL borate',S.FRAGMENTATION,((bulk_dna,R.TEMPLATE),(sonication,R.CONDITION)),(fragmented,))
 control,enriched=_capture(g,'bulk',fragmented,reserve='Reserve5µL of100µL input;95µL for enrichment',bead_volume=10)
 heat=g.input('bulk_heat','95C10min; immediately ice',K.CONDITION)
 for prefix,dna in [('bulk_input',control),('bulk_enriched',enriched)]:
  heated=g.product(prefix+'_heated','heat-treated DNA; reversible label removed',K.MOLECULAR,state=TaggedDNA('heat-deprotected DNA'))
  g.react(prefix+'_heat','Heat input or elute enriched DNA:95C10min',S.THERMAL_TREATMENT,((dna,R.TEMPLATE),(heat,R.CONDITION)),(heated,))
  _bulk_library(g,prefix,heated,'7–8' if prefix=='bulk_input' else '12–14')
 tn5=g.input('low_tn5','Tagmentation DNA Enzyme15027865',K.CATALYST)
 buffer=g.input('low_tn5_buffer','Buffer15027866;37C30min;1.5/2/5µL enzyme for1000/5000/10000cells',K.CONDITION)
 adapters=g.input('low_loaded_adapters','loaded adapter polymers of source-named Tn5 reagent; historical identity unassembled',K.ADAPTER,state=molecule('unassembled loaded adapters',(role('exact loaded Tn5 adapter polymers unresolved'),)))
 tagged=g.product('low_tagmented','Tn5-tagged biotinylated genomic DNA',K.MOLECULAR,state=TaggedDNA('biotin-clicked G adduct',end_architecture='Tn5-tagged; polymers unresolved'))
 g.react('low_tagment','Tagment biotinylated DNA before enrichment',S.TAGMENTATION,((low_dna,R.TEMPLATE),(adapters,R.STRUCTURAL),(tn5,R.CATALYST),(buffer,R.CONDITION)),(tagged,))
 purified=_cleanup(g,'low_tagment_cleanup',tagged,'DNA Clean & Concentrator D4043;55µL borate elution')
 low_control,low_enriched=_capture(g,'low',purified,reserve='Reserve5µL of55µL input;50µL for enrichment (not5%)',bead_volume=5)
 polymerase=g.input('low_q5','NEBNext UltraII Q5 M0544S',K.CATALYST)
 i5=g.input('low_i5','Illumina20027213 i5 index-primer member unassigned',K.OLIGO,state=molecule('i5 primer',(Segment('source-named i5 primer; exact member unresolved','[unassembled i5 primer]',placeholder=True,feature=feature('kas_low_i5','sample_index','unknown')),)))
 i7=g.input('low_i7','Illumina20027213 i7 index-primer member unassigned',K.OLIGO,state=molecule('i7 primer',(Segment('source-named i7 primer; exact member unresolved','[unassembled i7 primer]',placeholder=True,feature=feature('kas_low_i7','sample_index','unknown')),)))
 _low_library(g,'low_input',low_control,polymerase,i5,i7)
 _low_library(g,'low_enriched',low_enriched,polymerase,i5,i7)
 # Curator-written source inventory; do not derive coverage from graph keys.
 targets=('bulk_cells','low_cells','n3_kethoxal','dbco','rnase_a','bulk_labelled','bulk_label','bulk_extraction','bulk_extracted','bulk_extraction_waste','bulk_extract','bulk_clicked','bulk_click','bulk_rnase_treated','bulk_rnase','bulk_click_cleanup_condition','bulk_click_cleanup_dna','bulk_click_cleanup_waste','bulk_click_cleanup','low_labelled','low_label','low_extraction','low_extracted','low_extraction_waste','low_extract','low_clicked','low_click','low_rnase_treated','low_rnase','low_click_cleanup_condition','low_click_cleanup_dna','low_click_cleanup_waste','low_click_cleanup','bulk_sonication_condition','bulk_fragmented','bulk_fragment','bulk_split_condition','bulk_input','bulk_capture_input','bulk_split','bulk_beads','bulk_enriched','bulk_unretained','bulk_capture','bulk_heat','bulk_input_heated','bulk_input_heat','bulk_input_kit','bulk_input_library','bulk_input_library_prep','bulk_enriched_heated','bulk_enriched_heat','bulk_enriched_kit','bulk_enriched_library','bulk_enriched_library_prep','low_tn5','low_tn5_buffer','low_loaded_adapters','low_tagmented','low_tagment','low_tagment_cleanup_condition','low_tagment_cleanup_dna','low_tagment_cleanup_waste','low_tagment_cleanup','low_split_condition','low_input','low_capture_input','low_split','low_beads','low_enriched','low_unretained','low_capture','low_q5','low_i5','low_i7','low_input_extended','low_input_initial_extension','low_input_amplified','low_input_pcr','low_input_cleanup_condition','low_input_cleanup_dna','low_input_cleanup_waste','low_input_cleanup','low_enriched_extended','low_enriched_initial_extension','low_enriched_amplified','low_enriched_pcr','low_enriched_cleanup_condition','low_enriched_cleanup_dna','low_enriched_cleanup_waste','low_enriched_cleanup')
 g.declare_reviewed_scope(targets,g.authority)
 return g.require_publishable()


def render_page():
 from batch_page import render
 return render(sys.modules[__name__])


def selftest():
 from checks import Check
 c=Check()
 protocol_graph();render_page()
 # This model contains no source-exact oligo polymer claims. Graph construction
 # checks branch connectivity and explicit source coverage; source details are
 # source-qualified declarations, not checks that merely restate prose.
 c.report()

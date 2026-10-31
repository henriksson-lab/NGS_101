"""Source-qualified bulk and low-input KAS-seq; no imported kit recipe."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
from chemdraw import Construct, Scene, Segment, feature
from pcr import PCRGeometry
from protocol_model import molecule, role
from protocol_graph import InputRole as R, MaterialKind as K, ProtocolGraph, SequencingPlan, StepKind as S, Terminal, Disposition
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



def _take(g,key,stock,label):
 part=g.product(key+'_selected',label,K.MOLECULAR,state=stock.state)
 rest=g.product(key+'_remainder','Actual untransferred source DNA stock after selected contribution; amount is branch-qualified',K.MOLECULAR,state=stock.state)
 g.react(key,label,S.SPLIT,((stock,R.TEMPLATE),),(part,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Source DNA stock is not entirely consumed by the selected allocation'))
 return part

def _loading_portion(g,prefix,stock):
 lib=g.product(prefix+'_sequenced','Actual sequencing loading portion; volume unassigned',K.LIBRARY,state=stock.state)
 rest=g.product(prefix+'_library_stock','Untransferred prepared library stock; amount unassigned',K.LIBRARY,state=stock.state)
 g.react(prefix+'_load','Actual selected library loading allocation',S.SPLIT,((stock,R.TEMPLATE),),(lib,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'No entire-stock consumption or physical pool is invented'))
 return lib

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
 g.react(prefix+'_label','Label accessible guanines:5mM N3-kethoxal,10min37C',S.CHEMICAL_CONVERSION,((cells,R.TEMPLATE),(n3,R.STRUCTURAL),(g.input(prefix+'_label_conditions','Actual 5 mM labelling medium; 10 min at 37°C',K.CONDITION),R.CONDITION)),(labelled,))
 kit=g.input(prefix+'_extraction',extraction,K.CONDITION)
 extracted=g.product(prefix+'_extracted','isolated adduct-bearing DNA',K.MOLECULAR,state=labelled.state)
 waste=g.product(prefix+'_extraction_waste','DNA-extraction removed fractions',K.WASTE)
 g.react(prefix+'_extract','Isolate DNA into25mM borate pH7',S.EXTRACTION,((labelled,R.TEMPLATE),(kit,R.CONDITION)),(extracted,waste));g.discard(waste,'Source DNA-isolation kit removes unretained material; label-containing RNA contamination is not assumed absent.')
 extracted=_take(g,prefix+'_click_dna',extracted,'Actual 2 µg DNA selected into 85 µl borate for bulk click reaction'if prefix=='bulk'else'Actual 42.5 µl of 43 µl DNA stock for low-input click reaction')
 clicked=g.product(prefix+'_clicked','biotin-clicked DNA preparation with possible RNA carryover',K.MOLECULAR,state=TaggedDNA('biotin-clicked G adduct',rna_component='possible co-isolated RNA carrying RNaseA-sensitive phosphodiester bonds'))
 g.react(prefix+'_click',click_amount+'; DBCO-PEG4-biotin37C90min',S.CHEMICAL_CONVERSION,((extracted,R.TEMPLATE),(dbco,R.STRUCTURAL),(g.input(prefix+'_click_conditions','Actual click reaction in borate/PBS; 37°C for 90 min at 500 rpm',K.CONDITION),R.CONDITION)),(clicked,))
 treated=g.product(prefix+'_rnase_treated','RNaseA-treated DNA preparation',K.MOLECULAR,state=TaggedDNA('biotin-clicked G adduct',rna_component='RNA digestion fragments; completeness not quantified'))
 g.react(prefix+'_rnase','RNaseA37C for'+str(rnase_minutes)+'min',S.DIGESTION,((clicked,R.TEMPLATE),(rnase,R.CATALYST)),(treated,))
 return _cleanup(g,prefix+'_click_cleanup',treated,'DNA Clean & Concentrator D4043; borate elution',retained_state=TaggedDNA('biotin-clicked G adduct'))


def _capture(g,prefix,dna,*,reserve,bead_volume):
 aliquot=g.input(prefix+'_split_condition',reserve,K.CONDITION)
 control=g.product(prefix+'_input','reserved input DNA control',K.MOLECULAR,state=dna.state)
 capture=g.product(prefix+'_capture_input','DNA aliquot for streptavidin enrichment',K.MOLECULAR,state=dna.state)
 g.react(prefix+'_split',reserve,S.SPLIT,((dna,R.TEMPLATE),(aliquot,R.CONDITION)),(control,capture))
 raw=g.input(prefix+'_beads_raw',str(bead_volume)+' µl actual C1 bead stock',K.MOLECULAR)
 prep=g.input(prefix+'_bead_preparation','Source three prewashes with '+str(50 if prefix=='bulk'else 25)+' µl 1× binding/wash buffer, then resuspend in '+str(95 if prefix=='bulk'else 50)+' µl 2× buffer',K.CONDITION)
 beads=g.product(prefix+'_beads','Actual source prewashed/resuspended C1 capture stock',K.SELECTION)
 prewash=g.product(prefix+'_bead_prewash_waste','Actual three prewash fractions',K.WASTE)
 g.react(prefix+'_prepare_beads','Prepare actual C1 beads before capture',S.CLEANUP,((raw,R.TEMPLATE),(prep,R.CONDITION)),(beads,prewash));g.discard(prewash,'Actual source prewash fractions')
 enriched=g.product(prefix+'_enriched','bead-bound enriched DNA',K.MOLECULAR,state=dna.state)
 waste=g.product(prefix+'_unretained','capture supernatant and wash fractions',K.WASTE)
 g.react(prefix+'_capture','Capture biotin-bearing DNA; wash five times',S.AFFINITY_CAPTURE,((capture,R.TEMPLATE),(beads,R.SELECTION)),(enriched,waste));g.discard(waste,'Physical unretained capture/wash fractions; not a claim of perfect exclusion of untagged molecules.')
 return control,enriched


def _bulk_library(g,prefix,dna,cycles):
 kit=g.input(prefix+'_kit','Actual Swift Accel-NGS Methyl-seq 30024 aggregate supplies; stock packaging/formulation unassembled; no bisulfite recipe inferred',K.MOLECULAR)
 ol=g.product(prefix+'_kit_oligos','Actual source kit oligo inventory for library preparation/PCR; whole bodies/functional partition unassembled',K.OLIGO)
 enzyme=g.product(prefix+'_kit_enzyme','Actual PCR/library enzyme classes from kit supplies; identities/amounts unassigned',K.CATALYST)
 nts=g.product(prefix+'_kit_nt','Actual polymerization nucleotide class from kit supplies; identities/doses unassembled',K.MOLECULAR)
 other=g.product(prefix+'_kit_other','Other actual kit supplies; identities/packaging unassembled',K.MOLECULAR)
 g.react(prefix+'_kit_inventory','Bookkeeping actual kit supply classes, not a premix factory or new independent doses',S.GENERIC,((kit,R.TEMPLATE),),(ol,enzyme,nts,other),changes_molecule=False)
 library=g.product(prefix+'_library','bulk library; source '+cycles+'PCR cycles;200–500bp',K.LIBRARY,state=BULK_LIBRARY)
 g.react(prefix+'_library_prep','UNRESOLVED KIT BOUNDARY — Swift30024 library construction;'+cycles+'PCR cycles',S.MAGIC,((dna,R.TEMPLATE),(ol,R.STRUCTURAL),(enzyme,R.CATALYST),(nts,R.STRUCTURAL),(other,R.STRUCTURAL)),(library,),reason='Source names kit30024 and resulting sequencing library; exact historical adapter/end-processing/PCR geometry not reconstructed in this capped pass. No bisulfite or Adaptase recipe imported.')
 g.sequence(_loading_portion(g,prefix,library),SequencingPlan(library=BULK_LIBRARY,limitation='Bulk kit boundary unresolved; exact run primers and index mask unassigned. Associated accepted manuscript specifies NextSeq500 SE80; publisher protocol specifies single-end.',read_lengths={'Read 1':80},acquired_roles=('Read 1',)))


def _low_library(g,prefix,dna,i5,i7):
 stock=g.input(prefix+'_q5_stock','Actual 25 µl NEBNext Ultra II Q5 Master Mix M0544S allocation for this source PCR',K.MOLECULAR)
 polymerase=g.product(prefix+'_q5_enzyme','Polymerase component of the same actual 25 µl Q5 Master Mix',K.CATALYST)
 nts=g.product(prefix+'_q5_nt','Nucleotide component of the same actual 25 µl Q5 Master Mix; no fresh dNTP dose',K.MOLECULAR)
 buf=g.product(prefix+'_q5_buffer','Buffer component of the same actual 25 µl Q5 Master Mix',K.CONDITION)
 g.react(prefix+'_q5_components','Copresent components of the actual prepared stock; no physical separation or extra additions',S.GENERIC,((stock,R.TEMPLATE),),(polymerase,nts,buf),changes_molecule=False)
 filled=g.product(prefix+'_extended','initially extended and heated low-input template',K.MOLECULAR,state=dna.state)
 g.react(prefix+'_initial_extension','Initial72C5min extension then95C10min',S.EXTENSION,((dna,R.TEMPLATE),(polymerase,R.CATALYST),(nts,R.STRUCTURAL),(buf,R.CONDITION)),(filled,))
 amplified=g.product(prefix+'_amplified','low-input indexed PCR products',K.MOLECULAR,state=LOW_LIBRARY)
 g.react(prefix+'_pcr','i5/i7 PCR:15cycles98C10s/60C30s/72C1min',S.PCR,((filled,R.TEMPLATE),(i5,R.STRUCTURAL),(i7,R.STRUCTURAL),(polymerase,R.CATALYST),(nts,R.STRUCTURAL),(buf,R.CONDITION)),(amplified,),geometry=PCRGeometry.unavailable(molecule('low-input PCR substrate',(role('source-named Tn5-tagged genomic DNA; exact end polymers unresolved'),)),reason='Historical Tn5 adapter and source20027213 primer members remain unassembled; binding sites are not invented.'))
 library=_cleanup(g,prefix+'_cleanup',amplified,'Remove beads if enriched; MinElute PCR purification;200–600bp',output_kind=K.LIBRARY)
 g.sequence(_loading_portion(g,prefix,library),SequencingPlan(library=LOW_LIBRARY,limitation='Publisher low-input protocol specifies single-end; exact cycle lengths, Tn5 adapter polymers, i5/i7 members and loaded run primers remain unassigned.',acquired_roles=('Read 1',)))


def protocol_graph():
 g=ProtocolGraph(TITLE,authority='Published Supplementary Protocol, MOESM1 sections1–3; associated manuscript Results/Methods; correction2020')
 bulk=g.input('bulk_cells','bulk1–5million mammalian cells',K.BIOLOGICAL,state=TaggedDNA('accessible G'),primary=True)
 low=g.input('low_cells','alternative low-input1000/5000/10000cells',K.BIOLOGICAL,state=TaggedDNA('accessible G'),primary=True)
 n3=g.input('n3_kethoxal','Actual N3-kethoxal reagent stock; 500 mM stock diluted to 5 mM per source branch reaction, aggregate allocation unassigned',K.MOLECULAR)
 dbco=g.input('dbco','Actual 20 mM DBCO-PEG4-biotin reagent stock; 5 µl bulk/2.5 µl low-input additions in borate/PBS',K.MOLECULAR)
 rnase=g.input('rnase_a','RNaseA',K.CATALYST)
 bulk_dna=_front(g,'bulk',bulk,n3,dbco,rnase,extraction='PureLink K182001;25mM borate elution',rnase_minutes=15,click_amount='2µg DNA in100µL click mixture')
 low_dna=_front(g,'low',low,n3,dbco,rnase,extraction='Quick-DNA Microprep Plus D4074;43µL borate elution',rnase_minutes=5,click_amount='42.5µL DNA in50µL click mixture')
 bulk_dna=_take(g,'bulk_shear_dna',bulk_dna,'Actual 1 µg of biotinylated DNA stock selected for 100 µl sonication')
 sonication=g.input('bulk_sonication_condition','Bioruptor Pico30cycles30s-on/off;150–350bp',K.CONDITION)
 fragmented=g.product('bulk_fragmented','sonicated biotinylated genomic fragments',K.MOLECULAR,state=bulk_dna.state)
 g.react('bulk_fragment','Sonicate1µg tagged DNA in100µL borate',S.FRAGMENTATION,((bulk_dna,R.TEMPLATE),(sonication,R.CONDITION)),(fragmented,))
 control,enriched=_capture(g,'bulk',fragmented,reserve='Reserve5µL of100µL input;95µL for enrichment',bead_volume=10)
 heat=g.input('bulk_heat','95C10min; immediately ice',K.CONDITION)
 for prefix,dna in [('bulk_input',control),('bulk_enriched',enriched)]:
  heated=g.product(prefix+'_heated','heat-treated DNA; reversible label removed',K.MOLECULAR,state=TaggedDNA('heat-deprotected DNA'))
  outputs=(heated,)
  if prefix=='bulk_enriched':
   residual=g.product('bulk_residual_beads','Actual residual C1 beads/material after 15 µl water heat-elution and magnetic transfer',K.MOLECULAR)
   outputs+=(residual,)
  g.react(prefix+'_heat','Heat input or elute enriched DNA in 15 µl water at 95°C for 10 min; magnetically transfer actual eluate',S.THERMAL_TREATMENT,((dna,R.TEMPLATE),(heat,R.CONDITION)),outputs)
  if prefix=='bulk_enriched':g.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,'Residual source beads excluded from actual eluate; composition/fate unassigned'))
  _bulk_library(g,prefix,heated,'7–8' if prefix=='bulk_input' else '12–14')
 tn5=g.input('low_tn5','Tagmentation DNA Enzyme15027865',K.CATALYST)
 buffer=g.input('low_tn5_buffer','Buffer15027866;37C30min;1.5/2/5µL enzyme for1000/5000/10000cells',K.CONDITION)
 adapters=g.input('low_loaded_adapters','loaded adapter polymers of source-named Tn5 reagent; historical identity unassembled',K.ADAPTER,state=molecule('unassembled loaded adapters',(role('exact loaded Tn5 adapter polymers unresolved'),)))
 low_dna=_take(g,'low_tagment_dna',low_dna,'Actual 23.5/23/20 µl selected from 25 µl according to 1000/5000/10000-cell branch; selected branch unassigned')
 tagged=g.product('low_tagmented','Tn5-tagged biotinylated genomic DNA',K.MOLECULAR,state=TaggedDNA('biotin-clicked G adduct',end_architecture='Tn5-tagged; polymers unresolved'))
 g.react('low_tagment','Tagment biotinylated DNA before enrichment',S.TAGMENTATION,((low_dna,R.TEMPLATE),(adapters,R.STRUCTURAL),(tn5,R.CATALYST),(buffer,R.CONDITION)),(tagged,))
 purified=_cleanup(g,'low_tagment_cleanup',tagged,'DNA Clean & Concentrator D4043;55µL borate elution')
 low_control,low_enriched=_capture(g,'low',purified,reserve='Reserve5µL of55µL input;50µL for enrichment (not5%)',bead_volume=5)
 i5=g.input('low_i5','Illumina20027213 i5 index-primer member unassigned',K.OLIGO,state=molecule('i5 primer',(Segment('source-named i5 primer; exact member unresolved','[unassembled i5 primer]',placeholder=True,feature=feature('kas_low_i5','sample_index','unknown')),)))
 i7=g.input('low_i7','Illumina20027213 i7 index-primer member unassigned',K.OLIGO,state=molecule('i7 primer',(Segment('source-named i7 primer; exact member unresolved','[unassembled i7 primer]',placeholder=True,feature=feature('kas_low_i7','sample_index','unknown')),)))
 _low_library(g,'low_input',low_control,i5,i7)
 _low_library(g,'low_enriched',low_enriched,i5,i7)
 # Curator-written source inventory; do not derive coverage from graph keys.
 targets=('bulk_cells','low_cells','n3_kethoxal','dbco','rnase_a','bulk_labelled','bulk_label','bulk_extraction','bulk_extracted','bulk_extraction_waste','bulk_extract','bulk_clicked','bulk_click','bulk_rnase_treated','bulk_rnase','bulk_click_cleanup_condition','bulk_click_cleanup_dna','bulk_click_cleanup_waste','bulk_click_cleanup','low_labelled','low_label','low_extraction','low_extracted','low_extraction_waste','low_extract','low_clicked','low_click','low_rnase_treated','low_rnase','low_click_cleanup_condition','low_click_cleanup_dna','low_click_cleanup_waste','low_click_cleanup','bulk_sonication_condition','bulk_fragmented','bulk_fragment','bulk_split_condition','bulk_input','bulk_capture_input','bulk_split','bulk_beads','bulk_enriched','bulk_unretained','bulk_capture','bulk_heat','bulk_input_heated','bulk_input_heat','bulk_input_kit','bulk_input_library','bulk_input_library_prep','bulk_enriched_heated','bulk_enriched_heat','bulk_enriched_kit','bulk_enriched_library','bulk_enriched_library_prep','low_tn5','low_tn5_buffer','low_loaded_adapters','low_tagmented','low_tagment','low_tagment_cleanup_condition','low_tagment_cleanup_dna','low_tagment_cleanup_waste','low_tagment_cleanup','low_split_condition','low_input','low_capture_input','low_split','low_beads','low_enriched','low_unretained','low_capture','low_i5','low_i7','low_input_extended','low_input_initial_extension','low_input_amplified','low_input_pcr','low_input_cleanup_condition','low_input_cleanup_dna','low_input_cleanup_waste','low_input_cleanup','low_enriched_extended','low_enriched_initial_extension','low_enriched_amplified','low_enriched_pcr','low_enriched_cleanup_condition','low_enriched_cleanup_dna','low_enriched_cleanup_waste','low_enriched_cleanup')
 g.declare_reviewed_scope(targets,g.authority)
 g.declare_reviewed_scope(('bulk_label_conditions','bulk_click_conditions','low_label_conditions','low_click_conditions','bulk_click_dna', 'low_click_dna', 'bulk_shear_dna', 'low_tagment_dna', 'bulk_beads_raw', 'bulk_bead_preparation', 'bulk_prepare_beads', 'low_beads_raw', 'low_bead_preparation', 'low_prepare_beads', 'bulk_input_kit_inventory', 'bulk_enriched_kit_inventory', 'low_input_q5_stock', 'low_input_q5_components', 'low_enriched_q5_stock', 'low_enriched_q5_components', 'bulk_input_load', 'bulk_input_sequenced', 'bulk_enriched_load', 'bulk_enriched_sequenced', 'low_input_load', 'low_input_sequenced', 'low_enriched_load', 'low_enriched_sequenced'),'Published MOESM1 protocol sections 1/3: actual selected DNA contributions, prewashed C1 stock, 25 µl Q5 Master Mix, Swift PCR supply classes and actual sequencing loading fractions. Bookkeeping views preserve known same-stock components, not fresh doses; whole kit formulation unassembled.')
 return g.require_publishable()


def render_page():
 from page import head,caveat,table
 from chemdraw import panel
 from protocol_views import assemble_page,WorkflowSection
 g=protocol_graph()
 return assemble_page(g,(head(TITLE),'<main class="wrap" data-protocol-complete="true"><h1>'+TITLE+'</h1><p>'+SUMMARY+'</p><p>'+SOURCE+' · <a href="01_kas-seq.html">Research note and detailed source differences</a></p>',
 caveat('<strong>INFERRED — representative genomic site.</strong> Every guanine-bearing drawing is a possible accessible site, not a known genomic coordinate or complete labeling of all guanines. Reads report the genomic sequence of enriched fragments; they do not directly locate the adduct. Capture is not assumed exclusive to transcription or PolII. Exact bulk-kit and loaded-Tn5/index oligos remain unassembled.'),
 '<h2>Source-qualified chemical labeling</h2><ol><li>Label cells with 5 mM N3-kethoxal for 10 min at 37°C and recover DNA in borate buffer.</li><li>Couple DBCO-PEG4-biotin to the azide handle. RNase A treatment and cleanup remove RNA and unretained material without asserting perfect purity.</li><li>Bulk: select 1 µg biotinylated DNA for sonication to 150–350 bp; reserve 5 µl input and capture the other 95 µl on prepared C1 beads.</li><li>Low input: tagment first with the source-named Tn5 reagent, recover 55 µl, reserve 5 µl and capture 50 µl. The loaded adapter polymers remain unassembled.</li></ol>',
 panel(TaggedDNA('biotin-clicked G adduct').rows(),caption='INFERRED — representative accessible guanine carrying the source-described clicked label; actual genomic position and atom map unassigned'),
 '<h2>Distinct library branches</h2>',table(('Branch','Physical preparation and read setup'),(('Bulk input and enriched DNA','Heat both at 95°C for 10 min; enriched DNA is transferred as supernatant while residual C1 beads remain separate. Use Swift30024: 7–8 PCR cycles for input, 12–14 for enriched DNA; 200–500 bp library. Associated accepted-manuscript readout: NextSeq500, single-end 80 bp.'),('Low-input input and enriched DNA','Carry retained beads directly to the source PCR; do not insert bulk heat elution. Each reaction uses 25 µl UltraIIQ5 master mix with its co-present polymerase/nucleotide/buffer components, plus source i5/i7 primer allocations. Fifteen cycles follow the initial 72°C/95°C steps. Remove beads after PCR; MinElute cleanup gives a 200–600 bp library. Published protocol: single-end sequencing, cycle count unassigned.'))),
 '<p>Input controls and enriched libraries remain separate. Selected DNA transfers and sequencing portions retain explicit stock remainders; no physical pooling is invented.</p>',
 '<h2>Correction and source-version differences</h2><p>The published correction replaces 100 mM with 100 µM for a synthetic DNA oligo in an in-vitro characterization assay. It does not change the live-cell 5 mM labeling condition. The preparation below follows the published Supplementary Protocol; complete equivalence to the article-associated accepted manuscript is not asserted.</p>',
 table(('Item','Published Supplementary Protocol','Associated accepted manuscript'),(('Bulk click input / RNase A','2 µg /15 min','1 µg /5 min'),('Low-input extraction kit','D4074','D4068'),('Low-input Tn5 identifier','15027865; buffer15027866','10527865'),('Low-input input/capture allocation','5 µl of55 µl /50 µl','No separate volume assignment in the inspected paragraph'))),
 '<p>The Swift kit name does not imply bisulfite conversion. Historical finishing recipes, exact loaded adapters, primer/index members and sequencing-primer landing sites remain unassembled. The source’s analysis extension to 150 bp is not a 150-base sequencing read.</p>',
 '<details><summary>Full molecular workflow, physical fractions and unresolved chemistry</summary>',WorkflowSection('INFERRED representative-site views — source-qualified bulk and low-input workflows'),'</details></main>'))


def selftest():
 from checks import Check
 c=Check()
 g=protocol_graph();render_page()
 c('Labelling reagent molecular, not hidden condition',g.materials['n3_kethoxal'].kind,K.MOLECULAR)
 c('Click reagent molecular, not hidden condition',g.materials['dbco'].kind,K.MOLECULAR)
 c('Bulk click takes selected stock contribution',g.steps['bulk_click'].inputs[0].material,'bulk_click_dna_selected')
 c('Low click takes selected42.5ul portion',g.steps['low_click'].inputs[0].material,'low_click_dna_selected')
 c('Bulk sonication takes selected1ug portion',g.steps['bulk_fragment'].inputs[0].material,'bulk_shear_dna_selected')
 c('Low tagmentation leaves25ulstock remainder','low_tagment_dna_remainder'in g.terminals,True)
 c('Bulk residual beads excluded from library','bulk_residual_beads'in g.terminals,True)
 c('Low input does not inherit bulk heat',any(k.startswith('low')and k.endswith('_heat')for k in g.steps),False)
 for prefix in ('low_input','low_enriched'):
  c(prefix+' Q5 single stock provenance',g.steps[prefix+'_q5_components'].inputs[0].material,prefix+'_q5_stock')
  c(prefix+' PCR consumes same-stock NT',any(i.material==prefix+'_q5_nt'for i in g.steps[prefix+'_pcr'].inputs),True)
 for prefix in ('bulk_input','bulk_enriched','low_input','low_enriched'):
  c(prefix+' library stock remainder',prefix+'_library_stock'in g.terminals,True)
 c('Publication source gate',g.publication_problems(),[])
 # This model contains no source-exact oligo polymer claims. Graph construction
 # checks branch connectivity and explicit source coverage; source details are
 # source-qualified declarations, not checks that merely restate prose.
 c.report()

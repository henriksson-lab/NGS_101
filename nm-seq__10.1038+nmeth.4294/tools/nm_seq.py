"""Source-selected Nm-seq chemistry with qualified actual library finishing."""
from dataclasses import dataclass,replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SequencingPlan
from chemdraw import Scene,Segment
from adapter_ligation import AdapterLigationGeometry
SOURCE='Dai et al.2017 DOI10.1038/nmeth.4294; published SupplementaryProtocol1 steps1–19 controls the selected preparation; associated manuscript NIHMS922173 supplies qualified mechanism/sequencing context, not final procedural equivalence.'
PAGE_FILE='nm-seq.html'

@dataclass(frozen=True)
class EndPath:
    rna: str
    methylated: frozenset[int]=frozenset()
    end3: str='OH'
    released: tuple[str,...]=()
    def __post_init__(self):
        if not self.rna or set(self.rna)-set('ACGU'):raise ValueError('nonempty RNA required')
        if any(not 0<=p<len(self.rna) for p in self.methylated):raise ValueError('methylated site out of range')
        if self.end3 not in {'OH','dialdehyde','phosphate'}:raise ValueError('unknown terminal chemistry')
        if self.end3=='dialdehyde' and len(self.rna)-1 in self.methylated:raise ValueError('selected terminal Nm cannot have the vicinal-diol oxidation product')
    def oxidize(self):
        if self.end3!='OH':raise ValueError('vicinal-diol path requires dephosphorylated end')
        return self if len(self.rna)-1 in self.methylated else replace(self,end3='dialdehyde')
    def eliminate(self):
        if self.end3=='OH' and len(self.rna)-1 in self.methylated:return self
        if self.end3!='dialdehyde' or len(self.rna)==1:raise ValueError('selected elimination requires oxidized terminal ribose and surviving upstream RNA')
        return replace(self,rna=self.rna[:-1],end3='phosphate',released=self.released+(self.rna[-1],))
    def dephosphorylate(self):
        if self.end3=='dialdehyde':raise ValueError('eliminate before phosphatase')
        return replace(self,end3='OH')
    def cycle(self):return self.oxidize().eliminate().dephosphorylate()
    @property
    def ligatable3(self):return self.end3=='OH'
    def cleanup(self):return replace(self,released=())

def _clean(g,key,template,label):
    reagent=g.input(key+'_column',label,K.SELECTION)
    kept=g.product(key+'_rna','Source column-recovered RNA; purity/yield unassigned',K.MOLECULAR,chemical_state=template.chemical_state)
    waste=g.product(key+'_waste','Actual column reject/wash fractions',K.WASTE)
    g.react(key,label,S.CLEANUP,((template,R.TEMPLATE),(reagent,R.SELECTION)),(kept,waste))
    g.discard(waste,'Physical column separation, not removal by the preceding chemical reaction. Exact debris partition and recovery unassigned.')
    return kept


def end_classes():
    sc=Scene()
    sc.strand('nm',[Segment('source RNA ending in Nm','[RNA ending in Nm]',placeholder=True)],label='Selected source Nm-ended class',mod5='phosphate',mod3='OH')
    sc.strand('n',[Segment('susceptible RNA after final OE','[shortened RNA]',placeholder=True)],label='Selected susceptible end class',mod5='phosphate',mod3='phosphate')
    sc.note('nm','Source class declarations only; incomplete reaction background, site identity and ligation occupancy remain unassigned.')
    return sc


def protocol_graph():
    g=ProtocolGraph('Nm-seq — selected enriched-RNA preparation',authority=SOURCE,representation='published_instance')
    total=g.input('total_rna','Source10ugDNase-treated totalRNA',K.MOLECULAR,primary=True)
    fragmentation=g.input('fragmentation_supply','Source2uL10xRNAfragmentationreagent in20uL reaction;95C5min then2uLSTOP/ice',K.CONDITION)
    fragmented=g.product('fragmented','Source fragmented totalRNA; size/site distribution unassigned',K.MOLECULAR)
    g.react('fragment','Fragment source totalRNA',S.FRAGMENTATION,((total,R.TEMPLATE),(fragmentation,R.CONDITION)),(fragmented,))
    fragments=_clean(g,'fragment_cleanup',fragmented,'SourceRNA-Clean-and-Concentrator5; elute15uL')
    antarctic=g.input('antarctic','Source2uLAntarcticPhosphatase,10U',K.CATALYST)
    repair_buffer=g.input('repair_buffer','Source5uL10xreactionbuffer;RNA/water43uL,total50uL,37C30min',K.CONDITION)
    repaired=g.product('repaired','RNA3OH end class; internal Nm and conversionfractions unassigned',K.MOLECULAR,chemical_state={'RNA_ends3':('OH',),'scope':'sourceexpected endclass, not completeconversion guarantee'})
    g.react('repair3','Source3prime-end repair before OED',S.DEPHOSPHORYLATION,((fragments,R.TEMPLATE),(antarctic,R.CATALYST),(repair_buffer,R.CONDITION)),(repaired,))
    prior=_clean(g,'repair_cleanup',repaired,'PublishedProtocol1step6column; elute15uL')
    for n in (1,2,3,4,5,6,7,8,9):
        final=n==9
        periodate=g.input(f'periodate_{n}','Actual4uL200mMNaIO4,source20mM final in40uL',K.MOLECULAR)
        lysine=g.input(f'lysine_{n}','Actual4uL2MlysineHCl,pH8.5,source200mM final',K.MOLECULAR)
        condition=g.input(f'oe_program_{n}','RNA/water32uL;total40uL;37C30min withshaking/lightprotection',K.CONDITION)
        oe=g.product(f'oe_{n}','Nm-ended3OH class plus susceptible shortenedRNA3P and released residue products; mixture unquantified',K.MOLECULAR,chemical_state={'RNA_ends3':('OH','phosphate'),'Nm_path':'3OH; selected oxidation-resistant class','susceptible_path':'shortenedRNA3P plus released residue-derived products','yields':'unassigned'})
        g.react(f'oxidation_elimination_{n}','Final OE; no dephosphorylation follows' if final else f'Source OEDcycle{n}: oxidation/elimination',S.CHEMICAL_CONVERSION,((prior,R.TEMPLATE),(periodate,R.STRUCTURAL),(lysine,R.STRUCTURAL),(condition,R.CONDITION)),(oe,))
        glycol=g.input(f'ethylene_glycol_{n}','Actual2uLethylene-glycol quench reagent',K.MOLECULAR)
        quench_program=g.input(f'quench_program_{n}','SourceStep9mix/spin after ethylene-glycol addition; kinetics unassigned',K.CONDITION)
        quenched=g.product(f'quenched_{n}','Quenched source RNA/reagent/debris mixture',K.MOLECULAR,chemical_state=oe.chemical_state)
        g.react(f'quench_{n}','Quench before enzyme or column processing',S.QUENCHING,((oe,R.TEMPLATE),(glycol,R.STRUCTURAL),(quench_program,R.CONDITION)),(quenched,))
        prior=quenched
        if not final:
            sap=g.input(f'rsap_{n}','Source2uLrSAP,2U',K.CATALYST)
            sap_buffer=g.input(f'rsap_buffer_{n}','Source5uL10xCutSmart plus1uLwater;50uLtotal,37C30min',K.CONDITION)
            deph=g.product(f'deph_{n}','SourceexpectedRNA3OH class with copresent released material; completenessunassigned',K.MOLECULAR,chemical_state={'RNA_ends3':('OH',),'scope':'selected reactionclass, conversion/debrispartitionunassigned'})
            g.react(f'dephosphorylate_{n}','Dephosphorylate directly after quench; no intervening cleanup',S.DEPHOSPHORYLATION,((prior,R.TEMPLATE),(sap,R.CATALYST),(sap_buffer,R.CONDITION)),(deph,));prior=deph
        prior=_clean(g,f'cleanup_{n}',prior,'SourceStep15column,elute32uL; final3OH/3P distinction retained' if final else 'SourceStep12RNAcolumn,elute32uL before nextcycle')
    pnk=g.input('pnk','Source1uLT4PNK3primephosphatase-minus,10U',K.CATALYST)
    atp=g.input('atp','Source5uL10mMATP,1mM final',K.MOLECULAR)
    pnk_buffer=g.input('pnk_buffer','Source5uL10xPNKbuffer;RNA/water39uL,total50uL;37C30min',K.CONDITION)
    ready=g.product('prepared','Source5primephosphorylated RNA with final3OH/3P classes retained',K.MOLECULAR,chemical_state={'RNA_ends3':prior.chemical_state['RNA_ends3'],'end5_expected':'phosphate','scope':'selected endclasses, no completeconversion/purityclaim'})
    g.react('phosphorylate5','Source5primephosphorylation, not removal of final3primeP',S.PHOSPHORYLATION,((prior,R.TEMPLATE),(pnk,R.CATALYST),(atp,R.STRUCTURAL),(pnk_buffer,R.CONDITION)),(ready,))
    kit_rna=_clean(g,'post_pnk_cleanup',ready,'ActualpublishedStep18RNAcolumn;elute8uL forlibraryprep')
    kit=g.input('kit_supply','ActualNEBNextSmallRNALibraryPrepSetE7330S reagentcollection; multiple stocks, not asserted to be one premix',K.MOLECULAR)
    adapter=g.product('kit_adapter3','Sourcekit3prime-adapter class; wholepolymer/modifications/geometry unassembled',K.ADAPTER)
    ligase=g.product('kit_ligase3','Source3prime-ligation catalystclass; identity/formulation unassigned',K.CATALYST)
    remainder=g.product('kit_remaining','Remaining actualkit reagentcollection; stockcount/composition/ordering unassembled',K.MOLECULAR)
    g.react('kit_inventory','Expose source-positive ligation supplies within actualkit collection',S.MAGIC,((kit,R.TEMPLATE),),(adapter,ligase,remainder),changes_molecule=False,reason='Component views of one actual supplied kitcollection, not an invented premix or fresh duplicatekit. No synthetic fluorescent model-adapter is borrowed into sequencing prep.')
    ligation_condition=g.input('ligation3_condition','PublishedStep19override:3prime-adapterligation overnight16C',K.CONDITION)
    ligated=g.product('ligation3_mixture','Source3prime-ligated RNA classes plus nonproductive/copresent material; occupancy/yield unassigned',K.MOLECULAR)
    geometry=AdapterLigationGeometry('unavailable','unavailable',(end_classes(),),'Source3OH/3P substrateclasses are declared; fullkitadapter/occupiedjunction geometry unassembled. No unreviewed historicalkit oligo is substituted.')
    g.react('ligation3','Source overnight3prime-adapter ligation',S.LIGATION,((kit_rna,R.TEMPLATE),(adapter,R.STRUCTURAL),(ligase,R.CATALYST),(ligation_condition,R.CONDITION)),(ligated,),geometry=geometry)
    library=g.product('library','DNA-library material constructed by source NEBNextworkflow; wholepolymer and physicalpurity unassigned',K.LIBRARY)
    g.react('kit_finish','Source actualsmallRNA library-construction continuation',S.MAGIC,((ligated,R.TEMPLATE),(remainder,R.STRUCTURAL)),(library,),reason='PublishedStep19 explicitly constructs smallRNA library by manufacturer instructions. Remaining dual-ligation/directional conversion and PCR are reported in the qualified manuscript; precise adapter/RT/PCR/workup details, individualkitstocks and exactendproducts remain an opaque source-qualified interface; this is a DNA-library product, not sequencing the preparedRNA.')
    selected=g.product('selected_library','SelectedDNA-library contribution for source clusterpreparation; volume unassigned',K.LIBRARY)
    stored=g.product('library_remainder','Untransferred library stock; amount unassigned',K.LIBRARY)
    g.react('select_loading','Select actuallibrary-to-run contribution',S.SPLIT,((library,R.TEMPLATE),),(selected,stored),changes_molecule=False)
    g.mark_terminal(Terminal(stored.id,Disposition.UNRESOLVED,'No whole-library consumption, specific samplepool or untreated-input-control preparation is inferred.'))
    run=g.input('run_supply','Associated-manuscript sourceTruSeqRapidPECluster/SBS/MultiplexPrimerkit supplies; exact historical components unassembled',K.MOLECULAR)
    run_condition=g.input('run_condition','Associated-manuscript sourceHiSeq2500,10pMtemplatepersample; readlength/indexmask unassigned',K.CONDITION)
    loading=g.product('sequencing_input','Actual selectedDNA-library cluster/sequencing preparation interface',K.LIBRARY)
    g.react('cluster_prep','Prepare selected library for the source sequencing system',S.MAGIC,((selected,R.TEMPLATE),(run,R.STRUCTURAL),(run_condition,R.CONDITION)),(loading,),changes_molecule=False,reason='Qualified source sequencing context from associated manuscript; no completeAM/VOR equivalence or exactrun-primer geometry claimed.')
    g.sequence(loading,SequencingPlan(acquired_roles=('Read 1','Read 2'),limitation='Associated manuscript reports pairedreads/HiSeq2500 and10pMtemplate. Exact rawlengths, sampleindexes, loadedprimers and wholelibrary geometry unassigned; R2start is sourceanalysis of capturedRNA3ends, not a derived primerlanding here.'))
    control_rna=g.input('control_rna','Associated-manuscript actual untreatedinput RNAfragments; sourceendhandling unassembled, not equated with treatedRNA',K.MOLECULAR,primary=True)
    control_condition=g.input('control_condition','Associatedmanuscript sourcekitlibraryconstruction,3primeadapterligation overnight16C; earlier control endhandling unassigned',K.CONDITION)
    control_lib=g.product('control_library','DNA-library material from actualuntreated-input sourcecontrol; wholearchitecture unknown',K.LIBRARY)
    g.react('control_kit_prep','Construct separately reported untreated-input control library',S.MAGIC,
        ((control_rna,R.TEMPLATE),(adapter,R.STRUCTURAL),(ligase,R.CATALYST),(remainder,R.STRUCTURAL),(control_condition,R.CONDITION)),(control_lib,),
        reason='Associatedmanuscript explicitly constructs untreatedinput-fragment libraries with the kit and overnight3prime-ligation override. Opaqueactualcontrol interface: no invented shared OED/PNK/end-repair or finalprimer geometry. Kitports are aggregate allocationviews, not one reused physical aliquot.')
    control_aliquot=g.product('control_loading_aliquot','Selected untreatedcontrol library-to-run contribution; volumeunassigned',K.LIBRARY)
    control_rest=g.product('control_remainder','Untransferred untreatedcontrol library stock; amountunassigned',K.LIBRARY)
    g.react('control_select_loading','Select actualcontrol library contribution',S.SPLIT,((control_lib,R.TEMPLATE),),(control_aliquot,control_rest),changes_molecule=False)
    g.mark_terminal(Terminal(control_rest.id,Disposition.UNRESOLVED,'No completecontrol endhandling or whole-stock loading is inferred.'))
    control_loading=g.product('control_sequencing_input','Actualcontrol DNA-library clusterpreparation interface',K.LIBRARY)
    g.react('control_cluster_prep','Prepare separatelyreported control library for source sequencing',S.MAGIC,
        ((control_aliquot,R.TEMPLATE),(run,R.STRUCTURAL),(run_condition,R.CONDITION)),(control_loading,),changes_molecule=False,
        reason='Source untreatedcontrol librarycontribution; actualrun/poolgrouping remains unassigned, no invented joint treated/controlphysicalpool.')
    g.sequence(control_loading,SequencingPlan(acquired_roles=('Read 1','Read 2'),limitation='Associated-manuscript untreatedinput control: pairedreads andHiSeq2500 context. Earliercontrolendhandling, exactcycle/indexmask, finalstrand and actualrun/poolgrouping unassembled.'))
    g.declare_reviewed_scope(SOURCE_INVENTORY,SOURCE)
    return g.require_publishable()


# Curator-selected source counts: eight OED rounds plus one final OE, not graph enumeration.
SOURCE_INVENTORY=('control_rna','control_condition','control_kit_prep','control_select_loading','control_remainder','control_cluster_prep','control_sequencing_input','total_rna','fragmentation_supply','fragment','fragment_cleanup_column','fragment_cleanup','antarctic','repair_buffer','repair3','repair_cleanup_column','repair_cleanup',
 'pnk','atp','pnk_buffer','phosphorylate5','post_pnk_cleanup_column','post_pnk_cleanup','kit_supply','kit_inventory','ligation3_condition','ligation3','kit_finish','select_loading','library_remainder','run_supply','run_condition','cluster_prep','sequencing_input')+tuple(prefix+str(n)for n in (1,2,3,4,5,6,7,8,9)for prefix in ('periodate_','lysine_','oe_program_','oxidation_elimination_','ethylene_glycol_','quench_program_','quench_'))+tuple('cleanup_'+str(n)+suffix for n in (1,2,3,4,5,6,7,8,9)for suffix in ('','_column'))+tuple(prefix+str(n)for n in (1,2,3,4,5,6,7,8)for prefix in ('rsap_','rsap_buffer_','dephosphorylate_'))


def render_page():
    from page import head,info
    from chemdraw import panel
    from protocol_views import assemble_page,WorkflowSection
    return assemble_page(protocol_graph(),(head('Nm-seq — terminal enrichment of2prime-O-methylatedRNA'),'<div class="wrap" data-protocol-complete="true"><h1>Nm-seq</h1>',
      '<p class="research-notes"><a href="01_nm_seq.html">Research notes</a></p>',info('<a href="https://doi.org/10.1038/nmeth.4294">Published SupplementaryProtocol1 and associated Nm-seq manuscript</a>'),
      '<h2>Expose protected RNA ends, then build a DNA library</h2>',
      '<p>The selected preparation fragments10µg totalRNA and repairs3′ends. Eight oxidation–elimination–dephosphorylation cycles progressively remove susceptible terminal residues. A final oxidation/elimination omits dephosphorylation, leaving different3′end classes. PNK adds5′phosphate without removing the final3′block, followed by actual column cleanup, overnight3′adapterligation and the source smallRNA-kit library workflow.</p>',
      panel(end_classes().rows(),caption='Source end classes after finalOE/PNK: selected Nm-ended3′OH versus susceptible shortened3′P. This is not a literal genomic sequence or a prediction of perfect enrichment/ligation.'),
      '<h2>Source wording and version boundaries</h2>',
      '<p>The published protocol controls the preparation:20mMperiodate and30minPNK. Its section heading says5′phosphorylation, whileStep16calls the ATP/PNK reaction5′dephosphorylation. The model uses the heading and ATP-supported kinase direction; the printed wording difference is retained. The associated manuscript uses10mMperiodate in a synthetic-model assay and describes60minPNK for Nm-seq. These values are not imported into the selected supplementary route. Different experimental contexts and unreconciled versions are not declared errors in a reconciled final paper.</p>',
      '<h2>Known finishing stages, unassembled kit geometry</h2>',
      '<p>PublishedStep18purifies the PNKreaction and elutes8µL before Step19libraryconstruction with NEBNextE7330S. The3′ligation override is16°Covernight. The sequencing library is a source-positive DNA product; individual historical adapters, later RT/PCR/workups and whole-library geometry remain unassembled. The fluorescent synthetic ligation-control adapter is not substituted for sequencing reagents. Kit inventories are collections, not invented premixes.</p>',
      '<p>The associated manuscript reports HiSeq2500 pairedreads with10pMtemplate. Raw cycle/index masks and exact runprimers remain unassigned. The associated manuscript separately reports untreated input libraries. Their actual kit conversion and selected loading contribution are retained as a distinct opaque control branch; earlier end handling and run grouping remain unassigned, without borrowing the treated RNA conditions. Enrichment and a read-end pileup do not directly establish modification stoichiometry, site identity, specificity or yield.</p>',
      WorkflowSection('Full selected preparation and molecular limits'),'</div>'))

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();p=EndPath('ACGU',frozenset({1}));q=p.cycle().cycle()
    c('two cycles expose internal Nm and conserve released sequence',(q.rna,q.released,q.end3),('AC',('U','G'),'OH'))
    c('terminal Nm halts selected progressive shortening',q.cycle(),q)
    unmodified=EndPath('ACGU').oxidize().eliminate()
    c('final OE without phosphatase leaves unligatable RNA3′p',(unmodified.rna,unmodified.end3,unmodified.ligatable3),('ACG','phosphate',False))
    c('physical cleanup removes released material only',unmodified.cleanup().rna,unmodified.rna)
    try:unmodified.oxidize();bad=False
    except ValueError:bad=True
    c('3′phosphate cannot enter next diol oxidation without repair',bad,True)
    original=Source(sources_root()/'nm-seq__10.1038+nmeth.4294'/'nmeth.4294_41592_2017_BFnmeth4294_MOESM3_ESM.pdf.txt','publisher Supplementary Protocol1')
    if have(c,original,label='original protocol checks'):
        text=' '.join(original.path.read_text().split())
        c('source repeats OED seven additional times','Repeat steps 7-12 seven more times.' in text,True)
        c('source explicitly excludes final dephosphorylation','not perform dephosphorylation (steps 10-11)' in text,True)
        c('actual publishedfinishing explicitly constructs smallRNAlibrary','Use RNA from step 18 to construct a small RNA library using NEBNext' in text,True)
    g=protocol_graph()
    c('eight rounds repair and final ninthOE preserves blocking',sum(k.startswith('dephosphorylate_')for k in g.steps),8)
    c('source post-PNK purification feeds actual ligation',[p.material for p in g.steps['ligation3'].inputs if p.role is R.TEMPLATE],['post_pnk_cleanup_rna'])
    c('sourceATP is separate from PNKcatalyst',any(p.material=='atp' and p.role is R.STRUCTURAL for p in g.steps['phosphorylate5'].inputs),True)
    c('post-PNK cleanup preserves blocked3P RNAclass',g.materials['post_pnk_cleanup_rna'].chemical_state['RNA_ends3'],('OH','phosphate'))
    c('sourcekitfinishing produces actual libraryclass',g.materials['library'].kind,K.LIBRARY)
    c('sequencingdoesnot target preparedRNA',tuple(g.sequencing),('sequencing_input','control_sequencing_input'))
    c('untreatedcontrol doesnot pass through OED',g.steps['control_kit_prep'].inputs[0].material,'control_rna')
    c('untransferred library remains explicit',g.terminals['library_remainder'].disposition,Disposition.UNRESOLVED)
    c.report()
if __name__=='__main__':selftest()

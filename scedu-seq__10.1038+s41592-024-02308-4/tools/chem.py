"""Source-owned private partial scEdU preparation; no complete validation claim."""
TITLE='scEdU-seq selected single-pulse RPE-1 route'
AUTHORITY='Published PMC11239516 main54/150–188; reporting summary and peer-review administrative headers inspected; three-original cap closed. Workbook Table2/3 unreviewed.'
INITIAL='Source RPE-1 cells for single15min10uM EdU pulse; other labeling/recovery and FUCCI enrichment branches separate'
CELL_INPUT=True
READOUT='Published main188 literal: NextSeq500 or NextSeq2000, R1 100cycles (cell index/UMI), R2 100cycles (sample index). Full barcode widths, sequencer-index masks and molecular read landing remain unassembled; unusual source wording retained rather than silently repaired. UMI/cell/sample meanings known but exact historical adapter/primer inventories unassembled.'
DEFERRED='Published Table2/3 workbook originals deferred at3-original cap; historical EdUClick kit/copper formulation, streptavidin workup and full T7 adapter/PCR geometry not assembled. Main176 lists PNK reaction buffer twice while saying both enzymes and main54 names PNK: exact formulation unresolved.'
OLIGOS={
 'Source Klenow reconstitution primer':{'segments':lambda:[Segment('source adapter primer','ATGCCGGTAATACGACTCAC')],'chemistry':{'polymer':'DNA'}},
 'Source random-hexamer RT primer':{'segments':lambda:[Segment('source RT fixed arm','G'+revcomp(TRUSEQ_SMALL_RNA_RA3)),Segment('random priming hexamer','N'*6,placeholder=True,note='Random priming, not a source-assigned UMI')],'chemistry':{'polymer':'DNA'}}}

from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from chemdraw import Scene,Segment,feature,revcomp
from illumina import TRUSEQ_SMALL_RNA_RA3
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry

def scene(label):
    item=OLIGOS.get(label)
    if item is None:return None
    sc=Scene();sc.strand('oligo',item['segments'](),mod5=item.get('mod5',''),mod3=item.get('mod3',''));return sc

def chemistry(label):return OLIGOS.get(label,{}).get('chemistry')
def unknown_scene():
    sc=Scene();sc.strand('unassembled',[Segment('unassembled geometry','[source geometry unassembled]',placeholder=True)]);return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance')
    x=g.input('input',INITIAL,K.BIOLOGICAL if CELL_INPUT else K.MOLECULAR,primary=True)
    def stock(key,label):return g.input(key,label,K.OLIGO,state=scene(label),chemical_state=chemistry(label))
    def op(key,label,kind,reagent,product,oligos=(),symmetric=False):
        nonlocal x
        enzyme=kind in {S.DIGESTION,S.END_REPAIR,S.EXTENSION,S.TAILING,S.LIGATION,S.PCR,S.DEPROTEINIZATION,S.TAGMENTATION,S.TRANSCRIPTION}
        c=g.input(key+'_stock',reagent,K.CATALYST if enzyme else K.CONDITION)
        ins=[(x,R.TEMPLATE),(c,R.CATALYST if enzyme else R.CONDITION)]
        for i,o in enumerate(oligos):ins.append((stock(key+'_oligo_'+str(i),o),R.STRUCTURAL))
        kw={}
        if kind==S.MAGIC:kw['reason']='Source-specific compound transition or molecular geometry remains unassembled in this selected model.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        if kind==S.LIGATION:kw['geometry']=AdapterLigationGeometry('unavailable','unavailable',(unknown_scene(),),'Actual source ligation; full source-specific substrate/junction geometry unassembled, not a source absence claim.')
        if kind==S.PCR:kw['geometry']=PCRGeometry.unavailable(unknown_scene(),reason='Source amplification known; whole source-specific primer landing and product geometry unassembled.',symmetric_species=symmetric)
        y=g.product(key+'_product',product,K.MOLECULAR);g.react(key,label,kind,ins,(y,),**kw);x=y
    def clean(key,label):
        nonlocal x
        c=g.input(key+'_condition',label,K.CONDITION);y=g.product(key+'_kept','Physically recovered source preparation; exact recovery/carryover unassigned',K.MOLECULAR);w=g.product(key+'_waste','Physically separated source workup fractions',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(c,R.CONDITION)),(y,w));g.discard(w,'Source physical workup, not perfect molecular exclusion.');x=y
    def split(key,label,rest_label,disposition=Disposition.UNRESOLVED):
        nonlocal x
        y=g.product(key+'_used',label,K.MOLECULAR);r=g.product(key+'_rest',rest_label,K.QC if disposition==Disposition.QC else K.MOLECULAR);g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,r));g.mark_terminal(Terminal(r.id,disposition,rest_label));x=y
    def qc(key,label):split(key,'Production preparation retained after separate analytical aliquot',label,Disposition.QC)
    def external_anneal(key,primer,label,enzyme,programme,product,rt=False):
        nonlocal x
        p=stock(key+'_primer',primer);c=g.input(key+'_anneal_condition',label,K.CONDITION)
        t=g.product(key+'_template','Same-tube template component after annealing; external primer tracked separately, not an isolated fraction',K.MOLECULAR)
        retained=g.product(key+'_retained_primer','Same-tube co-retained external primer component, not a second stock addition or intramolecular geometry',K.OLIGO,state=scene(primer),chemical_state=chemistry(primer))
        g.react(key+'_anneal',label,S.THERMAL_TREATMENT,((x,R.TEMPLATE),(p,R.STRUCTURAL),(c,R.CONDITION)),(t,retained))
        e=g.input(key+'_enzyme',enzyme,K.CATALYST);y=g.product(key+'_product',product,K.MOLECULAR)
        g.react(key,programme,S.REVERSE_TRANSCRIPTION if rt else S.EXTENSION,((t,R.TEMPLATE),(retained,R.STRUCTURAL),(e,R.CATALYST)),(y,));x=y

    op('edu','Source10uM EdU15min',S.MAGIC,'Source in-cell nucleotide incorporation; replication chemistry unassembled','Source pulse-labeled cells; label occupancy and repair incorporation unassigned')
    clean('pulse_wash','Source trypsinization/PBS wash')
    op('fix','Source fixation75percent ice-cold ethanol; selected short storage -20C<=24h',S.GENERIC,'Source fixation/storage','Source fixed EdU-labeled cells; longer storage buffer branch separate')
    clean('block','Source BSA-precoated ProteinLoBind tubes/pellet600g3min;0.25percentBSA/PBS4C30min block/pellet')
    azide=g.input('azide_biotin','Source2mM azidePEG3biotin covalent labeling reagent replacingazide647',K.MOLECULAR)
    copper=g.input('click_catalyst','Source copper(I)-catalyzed EdUClick647 kit with6mMTHPTA; exact historical copper formulation unreviewed',K.CATALYST)
    click=g.product('click_product','Source covalently biotin-linked EdU-containing DNA in fixed cells; conversion yield unassigned',K.MOLECULAR)
    g.react('click','Source50ul in-situ EdUClick biotin substitution',S.CHEMICAL_CONVERSION,((x,R.TEMPLATE),(azide,R.STRUCTURAL),(copper,R.CATALYST)),(click,));x=click
    clean('click_wash','Source PBSwash/resuspend0.25percentBSA/PBS10ugmlDAPI/20um mesh')
    op('sort','Source index-sort single cells into384wells each5ul mineral oil',S.MAGIC,'Source BD FACS Influx gating excludes debris/doublets; physical partition occupancy unassembled','Source selected single-cell well classes with unselected cell fractions unassigned')
    split('chosen_cell','Chosen source single-cell well','Other sorted wells/uncollected cells remain outside selected route; not implicitly discarded')
    op('lysis','Source100nllysis55C2h/80C20min',S.DEPROTEINIZATION,'Source10nlproteinaseK/10nlCutSmart/80nlwater','Source cell DNA/RNA and proteolysis products; no perfect nuclear purification')
    op('digest','Source100nlNlaIII mix37C4h/65C30min',S.DIGESTION,'Source NlaIII/declaredCutSmart formulation','Source NlaIII-cut DNA and residual substrates; exact end classes unassembled')
    op('endprep','Source100nl endrepair37C30min/75C20min; apparent PNK formulation discrepancy retained',S.END_REPAIR,'Source Klenow large fragment; PNK named inmain54 butmain176 lists buffer twice, exactenzymeaddition unresolved; dNTP/ATP/Mg/PEG/BSA','Source end-repaired fragment preparation; exact PNK composition/end distributions unassembled')
    op('atail','Source100nl Atailing72C15min',S.TAILING,'Source AmpliTaq360/dATP/KCl/PEG/BSA','Source A-tailed fragment classes plus unprocessed/offproducts')
    op('ligate','Source50nl5mM T7adapters+150nl ligasemix;4C20min16C16h65C20min',S.LIGATION,'Source T4 DNA ligase/Mg/Tris/DTT/ATP/PEG/BSA','Source T7/cellbarcode/UMI-bearing ligation classes; all physical molecules retained beforeworkup',('Source T7promoter T-overhang adapter with cellbarcode/UMI; Table2unreviewed',))
    op('pool_cells','Source plate contents centrifuged300g1min to oil-coated reservoirs; aggregate cell pool',S.GENERIC,'Physical pooling of indexed wells; individual members unassembled','Source pooled aqueous/oil preparation')
    clean('oil_remove','Source separate aqueous phase from residual mineral oil by centrifugation')
    c=g.input('streptavidin','Source MyOne StreptavidinC1 beads; historicalworkup unreviewed',K.SELECTION)
    bound=g.product('edu_bound','Physical source bead-associated EdU-biotin DNA preparation; nonspecific carryover unassigned',K.MOLECULAR);unbound=g.product('edu_unbound','Source noncaptured DNA/reagents; exact fate unassigned in selected main',K.MOLECULAR)
    g.react('edu_capture','Source affinity purification of EdU-biotin DNA',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(c,R.SELECTION)),(bound,unbound));g.mark_terminal(Terminal(unbound.id,Disposition.UNRESOLVED,'No ideal perfect biotin capture or unreported discard asserted.'))
    heat=g.input('release_heat','Source95C denaturation',K.CONDITION);free=g.product('complement_strand','Source retrieved non-EdU complementary DNA strand preparation; not the biotinylated strand',K.MOLECULAR);beads=g.product('released_beads','Co-separated bead-associated EdU-biotin strands and residualDNA/reagents; exact laterfate unassigned',K.MOLECULAR)
    g.react('complement_release','Source heat release complementary strand95C; physical fraction boundary explicit',S.STRAND_SEPARATION,((bound,R.TEMPLATE),(heat,R.CONDITION)),(free,beads));g.mark_terminal(Terminal(beads.id,Disposition.UNRESOLVED,'Bound EdU strands are not silently sequenced.'));x=free
    external_anneal('reconstitute','Source Klenow reconstitution primer','Source ramp95Cto20C0.1Csec;20mMTrispH8/1mMMg/100mMNaCl','Source0.5UKlenowlargefragment/NEBuffer2/50mMdNTP literal','Source25C45min75C20min','Source rebuilt complementary DNA duplex classes; exact ends/primer geometry unassembled')
    clean('dna_clean','Source1:1 sample:AMPure/elute7ulwater')
    op('ivt','Source MEGAscriptT7 IVT37C12h',S.TRANSCRIPTION,'Source MEGAscriptT7 kit','Source amplified RNA plus templateDNA/reagents; no exact yield assumption')
    op('dnase','Source2ulTurboDNase37C15min',S.DIGESTION,'Source TurboDNase','Source RNA and DNA digestion products; no physical DNA disappearance')
    op('rna_fragment','Source aRNA94C2min inTrisacetate/KOAc/MgOAc;icecool',S.FRAGMENTATION,'Source fragmentation buffer','Source fragmented RNA with DNA/reagent remnants')
    op('fragment_stop','Source50mMEDTA',S.QUENCHING,'EDTA','Source quenched RNA mixture')
    clean('rna_clean','Source1:1RNAcleanXP/elute12ulwater')
    split('rt_aliquot','Source5ulof12ulRNA for RT','Untransferred RNA remainder; not implicitly discarded')
    external_anneal('rt','Source random-hexamer RT primer','Source0.5ul10mMdNTP/1ul20uMprimer65C5min/ice','Source0.5ulSuperscriptII/2ulFirstStrand/1ul0.1MDTT/0.5ulRNaseOUT','Source25C10min42C60min70C20min','Source cDNA/RNA hybrids plus unextended RNA/primer/reagents',rt=True)
    op('rnase','Source0.5ulRNaseA37C30min',S.DIGESTION,'RNaseA','Source ss-cDNA and RNA fragments; digestion distinct from physical removal')
    op('pcr','Source NEBNextUltraIIQ5 PCR adds smallRNAhandles/samplebarcodes; actual programme unassembled',S.PCR,'Source25ulMasterMix/11ulwater;listed2ulRP1/RPIx10uM allocationunassigned','Source amplified Illumina library classes',('SourceRP1 actual Table3unreviewed','SourceRPIx actual Table3unreviewed'))
    clean('final_clean1','SourceAMPure sample:beads0.8:1 first cleanup')
    clean('final_clean2','SourceAMPure sample:beads0.8:1 second cleanup/elute10ulwater')
    qc('final_qc','Separate source Qubit/Bioanalyzer analytical material')
    library=g.product('library','Source sequencing preparation; full source-specific library geometry unassembled',K.LIBRARY)
    g.react('loading','Physical source loading preparation',S.GENERIC,((x,R.TEMPLATE),),(library,),changes_molecule=False)
    g.sequence(library,SequencingPlan(limitation=READOUT))
    g.add_source_item(SourceItem('pending',DEFERRED,'Selected source scope only; incomplete linked-source/whole molecular review, no accepted absence claim.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Click precedes single-cell sorting',g.steps['click_wash'].inputs[0].material,'click_product')
    c('Click reagent structurally incorporated',[(p.material,p.role)for p in g.steps['click'].inputs if p.material=='azide_biotin'],[('azide_biotin',R.STRUCTURAL)])
    c('Complement released from captured material',g.steps['complement_release'].inputs[0].material,'edu_bound')
    c('Bead-bound EdU strand not silently sequenced',g.terminals['released_beads'].disposition,Disposition.UNRESOLVED)
    c('Reconstituted DNA precedes T7 IVT',g.steps['ivt'].inputs[0].material,'dna_clean_kept')
    c('DNase precedes RNA fragmentation',g.steps['rna_fragment'].inputs[0].material,'dnase_product')
    c('Only allocated RNA enters RT',g.steps['rt_anneal'].inputs[0].material,'rt_aliquot_used')
    c('RNA remainder explicit',g.terminals['rt_aliquot_rest'].disposition,Disposition.UNRESOLVED)
    for key in ['reconstitute','rt']:
        c(key+' uses produced external primer as STRUCTURAL',[(p.material,p.role)for p in g.steps[key].inputs if p.role in {R.STRUCTURAL,R.SELF_PRIMER}],[(key+'_retained_primer',R.STRUCTURAL)])
        c(key+' stock added once',sum(p.material==key+'_primer'for st in g.steps.values()for p in st.inputs),1)
    c('Reconstitution source primer transcription',''.join(s.top for s in scene('Source Klenow reconstitution primer').strands['oligo'].segs),'ATGCCGGTAATACGACTCAC')
    c('RT full primer pattern transcription',''.join(s.top for s in scene('Source random-hexamer RT primer').strands['oligo'].segs),'GCCTTGGCACCCGAGAATTCCANNNNNN')
    c('Private partial graph remains incomplete',bool(g.publication_problems()),True)
    c.report()
if __name__=='__main__':selftest()

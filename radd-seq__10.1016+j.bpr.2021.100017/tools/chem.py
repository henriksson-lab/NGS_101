"""Protocol-owned RADD oxidation damage preparation; private partial model."""
TITLE='RADD-seq selected oxidation-damage route'
AUTHORITY='Published PMC8651515 main76–108/130–141; DocumentS1 FigS1–S5captions only;2originals. HistoricalTruSeqDNA-Nanokit and exact libraryends unreviewed.'
INITIAL='Source2.2ug extracted U2OSgenomicDNA after50mMKBrO3 1h andPBSwash; damage sample selected, independent1hrepairedsample and negativecontrols separate'
CELL_INPUT=False
READOUT='Source HiSeq2500 50base single-end reads; no index/primer/wholeDNAjunction geometry assigned from unreviewed historical TruSeqDNANano kit. Biotin-dUTP repairincorporation is chemicalmarking, not a source UMI/sampleindex or absolute lesioncount from reads. Computationalbackgroundmasking/repairratios are not physical DNA removal.'
DEFERRED='Historical TruSeqDNA-Nano/GenElute/QIAquick kit originals and exact source-specific DNArepair end/adapter/read states unassembled. Main30ulrepair→120ulhydratedsonicationtransfer→500ulIP-input volume/mass bridge unassigned, not a fabricated inaccessible disclosure document. No confirmed missing claim.'
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance')
    x=g.input('input',INITIAL,K.BIOLOGICAL if CELL_INPUT else K.MOLECULAR,primary=True)
    def op(key,label,kind,reagent,product):
        nonlocal x
        enzyme=kind in {S.DIGESTION,S.END_REPAIR,S.EXTENSION,S.TAILING,S.LIGATION,S.PCR,S.DEPROTEINIZATION,S.TAGMENTATION,S.TRANSCRIPTION}
        c=g.input(key+'_stock',reagent,K.CATALYST if enzyme else K.CONDITION)
        ins=[(x,R.TEMPLATE),(c,R.CATALYST if enzyme else R.CONDITION)]
        kw={}
        if kind==S.MAGIC:kw['reason']='Source-specific compound transition or molecular geometry remains unassembled in this selected model.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        y=g.product(key+'_product',product,K.MOLECULAR);g.react(key,label,kind,ins,(y,),**kw);x=y
    def clean(key,label):
        nonlocal x
        c=g.input(key+'_condition',label,K.CONDITION);y=g.product(key+'_kept','Physically recovered source preparation; exact recovery/carryover unassigned',K.MOLECULAR);w=g.product(key+'_waste','Physically separated source workup fractions',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(c,R.CONDITION)),(y,w));g.discard(w,'Source physical workup, not perfect molecular exclusion.');x=y
    def split(key,label,rest_label,disposition=Disposition.UNRESOLVED):
        nonlocal x
        y=g.product(key+'_used',label,K.MOLECULAR);r=g.product(key+'_rest',rest_label,K.QC if disposition==Disposition.QC else K.MOLECULAR);g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,r));g.mark_terminal(Terminal(r.id,disposition,rest_label));x=y
    op('ogg','Source20ul:2.2ugDNA/2ul10xNEBuffer4/2ulBSA1mgml/1ulhOGG10.5mgml,37C30min',S.DIGESTION,'Source hOGG1','Source excision-generated AP-sites, released lesion bases and incompletely processedDNA; no instantaneous lesionremoval fromphysicaltube')
    op('endoiv','Source add1ulEndoIV10000Uml,37C30min',S.DIGESTION,'Source EndonucleaseIV','Source AP-incised sites/source-described gaps plus uncleavedAP/lesionDNA; exact sugar-remnant/end chemistry unassembled')
    bs=g.input('bst','Source0.5ulBstDNApolymeraseLargeFragment8000Uml',K.CATALYST);ta=g.input('taq','Source0.5ulTaqDNApolymerase5000Uml',K.CATALYST);li=g.input('taq_ligase','Source0.4ulTaqDNA Ligase40000Uml',K.CATALYST)
    nts=g.input('repair_nts','Source0.3ul1mMdNTPmixture A/G/C plusbiotinylateddUTP (no source dTTP substitution)',K.MOLECULAR,chemical_state={'nucleotides':('dATP','dGTP','dCTP','biotin-dUTP'),'biotinylated_base':'U'})
    rc=g.input('repair_condition','Source total30ul,add0.5ul10xbuffer4/0.5ul1mgmlBSA/1.2ul50mMNAD+/water;65C30min',K.CONDITION)
    marked=g.product('repair_product','Source Bst-displacement/Taq-exonucleolysis-extension/Taq-ligation mixture with biotin-dU-incorporated repairedDNA, unprocessedlesion/AP/nickedDNA and free reagents; no exact completedrepair guarantee',K.MOLECULAR)
    g.react('repair','Source simultaneous Bst/Taq/Taqligase thirdreaction, not invented serialenzymeadditions',S.EXTENSION,((x,R.TEMPLATE),(bs,R.CATALYST),(ta,R.CATALYST),(li,R.CATALYST),(nts,R.STRUCTURAL),(rc,R.CONDITION)),(marked,));x=marked
    op('hydrate','Source hydratedDNA preparation before120ultransfer;30→120ulvolume/mass bridge unassigned',S.GENERIC,'Source aqueous preparation; hydrationamount/wholemass relationunassigned','Source hydrated DNA preparation; no added DNA molecules assumed')
    split('sonication_aliquot','Source120ulhydratedDNA transferredtoCovarismicrotube','Untransferred preparation mass/fateunassigned')
    op('sonicate','SourceCovarisS22010min175W20percentduty200cycles/burst,prechilled4C;target~150bp',S.FRAGMENTATION,'Source sonication','Source fragmented biotin-marked/unmarked DNA and free nucleotides/reagents; nominal150bp not exactgeometry')
    clean('post_shear_clean','SourceQIAquickPCRPurification column removes excessnucleotides aftersonication')
    op('ip_volume_bridge','Source DNApreparation for500ulIP; columneluate-to500ulvolume/mass bridge unassigned',S.GENERIC,'Source aqueousIPpreparation; exactmass/hydrationunassigned','Source500ulscaleIP-input preparation')
    split('ip_aliquot','Source500ulfragmentedDNA selectedforpreclear','Any untransferred DNA stock; amount/fateunassigned')
    # Independent preparation of source70ul ProteinG, with real50/20ul physical subdivision.
    pg=g.input('protein_g_stock','Source70ulProteinGbeads stock',K.MOLECULAR);pw=g.input('protein_g_wash_condition','Source twice70ulIPbuffer10mMTrispH7/1mMEDTA/150mMNaCl;resuspend70ul',K.CONDITION)
    washed=g.product('washed_protein_g','Source washed70ulProteinG stock',K.MOLECULAR);washw=g.product('protein_g_waste','Separated source beadwash fractions',K.WASTE)
    g.react('protein_g_wash','Source independent ProteinG washing',S.CLEANUP,((pg,R.TEMPLATE),(pw,R.CONDITION)),(washed,washw));g.discard(washw,'Source physical beadwash fractions.')
    pg50=g.product('protein_g_50','Source50ulProteinGportion foranti-biotinloading',K.MOLECULAR);pg20=g.product('protein_g_20','Source20ulProteinGportion forDNApreclear',K.SELECTION)
    g.react('protein_g_split','Physical source70ul→50ulantibodytube1+20ulprecleartube2',S.SPLIT,((washed,R.TEMPLATE),),(pg50,pg20))
    ab=g.input('antibiotin','Source5ug anti-biotinantibody from1mgmlstock',K.SELECTION);ac=g.input('antibody_condition','Source tube1incubation2h4Crotation/vibration',K.CONDITION);abmix=g.product('antibody_loading_mix','Source antibody/ProteinG binding mixture; includesunboundantibody',K.MOLECULAR)
    g.react('antibody_load','Source independent anti-biotin loadingon50ulProteinG',S.BINDING,((pg50,R.TEMPLATE),(ab,R.SELECTION),(ac,R.CONDITION)),(abmix,))
    mag=g.input('antibody_magnet_condition','Source tube1magneticpulldown',K.CONDITION);loaded=g.product('antibody_beads','Source prepared anti-biotinProteinGbeads fromtube1',K.SELECTION);super1=g.product('antibody_supernatant','Source tube1supernatant afterloading',K.WASTE)
    g.react('antibody_magnet','Source removeantibodyloading supernatant',S.CLEANUP,((abmix,R.TEMPLATE),(mag,R.CONDITION)),(loaded,super1));g.discard(super1,'Source explicitlydiscardstube1supernatant beforeDNAtransfer.')
    pc=g.input('preclear_condition','Source tube2 DNA+20ulProteinG2h4Crotation/vibration;magnet',K.CONDITION);free=g.product('precleared_dna','Source transferredtube2 supernatant DNA',K.MOLECULAR);prebound=g.product('preclear_beads','Source co-separated20ulProteinG/nonspecificboundDNA; laterfateunassigned',K.MOLECULAR)
    g.react('preclear','Source DNApreclear independentof antibodyloading',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(pg20,R.SELECTION),(pc,R.CONDITION)),(free,prebound));g.mark_terminal(Terminal(prebound.id,Disposition.UNRESOLVED,'Nonspecificboundfraction is physically separated, not computationally removed.'))
    bc=g.input('capture_condition','Source transferDNA-supernatant2toantibodybeads1;overnight4Crotation/vibration',K.CONDITION);bound=g.product('ip_bound','Source anti-biotin-capturedDNA/beads with imperfectrecovery/specificity',K.MOLECULAR);unbound=g.product('ip_unbound','Source noncapturedDNA fraction afterIP; laterfateunassigned',K.MOLECULAR)
    g.react('capture','Source captureonpreparedanti-biotinProteinG, no fresh antibody/bead addition',S.AFFINITY_CAPTURE,((free,R.TEMPLATE),(loaded,R.SELECTION),(bc,R.CONDITION)),(bound,unbound));g.mark_terminal(Terminal(unbound.id,Disposition.UNRESOLVED,'Source workup transfersretainedbeads; exact noncaptured fractionfate unassigned.'));x=bound
    clean('ip_wash','Source seven700ulIPbufferwashes5min4Crotation/vibration; capturedbeadsretained')
    op('elute_digest','Source40ulelutionbuffer10mMTE/1mgmlProteinaseK/0.5percentSDS50C3h900rpm',S.DEPROTEINIZATION,'SourceProteinaseK','Source releasedDNA/proteinfragments/SDS andbeads in same tube beforemagneticseparation')
    ec=g.input('elution_magnet','Source postelutionmagneticpulldown',K.CONDITION);eluate=g.product('eluate','Source recoveredDNA-containing eluate',K.MOLECULAR);eb=g.product('elution_beads','Source co-separatedspentbeads/residualDNA; laterfateunassigned',K.MOLECULAR)
    g.react('elution_separate','Source separateDNAeluatefrombeads',S.CLEANUP,((x,R.TEMPLATE),(ec,R.CONDITION)),(eluate,eb));g.mark_terminal(Terminal(eb.id,Disposition.UNRESOLVED,'No assumed completeDNArelease or beadrecycling.'));x=eluate
    clean('ip_clean','SourceQIAquickPCRPurification column afterIPelution')
    op('library_kit','SourceTruSeqDNA-Nano prep WITHOUTadditionalfragmentation orsize selection',S.MAGIC,'HistoricalTruSeqDNA-Nano kit oligos/stages/ends/indexes unassembled','Source amplified Illumina library; no new fragmentation/size-selection imposed from generic kit knowledge')
    library=g.product('library','Source sequencing preparation; full source-specific library geometry unassembled',K.LIBRARY)
    g.react('loading','Physical source loading preparation',S.GENERIC,((x,R.TEMPLATE),),(library,),changes_molecule=False)
    g.sequence(library,SequencingPlan(limitation=READOUT))
    g.add_source_item(SourceItem('pending',DEFERRED,'Selected source scope only; incomplete linked-source/whole molecular review, no accepted absence claim.'))
    return g

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('EndoIV follows actual hOGG1stage',g.steps['endoiv'].inputs[0].material,'ogg_product')
    c('Repair consumes incised preparation',g.steps['repair'].inputs[0].material,'endoiv_product')
    c('Bst/Taq/Taqligase added in one reaction',[p.material for p in g.steps['repair'].inputs if p.role==R.CATALYST],['bst','taq','taq_ligase'])
    c('Source nucleotide pool retains biotin-dU',g.materials['repair_nts'].chemical_state['nucleotides'],('dATP','dGTP','dCTP','biotin-dUTP'))
    c('Physical free-nucleotide cleanup after sonication',g.steps['post_shear_clean'].inputs[0].material,'sonicate_product')
    c('ProteinG50/20 from same washedstock',g.steps['protein_g_split'].inputs[0].material,'washed_protein_g')
    c('Preclear uses only20ulportion',[p.material for p in g.steps['preclear'].inputs if p.role==R.SELECTION],['protein_g_20'])
    c('Antibody50ulportion prepared independently',g.steps['antibody_load'].inputs[0].material,'protein_g_50')
    c('Antibody supernatant explicitly discarded',g.terminals['antibody_supernatant'].disposition,Disposition.DISCARD)
    c('IP receives transferred preclear supernatant',g.steps['capture'].inputs[0].material,'precleared_dna')
    c('IP consumes produced antibodybeads',[p.material for p in g.steps['capture'].inputs if p.role==R.SELECTION],['antibody_beads'])
    c('No fresh antibody added during DNAcapture',sum(p.material=='antibiotin'for st in g.steps.values()for p in st.inputs),1)
    c('Nonspecific bead-bound DNA fate explicit',g.terminals['preclear_beads'].disposition,Disposition.UNRESOLVED)
    c('Proteinase precedes physicalelution',g.steps['elution_separate'].inputs[0].material,'elute_digest_product')
    c('TruSeq receives physically purified eluate',g.steps['library_kit'].inputs[0].material,'ip_clean_kept')
    c('No second fragmentation or size-selection transition',sum(st.kind==S.FRAGMENTATION for st in g.steps.values()),1)
    c('Private partial graph remains incomplete',bool(g.publication_problems()),True)
    c.report()
if __name__=='__main__':selftest()

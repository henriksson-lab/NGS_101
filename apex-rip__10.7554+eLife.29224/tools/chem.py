"""Protocol-owned source-qualified partial population graph, no finished marker."""
from dataclasses import dataclass,replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition

def clean(g,key,x,procedure,retained='source-recovered preparation; exact molecular partition/carryover unquantified'):
    condition=g.input(key+'_condition',procedure,K.CONDITION)
    out=g.product(key+'_kept',retained,K.MOLECULAR)
    waste=g.product(key+'_separated','Source physically separated fractions; no perfect molecular recovery/removal claimed',K.WASTE)
    g.react(key,procedure,S.CLEANUP,((x,R.TEMPLATE),(condition,R.CONDITION)),(out,waste));g.discard(waste,'Actual source physical separation only.');return out

def op(g,key,x,label,kind,*,reagent=None,catalyst=False,product=None,reason='',changes_molecule=True):
    inputs=[(x,R.TEMPLATE)]
    if reagent:
        stock=g.input(key+'_stock',reagent,K.CATALYST if catalyst else K.CONDITION);inputs.append((stock,R.CATALYST if catalyst else R.CONDITION))
    out=g.product(key+'_product',product or label+'; reaction products/remnants/carryover remain unassigned',K.MOLECULAR)
    g.react(key,label,kind,inputs,(out,),reason=reason,changes_molecule=changes_molecule);return out
SOURCE='Kaewsapsak et al.2017 DOI10.7554/eLife.29224; publisher v1VOR MethodsPartsI–III516–540; SuppFile5construct/qPCRinventory. Selected optimizedBP-firstquench-crosslink route.'
@dataclass(frozen=True)
class CrosslinkedRNP:
    rna_id:str
    protein_ids:frozenset[str]
    biotinylated_proteins:frozenset[str]
    crosslinked_proteins:frozenset[str]
    def __post_init__(self):
        if not self.rna_id or not self.protein_ids or not self.biotinylated_proteins<=self.protein_ids or not self.crosslinked_proteins<=self.protein_ids:raise ValueError('explicit RNA/protein identities and source-qualified selected relationships required')
    def selected_capture_path(self):
        return bool(self.biotinylated_proteins&self.crosslinked_proteins)
    def proteolyzed(self):
        return {'rna_id':self.rna_id,'direct_rna_biotinylation':False,'possible_peptide_remnants':True,'original_protein_ids':self.protein_ids}

def protocol_graph():
    g=ProtocolGraph('APEX-RIP — optimized published route',authority=SOURCE,representation='published_instance')
    current=g.input('cells','Source stable targetedAPEX2expressingHEK293T cells; alternatives/compartments not physically pooled',K.MOLECULAR,primary=True)
    current=op(g,'bp',current,'BP500µM30min37°C',S.MAGIC,reagent='source BP/media',reason='Source exposure, distribution unassembled.')
    current=op(g,'protein_label',current,'H2O2final1mM1min initiates source local proteinbiotinylation',S.CHEMICAL_CONVERSION,reagent='source H2O2 plus already expressedAPEX2',product='Selected biotinylated proteins plus RNA/other cell material; not a declaration that RNA is directly biotinylated')
    current=op(g,'radical_quench',current,'Source azide-free Trolox/ascorbatequench1min',S.QUENCHING,reagent='ascorbate10mM/Trolox5mM')
    current=op(g,'crosslink',current,'0.1%formaldehyde crosslinkquench1+9min',S.CROSS_LINKING,reagent='source crosslinkquenchsolution',product='Source formaldehyde RNA–protein complexes, labeled/unlabeled proteins and background; no exclusively direct-binding or complete crosslink guarantee')
    current=op(g,'glycine',current,'Glycinefinal125mM5min',S.QUENCHING,reagent='source glycine')
    current=clean(g,'harvest',current,'Source2PBSwashes/scrape/pellet')
    current=op(g,'lyse',current,'RIPAlysis5min4°C and source sonication3x30s',S.MAGIC,reagent='source RIPA/quenchers',reason='Source extraction/sonication envelope, actual break distributions unassembled.')
    current=clean(g,'clarify',current,'Source15,000g5min4°C clarifiedlysate transfer')
    current=op(g,'dilute',current,'Source1mlNLBdilution',S.GENERIC,reagent='source NLB/quenchers',reason='Buffer handling, no new covalent transformation assigned.',changes_molecule=False)
    control=g.product('input_control','Source20%inputaliquot; separately processed control route unassembled',K.MOLECULAR);rest=g.product('ip_input','Remaining sourcepreparation',K.MOLECULAR)
    g.react('split','Allocate20%input separately',S.SPLIT,((current,R.TEMPLATE),),(control,rest));g.mark_terminal(Terminal(control.id,Disposition.UNRESOLVED,'Independent source inputcontrol downstream preparation not fully assembled; not pooled by normalization.'))
    bead=g.input('beads','Source50µlstreptavidinbeadslurry2h4°C;8specifiedwashes',K.SELECTION)
    bound=g.product('bound','Captured protein-biotin/RNAcomplexes plus background',K.MOLECULAR);waste=g.product('ip_waste','Actual source captureflowthrough/wash fractions',K.WASTE)
    g.react('capture','Source streptavidincapture/RIPA-highsalt-urea washes',S.AFFINITY_CAPTURE,((rest,R.TEMPLATE),(bead,R.SELECTION)),(bound,waste));g.discard(waste,'Actual source beadwashing fractions; no perfect specificity.')
    current=op(g,'proteolyze',bound,'ProteinaseK200µg42°C1h/55°C1h in sarcosidebuffer',S.DEPROTEINIZATION,reagent='source ProteinaseK',catalyst=True,product='Released RNA-containing mixture with possible crosslinked peptide remnants/protein fragments/beads; no direct transfer of proteinbiotin to RNA asserted')
    current=clean(g,'rna_clean1',current,'Source RNACleanXP purification')
    current=op(g,'dnase',current,'RQ1RNase-freeDNaseI5U37°C30min',S.DIGESTION,reagent='source DNaseI',catalyst=True,product='RNApreparation with DNA-digestion fragments retained until cleanup')
    current=clean(g,'rna_clean2',current,'Source RNACleanXP cleanup after DNase')
    current=op(g,'ribo_mix',current,'RiboZeroGoldremoval solution68°C10min/25°C5min',S.MAGIC,reagent='source mass-dependent1or2µlremovalsolution/20µlreaction',reason='Explicit incubation known, internal molecular recognition unassembled.')
    bead=g.input('ribo_beads','Source magneticbeadbindingRT5min/50°C5min',K.SELECTION)
    sup=g.product('ribodepleted','Source RNA-containingSUPERNATANT; exact recovery/carryover unquantified',K.MOLECULAR);pellet=g.product('ribo_pellet','Source companion beadfraction; exact RNA/remnant composition unresolved',K.MOLECULAR)
    g.react('ribo_partition','Source recoverSUPERNATANT after magneticbead incubation',S.SPLIT,((current,R.TEMPLATE),(bead,R.SELECTION)),(sup,pellet));g.mark_terminal(Terminal(pellet.id,Disposition.UNRESOLVED,'Source explicit supernatant recovery, companion beadfraction fate not asserted.'))
    current=clean(g,'ribo_cleanup',sup,'Source RNACleanXP, elute19.5µlTruSeqv2ElutePrimeFragmentmix')
    kit=g.input('truseq','Source TruSeqRNAsampleprepv2 fromIncubateRFPstep; unique index assigned, exact member/oligos/steps unreviewed',K.CONDITION)
    lib=g.product('library','Source uniquely indexedRNA-derived library; exact strand/end/read geometry unassembled',K.LIBRARY)
    g.react('library_prep','Source historicalTruSeqv2kit handoff',S.MAGIC,((current,R.TEMPLATE),(kit,R.CONDITION)),(lib,),reason='Namedkit/version/startpoint known, no internal chemistry or modern substitute inferred.')
    others=g.input('other_libraries','Other source separately prepared indexedlibraries; equimolar≤10perpool',K.LIBRARY)
    pool=g.product('pool','Source equimolar pooledlibrary preparation',K.LIBRARY);g.react('pool_libraries','Source librarypooling',S.POOL,((lib,R.TEMPLATE),(others,R.TEMPLATE)),(pool,))
    g.sequence(pool,SequencingPlan(acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':50,'Read 2':50},limitation='Source HiSeq2500rapidPE50thenindex; exact indexcycles/architecture and full primergeometry unassembled.'))
    g.add_source_item(SourceItem('pending_library','Historical TruSeqv2/RiboZerokit originals and full polymer/read assembly','Published Methods536/540'))
    return g

def selftest():
    from checks import Check
    c=Check();x=CrosslinkedRNP('rna1',frozenset({'p1','p2'}),frozenset({'p1'}),frozenset({'p1'}));p=x.proteolyzed()
    c('selected protein-label/crosslink path can recover linked RNA',x.selected_capture_path(),True)
    c('uncrosslinked labeled protein does not supply declared RNAcapture path',replace(x,crosslinked_proteins=frozenset({'p2'})).selected_capture_path(),False)
    c('protein digestion preserves RNAidentity',p['rna_id'],'rna1')
    c('proteinborne biotin is not silently transferred onto RNA',p['direct_rna_biotinylation'],False)
    c('peptide/remnantfates remain possible after proteolysis',p['possible_peptide_remnants'],True)
    c('full sourceknown route retains publicationboundary',bool(protocol_graph().publication_problems()),True);c.report()
if __name__=='__main__':selftest()

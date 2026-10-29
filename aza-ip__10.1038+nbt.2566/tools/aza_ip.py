"""Aza-IP selected covalent trap and inferred ring-open copying explanation."""
from dataclasses import dataclass,replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem
SOURCE='Khoddami/Cairns2013 DOI10.1038/nbt.2566; associated accepted manuscript PMC3791587 mechanism27–42, Methods90–100, Fig2caption372; VOR applicability unverified'

@dataclass(frozen=True)
class Site:
    origin: str='RNA_5azaC'
    stage: str='enzyme_bound'
    inferred: bool=False
    first_dna: str|None=None
    def __post_init__(self):
        if self.origin!='RNA_5azaC' or self.stage not in {'enzyme_bound','processed_unassigned','ring_open_hypothesis'}:raise ValueError('RNA5azaC-selected lineage only; not DNA5aza-dC or genericm5C')
        if self.stage=='ring_open_hypothesis' and not self.inferred:raise ValueError('proposed RNA ring opening must retain inference qualification')
        if self.first_dna is not None and (self.stage!='ring_open_hypothesis' or self.first_dna!='C'):raise ValueError('only selected hypothesized C-pairing explanation represented; other outcomes unassigned')
    def processed(self):return replace(self,stage='processed_unassigned')
    def proposed_ring_open(self):return replace(self,stage='ring_open_hypothesis',inferred=True)
    def copy_selected(self):
        if self.stage!='ring_open_hypothesis':raise ValueError('cannot assert selected copying identity of unknown/closed processed substrate')
        return replace(self,first_dna='C')
    @property
    def reported(self):
        if self.first_dna is None:raise ValueError('DNA/reference report requires copying first')
        return 'G'

def protocol_graph():
    g=ProtocolGraph('Aza-IP — fragmentation and RNA recovery prefix',authority=SOURCE,representation='published_instance')
    prior=g.input('bound_rnp','Source washed antiV5-IP beads with covalent RNA5azaC–RMT complexes; input/bound/free proportions unassigned',K.MOLECULAR,primary=True,chemical_state=Site())
    fragmentation=g.input('fragmentation_reagent','Ambion fragmentation reagent94°C5min then ice1min',K.CONDITION)
    mixture=g.product('fragmented','Fragmented/processed RNA, beads/RMT and released molecular material; site bond state and release fractions unassigned',K.MOLECULAR,chemical_state=Site().processed())
    g.react('fragment','Source bead-associated RNA fragmentation; no universal complete adduct release inferred',S.FRAGMENTATION,((prior,R.TEMPLATE),(fragmentation,R.CONDITION)),(mixture,))
    stop=g.input('stop_solution','Source diluted Ambion fragmentation stop solution',K.CONDITION)
    stopped=g.product('stopped','Stopped fragmentation mixture; all material retained before separation',K.MOLECULAR,chemical_state=Site().processed())
    g.react('stop','Source stop addition after cooling',S.QUENCHING,((mixture,R.TEMPLATE),(stop,R.CONDITION)),(stopped,))
    collect=g.input('supernatant_collection','Source transfer of soluble fractions from IP samples within each replicate',K.SELECTION)
    supernatant=g.product('supernatant','Collected RNA-containing soluble population; exact site/adduct chemistry unassigned',K.MOLECULAR,chemical_state=Site().processed())
    remaining=g.product('remaining_beads','Untransferred bead/RMT-associated material; exact remaining RNA fraction and final fate unassigned',K.MOLECULAR)
    g.react('collect','Physical source supernatant transfer; do not pool biological replicates in this model',S.SELECTION,((stopped,R.TEMPLATE),(collect,R.SELECTION)),(supernatant,remaining))
    g.mark_terminal(Terminal(remaining.id,Disposition.UNRESOLVED,'Final fate of remaining bead-associated material is not assigned by this scoped model.'))
    precipitation=g.input('precipitation_reagents','Source GlycoBlue/sodiumacetate/ethanol precipitation and70%ethanol wash',K.SELECTION)
    pellet=g.product('pellet','Recovered RNA-containing pellet with possible carrier/adduct/reagent carryover unquantified',K.MOLECULAR,chemical_state=Site().processed())
    waste=g.product('ethanol_waste','Physical precipitation supernatants/washes',K.WASTE)
    g.react('precipitate','Source precipitation before Trizol extraction',S.CLEANUP,((supernatant,R.TEMPLATE),(precipitation,R.SELECTION)),(pellet,waste));g.discard(waste,'Physical precipitation/wash fractions')
    trizol=g.input('trizol','Source Trizol RNA extraction procedure',K.SELECTION)
    rna=g.product('recovered_rna','Physically extracted RNA-containing material; residual adduct chemistry/purity unassigned',K.MOLECULAR,chemical_state=Site().processed())
    extract_waste=g.product('extraction_waste','Physical source extraction rejects',K.WASTE)
    g.react('extract','Source RNA extraction after fragmentation/precipitation',S.EXTRACTION,((pellet,R.TEMPLATE),(trizol,R.SELECTION)),(rna,extract_waste));g.discard(extract_waste,'Physical source extraction rejects; no perfect depletion asserted')
    g.mark_terminal(Terminal(rna.id,Disposition.UNRESOLVED,'Exact historical directional RNA-adapter/RT/PCR/run geometry, actual processed-site chemistry and final VOR applicability remain unassembled. Hypothesized ring-opening copy path is a separately qualified explanation, not a graph-proven reaction.'))
    g.declare_reviewed_scope(['bound_rnp','fragmentation_reagent','fragment','stop_solution','stop','supernatant_collection','collect','precipitation_reagents','precipitate','trizol','extract'],SOURCE)
    g.add_source_item(SourceItem('unreviewed_full_library','Final VOR, complete directional mRNA library polymers/read plan and proposed ring-open mechanism applicability','DOI10.1038/nbt.2566; cited earlier ring-opening mechanism and detailed2014 protocol'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();s=Site().processed();p=s.proposed_ring_open().copy_selected()
    c('physical RNA origin remains distinct from DNA/referenceG signal',(p.origin,p.stage,p.first_dna,p.reported,p.inferred),('RNA_5azaC','ring_open_hypothesis','C','G',True))
    try:s.copy_selected();bad=False
    except ValueError:bad=True
    c('unassigned processed RNA cannot silently become proposed ring-open copy',bad,True)
    try:Site(stage='ring_open_hypothesis',inferred=False);bad=False
    except ValueError:bad=True
    c('inference qualification enforced by molecular constructor',bad,True)
    c('remaining bead material is not silently discarded',protocol_graph().terminals['remaining_beads'].disposition,Disposition.UNRESOLVED)
    original=Source(sources_root()/'aza-ip__10.1038+nbt.2566'/'nbt.2566_PMC3791587.xml.txt','associated manuscript Methods96/100')
    if have(c,original,label='source recovery/readout'):
        text=' '.join(original.path.read_text().split());c('source library route follows recovered RNA','50-cycle single-end' in text,True)
    c.report()
if __name__=='__main__':selftest()

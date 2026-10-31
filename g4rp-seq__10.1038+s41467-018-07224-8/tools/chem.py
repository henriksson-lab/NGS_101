"""Protocol-owned source-qualified partial model; no completed-publication claim."""
from dataclasses import dataclass,replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from pcr import PCRGeometry

def state(label):
    sc=Scene();sc.strand('population',[Segment(label,'['+label+']',placeholder=True)]);return sc

def clean(g,key,current,procedure):
    reagent=g.input(key+'_reagent','Actual TRIZOL molecular extraction reagent; quantity/partition chemistry unassembled',K.MOLECULAR)
    kept=g.product(key+'_kept','Recovered '+procedure+'; exact molecular recovery/carryover unquantified',K.MOLECULAR)
    waste=g.product(key+'_waste','Physical separation fractions; no complete target specificity assumed',K.WASTE)
    g.react(key,procedure,S.CLEANUP,((current,R.TEMPLATE),(reagent,R.STRUCTURAL)),(kept,waste));g.discard(waste,'Actual source physical separation only.');return kept
SOURCE='Yang et al.2018 DOI10.1038/s41467-018-07224-8; published non-manuscript PMC6226477 Methods138–148, limitations116 and SI Fig6/7. Exact historical kits unreviewed.'
@dataclass(frozen=True)
class SelectedPartition:
    inputs:tuple[str,...]
    retained:frozenset[str]
    def __post_init__(self):
        if not self.inputs or len(set(self.inputs))!=len(self.inputs) or not self.retained<=set(self.inputs):raise ValueError('explicit unique molecular identities and selected binding outcome required')
    @property
    def flowthrough(self):return frozenset(self.inputs)-self.retained
    def recovered(self,retained_out):
        if not retained_out<=self.retained:raise ValueError('release cannot invent molecular identities')
        return retained_out,self.retained-retained_out

def protocol_graph():
    g=ProtocolGraph('G4RP-seq — source baseline MCF7 route',authority=SOURCE,representation='published_instance')
    current=g.input('cells','Source baseline MCF7 cells; optional drug-treated experiments remain separate inputs, no ligand exposure inferred',K.MOLECULAR,primary=True)
    fa=g.input('formaldehyde','Actual source 1% formaldehyde/PBS molecular crosslinking contribution',K.MOLECULAR)
    fixation_environment=g.input('fixation_environment','Source 5min25°C fixation environment',K.CONDITION)
    out=g.product('crosslinked','Crosslinked cellular RNA/protein/nucleic acid complexes; no guarantee all folds trapped or direct interactions only',K.MOLECULAR)
    g.react('crosslink','Source formaldehyde crosslink before BioTASQ exposure',S.CROSS_LINKING,((current,R.TEMPLATE),(fa,R.STRUCTURAL),(fixation_environment,R.CONDITION)),(out,))
    gly=g.input('glycine','Actual source0.125M glycine molecular quencher',K.MOLECULAR);quench_environment=g.input('quench_environment','Source glycine quench5min; exact stoichiometry unassembled',K.CONDITION);quenched=g.product('quenched','Quenched crosslinked material plus residual reagents/products',K.MOLECULAR)
    g.react('quench','Source glycine quench',S.QUENCHING,((out,R.TEMPLATE),(gly,R.STRUCTURAL),(quench_environment,R.CONDITION)),(quenched,),reason='Reaction product identities and exact quench stoichiometry unassembled.')
    buffer=g.input('lysis_buffer','Actual source G4RP working-buffer aliquot for scraped/resuspended cells; totalvolume unassembled',K.MOLECULAR)
    components=[]
    for key,label in (('kcl','150mM KCl'),('tris','25mM Tris pH7.4'),('edta','5mM EDTA'),('dtt','0.5mM DTT'),('np40','0.5% NP40'),('rnase_inhibitor','Roche RNase inhibitor; dose unassembled'),('protease_inhibitors','Homebrew protease inhibitor cocktail; components/dose unassembled')):
        components.append(g.product('g4rp_'+key,'Co-present '+label+' from actual G4RP buffer, not a fresh addition',K.MOLECULAR))
    g.react('buffer_components','Expose source working-buffer co-components; no extra doses or invented cocktail',S.GENERIC,((buffer,R.STRUCTURAL),),tuple(components),changes_molecule=False)
    lysate=g.product('lysate','Sonicated crosslinked cell extract with RNA, DNA/protein/debris retained',K.MOLECULAR)
    g.react('sonicate','Source Covarism22010%duty2min after scraping/resuspension',S.MAGIC,((quenched,R.TEMPLATE),)+tuple((z,R.STRUCTURAL) for z in components),(lysate,),reason='Source extraction/sonication operation known; molecular break/retention distributions unassembled.')
    probe=g.input('biotasq','Source100µMBioTASQ overnight4°C, biotinylated G4-binding probe; noncovalent association, not RNA-biotin covalent labeling',K.MOLECULAR)
    association=g.product('probe_associated','BioTASQ-associated crosslinked RNA complexes plus free probe and nonspecific/background material',K.MOLECULAR)
    g.react('probe_binding','Source noncovalent BioTASQ association; source-specific folding/binding not automatically predicted',S.BINDING,((lysate,R.TEMPLATE),(probe,R.STRUCTURAL)),(association,))
    input_control=g.product('input_control','Source5%sonicate input-control aliquot; exact molecular partition unassigned',K.MOLECULAR)
    remaining=g.product('remaining','Remaining preparation for source affinity capture',K.MOLECULAR)
    g.react('aliquot','Source5%sonicate collected as input control; molecular selectivity not inferred',S.SPLIT,((association,R.TEMPLATE),),(input_control,remaining))
    g.mark_terminal(Terminal(input_control.id,Disposition.UNRESOLVED,'Independent source input-control RNA/library route not assembled; not physically pooled into pull-down.'))
    beads=g.input('beads','Source10µgstreptavidin magneticbeads2h4°C;4G4RP-buffer washes5min',K.SELECTION)
    bound=g.product('captured','Bead-associated BioTASQ complexes plus possible nonspecific/carryover material',K.MOLECULAR)
    flow=g.product('flowthrough','Actual capture flowthrough/wash fractions; no perfect G4 specificity asserted',K.WASTE)
    g.react('capture','Source streptavidin capture and washes',S.AFFINITY_CAPTURE,((remaining,R.TEMPLATE),(beads,R.SELECTION)),(bound,flow));g.discard(flow,'Physical bead isolation/washes, exact bound/unbound RNA distribution unquantified.')
    heat=g.input('release_heat','Source70°C1h reverse-crosslink incubation',K.CONDITION)
    released=g.product('released','Heat-treated capture material, released RNA plus beads/probe/protein/nucleic acid remnants; no complete crosslink reversal assumed',K.MOLECULAR)
    g.react('reverse_crosslink','Source heated reversal on beads; no physical removal at heat alone',S.THERMAL_TREATMENT,((bound,R.TEMPLATE),(heat,R.CONDITION)),(released,))
    current=clean(g,'trizol',released,'Source TRIZOL RNA extraction from beads per manufacturer')
    dnase=g.input('dnase','Source DNase treatment, exact preparation/workup unreviewed',K.CATALYST)
    treated=g.product('dnase_product','Source RNA-containing DNase-treated mixture with DNA digestion debris retained',K.MOLECULAR)
    g.react('dnase_digest','Source DNAase treatment before thermal fragmentation',S.DIGESTION,((current,R.TEMPLATE),(dnase,R.CATALYST)),(treated,))
    frag=g.input('fragment_condition','Source brief thermal fragmentation; exact conditions unreviewed',K.CONDITION)
    fragments=g.product('fragmented','RNA fragments and reaction products retained; no exact ends/size distribution assigned',K.MOLECULAR)
    g.react('fragment','Source thermal RNA fragmentation',S.FRAGMENTATION,((treated,R.TEMPLATE),(frag,R.CONDITION)),(fragments,))
    ribozero=g.input('ribozero','Source IlluminaRibo-Zero rRNA removal kit; historical edition/internal mechanism unreviewed',K.CONDITION)
    depleted=g.product('ribodepleted','Source kit-processed preparation containing ribodepleted RNA and unresolved RNA/DNA/remnant/carryover fates; no physical separation or disposal asserted',K.MOLECULAR)
    g.react('ribodeplete','Opaque source Ribo-Zero processing after fragmentation',S.MAGIC,((fragments,R.TEMPLATE),(ribozero,R.CONDITION)),(depleted,),reason='Main names historical kit only; no affinity mechanism, internal reactions or physical removal/discard inferred. Molecular remnant/partition fates remain unresolved.')
    kit=g.input('truseq','Source Illumina TruSeq RNA Library Prep kit; exact edition/strand/oligos/index member unreviewed',K.CONDITION)
    library=g.product('library','Source RNA-derived DNA library; exact physical template/remnant/end geometry unassembled',K.LIBRARY)
    g.react('library_prep','Source named TruSeq library preparation handoff',S.MAGIC,((depleted,R.TEMPLATE),(kit,R.CONDITION)),(library,),reason='Historical kit internal transitions/read assembly unreviewed; no modern kit substitute.')
    g.sequence(library,SequencingPlan(acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':75,'Read 2':75},limitation='Source NextSeq500PE75; full adapter/strand/run-primer/index geometry unassembled.'))
    g.add_source_item(SourceItem('unreviewed_kit','Named historical RiboZero/TruSeq/TRIZOL documents and full read geometry','Published mainMethods138/148'))
    return g

def selftest():
    from checks import Check
    c=Check();x=SelectedPartition(('foldedRNA','indirectRNA','backgroundDNA'),frozenset({'foldedRNA','indirectRNA'}))
    c('selected capture preserves molecular input identities across fractions',x.retained|x.flowthrough,frozenset(x.inputs))
    c('selected capture does not assert direct-G4-only specificity',x.retained,frozenset({'foldedRNA','indirectRNA'}))
    c('release retains unreleased material explicitly',x.recovered(frozenset({'foldedRNA'})),(frozenset({'foldedRNA'}),frozenset({'indirectRNA'})))
    try:x.recovered(frozenset({'inventedRNA'}));bad=False
    except ValueError:bad=True
    c('release cannot invent RNA identities',bad,True)
    g=protocol_graph()
    c('Actual formaldehyde and glycine are molecular contributors',(g.materials['formaldehyde'].kind,g.materials['glycine'].kind),(K.MOLECULAR,K.MOLECULAR))
    c('Source glycine has typed quenching operation',g.steps['quench'].kind,S.QUENCHING)
    c('Known G4RP co-components derive from one working buffer',all(g.producer_of(x).id=='buffer_components' for x in ('g4rp_kcl','g4rp_tris','g4rp_edta','g4rp_dtt','g4rp_np40','g4rp_rnase_inhibitor','g4rp_protease_inhibitors')),True)
    c('TRIZOL extraction reagent is not a selection matrix',(g.materials['trizol_reagent'].kind,g.steps['trizol'].inputs[1].role),(K.MOLECULAR,R.STRUCTURAL))
    c('Actual historical-kit/fullread item remains uncovered',any("'unreviewed_kit'" in p for p in g.publication_problems()),True)
    c('known downstream route remains publication-limited',bool(protocol_graph().publication_problems()),True);c.report()
if __name__=='__main__':selftest()

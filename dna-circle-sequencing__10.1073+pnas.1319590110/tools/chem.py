"""Source-local 2013 DNA circle-sequencing preparation, with explicit unresolved SI."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
from chemdraw import Construct,Segment,Scene
from circular import circularize_ssdna,rolling_circle,RCAProduct,CircularizationGeometry
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan
SOURCE='Printed finalPNAS2013 DOI10.1073/pnas.1319590110, UC repository original: Results128–177/215–260,Methods499–513. SI MaterialsMethods unavailable via printed index403; exact wet-workup/adapter details unassembled.'

@dataclass(frozen=True)
class SelectedRCAReference:
    """Caller-selected exact pairing/traversal, not a genomic product or yield prediction."""
    product:RCAProduct
    inferred:bool=True
    def scene(self):
        s=Scene();s.strand('caller-selected ideal copy',[Segment('selected traversal '+str(i+1),self.product.unit,inferred=self.inferred)for i in range(self.product.copies)]);return s

def selected_reference(template,primer,copies,*,five_prime_phosphate,three_prime_oh):
    if not template or set(template+primer)-set('ACGT') or type(copies)is not int or copies<1:raise ValueError('Explicit DNA, primer and positive integer traversal count required')
    c=Construct([Segment('caller-selected strand',template,inferred=True)],'Caller-selected hypothetical circular substrate')
    circle=circularize_ssdna(c,five_prime_phosphate=five_prime_phosphate,three_prime_oh=three_prime_oh)
    return SelectedRCAReference(rolling_circle(circle,primer,copies=copies))

def opaque(g,key,x,label):
    c=g.input(key+'_condition','Source conditions/handling not fully assembled',K.CONDITION);y=g.product(key+'_product',label+'; complete molecular mixture and residual classes unassembled',K.MOLECULAR)
    g.react(key,label,S.MAGIC,((x,R.TEMPLATE),(c,R.CONDITION)),(y,),reason='Main establishes selected operation; detailed SI/geometry/handling unassembled.');return y

def protocol_graph():
    g=ProtocolGraph('DNA circle sequencing — finalmain glycosylase-containing route',authority=SOURCE,representation='published_instance')
    x=g.input('genomic','Source extractedS288C yeastgDNA; reference sequence/cut coordinates unassigned',K.MOLECULAR,primary=True)
    shear=g.input('shearing','Source shearing conditions unassigned in inspected main',K.CONDITION);y=g.product('sheared','Fragmented genomicDNA plus all fragments; no exact boundaries predicted',K.MOLECULAR);g.react('shear','Source genomicDNA shearing',S.FRAGMENTATION,((x,R.TEMPLATE),(shear,R.CONDITION)),(y,))
    gel=g.input('gel','Source1.5percentlowmeltingagarose/narrow150bpslice; full extraction workup unreviewed',K.SELECTION);x=g.product('selected','Source physically selected DNA fraction; no perfect length cutoff or purity assertion',K.MOLECULAR);w=g.product('outside_slice','Actual material outside selected gel slice',K.WASTE);g.react('gel_select','Source narrow150bp nominal slice',S.SELECTION,((y,R.TEMPLATE),(gel,R.SELECTION)),(x,w));g.discard(w,'Physical gel selection, not computationally deleted reads.')
    e=g.input('phosphorylation_stock','Source phosphorylation catalyst/buffer; identity and exact end populations unreviewed',K.CATALYST);y=g.product('phosphorylated','Source phosphorylated DNA and reagent/partial-conversion classes; full end inventory unassembled',K.MOLECULAR);g.react('phosphorylate','Source DNA phosphorylation',S.PHOSPHORYLATION,((x,R.TEMPLATE),(e,R.CATALYST)),(y,))
    x=opaque(g,'denature',y,'Source DNA denaturation; both original strand classes and residualduplex/other material retained in same preparation')
    ligase=g.input('circligase','SourceCircLigaseIIssDNAligase;300ngDNA/~3pmol/20ulreaction; full formulation unreviewed',K.CATALYST);y=g.product('closure_mixture','Source ssDNAcircle/uncircularizedDNA/reagent populations; no perfect closure or end product count asserted',K.MOLECULAR)
    g.react('circle','Source single-strand intramolecular circularization',S.CIRCULARIZATION,((x,R.TEMPLATE),(ligase,R.CATALYST)),(y,),geometry=CircularizationGeometry.known('single-strand intramolecular ligation',statement='MainFig1B/Methods explicitly denaturedssDNA andCircLigaseII; exact source base/end/junction inventory unassembled.'))
    exo=g.input('exonuclease','Source exonuclease for uncircularized DNA; enzyme/end specificity unassigned in main',K.CATALYST);x=g.product('post_exo','Circular/uncut populations plus linearDNAhydrolysis products and retained enzyme/reagents; no undocumented physical cleanup',K.MOLECULAR);g.react('digest_linear','Source exonucleolytic processing of uncircularizedDNA',S.DIGESTION,((y,R.TEMPLATE),(exo,R.CATALYST)),(x,))
    primer=g.input('random_primers','Source exonuclease-resistant random DNA primers; length/bond modification/whole sequence unassigned',K.OLIGO);c=g.input('anneal_conditions','Source annealing conditions; exact temp/time/workup in inaccessibleSI',K.CONDITION);y=g.product('primed_mixture','Source primed/unprimed circles, remaining substrates and primers/reagents',K.MOLECULAR);g.react('anneal','Source random primers annealed to variable circle amounts',S.HYBRIDIZATION,((x,R.TEMPLATE),(primer,R.STRUCTURAL),(c,R.CONDITION)),(y,))
    pol=g.input('phi29','SourcePhi29DNApolymerase',K.CATALYST);ppase=g.input('pyrophosphatase','Sourceinorganicpyrophosphatase',K.CATALYST);glyco=g.input('glycosylases','MainselectedUDG+Fpg; earlier/noenzyme comparator remains separate; excision is not faithful base repair',K.CATALYST);mix=g.input('rca_mix','Source reactionbuffer/dNTPs/BSA; exact composition/time unreviewed',K.CONDITION)
    x=g.product('rca_mixture','Source tandemcopy ss/dsDNA products, original circles/otherDNA, primers and reagent/lesion-processing remnants; no copy-number/yield/independent-error guarantee',K.MOLECULAR)
    g.react('rca','SourcePhi29 strand displacement and secondary random priming; UDG/Fpg co-treatment',S.EXTENSION,((y,R.TEMPLATE),(pol,R.CATALYST),(ppase,R.CATALYST),(glyco,R.CATALYST),(mix,R.CONDITION)),(x,))
    # Main gives no complete downstream template conversion/adapter sequence recipe.
    stock=g.input('downstream_stock','Historical source-compatible Illumina library conversion components; actual workup/end/adapter chemistry unassembled',K.CONDITION);lib=g.product('library','Source MiSeq sequencing preparation; full concatemer-to-library molecule and companions unassembled',K.LIBRARY)
    g.react('ngs_boundary','Unassembled downstream workup/adapter/library/acquisition boundary',S.MAGIC,((x,R.TEMPLATE),(stock,R.CONDITION)),(lib,),reason='Definingmain says protocol precedes standardlibrarypreparation; actualsourceSIunavailable. No modern kit or Duplex comparator constructor imported.')
    g.sequence(lib,SequencingPlan(limitation='SourceMiSeq PE150/PE250 reported; selected run/cycles/indices/full primer geometry unassigned. Consensus/periodicity/reference filters computational. No proven error-rate, molecular independence or perfectdamage removal claim.'))
    g.add_source_item(SourceItem('pending','ActualSI MaterialsMethods/FigS1 detailed preparation and end/workup/library inventories','Printed supplement index returned403; main repositoryoriginalread. Randomprimer resistance not silently translated to PS/hexamer.2014comment/reply identities noted,not inspected; no performance claim.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();ref=selected_reference('ACGTTGCA','AAC',2,five_prime_phosphate=True,three_prime_oh=True)
    c('Selectedcopycandidate inferredbydefault',ref.inferred,True);c('Selectedcopy scene inherits uncertainty',bool(ref.scene()),True)
    try:selected_reference('ACGTTGCA','AAC',2,five_prime_phosphate=False,three_prime_oh=True);bad=False
    except ValueError:bad=True
    c('No ideal closure without explicitphosphate',bad,True)
    g=protocol_graph();c('Main exonuclease not drawn as invented beadcleanup',g.steps['anneal'].inputs[0].material,'post_exo');c('RCA uses already annealed primers without a second addition',not any(i.material=='random_primers'for i in g.steps['rca'].inputs),True);c('Glycosylases remain source co-treatment inputs',any(i.material=='glycosylases'for i in g.steps['rca'].inputs),True);c('Downstream comparator library not silently inherited',bool(g.publication_problems()),True)
    src=Source(sources_root()/'dna-circle-sequencing__10.1073+pnas.1319590110'/'source_pdf_25b44e781411953c4d40.pdf.txt','UC published2013original,Methods499–505')
    if have(c,src,label='PublishedMain'):
        t=src.path.read_text();c('Printed definingtitle/DOI identity','High-throughput DNA sequencing errors are reduced'in t and'10.1073/pnas.1319590110'in t,True);c('ActualCircleMethods catalyst inventory','CircLigase II ssDNA ligase'in t and'inorganic pyrophosphatase, uracil-DNA glycosylase'in t,True)
    c.report()
if __name__=='__main__':selftest()

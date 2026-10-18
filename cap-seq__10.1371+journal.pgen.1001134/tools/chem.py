"""Source-owned2010 CAP-seq affinity chemistry and bounded predecessor stages."""
from chemdraw import Construct,Scene,Segment,complement_segments
from protocol_model import role
from protocol_graph import ProtocolGraph,MaterialKind as M,InputRole as I,StepKind as K,SequencingPlan
from adapter_ligation import AdapterLigationGeometry
import seqprimers as sp
TITLE="CAP-seq — native CXXC affinity purification"
NOTES="01_cap_seq.html"
SOURCE='<a href="https://doi.org/10.1371/journal.pgen.1001134">Illingworth et al.,PLOS Genetics2010</a>; Methods DNAchromatography and Solexasequencing.'
SUMMARY="Solexa adapters precede CXXC enrichment of native genomic DNA. Salt washes and elution enrich clustered unmethylated CpGs without converting DNA bases."
CAVEAT="Selected human600mM wash route; mouse560mM is a separate experimental setting. Defining2010Methods delegates input preparation, adapters and final library construction to Skene2010MAP-seq. The directly cited Skene author-manuscript Methods points to SupplementalData, which was not inspected in this capped pass. Native fragment ends, adapter/PCR/loadedprimer identities and cycles remain unassigned; no canonical Solexa shell or postcapture PCR order is supplied."
CAP_HUMAN_WASH_MM=600
CAP_MOUSE_WASH_MM=560
SEQ_PRIMERS=()
FINAL_LIBRARY=Construct([role('native Solexa endA'),role('native CXXC-enriched genomicDNA; methylation unchanged'),role('native Solexa endB')],name='CAP native library componentreference')
FINAL_CAPTION="Native adapter roles; complete predecessor strand assembly unassigned."
SEQUENCING_ENDING="Defining source specifies single-end Solexa; acquired read length and index acquisition unassigned. Mapping score≥30 is a computational filter."
READ_SPECS=(sp.read_spec('Read 1','unassigned','genomicDNA; single-end Solexa','include'),sp.read_spec('Index 1 (i7)','unassigned','nativeindex identity/acquisition unassigned','undetermined'),sp.read_spec('Index 2 (i5)','unassigned','nativeindex identity/acquisition unassigned','undetermined'))
def dna_state(adapted=False,enriched=False):
 sc=Scene();parts=[role('native genomicDNA; CpGmethylation retained; enrichment is not conversion')]
 if adapted:parts=[role('native Solexa adapterA; predecessorgeometry unassigned'),*parts,role('native Solexa adapterB; predecessorgeometry unassigned')]
 sc.strand('DNA',parts,label='CXXC-enriched population; not a pureunmethylated class' if enriched else 'original genomicDNA')
 sc.strand('opposite',[role('opposite original genomicDNA strand; native ends/adapterconnectivity unassigned')],label='independent opposite strand')
 sc.strand('residuals',[role('possible co-carried inputpolymers and off-affinity DNA; class-specific recovery unassigned')],label='independent residualpopulations')
 if enriched:sc.strand('reagent_residuals',[role('possible CXXC-preparation bacterialDNA/DNase-generatedfragments or residualDNaseI; physicalcarryover unassigned')],label='independent reagent-origin residuals')
 return sc

def protocol_graph():
 source='https://doi.org/10.1371/journal.pgen.1001134; Methods218/226; Skene2010PMC4338610Methods113/121'
 g=ProtocolGraph(TITLE,authority=source)
 raw=g.input('genomic_dna','original genomicDNA sample',M.BIOLOGICAL,state=dna_state(),primary=True)
 prep=g.input('map_prep','cited Skene2010MAP inputDNApreparation',M.CONDITION)
 prepared=g.product('prepared','native preparedDNA; fragment ends/geometry unassigned',M.MOLECULAR,state=dna_state())
 g.react('prepare','Native MAP-seq inputDNApreparation boundary',K.MAGIC,((raw,I.TEMPLATE),(prep,I.CONDITION)),(prepared,),reason='DefiningMethods delegates preparation to Skene2010. Inspected predecessor mainMethods refers to SupplementalData; those detailed chemistry sections were not inspected. No fragmenting/repair defaults supplied.')
 adapters=g.input('adapters','native Solexa adapters from citedMAP preparation',M.ADAPTER)
 adapted=g.product('adapted','adapter-ligated genomicDNA before CAP',M.MOLECULAR,state=dna_state(True))
 g.react('ligate','Sourcepositive Solexa adapter ligation beforeCXXCpurification',K.LIGATION,((prepared,I.TEMPLATE),(adapters,I.STRUCTURAL)),(adapted,),geometry=AdapterLigationGeometry.unavailable('PublishedMethods226 gives operation order; specific adapterstrands/endgeometry/ligationreagents delegated to unread predecessorSupplementalData.'))
 protein_state=Scene();protein_state.strand('protein',[role('purified recombinantCXXC protein')]);protein_state.strand('bacterialDNA',[role('possible contaminatingbacterialDNA targetedbyDNaseI')])
 protein=g.input('cxxc_protein','CXXC reagentpreparation withpossiblebacterialDNA;100mg source amount',M.MOLECULAR,state=protein_state)
 dnase=g.input('dnase','DNaseI700U;roomtemperature90min; acts on CXXC preparation, not sampleDNA',M.CATALYST)
 treated_state=Scene();treated_state.strand('protein',[role('originalCXXC protein; not aDNaseIsubstrate')]);treated_state.strand('DNAproducts',[role('bacterialDNA/DNaseIdigestionproducts; hydrolysisnotphysicalremoval')])
 treated=g.product('dnase_treated_cxxc','CXXC preparation after bacterialDNA digestion',M.MOLECULAR,state=treated_state)
 g.react('clean_cxxc','Digest bacterialDNA in CXXC reagent preparation',K.DIGESTION,((protein,I.TEMPLATE),(dnase,I.CATALYST)),(treated,))
 ni=g.input('nickel','nickel-charged Sepharose;50mgCXXC/mlbeadvolume',M.SELECTION)
 immobilized=g.product('immobilized_cxxc','DNaseI-treated CXXC on nickelSepharose; reagentDNAfragmentretentionunassigned',M.SELECTION,state=treated_state)
 g.react('immobilize','Immobilize recombinantCXXC on nickelSepharose',K.BINDING,((treated,I.TEMPLATE),(ni,I.SELECTION)),(immobilized,))
 buffer=g.input('binding_buffer','columnbuffer0.1MNaCl',M.CONDITION)
 bound=g.product('bound','CXXC-bound adapterDNA population',M.MOLECULAR,state=dna_state(True,True));flow=g.product('flowthrough','physically unbound columnfraction',M.WASTE)
 g.react('bind','Bind adapterDNA toCXXCmatrix in0.1MNaCl',K.AFFINITY_CAPTURE,((adapted,I.TEMPLATE),(immobilized,I.SELECTION),(buffer,I.CONDITION)),(bound,flow));g.discard(flow,'Columnflowthrough physically separated; no complete methylationclass specificity asserted.')
 wash=g.input('wash_buffer','selectedhuman600mMNaCl; mouse560mM is separate',M.CONDITION)
 washed=g.product('washed','salt-washed retainedDNA',M.MOLECULAR,state=bound.state);waste=g.product('wash_waste','physically removed washfractions',M.WASTE)
 g.react('wash','Wash CXXCcolumn with600mMNaCl',K.CLEANUP,((bound,I.TEMPLATE),(wash,I.CONDITION)),(washed,waste));g.discard(waste,'Actual columnwash fractions; no CpGbase is chemically changed.')
 salt=g.input('elution','1MNaCl elutionbuffer',M.CONDITION)
 eluate=g.product('eluate','pooled1MNaCl eluate fractions',M.MOLECULAR,state=washed.state)
 g.react('elute','Elute with1MNaCl and poolfractions',K.EXTRACTION,((washed,I.TEMPLATE),(salt,I.CONDITION)),(eluate,))
 precip=g.input('concentration_precipitation','concentrate andprecipitate as citedMAP; detailedconditionsunassigned',M.SELECTION)
 recovered=g.product('recovered','recovered CAPDNA',M.MOLECULAR,state=eluate.state);sup=g.product('precip_supernatant','physically unretained precipitationfraction',M.WASTE)
 g.react('precipitate','Concentrate andprecipitate pooledeluates',K.CLEANUP,((eluate,I.TEMPLATE),(precip,I.SELECTION)),(recovered,sup));g.discard(sup,'Precipitation supernatant; classspecificpurity/recovery unassigned.')
 native=g.input('native_library_reagents','citedSkeneMAP libraryconstruction; detailedSupplementalDataunread',M.CONDITION)
 lib=g.product('library','native CAP Solexalibrary; exactpostenrichmentoperationsunassigned',M.LIBRARY,state=recovered.state)
 g.react('native_library_boundary','PostCAP nativeSolexalibrary handoff',K.MAGIC,((recovered,I.TEMPLATE),(native,I.CONDITION)),(lib,),reason='DefiningMethods states CAP libraries prepared as previouslydescribed forMAP. Genuine unread predecessorSupplementalData boundary: no PCR/gel/ligationrounds invented aftercapture.')
 g.sequence(lib,SequencingPlan(acquired_roles=('Read 1',),limitation=SEQUENCING_ENDING))
 g.declare_reviewed_scope(('genomic_dna','map_prep','prepare','adapters','ligate','cxxc_protein','dnase','clean_cxxc','nickel','immobilize','binding_buffer','bind','wash_buffer','wash','elution','elute','concentration_precipitation','precipitate','native_library_reagents','native_library_boundary','library'),source)
 return g.require_publishable()

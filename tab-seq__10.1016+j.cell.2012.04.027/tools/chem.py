"""Yu2012 TAB-seq: source-owned protection/oxidation and native WGBS boundary."""
from types import SimpleNamespace
from chemdraw import Construct, Scene
from protocol_model import role
from base_conversion import BasePath
from pcr import PCRGeometry
from protocol_graph import Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
import seqprimers as sp

def tab_path(base):
 paths={'5hmC':('5hmC','5gmC','5gmC','C'),'5mC':('5mC','5caC','U','T'),'C':('C','C','U','T')}
 if base not in paths:raise ValueError('TAB-seq input must be C,5mC or5hmC')
 return BasePath(paths[base])
GLUCOSYLATION_MIN=60
OXIDATION_2012_MIN=90
TITLE='TAB-seq — TET-assisted bisulfite sequencing'
NOTES='01_tab_seq.html'
SOURCE='<a href="https://doi.org/10.1016/j.cell.2012.04.027">Yu et al., Cell2012</a>, ExperimentalProcedures144–146; authors\' later published2012 detailed protocol PMC3641661 is a separately scoped companion reference.'
SUMMARY='βGT protects5hmC as5gmC; mTet1 oxidizes5mC to5caC. After bisulfite conversion and DNA copying, original5hmC positions preferentially reportC while C/5mC positions reportT.'
CAVEAT='Physical5gmC-bearing originals are distinct from normalC daughter-DNA readout. Conversion is not assumed complete. Native WGBS adapter/index/primer assembly and acquisition remain unassigned; the2012 whole-genome protocol points to Illumina guide15021861A, which was not inspected in this capped pass. No canonical TruSeq shell or fixed index width is substituted. The2012 locus-specific PfuTurboCx assay is not silently assigned to whole-genome libraries.'
FINAL_LIBRARY=Construct([role('native WGBS endA; sequence/geometry unassigned'),role('daughter DNA:original5hmC can readC; originalC/5mC can readT'),role('native opposite/index end; identity/count/width unassigned')],name='TAB native WGBS readout reference')
FINAL_CAPTION='Native endpoint roles; original modification is not copied as a glucose group onto daughter DNA.'
SEQ_PRIMERS=()
SEQUENCING_INTRO='The authors\'2012 whole-genome companion protocol names HiSeq2000; defining2012 study cycle allocation and native read/index roles remain unassigned here.'
SEQUENCING_UNAVAILABLE='Selected historical WGBS native adapters, primer landing and index/read acquisition remain unassigned.'
SEQUENCING_ENDING=SEQUENCING_UNAVAILABLE
SHOW_READ_CONTENT=False
READ_SPECS=tuple(sp.read_spec(r,'unassigned','Native whole-genome read/index acquisition not assigned by this scoped reconstruction','undetermined') for r in ('Read 1','Read 2','Index 1 (i7)','Index 2 (i5)'))

def sites(stage,context=False,include_controls=True):
 states={'initial':('C','5mC','5hmC'),'protected':('C','5mC','5gmC'),'oxidized':('C','5caC','5gmC'),'converted':('U','U','5gmC')}
 sc=Scene()
 for population in (('genomic DNA','independent spike-in control DNA') if include_controls else ('genomic DNA',)):
  for origin,state in zip(('C','5mC','5hmC'),states[stage]):
   parts=[role(population+' site originally'+origin+': physical'+state)]
   if context:parts=[role('native WGBS template context/endA unassigned'),*parts,role('native opposite template context/endB unassigned')]
   sc.strand(population+origin,parts,label=('physical converted single-strand view' if stage=='converted' else 'representative native DNA strand view; opposite strand not drawn')+'; original'+origin)
 return sc

def clean(g,key,template,label):
 reagent=g.input(key+'_reagent',label,MaterialKind.SELECTION)
 kept=g.product(key+'_retained','recovered '+template.label,MaterialKind.MOLECULAR,state=template.state)
 waste=g.product(key+'_waste','source physical removed fraction; composition unassigned',MaterialKind.WASTE)
 g.react(key,label,StepKind.CLEANUP,((template,InputRole.TEMPLATE),(reagent,InputRole.SELECTION)),(kept,waste),changes_molecule=False)
 g.discard(waste,'Source physical column/wash fractions; modification-specific recovery not assumed complete.')
 return kept

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 dna=g.input('dna','genomic DNA',MaterialKind.BIOLOGICAL,state=sites('initial',include_controls=False),primary=True)
 shear=g.input('shear','source sonication; native fragment-size distribution unassigned',MaterialKind.CONDITION)
 fragmented=g.product('fragmented','sonicated genomic DNA',MaterialKind.MOLECULAR,state=sites('initial',include_controls=False))
 g.react('sonicate','Sonicate genomic DNA before treatment',StepKind.FRAGMENTATION,((dna,InputRole.TEMPLATE),(shear,InputRole.CONDITION)),(fragmented,))
 controls=g.input('controls','source C/5mC/5hmC spike-in DNA; selected constructs unassigned',MaterialKind.MOLECULAR)
 bgt=g.input('beta_gt',f'wild-type βGT;2012 source1µM,37°C{GLUCOSYLATION_MIN}min',MaterialKind.CATALYST)
 udp=g.input('udp_glucose','UDP-glucose200µM in source HEPES/MgCl2 reaction',MaterialKind.CONDITION)
 protected=g.product('protected','glucosylated genomic/control DNA',MaterialKind.MOLECULAR,state=sites('protected'))
 g.react('glucosylate','Protect original5hmC as5gmC',StepKind.CHEMICAL_CONVERSION,((fragmented,InputRole.TEMPLATE),(controls,InputRole.SPIKE_IN),(bgt,InputRole.CATALYST),(udp,InputRole.CONDITION)),(protected,))
 protected=clean(g,'glucose_cleanup',protected,'2012 QIAquick Nucleotide Removal Kit')
 tet=g.input('mtet1',f'recombinant mTet1;2012 source3µM,37°C{OXIDATION_2012_MIN}min',MaterialKind.CATALYST)
 cofactors=g.input('tet_cofactors','HEPES,FeII,αKG,ascorbate,DTT,NaCl,ATP per2012 reaction',MaterialKind.CONDITION)
 oxidized=g.product('oxidized','protected DNA after5mC oxidation to5caC',MaterialKind.MOLECULAR,state=sites('oxidized'))
 g.react('oxidize','Oxidize5mC while protected5gmC retained',StepKind.CHEMICAL_CONVERSION,((protected,InputRole.TEMPLATE),(tet,InputRole.CATALYST),(cofactors,InputRole.CONDITION)),(oxidized,))
 protease=g.input('proteinase','source proteinaseK treatment; exact2012 conditions unassigned',MaterialKind.CATALYST)
 deprot=g.product('deproteinized','deproteinized oxidized DNA',MaterialKind.MOLECULAR,state=oxidized.state)
 g.react('proteolyse','ProteinaseK treatment after oxidation',StepKind.DEPROTEINIZATION,((oxidized,InputRole.TEMPLATE),(protease,InputRole.CATALYST)),(deprot,))
 purified=clean(g,'biospin',deprot,'MicroBioSpin30 columns')
 purified=clean(g,'qiagen',purified,'2012 QIAquick PCRPurification;2012 companion specifies NucleotideRemoval instead')
 native=g.input('native_wgbs_reagents','native whole-genome library preparation context;2012 companion points to Illumina15021861A',MaterialKind.CONDITION)
 prepared=g.product('native_template','native WGBS preconversion context; exact end/adapter order unassigned',MaterialKind.MOLECULAR,state=sites('oxidized',True))
 g.react('wgbs_context','Bounded native WGBS template preparation before conversion',StepKind.MAGIC,((purified,InputRole.TEMPLATE),(native,InputRole.CONDITION)),(prepared,),reason='Detailed2012 whole-genome Step25B delegates library construction to historicalIllumina15021861A. Its endprep/adapter/index order and native strands were not inspected; no generic ligation recipe is supplied. This node records unassigned preconversion template context, not proven adapter ligation.')
 bisulfite=g.input('bisulfite','TAB bisulfite conversion; whole-genome native implementation unassigned',MaterialKind.CONDITION)
 converted=g.product('converted','physicalU/U/5gmC converted template strands',MaterialKind.MOLECULAR,state=sites('converted',True))
 g.react('bisulfite_convert','Bisulfite-convert treated DNA: finished unprotected sitesU, protected5gmC retained',StepKind.CHEMICAL_CONVERSION,((prepared,InputRole.TEMPLATE),(bisulfite,InputRole.CONDITION)),(converted,))
 f=g.input('forward','native whole-genome library primerA; strand unassigned',MaterialKind.OLIGO)
 r=g.input('reverse','native whole-genome library primerB; index member/width unassigned',MaterialKind.OLIGO)
 pol=g.input('polymerase','native whole-genome amplification enzyme; not inferred from locus assay',MaterialKind.CATALYST)
 library=g.product('library',FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 originals=g.product('originals','modified original genomic/control strands and byproducts; recovery unassigned',MaterialKind.MOLECULAR,state=converted.state)
 g.react('library_pcr','Copy/amplify native converted templates',StepKind.PCR,((converted,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(library,originals),geometry=PCRGeometry.unavailable(converted.state,reason='Actual physical convertedU/5gmC native template strands; expected normalC/T daughter endpoint is not the PCR substrate.'))
 g.mark_terminal(Terminal(originals.id,Disposition.UNRESOLVED,'PCR readout doesnot establish destruction/removal of modified original strands or spike controls.'))
 g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,acquisition='unspecified',limitation=SEQUENCING_UNAVAILABLE))
 g.declare_reviewed_scope(REVIEWED,SOURCE)
 return g.require_publishable()

REVIEWED=('dna','shear','sonicate','controls','beta_gt','udp_glucose','glucosylate','glucose_cleanup_reagent','glucose_cleanup','mtet1','tet_cofactors','oxidize','proteinase','proteolyse','biospin_reagent','biospin','qiagen_reagent','qiagen','native_wgbs_reagents','wgbs_context','bisulfite','bisulfite_convert','forward','reverse','polymerase','library_pcr','library')
SPEC=SimpleNamespace(**{k:v for k,v in globals().copy().items() if k.isupper()},protocol_graph=protocol_graph,oligos=lambda:())

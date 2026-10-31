"""Source-qualified public preparation with partial molecular detail; no full source/absence audit."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SequencingPlan,SourceItem,SourceCoverage
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from chemdraw import Scene,Segment

def unknown_scene():
 s=Scene();s.strand('unassembled',[Segment('unknown_template','[unassembled]',placeholder=True,inferred=True)]);return s
class T:
 def __init__(self,g,x):self.g=g;self.x=x
 def inp(self,k,label,kind=K.CONDITION,role=R.CONDITION,**kw):return(self.g.input(k,label,kind,**kw),role)
 def op(self,k,label,ports=(),kind=S.MAGIC,chemical_state=None,**kw):
  y=self.g.product(k+'_product','Source preparation: '+label+'; whole occupied molecules/yield/sideproducts unassembled',K.MOLECULAR,chemical_state=chemical_state)
  if kind==S.MAGIC:kw.setdefault('reason','Positive source operation; detailed molecular products/occupancy unassembled, not source-absence evidence.')
  if kind==S.GENERIC:kw['changes_molecule']=False
  self.g.react(k,label,kind,((self.x,R.TEMPLATE),)+tuple(ports),(y,),**kw);self.x=y;return y
 def clean(self,k,label,ports=()):
  y=self.g.product(k+'_kept','Actual retained source workup fraction',K.MOLECULAR);w=self.g.product(k+'_excluded','Actual excluded source fraction; purity not asserted',K.WASTE)
  self.g.react(k,label,S.CLEANUP,((self.x,R.TEMPLATE),)+tuple(ports)+(self.inp(k+'_condition',label),),(y,w));self.g.discard(w,'Actual excluded workup fraction; no invented molecular partition or complete purification');self.x=y;return y
 def aliquot(self,k,label,qc=False):
  y=self.g.product(k+'_selected','Actual allocation: '+label,K.QC if qc else K.MOLECULAR);r=self.g.product(k+'_remainder','Untransferred source stock; actual amount/later fate unassigned',K.MOLECULAR)
  self.g.react(k,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r))
  self.g.mark_terminal(Terminal(y.id if qc else r.id,Disposition.QC if qc else Disposition.UNRESOLVED,label if qc else 'Untransferred stock not silently consumed/discarded'))
  self.x=r if qc else y;return y
 def library(self,k,label,platform,reads):
  y=self.g.product(k+'_library',label+'; whole library ends/primerlanding unassembled',K.LIBRARY)
  self.g.react(k+'_loading','Actual sequencing loading contribution',S.GENERIC,((self.x,R.TEMPLATE),),(y,),changes_molecule=False)
  self.g.sequence(y,SequencingPlan(read_lengths=reads,limitation=platform+'; source run lengths only; whole adapter/index/primerlanding geometry unassembled.'));return y

AUTHORITY='Published OUP/NAR2024 PDF Methods77–122/Resultsmethodologyp4/actualFig1p5; associatedSIrolemap p1 andTableS8rawXML validationroles. AMXMLidentityheaderonly excluded; conservative4identity exposure/capexception recorded. No vendor/predecessor recipe inheritance.'

def _finish(g,x,prefix,ng):
 b=T(g,x);b.aliquot(prefix+'_rna_input',str(ng)+'ng RNA selected for separate Collibri library')
 stock=g.input(prefix+'_kit_inventory','Actual Collibri Stranded RNA library preparation aggregate supplies; number/packaging/addition schedule unassembled, not a proven premix',K.MOLECULAR)
 ol=g.product(prefix+'_kit_oligos','Actual kit oligo contributors; bodies/modifications/functional partition unassembled',K.OLIGO)
 rt=g.product(prefix+'_kit_rt','Actual reverse-transcription catalyst class from RNA-to-cDNA preparation supplies; identity/dose unassigned',K.CATALYST)
 nt=g.product(prefix+'_kit_nt','Actual polymerization nucleotide supply class; identities/doses/sugar mixtures unassembled',K.MOLECULAR)
 rest=g.product(prefix+'_kit_other','Other actual kit molecular supplies; identities/packaging unassembled',K.MOLECULAR)
 g.react(prefix+'_kit_components','Bookkeeping inventory view of actual kit supplies; not fresh independent doses',S.GENERIC,((stock,R.TEMPLATE),),(ol,rt,nt,rest),changes_molecule=False)
 b.op(prefix+'_library_finish','Source Collibri stranded cDNA library preparation; whole sequence/internal order unassembled',((ol,R.STRUCTURAL),(rt,R.CATALYST),(nt,R.STRUCTURAL),(rest,R.STRUCTURAL)))
 b.aliquot(prefix+'_library_qc','Actual TapeStation4200 DNA ScreenTape QC allocation; amount unassigned',True)
 b.aliquot(prefix+'_loading_allocation','Actual source loading allocation; amount unassigned')
 b.library(prefix,'Separate '+prefix+' library','Illumina NextSeq550',{'Read1':150})

def protocol_graph():
 g=ProtocolGraph('CLiPPR-seq selected betaTC6 capture/RNaseR and paired CLiPP control',authority=AUTHORITY,representation='published_instance')
 cells=g.input('cells','Actual betaTC6 cells; source cell number/culture allocation unassigned',K.BIOLOGICAL,primary=True);b=T(g,cells)
 g.add_source_item(SourceItem('terminology','Abstract/Methodsheading says proximity ligation; actualFig1p5 and overview describe AMT365nm RNA-duplex crosslinking, polyA capture,254nm photoreversal.','No phosphodiester RNA-RNA partner join or physically joined partner-read inferred; source BLAST comparison and validation are computational/assay stages outside molecular production.'))
 g.add_source_item(SourceItem('wavelength','Methodsprints254nM; actualFigure1p5/Results says254nm UV.','Model adopts originalFigure1 wavelength interpretation, preserves literalMethods notation; no publication-error assertion or normalized source text.'))
 g.add_source_item(SourceItem('bounded_exception','3selected chemistry originals plus distinctAMXMLidentityheaderonly exposure =4 conservative identities.','AMchemistryexcluded; capoverrunrecorded, no more investigation. This partial source-positive model is not a cap-compliance or absence claim.'))
 amt=g.input('amt','Actual prepared AMT20ugml irradiation medium;1mgmlstock diluted inDMEM/2xPBS',K.MOLECULAR,chemical_state={'class':'PsoralenAMT duplex crosslinking reagent','RNA_relation':'Source365nm induced RNA-duplex crosslink class, not phosphodiesterligation; atommap/occupancy unassembled'})
 b.op('crosslink','37C30min AMT incubation then two7min365nm exposures onice5cm',((amt,R.STRUCTURAL),b.inp('uv_conditions','Actual365nm programme')),kind=S.CROSS_LINKING,chemical_state={'RNA':'AMT-crosslinked duplex-associated and unreacted RNA classes; not perfect capture of all direct interactions'})
 b.clean('rna_isolation','Selected conventional TRIzol isolation; source aqueous RNA precipitation, full workup recipe unassembled',ports=(b.inp('trizol_supplies','Actual TRIzol/precipitation molecular supplies, packaging/doses unassembled',K.MOLECULAR,R.STRUCTURAL),))
 b.op('dnase','Actual source DNase treatment on precipitation-isolated RNA',(b.inp('dnase_enzyme','Actual source DNase; dose/subtype unassigned',K.CATALYST,R.CATALYST),),kind=S.DIGESTION)
 b.aliquot('total_rna_qc','Actual RNA quality/quantity assay allocations; QubitRNAHS and spectrophotometry/TapeStation alternatives not normalized to fixed doses',True)
 b.aliquot('capture_rna_input','Actual10ug total RNA allocation for NEBNextpolyA magnetic isolation')
 module=g.input('polya_module','Actual NEBNextpolyA magnetic isolation aggregate supplies; stock packaging/doses unassembled',K.MOLECULAR)
 cap=g.product('polya_support','Actual source oligo-dT-bearing magnetic capture support; occupancy/oligo length unassembled',K.SELECTION,chemical_state={'capture_oligo':'Source oligo-dT/polyA hybridization class; whole oligo length/modifications unassembled'})
 ol=g.product('polya_oligo','Same support-associated oligo-dT molecular inventory; not a second added free oligo',K.OLIGO,chemical_state={'attachment':'Copresent on capture support; whole body/length unassembled'})
 work=g.product('polya_other','Other actual isolation module supplies, detailed internal workup unassembled',K.MOLECULAR)
 g.react('capture_inventory','Actual module component inventory, not independent fresh additions',S.GENERIC,((module,R.TEMPLATE),),(cap,ol,work),changes_molecule=False)
 selected=g.product('captured_rna','Actual source polyA/crosslinked-associated RNA capture preparation',K.MOLECULAR);unbound=g.product('capture_unbound','Actual unbound capture fraction',K.WASTE)
 g.react('polya_capture','Source10ugRNA oligo-dT magnetic isolation; exact whole duplex/end geometry unassembled',S.AFFINITY_CAPTURE,((b.x,R.TEMPLATE),(cap,R.SELECTION),(ol,R.STRUCTURAL),(work,R.STRUCTURAL)),(selected,unbound));g.discard(unbound,'Actual nonselected fraction, no perfect polyA or partner-specificity claim');b.x=selected
 elu=g.product('polya_eluate','Actual captured RNA eluate; volume unassigned for NGS branch',K.MOLECULAR);res=g.product('polya_residual','Residual support/material after actual isolation',K.MOLECULAR)
 g.react('polya_elution','Source recovered mRNA/associatedRNA for QC and photoreversal; exact kit operations unassembled',S.MAGIC,((b.x,R.TEMPLATE),),(elu,res),reason='Positive RNA recovery; validation20ulelution/12ulRT aliquot not imported into production.');g.mark_terminal(Terminal(res.id,Disposition.UNRESOLVED,'Residual support fate unassigned'));b.x=elu
 b.aliquot('polya_qc','Actual capturedRNA quantity/quality assay allocations; amounts unassigned',True)
 b.op('photoreverse','Actual5min254nmUV photoreversal; source Methods literal254nM preserved',ports=(b.inp('photoreversal_conditions','SourceFigure1/Results254nmUV5min'),),chemical_state={'RNA':'Intended photoreversed crosslinked RNA classes; damage/unreversed species/yield unassembled'})
 treated=g.product('clippr_half','Actual one-half de-crosslinked RNA allocation for RNaseR',K.MOLECULAR);control=g.product('clipp_half','Actual other-half untreated CLiPP control',K.MOLECULAR)
 g.react('branch_halves','Actual source half/half physical split before RNaseR',S.SPLIT,((b.x,R.TEMPLATE),),(treated,control))
 b.x=treated;b.op('rnase_r','Source RNaseR enrichment; detailed citedref20 conditions unreviewed',(b.inp('rnase_r_enzyme','Actual RNaseR exonuclease class; dose unassigned',K.CATALYST,R.CATALYST),),kind=S.DIGESTION,chemical_state={'RNA':'CircularRNA-enriched heterogeneous preparation with copresent degradation products; no instantaneous separation or perfect specificity'})
 _finish(g,b.x,'clippr',10);_finish(g,control,'clipp',25)
 return _publication_coverage(g)

def selftest():
 from checks import Check,Source,have
 from storage import sources_root
 c=Check();g=protocol_graph()
 c('AMT actual molecular contributor',g.materials['amt'].kind,K.MOLECULAR)
 c('No phosphodiester partnerligation',any(s.kind==S.LIGATION for s in g.steps.values()),False)
 c('DNase applies to isolated RNA',g.steps['dnase'].inputs[0].material,'rna_isolation_kept')
 c('Capture uses actual10ug allocation',g.steps['polya_capture'].inputs[0].material,'capture_rna_input_selected')
 c('Attached oligo is produced not fresh', 'polya_oligo' in g.steps['capture_inventory'].outputs,True)
 c('Photoreversal follows captureQC stock',g.steps['photoreverse'].inputs[0].material,'polya_qc_remainder')
 c('RNaseR only treated half',g.steps['rnase_r'].inputs[0].material,'clippr_half')
 c('Control independently takes untreated half',g.steps['clipp_rna_input'].inputs[0].material,'clipp_half')
 c('Two independent NGS outputs',len(g.sequencing),2)
 c('No physical pool',any(s.kind==S.POOL for s in g.steps.values()),False)
 c('Support remainder explicit','polya_residual' in g.terminals,True)
 c('No validation ASO in production',any('aso' in k.lower()for k in g.materials),False)
 for p in ['clippr','clipp']:
  c(p+' library consumes produced oligos',any(i.material==p+'_kit_oligos'for i in g.steps[p+'_library_finish'].inputs),True)
  c(p+' QC distinct',g.terminals[p+'_library_qc_selected'].disposition,Disposition.QC)
  c(p+' loading remainder explicit',p+'_loading_allocation_remainder'in g.terminals,True)
 p=sources_root()/HERE.parent.name/'nar_gkae058_PMC11014417.pdf.txt'
 if have(c,Source(p,'https://doi.org/10.1093/nar/gkae058')):
  t=p.read_text();c('Printed published identity','Nucleic Acids Research, 2024'in t,True);c('Source single-end150','150 bp single-end sequencing'in t,True)
 c('Publication graph validates with declared limitations',g.publication_problems(),[]);c.report()


def _cover_reviewed_targets(g, groups):
 """Only the literal curator-reviewed targets below receive these source assignments."""
 for source,targets in groups:
  for target in targets:
   label=g.materials[target].label if target in g.materials else g.steps[target].label
   key='reviewed_'+target
   g.add_source_item(SourceItem(key,label,source))
   g.cover_source_item(SourceCoverage(key,'modelled',target))

PAGE_FILE='clippr-seq.html'
SOURCE_TARGETS=(
 ('Printed published OUP/NAR2024 Methods77–108: selected betaTC6 AMT/TRIzol/DNase and10ugpolyA isolation/QC; inventory views are bookkeeping, not separate fresh doses.',
  ('cells','amt','uv_conditions','crosslink','trizol_supplies','rna_isolation_condition','rna_isolation','dnase_enzyme','dnase','total_rna_qc','capture_rna_input','polya_module','capture_inventory','polya_capture','polya_elution','polya_qc')),
 ('Same published Methods109–114 and Results methodologyPDFp4/Fig1p5: photoreversal, physicalhalves andRNaseR enrichment; detailed ref20 conditions not inherited.',
  ('photoreversal_conditions','photoreverse','branch_halves','rnase_r_enzyme','rnase_r')),
 ('Same published Methods115–122/Fig1:10ngtreated and25nguntreatedCollibrilibraries/DNAScreenTapeQC/NextSeq550SE150. Componentviews are class-level inventories; wholekitrecipe/bodies unassembled.',
  ('clippr_rna_input','clippr_kit_inventory','clippr_kit_components','clippr_library_finish','clippr_library_qc','clippr_loading_allocation','clippr_loading','clippr_library',
   'clipp_rna_input','clipp_kit_inventory','clipp_kit_components','clipp_library_finish','clipp_library_qc','clipp_loading_allocation','clipp_loading','clipp_library')),
)

def _publication_coverage(g):
 g.cover_source_item(SourceCoverage('terminology','modelled','crosslink'))
 g.cover_source_item(SourceCoverage('wavelength','modelled','photoreverse'))
 g.cover_source_item(SourceCoverage('bounded_exception','excluded',reason='Historical source-identity/cap record, not a molecular operation; recorded in the source review and publicly qualified source scope.'))
 _cover_reviewed_targets(g,SOURCE_TARGETS)
 return g.require_publishable()

def render_page():
 from page import head,caveat,table
 from protocol_views import assemble_page,WorkflowSection
 g=protocol_graph()
 return assemble_page(g,[head('CLiPPR-seq — capture of mRNA-associated circular RNA'),
 '<main class="wrap" data-protocol-complete="true"><h1>CLiPPR-seq</h1><p>AMT crosslinking preserves RNA-duplex associations before poly(A) capture. Photoreversal releases the associated RNAs; RNase R enriches circular RNA for sequencing. A separate untreated CLiPP library supplies the mRNA comparison.</p><p><a href="https://doi.org/10.1093/nar/gkae058">Singh et al., Nucleic Acids Research (2024)</a> · <a href="01_clippr-seq.html">Research note and evidence scope</a></p>',
 caveat('<strong>Source-qualified preparation, with partial molecular detail.</strong> The published PDF supports the betaTC6 workflow below. RNA sequences, occupied crosslinks, capture-oligo length and complete Collibri library molecules remain unassembled. The method enriches candidate associations; it does not establish exclusive direct binding or a covalently joined partner read.'),
 '<h2>Preparation</h2><ol><li>Incubate cells with AMT at 20 µg/ml for 30 min at 37°C, then give two 7 min 365-nm UV exposures on ice.</li><li>Isolate RNA using the selected TRIzol route, treat with DNase and allocate 10 µg for oligo-dT magnetic capture.</li><li>Photoreverse for 5 min and split the recovered RNA into halves. Treat one half with RNase R; retain the other as CLiPP control.</li><li>Prepare separate Collibri stranded libraries from 10 ng treated RNA and 25 ng untreated RNA, then allocate DNA for TapeStation QC and sequencing.</li></ol>',
 '<h2>Known molecular parts</h2>',table(('Part','Role and limit'),(('AMT plus 365-nm UV','RNA-duplex crosslinking class; occupied sites and atom maps unassembled.'),('Oligo-dT magnetic support','Captures poly(A) RNA and associated crosslinked RNA. Whole support oligo and hybrid occupancy remain unresolved.'),('RNase R','Circular-RNA enrichment; linear fragments and incompletely digested material are not assumed absent.'),('Collibri stranded RNA preparation','Actual kit oligo, catalyst and nucleotide inventory classes are retained; exact bodies, packaging and finishing order are unresolved.'))),
 '<h2>Read setup</h2><p>Illumina NextSeq 550, single-end 150 bases for each library. No physical pooling is asserted. Whole read-primer landing sites and index layout remain unassembled; candidate partners are inferred by comparing circRNA sequences with the CLiPP mRNA library and by validation.</p>',
 '<h2>Source wording and scope</h2><p>The abstract and Methods heading say “proximity ligation”; Figure 1 and the executed overview describe AMT crosslinking and photoreversal. This page models that physical capture without inventing a phosphodiester RNA–RNA join. Methods prints “254 nM”; Figure 1 and Results identify 254 nm UV. The page uses the figure wavelength interpretation and preserves both source values.</p><p>The PMC-hosted PDF visibly uses the published OUP/NAR2024 layout. Associated supplement index and TableS8 were checked for reagent roles; their validation, knockdown and reporter oligos are not assigned as sequencing adapters. The associated author manuscript is not used as a substitute for the inspected published chemistry. Detailed RNase R and kit originals remain outside this source pass.</p>',
 '<details><summary>Full molecular workflow, physical fractions and unresolved chemistry</summary>',WorkflowSection('Source-qualified molecular workflow'),'</details></main>'])

if __name__=='__main__':selftest()

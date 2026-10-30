from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from dataclasses import dataclass,replace
from chemdraw import Scene,Segment,Construct,feature,workflow_from_protocol_graph,workflow_panel
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from page import head,info
PAGE_FILE='chem.html'

class Operations:
    def __init__(self,g):self.g=g
    def op(self,key,x,label,reagent=None,state=None,chemical_state=None):
        g=self.g;inputs=[(x,R.TEMPLATE)]
        if reagent:
            y=g.input(key+'_condition',reagent,K.CONDITION);inputs.append((y,R.CONDITION))
        y=g.product(key+'_out',label+'; preparation includes unremoved companion species',K.MOLECULAR,state=state,chemical_state=chemical_state)
        g.react(key,label,S.GENERIC,inputs,(y,));return y
    def join(self,key,inputs,label,state=None,chemical_state=None):
        g=self.g;y=g.product(key+'_out',label+'; reaction components remain until an explicit separation',K.MOLECULAR,state=state,chemical_state=chemical_state)
        g.react(key,label,S.GENERIC,tuple((x,R.TEMPLATE if i==0 else R.STRUCTURAL) for i,x in enumerate(inputs)),(y,));return y
    def clean(self,key,x,label):
        g=self.g;z=g.input(key+'_selection',label,K.SELECTION)
        y=g.product(key+'_kept','Source physically retained fraction: '+label,K.MOLECULAR)
        w=g.product(key+'_removed','Source physically removed fractions; exact composition/purity unassigned',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(z,R.SELECTION)),(y,w));g.discard(w,'Only the source physical separation; no perfect selectivity asserted.');return y
    def split(self,key,x,label,other_label,other=Disposition.UNRESOLVED):
        g=self.g;y=g.product(key+'_production',label,K.MOLECULAR)
        z=g.product(key+'_other',other_label,K.MOLECULAR)
        g.react(key,'Source physical allocation: '+label,S.SPLIT,((x,R.TEMPLATE),),(y,z));g.terminals[z.id]=Terminal(z.id,other,other_label);return y,z
    def qc(self,key,x,label):
        g=self.g;y=g.product(key+'_production','Production preparation retained independently from analytical aliquot',K.MOLECULAR)
        z=g.product(key+'_assay',label,K.QC)
        g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,z));g.terminals[z.id]=Terminal(z.id,Disposition.QC,'Analytical aliquot does not become the sequencing substrate.');return y
    def opaque(self,key,x,label):
        g=self.g;y=g.product(key+'_out',label,K.MOLECULAR)
        g.react(key,label,S.MAGIC,((x,R.TEMPLATE),),(y,),reason='Selected originals establish this operational continuation; full molecular inventory/geometry remains unassembled.');return y

def finish(g,x,readout):
    g.sequence(x,SequencingPlan(acquisition='unspecified',limitation=readout))
    g.declare_reviewed_scope(tuple(g.materials)+tuple(g.steps),SOURCE)
    g.add_source_item(SourceItem('unassembled','Unassembled source inputs, molecular geometry and directly linked targets',LIMIT));return g

def render_page():
    return head(TITLE)+'<main class="wrap"><h1>'+TITLE+'</h1>'+info(LIMIT)+'<p><a href="01_chem.html">Research notes</a></p>'+workflow_panel(workflow_from_protocol_graph(protocol_graph(),ENDPOINT))+'</main>'

TITLE='PLATO — ribosome display to recovered-RNA sequencing'
SOURCE='Published Zhu et al.2013 nbt.2539 first-party main; published MOESM14 OnlineMethods Ribosome display/ORFeome precipitation/mRNA preparation/Analysis'
LIMIT='Private incomplete operational graph. Noncovalent ribosome/mRNA/protein association is separate from subsequently copied DNA. Exact Table4 primers, full product geometry, detailed linked ribosome-display backbone and kit revisions remain unreviewed. Source polyadenylates cDNA,not the original RNA. No ideal binding or recovery predicted.'
ENDPOINT='pool_out'
@dataclass(frozen=True)
class SampleIndexInventory:
    width:int=7
    semantic:object=feature('plato.sample','sample_index','unknown',note='Source7nt samplebarcode;actualmember/position unassembled')
@dataclass(frozen=True)
class DisplayPopulation:
    encoding_polymer:str='RNA'
    linkage:str='noncovalent ribosome association'
    released:bool=False
    def disrupt(self):
        if self.released:raise ValueError('Already released population')
        return replace(self,released=True)

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE);o=Operations(g)
    x=g.input('entry_pool','Source normalized humanORFeomev5.1 entry clone super-pools; full backbone/insert inventory unassembled',K.MOLECULAR,primary=True)
    x=o.op('lr',x,'Source Gateway LR into pRD-DEST; recovered source constructs/physical bacterial workup unassembled','Gateway LR system; exact backbone version unreviewed')
    t7=g.input('t7b','Source T7B primer,Table4 available but uninspected',K.OLIGO);tol=g.input('tolak','Source TolAk primer,Table4 uninspected',K.OLIGO)
    x=o.join('template_pcr',(x,t7,tol),'Source PrimeSTAR10cycles amplifies linear stop-codon-free RD templates from18super-pools')
    x=o.clean('dna_column',x,'QIAquickPCRpurification of each source super-pool')
    x=o.op('superpool_mix',x,'Equal DNA weight source18super-pool combination; individual constituents/amounts unassembled')
    x=o.op('ivt',x,'T7Ribomax creates new RNA transcripts; original DNA/reagents retained','T7RibomaxLarge transcription kit')
    x=o.clean('rna_column',x,'RNeasyRNAproduct purification')
    x=o.op('translation',x,'RTS100E.coliHY30°C15min creates source ribosome-displayed ORFs; genotype remains RNA','7.5µgRNA/50µL withRNAseOUT',chemical_state=DisplayPopulation())
    x,rest=o.split('translation_allocate',x,'Source12.5µLtranslation aliquot into ice-coldRDselectionbuffer','Other translation preparation aliquots; later allocation unassembled')
    x=o.clean('clarify',x,'14000g5min4°C; move source supernatant to new ice-coldtube; sediment removed')
    bait=g.input('bait','Prepared blocked bait-bearing beads; source GST-LYN/GSTcontrols selected representative route, bait manufacture unassembled',K.MOLECULAR)
    x=o.join('bind',(x,bait),'Source100µLdisplayedRNA preparation binds blocked beads4°C6h')
    x=o.clean('wash',x,'Six500µLicecoldRDwashbuffer washes; retain source beads with associated population')
    x=o.op('disrupt',x,'EB20/20mMEDTA37°C10min disrupts ribosomal complexes and releases encoding RNA; beads/protein/RNA not automatically separated','50µLEB20+RNAseOUT',chemical_state=DisplayPopulation().disrupt())
    x=o.clean('recover_rna',x,'Source RNeasy purification of elutedRNA; exact bead-transfer geometry unassembled')
    x=o.op('fragment',x,'NEBNextmagnesium fragmentation94°C3min; RNA fragments are actual substrates','NEBNextMagnesiumRNAFragmentationModule')
    x=o.clean('spin50',x,'Spin-50mini-column cleanup of fragmentedRNA')
    rt=g.input('tolart','Source TolA100RT primer,Table4 uninspected',K.OLIGO)
    x=o.join('rt',(x,rt),'SuperScriptIII makes new cDNA from retained ORF3prime-fragment class; originalRNA and unusedprimer remain')
    x=o.op('exoi',x,'ExonucleaseI37°C30min digests susceptible free ssDNA primer; not physical removal or guaranteed universal cDNA resistance','ExonucleaseI')
    x=o.op('rnaseh',x,'RNaseH37°C30min digests RNA in RNA–DNA hybrids; digestion remnants remain until column','RNaseH')
    x=o.clean('cdna_column',x,'QIAquickPCRpurification retains source cDNA preparation')
    x=o.op('tail',x,'TdT polyadenylates cDNA3prime ends; no RNA-tail or UMI inference','InvitrogenTdTreactionkit')
    x=o.clean('tail_column',x,'QIAquickPCRpurification afterTdT')
    dt=g.input('adaptor_dt','Source Adaptor-(dT)24 primer; wholeTable4body uninspected',K.OLIGO)
    x=o.join('first_extension',(x,dt),'Source single95°C2min/50°C1min/72°C7min extension before forward primer addition')
    f=g.input('p5_prdrev','Source P5-PRDREV primer; wholeTable4body uninspected',K.OLIGO)
    x=o.join('pcr1',(x,f),'Add forwardprimer/additional25µLmixture;30cyclesHerculaseII creates new dsDNA')
    x=o.clean('pcr1_column',x,'QIAquickPCRpurification of firstPCR')
    index=g.input('index','Source INDEX primer;7nt sample barcode,not moleculeUMI; exact member/Table4body uninspected',K.OLIGO,chemical_state=SampleIndexInventory())
    x=o.join('pcr2',(x,index),'Source100ngfirstPCR+P5-PRDREV/INDEX,HerculaseII10cycles')
    x=o.clean('pcr2_column',x,'QIAquickPCRpurification of secondPCR')
    x=o.qc('quantify',x,'Separate source7500Fast qPCR analytical aliquot,P5/P7_2; productionDNA retained')
    others=g.input('other_libraries','Other source independently prepared/quantified secondPCR libraries; actual batch identities unassembled',K.MOLECULAR)
    x=o.join('pool',(x,others),'Equal source secondPCR product amounts combined; samplebarcode,not moleculardeduplication')
    return finish(g,x,'PublishedHiSeq2000,multiplexSE50,customPRDREV-attB2-SP; exactTable4primer/fullstrandlanding not assembled.7ntsamplebarcodes assigned computationally,not additional chemical steps.')

def selftest():
    from checks import Check
    c=Check();g=protocol_graph();d=DisplayPopulation()
    c('Disruption releases RNA without relabeling genotypeDNA',(d.disrupt().encoding_polymer,d.disrupt().released),('RNA',True))
    c('Original population retained as immutable source class',d.released,False)
    c('Translation aliquot remainder has explicit unresolved fate',g.terminals['translation_allocate_other'].disposition,Disposition.UNRESOLVED)
    c('EDTA release precedes physicalRNApurification',g.steps['recover_rna'].inputs[0].material,'disrupt_out')
    c('RNaseH receives unremovedExoI reaction',g.steps['rnaseh'].inputs[0].material,'exoi_out')
    c('TdT follows physical cDNA purification',g.steps['tail'].inputs[0].material,'cdna_column_kept')
    c('Forward primer added after firstsingleextension',g.steps['pcr1'].inputs[0].material,'first_extension_out')
    c('Sequencing pool excludes analytical qPCRaliquot',g.steps['pool'].inputs[0].material,'quantify_production')
    c('Fullsourcegeometry incomplete',g.publishable,False);c.report()

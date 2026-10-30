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

from illumina import P5,P7
TITLE='CFMS mRNA display — single-round fusion selection and NGS'
SOURCE='Published2012main+SupportingInformation combined institutionalPDF source_pdf_3f348b50ea6f42f28215.pdf; DetailedMethods950–1054,TableS1 renderedPDFp21(S15),Fig1/FigS1'
LIMIT='Private incomplete source graph. Published combinedPDF contains both final article and SI,despite SI-only source-index description. Puromycin/mRNA/protein linkage is distinct from newly copied cDNA. Exact pF30P/splint chemistry,cassette end geometry/full display-library strands and prior2011procedure remain unassembled; no ideal purity,binding or uniqueclone count predicted. SI Sequencing reportsPE42 but its read1expected pattern totals44bases (A+19random+24scaffold); model retainsPE42 and leaves this printed pattern discrepancy unresolved.'
ENDPOINT='selected_pcr_out'
PRIMERS={'FnRTDS':P5+'TGTGCGGTAGTTGATGGAGATCG','Oligo C':P5[:21],'Oligo D':'G'+P7[:22], 'FNBCSEQ':'CCAGTTTCACCGTAGGTGATGCGGTAGTAGCGA','FNFGSEQ':'TGTCGACTATACCATCACGGTGTACGCCGTCAC'}
def primer(name):return Construct([Segment(name,PRIMERS[name])])
@dataclass(frozen=True)
class FusionClass:
    encoding_polymer:str='RNA'
    peptide_link:str='puromycin-mediated covalent peptide linkage through source linker'
    cdna:bool=False
    def reverse_transcribe(self):return replace(self,cdna=True)

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE);o=Operations(g)
    bc=g.input('bc_oligos','Source FnOligo3C25/FnOligo2DS/FnOligo4 thenFnOligo1DS/FnOligo5 BCcassette reagents; assembly unassembled',K.OLIGO,primary=True)
    fg=g.input('fg_oligos','Source FnOligo7/FnOligo6/FnOligo9 FGcassette reagents; fullassembly unassembled',K.OLIGO,primary=True)
    bc=o.op('bc_pcr',bc,'KODsourceBCcassette stagedextension/amplification25fold then8fold; notguaranteedcloneuniformity','KODpolymerase system')
    fg=o.op('fg_pcr',fg,'KODsourceFGcassette extension/amplification20fold','KODpolymerase system')
    bc=o.clean('bc_cassette_column',bc,'Source spincolumnpurification ofBCcassette beforeBsaI')
    fg=o.clean('cassette_column',fg,'Source spincolumnpurification beforeBsaI; individualcassette allocation unassembled')
    bc=o.op('bc_digest',bc,'BsaI digestion of sourceBCcassette; exact cohesive ends unassembled','BsaI')
    fg=o.op('fg_digest',fg,'BsaI digestion of sourceFGcassette; exact cohesive ends unassembled','BsaI')
    bc=o.clean('bc_digest_column',bc,'Spincolumnpurification ofBsaIdigestedBCcassette')
    fg=o.clean('fg_digest_column',fg,'Spincolumnpurification ofBsaIdigestedFGcassette')
    x=o.join('cassette_join',(bc,fg),'SourceT4DNAligase60µLjoinsBsaI-digestedcassettes; exactjunctions/endinventory unassembled')
    x=o.clean('library_gel',x,'Source2percentagarose extraction andQiagenspincolumnpurification')
    x=o.qc('initial_measure',x,'Separate source spectrophotometric concentration measurement; no invented success/uniqueness guarantee')
    x,remaining=o.split('library_allocate',x,'Serially diluted2.5e9molecule source aliquot for thisselection experiment','Other libraryDNA aliquots/fates unassigned')
    inputdna=g.product('input_dna_ngs','Separate source amplified inputDNA library sequenced for composition; detailed allocation/amplification unassembled',K.MOLECULAR)
    main=g.product('fusion_dna','Source productionDNA allocated to fusion preparation',K.MOLECULAR)
    g.react('input_reference','Separate source inputDNA composition-reference library',S.SPLIT,((x,R.TEMPLATE),),(main,inputdna))
    g.sequence(inputdna,SequencingPlan(acquisition='unspecified',limitation='PublishedGAIIxPE42/customBCandFGprimers,inputDNAreference; exactwholelibrary/clustergeometry unassembled.'))
    x=o.op('scale_pcr',main,'FnOligo1DS/FnOligo9 source19cycleKODPCR in5mL','Full primers remain in sourceTableS1; product unassembled')
    x=o.clean('dna_extract',x,'Phenol/chloroform extraction andethanolprecipitation of sourceDNA')
    x=o.op('ivt',x,'Entire sourceDNA used for5mLT7runofftranscription; originalDNA/reagents remain','AmbionT7RNApolymerase')
    x=o.clean('rna_extract',x,'Phenol/chloroform extraction andethanolprecipitation of sourceRNA')
    linker=g.input('pf30p','SourcepF30P puromycinlinker; wholepolymer/modifications unreviewed in directlycited predecessor',K.OLIGO)
    splint=g.input('splint','Source splintoligo forT4DNAligase association; wholepolymerunreviewed',K.OLIGO)
    x=o.join('linker_ligation',(x,linker,splint),'Source33.6nmolRNA/40nmolpF30P/37nmolsplint,T4DNAligase3mL; ligationendgeometry unassembled')
    x=o.clean('rna_linker_gel',x,'4.5percentdenaturingureaPAGE/electroelution thenethanolprecipitation retains sourceRNA–puromycintemplate')
    x,unused=o.split('template_allocate',x,'Source800pmolRNA–puromycintemplate for twotarget preparation','Other RNA–puromycin preparation sufficient for otherexperiments; laterallocationunassigned')
    x=o.op('translation',x,'Rabbitreticulocytelysate30°C1htranslation','2mLtranslationreaction',chemical_state=FusionClass())
    x=o.op('fusion',x,'KCl500mM/MgCl2 60mM thenroomtemperature30minpromotespuromycin-mediatedfusion; unreactedclasses remain','KCl/MgCl2',chemical_state=FusionClass())
    x=o.clean('dt_capture',x,'OligodTcellulose1hbinding then2bindingbuffer+2icecoldTBSwashes; retainassociatedpopulation')
    x=o.op('dt_release',x,'Three500µL5mMTrispH8roomtemperature elutions; matrixnotautomaticallysequenced')
    x=o.clean('dt_collect',x,'Source collect sourceoligodT-cellulose eluates apartfrom retainedmatrix')
    rt=g.input('fnrtds','SourceFnRTDS reverse-transcriptionprimer',K.OLIGO,state=primer('FnRTDS'))
    x=o.join('rt',(x,rt),'SuperScriptII/FnRTDS2µM creates DNAcopy; covalentpeptide–RNA linkage andoriginalRNA remain',chemical_state=FusionClass().reverse_transcribe())
    x=o.op('edta',x,'AddEDTA afterRT; no inventedRNAremoval','EDTA')
    x=o.clean('flag_capture',x,'M2resin1h4°Ccapture then4TBS/0.05percentTweenwashes; sourceintendeddepletionnotperfectfunctionalfilter')
    x=o.op('flag_release',x,'800µL0.15mg/mL3XFLAGpeptide competitivelyelutesfusions')
    x=o.clean('flag_collect',x,'Collect sourceFLAGeluate separatelyfromM2resin')
    x=o.op('bsa',x,'AddBSA1mg/mL prioraffinityenrichment','BSA')
    enriched,ref=o.split('fusion_reference',x,'Source99percentfusion preparation continuestoaffinityenrichment','Source1percentpurifiedfusionpool reservedforsequencing',Disposition.UNRESOLVED)
    del g.terminals[ref.id]
    g.sequence(ref,SequencingPlan(acquisition='unspecified',limitation='Publishedpurifiedinputfusion-cDNA compositionreference sequencedGAIIxPE42; exactreferencePCR/allocation unassembled.'))
    x=enriched
    for n in range(1,5):x=o.clean('preclear'+str(n),x,'Source30minneutravidin-bead preclear '+str(n)+'; retainedunboundfusionfraction,beadsphysicallyremoved')
    x,other=o.split('target_allocate',x,'One equal sourcepreclearedfraction forthe selected targetbead experiment','Other equal sourcefraction forthe independent secondtarget experiment; no forceddiscard')
    bait=g.input('target_beads','SourcebiotinylatedMBPorIgG(Fc) onneutravidinM270beads; selectedtargetidentity notforced; upstreamtargetmodification unassembled',K.MOLECULAR)
    x=o.join('target_bind',(x,bait),'Source1hroomtemperature targetbinding')
    x=o.clean('cfms',x,'CFMS30mL/hwithdrawal trapsbeads,3.5minbufferwash,evacuatebuffer; retaintrappedbead-associatedpopulation')
    x=o.op('bead_resuspend',x,'Remove magnets; resuspendbeads100µLPCRmixture; beadsstayintubereaction')
    x=o.op('release',x,'95°C2min dissociatesencodingcDNA; no physicalbeadremovalasserted')
    c=g.input('oligo_c','SourceOligoC',K.OLIGO,state=primer('Oligo C'));d=g.input('oligo_d','SourceOligoD',K.OLIGO,state=primer('Oligo D'))
    x=o.join('selected_pcr',(x,c,d),'KODsource17cyclesPCR makesnewDNAcopies; originalfusion/beadmixtureremainsuntilunreviewedloadingworkup')
    return finish(g,x,'PublishedGAIIx6pM,PE42; FNBCSEQthenFNFGSEQ sourcecustomprimers. RandomBC/FGbindingloops are proteinvariants,notUMIs. ExactpostPCRloadingworkup/wholelibrarylanding unassembled.')

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();f=FusionClass()
    c('RT creates cDNA without deletingencodingRNA',(f.reverse_transcribe().encoding_polymer,f.reverse_transcribe().cdna),('RNA',True))
    c('Translation startswithphysicallypurifiedRNA–puromycinaliquot',g.steps['translation'].inputs[0].material,'template_allocate_production')
    c('RT followscollectedoligodTeluate',g.steps['rt'].inputs[0].material,'dt_collect_kept')
    c('Affinityinputexcludesreservedreference',g.steps['preclear1'].inputs[0].material,'fusion_reference_production')
    c('Secondtargetfractionnotdiscarded',g.terminals['target_allocate_other'].disposition,Disposition.UNRESOLVED)
    c('PCRtemplateheatstillcontainsbeads',g.steps['selected_pcr'].inputs[0].material,'release_out')
    c('Resinworkupdoesnotclaimperfectfunctionalfilter',g.steps['flag_capture'].outputs[0],'flag_capture_kept')
    c('QCmeasurementdoesnotfeedfusionPCR',g.steps['scale_pcr'].inputs[0].material,'fusion_dna')
    src=Source(sources_root()/'cfms-mrna-display__10.1002+anie.201207005'/'source_pdf_3f348b50ea6f42f28215.pdf.txt','Published combinedmain+SI institutionalPDF,TableS1p21(S15),original rendered /tmp/cfms18-table21.png')
    if have(c,src,label='PublishedCFMSTableS1primertranscription'):
        txt=src.path.read_text();table=txt[txt.index('Library construction and sequencing oligos'):txt.index('Clone reconstruction oligos')]
        c('FiveprimerbodiesoccurverbatiminoriginalTableS1',all(s in table for s in PRIMERS.values()),True)
    c('Unassembledgeometrypreventsfinishedclaim',g.publishable,False);c.report()

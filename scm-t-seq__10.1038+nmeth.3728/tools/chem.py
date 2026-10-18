"""scM&T-seq: local source-qualified physical separation and native branch inputs."""
from chemdraw import Construct, Scene, Segment, feature
from protocol_model import role
from protocol_graph import (Disposition, Terminal, Highlight, InputRole, MaterialKind,
                            ProtocolGraph, SequencingPlan, StepKind)
from pcr import PCRGeometry
import illumina as il
import seqprimers as sp

TITLE='scM&T-seq — transcriptome and DNA methylome'
NOTES='01_scm-t-seq.html'
SOURCE='<a href="https://doi.org/10.1038/nmeth.3728">Angermueller et al., Nature Methods2016</a>, OnlineMethods82–88; explicitly cited2015G&T-seq and2014scBS-seq component scopes are named locally.'
SUMMARY='G&T-style bead capture physically separates poly(A) RNA from genomic DNA. The RNA route uses inherited template-switch cDNA amplification and scaled Nextera XT; the DNA route modifies scBS-seq with random hexamers and omission of first-strand biotin capture.'
CAVEAT='The native graph retains source-specific intermediates and both125bp paired-end outputs. Historical Nextera XT end/primer geometry, selected index members and acquired index cycles remain unassigned. The cited2015G&T TSO has AC+rGrG+G; this is an explicitly source-qualified inherited reference, not a laterACAT substitution. The exact2016TSO/linker assignment remains inferred. Primer-derived sequenced DNA is not an experimentally discarded fraction; exact concatemer assembly remains unresolved.'

# Exact components read in the cited2015G&T article, pp.13–14.
SMART_HANDLE='AAGCAGTGGTATCAACGCAGAGT'
CAPTURE_LINKER='AC'
G_AND_T_TSO_LINKER='AC'
CAPTURE=SMART_HANDLE+CAPTURE_LINKER+'T'*30+'VN'
ISPCR=SMART_HANDLE
TSO_WRITTEN=SMART_HANDLE+G_AND_T_TSO_LINKER+'rGrG+G'
# Smallwood2014 Methods prints these handles and PE1;2016 explicitly replaces N9 with N6.
DNA_H1='CTACACGACGCTCTTCCGATCT'
DNA_H2='TGCTGAACCGCTCTTCCGATCT'
DNA_PE1=il.P5+il.TRUSEQ_READ1[4:]
N_RANDOM=6
RNA_INDEX=feature('scmandt_rna_sample_index','sample_index','unknown',group='paired_cell',part='RNA well/library identity; member and native position unassigned')
DNA_INDEX=feature('scmandt_dna_ipcrtag','sample_index','unknown',group='paired_cell',part='matched DNA well/library identity; member and width unassigned')
RNA=Construct([role('native XT end regions; arrangement and index width unassigned',inferred=True,molecular_feature=RNA_INDEX),role('transcript-derived cDNA fragment'),role('opposite native XT end region; architecture unassigned',inferred=True)],name='scM&T RNA library')
DNA=Construct([Segment('source PE1.0',DNA_PE1),Segment('oligo1 random hexamer','N'*N_RANDOM,placeholder=True),role('DNA-copy sequence reflecting bisulfite conversion'),Segment('opposite random-priming hexamer','N'*N_RANDOM,placeholder=True),role('native iPCRTag-associated end; strand/index architecture unassigned',inferred=True,molecular_feature=DNA_INDEX)],name='scM&T methylome library')
RP=();DP=()
FINAL_LIBRARIES=(('RNA library',RNA,RP,'INFERRED native XT end roles — transcriptome endpoint; no default completed XT shell.','Published125bp paired-end reads.'),('DNA methylome library',DNA,DP,'INFERRED native opposite end — PE1 and random-hexamer components are source-qualified; native iPCRTag assembly unassigned.','Published125bp paired-end reads; index acquisition unassigned.'))
SEQUENCING_UNAVAILABLE='Physical native sequencing-primer sites and selected iPCRTag/index species remain unassigned within this capped reconstruction. Canonical arm compatibility is not proof of acquired index reads.'
READ_SPECS={name:(sp.read_spec('Read 1',125,'transcript cDNA' if name=='RNA library' else 'DNA-copy methylation signal; random-hexamer read prefix','include'),sp.read_spec('Read 2',125,'opposite insert end','include'),sp.read_spec('Index 1 (i7)','unassigned','Index acquisition and cycles unassigned','undetermined'),sp.read_spec('Index 2 (i5)','unassigned','Index acquisition and cycles unassigned','undetermined')) for name,*_ in FINAL_LIBRARIES}


def oligos():
    return (f'G&T2015 inherited capture primer:5′-biotin-triethyleneglycol-{CAPTURE}-3′',f'G&T2015 TSO reference:5′-{SMART_HANDLE}<inf>{G_AND_T_TSO_LINKER}</inf>rGrG+G-3′',f'G&T2015 ISPCR:5′-{ISPCR}-3′',f'scM&T oligo1:5′-{DNA_H1}NNNNNN-3′; native5′modification assignment unassigned',f'scM&T oligo2:5′-{DNA_H2}NNNNNN-3′',f'inherited PE1.0:5′-{DNA_PE1}-3′')


def populations(*names):
    sc=Scene()
    for i,name in enumerate(names):sc.strand('species_'+str(i),[role(name)],label=name)
    return sc


def capture_primer():
    sc=Scene();sc.strand('capture',[Segment('SMART handle',SMART_HANDLE),Segment('capture linker',CAPTURE_LINKER),Segment('dT30VN','T'*30+'VN',placeholder=True)],label='G&T2015 source capture oligo',mod5='biotin-triethyleneglycol')
    return sc


def rna_state(switched=False,amplified=False):
    sc=Scene()
    if amplified:
        body=[role('amplified full-length cDNA; native end-copy extent unassigned',inferred=True)]
        sc=Scene.duplex(body,unpaired=(body[0].name+"'",),label='amplified cDNA')
    else:
        parts=[Segment('capture-derived SMART handle',SMART_HANDLE),Segment('capture linker',CAPTURE_LINKER),Segment('dT30VN','T'*30+'VN',placeholder=True),role('antisense transcript cDNA')]
        if switched:parts.append(role('copied G&T TSO boundary; native extension extent unassigned',inferred=True))
        sc.strand('cDNA',parts,label='single-strand cDNA; complete first-copy geometry unassigned',mod5='biotin-triethyleneglycol')
    sc.strand('RNA',[role('independent RNA template population; later recovery/degradation unassigned')])
    return sc


def first_dna_copies(second=False):
    sc=Scene();sc.strand('first copies',[Segment('oligo1 handle',DNA_H1),Segment('N6','N'*N_RANDOM,placeholder=True),role('first-primer DNA extension copies; native extent unassigned')],label='oligo1-derived DNA strand classes',mod5='native5′modification unassigned')
    sc.strand('original DNA',[role('carried converted DNA originals with U and retained methylcytosine; partition unassigned')])
    if second:sc.strand('second copies',[Segment('oligo2 handle',DNA_H2),Segment('N6','N'*N_RANDOM,placeholder=True),role('complementary DNA-copy classes; native complete assembly unassigned')],label='separate oligo2-derived copies')
    return sc


def cleanup(g,key,template,reagent,library=False):
    select=g.input(key+'_selection',reagent,MaterialKind.SELECTION)
    kept=g.product(key+'_retained','source-recovered '+template.label,MaterialKind.LIBRARY if library else MaterialKind.MOLECULAR,state=template.state)
    waste=g.product(key+'_waste','physically removed bead supernatant and washes',MaterialKind.WASTE)
    g.react(key,'Recover source bead-selected fraction',StepKind.CLEANUP,((template,InputRole.TEMPLATE),(select,InputRole.SELECTION)),(kept,waste),changes_molecule=False)
    g.discard(waste,'Physically removed source fractions; individual RNA/primer/hydrolysate partition is not assumed complete.')
    return kept


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    lysate=g.input('lysate','single-cell RLTPlus/SUPERase-In lysate',MaterialKind.BIOLOGICAL,state=populations('genomic DNA','poly(A) RNA'),primary=True)
    beads=g.input('bead_capture','G&T2015 biotin-triethyleneglycol oligo(dT) streptavidin beads',MaterialKind.SELECTION,state=capture_primer())
    rna=g.product('rna_fraction','bead-bound poly(A) RNA fraction',MaterialKind.MOLECULAR,state=populations('bead-bound poly(A) RNA; native pairing unassigned'))
    dna=g.product('dna_fraction','genomic DNA supernatant plus retained washes',MaterialKind.MOLECULAR,state=populations('genomic DNA'))
    g.react('separate','Capture RNA; physically transfer DNA supernatant and washes',StepKind.AFFINITY_CAPTURE,((lysate,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(rna,dna),changes_molecule=False)
    rt_primer=g.input('rt_primer','inherited G&T capture primer as RT primer',MaterialKind.OLIGO,state=capture_primer())
    rt=g.input('rt','G&T2015 SuperScriptII reference; study-specific formulation unassigned',MaterialKind.CATALYST)
    first=g.product('rna_first','first-strand cDNA beside RNA template',MaterialKind.MOLECULAR,state=rna_state())
    g.react('rna_rt','Reverse-transcribe bead-captured RNA',StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(first,))
    tso_state=Scene();tso_state.strand('TSO',[Segment('SMART handle',SMART_HANDLE),Segment('inherited AC linker',G_AND_T_TSO_LINKER,inferred=True),Segment('rGrG+G','GGG',note='Source first two residues are RNA G; terminal residue is LNA G')],label='G&T2015 inherited reference, native2016 linker assignment inferred')
    tso=g.input('tso','G&T2015 TSO; native linker assignment inferred',MaterialKind.OLIGO,state=tso_state,inferred=True)
    switched=g.product('rna_switched','template-switched single-strand cDNA mixture',MaterialKind.MOLECULAR,state=rna_state(True),inferred=True)
    g.react('rna_switch','INFERRED inherited TSO choice — switch nascent cDNA',StepKind.TEMPLATE_SWITCHING,((first,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(switched,))
    ispcr=g.input('ispcr','G&T2015 ISPCR species; two complete native sites are not assumed on ss input',MaterialKind.OLIGO,state=ISPCR)
    pol=g.input('rna_preamplification_polymerase','G&T2015 KAPAHiFi reference; individual2016 cycle choice unassigned',MaterialKind.CATALYST)
    preamp=g.product('rna_preamp','amplified full-length cDNA mixture',MaterialKind.MOLECULAR,state=rna_state(True,True))
    g.react('rna_preamp_pcr','Amplify native template-switched cDNA',StepKind.PCR,((switched,InputRole.TEMPLATE),(ispcr,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(preamp,),geometry=PCRGeometry.unavailable(switched.state,symmetric_species=True,reason='Actual ss cDNA mixture, not completed duplex, enters one-species ISPCR. Native first-copy/symmetric site geometry unassigned.'))
    preamp=cleanup(g,'rna_preamp_cleanup',preamp,'G&T2015 reference1:1AMPure recovery; exact2016 elution choice unassigned')
    trans=g.input('xt_transposome','historical NexteraXT native loaded adapters; strands/gaps unassigned',MaterialKind.ADAPTER)
    tn5=g.input('xt_tn5','source NexteraXT one-fifth-volume preparation',MaterialKind.CATALYST)
    tagged=g.product('rna_tagged','actual native XT-tagged cDNA strand classes',MaterialKind.MOLECULAR,state=populations('native XT-tagged cDNA; transferred strands/gaps unassigned','carried RNA/other reaction species; native partition unassigned'))
    other=g.product('rna_other_ends','other XT insertion/end classes',MaterialKind.MOLECULAR)
    g.react('rna_tagment','Tagment recovered cDNA',StepKind.TAGMENTATION,((preamp,InputRole.TEMPLATE),(trans,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,other))
    g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Nonamplification does not prove physical removal; source cleanup recovery unassigned.'))
    xf=g.input('xt_forward','native XT forward indexing primer; selected species unassigned',MaterialKind.OLIGO)
    xr=g.input('xt_reverse','native XT reverse indexing primer; selected species unassigned',MaterialKind.OLIGO)
    xp=g.input('xt_polymerase','XT PCR polymerase/cycle choice unassigned',MaterialKind.CATALYST)
    rna_lib=g.product('rna_library','scM&T RNA library',MaterialKind.LIBRARY,state=RNA,inferred=True)
    g.react('rna_index_pcr','Complete source XT RNA library',StepKind.PCR,((tagged,InputRole.TEMPLATE),(xf,InputRole.STRUCTURAL),(xr,InputRole.STRUCTURAL),(xp,InputRole.CATALYST)),(rna_lib,),geometry=PCRGeometry.unavailable(tagged.state,reason='Actual transferred/gapped native input retained. No final XT duplex, default index width or kit primer geometry supplies the precursor.'))
    g.sequence(rna_lib,SequencingPlan(library=RNA,acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':125,'Read 2':125},limitation=SEQUENCING_UNAVAILABLE))
    dna=cleanup(g,'dna_recovery',dna,'G&T-style AMPure recovery; defining2016 DNA elution10µlH2O')
    reagent=g.input('bisulfite','EZMethylationDirect reagent with AMPureXP present; magnetic desulphonation/purification',MaterialKind.CONDITION)
    converted=g.product('dna_converted','recovered converted genomic DNA originals',MaterialKind.MOLECULAR,state=populations('native bisulfite-converted DNA: C deamination to U before copying; protected methylcytosines retained'))
    g.react('dna_bisulfite','Convert purified genomic DNA and recover into first-strand mixture',StepKind.CHEMICAL_CONVERSION,((dna,InputRole.TEMPLATE),(reagent,InputRole.CONDITION)),(converted,))
    one_state=Scene();one_state.strand('oligo1',[Segment('inherited handle',DNA_H1),Segment('random hexamer','N'*N_RANDOM,placeholder=True)],label='source-qualified scBS handle,2016 N6 substitution;5′modification unassigned',mod5='native5′chemistry unassigned')
    one=g.input('dna_oligo1','oligo1 with source N6 substitution; biotin capture omitted',MaterialKind.OLIGO,state=one_state)
    kl=g.input('dna_klenow','inherited scBS Klenowexo− strand synthesis',MaterialKind.CATALYST)
    copies=g.product('dna_first_copies','first-primer extension products plus carried DNA',MaterialKind.MOLECULAR,state=first_dna_copies())
    g.react('dna_first_synthesis','Repeated scBS-style first-strand synthesis with N6 oligo1',StepKind.EXTENSION,((converted,InputRole.TEMPLATE),(one,InputRole.STRUCTURAL),(kl,InputRole.CATALYST)),(copies,))
    exo=g.input('dna_exoi','inherited scBS exonucleaseI excess-primer digestion',MaterialKind.CATALYST)
    digested=g.product('dna_after_exoi','first-strand mixture after excess-primer digestion',MaterialKind.MOLECULAR,state=copies.state)
    hydro=g.product('dna_primer_hydrolysate','primer digestion products',MaterialKind.MOLECULAR)
    g.react('dna_exonuclease','Digest excess first-strand primers',StepKind.DIGESTION,((copies,InputRole.TEMPLATE),(exo,InputRole.CATALYST)),(digested,hydro))
    g.mark_terminal(Terminal(hydro.id,Disposition.UNRESOLVED,'Catalysis is not separation; exact recovery through subsequent beads unassigned.'))
    digested=cleanup(g,'dna_first_cleanup',digested,'inherited scBS0.8×AMPure recovery; no streptavidin capture')
    two=g.input('dna_oligo2','oligo2 with source N6 substitution',MaterialKind.OLIGO,state=Construct([Segment('oligo2 handle',DNA_H2),Segment('random hexamer','N'*N_RANDOM,placeholder=True)]))
    completed=g.product('dna_second_copies','oligo1/oligo2-derived DNA copy population',MaterialKind.MOLECULAR,state=first_dna_copies(True))
    g.react('dna_second_synthesis','Copy complementary strands with oligo2',StepKind.EXTENSION,((digested,InputRole.TEMPLATE),(two,InputRole.STRUCTURAL),(kl,InputRole.CATALYST)),(completed,))
    completed=cleanup(g,'dna_extra_cleanup',completed,'defining2016 extra0.8×AMPureXP between second synthesis and PCR')
    f=g.input('dna_pe1','inherited PE1.0 forward primer',MaterialKind.OLIGO,state=DNA_PE1)
    r=g.input('dna_ipcrtag','inherited indexed iPCRTag family; native selected strand/member/width unassigned',MaterialKind.OLIGO)
    dp=g.input('dna_pcr_polymerase','inherited scBS KAPAHiFi reference; selected2016 cycle choice unassigned',MaterialKind.CATALYST)
    dna_lib=g.product('dna_library','scM&T methylome library',MaterialKind.LIBRARY,state=DNA,inferred=True)
    contam=g.product('primer_derived_library','primer-derived non-genomic DNA represented among source reads; assembly/formation stage unassigned',MaterialKind.LIBRARY,state=populations('primer-derived non-genomic DNA: exact native molecular assembly unassigned'))
    g.react('dna_pcr','Amplify native oligo-derived copies with PE1.0/iPCRTag',StepKind.PCR,((completed,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(dp,InputRole.CATALYST)),(dna_lib,contam),geometry=PCRGeometry.unavailable(completed.state,reason='Actual oligo1/oligo2 copied mixture enters PCR. Selected iPCRTag/native complete junctions remain unassigned; final endpoint is not a template.'))
    final_beads=g.input('dna_final_beads','inherited scBS0.8×AMPure post-PCR recovery; exact2016 native-class recovery unassigned',MaterialKind.SELECTION)
    final_dna=g.product('dna_final_library','source recovered methylome library pool target class',MaterialKind.LIBRARY,state=DNA)
    final_contam=g.product('dna_final_contamination','co-recovered sequenced primer-derived class; no asserted complete removal',MaterialKind.LIBRARY,state=contam.state)
    final_waste=g.product('dna_final_waste','physically removed post-PCR bead supernatant/washes',MaterialKind.WASTE)
    g.react('dna_final_cleanup','Recover coexisting library classes through source bead purification',StepKind.CLEANUP,((dna_lib,InputRole.TEMPLATE),(contam,InputRole.TEMPLATE),(final_beads,InputRole.SELECTION)),(final_dna,final_contam,final_waste),changes_molecule=False)
    g.discard(final_waste,'Bead supernatant/washes physically removed; this is not a claim that all primer-derived species were excluded.')
    g.sequence(final_dna,SequencingPlan(library=DNA,acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':125,'Read 2':125},limitation=SEQUENCING_UNAVAILABLE))
    g.sequence(final_contam,SequencingPlan(limitation='Source reports high primer contamination among sequenced DNA reads. Exact concatemer assembly/entry sites are unassigned; this is not an experimental discard.'))
    g.annotate(Highlight('dna_pcr','Sequenced primer-derived products are not waste','The defining article reports high primer contamination and refers to SupplementaryFig3. Earlier notes describe first-primer concatemers from that figure; its exact assembly was not re-read here. The graph preserves sequenced non-genomic products without claiming a formation mechanism or physical removal.'))
    g.annotate(Highlight('rna_switch','Source/version scope for inherited RNA oligos','The cited2015G&T original prints the AC+rGrG+G form and biotin-triethyleneglycol capture primer. Its common sequences are retained locally; native2016 linker assignment remains inferred. A laterACAT reference does not override this selected cited reference.'))
    g.annotate(Highlight('dna_first_synthesis','scBS inheritance and local changes','The2016 method explicitly inherits Smallwood2014 with N6 in both primers, no first-strand biotin capture, and extra0.8×cleanup. The original5′biotin primer modification is not assumed removed merely because capture was omitted. The named iPCRTag family is not replaced by a guessed8nt canonical index.'))
    g.declare_reviewed_scope(REVIEWED,SOURCE)
    return g.require_publishable()

REVIEWED=('lysate','bead_capture','separate','rt_primer','rt','rna_rt','tso','rna_switch','ispcr','rna_preamplification_polymerase','rna_preamp_pcr','rna_preamp_cleanup_selection','rna_preamp_cleanup','xt_transposome','xt_tn5','rna_tagment','xt_forward','xt_reverse','xt_polymerase','rna_index_pcr','rna_library','dna_recovery_selection','dna_recovery','bisulfite','dna_bisulfite','dna_oligo1','dna_klenow','dna_first_synthesis','dna_exoi','dna_exonuclease','dna_first_cleanup_selection','dna_first_cleanup','dna_oligo2','dna_second_synthesis','dna_extra_cleanup_selection','dna_extra_cleanup','dna_pe1','dna_ipcrtag','dna_pcr_polymerase','dna_pcr','dna_library','primer_derived_library','dna_final_beads','dna_final_cleanup','dna_final_library','dna_final_contamination')

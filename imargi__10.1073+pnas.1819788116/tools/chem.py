"""Source-owned iMARGI linker, strand transfer and circle relinearization."""
from chemdraw import Construct,Scene,Segment,MolecularState,feature,revcomp,complement_segments,circle_rows
from protocol_model import role
from protocol_graph import InputRole,MaterialKind,ProtocolGraph,SequencingPlan,StepKind,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from circular import CircularizationGeometry
from pcr import PCRGeometry
import seqprimers as sp
TITLE="iMARGI — in situ RNA–genome interaction mapping"
NOTES="01_imargi.html"
SOURCE='<a href="https://doi.org/10.1073/pnas.1819788116">Yan et al.,PNAS2019</a>; directly linked <a href="https://doi.org/10.1038/s41596-019-0229-4">Wu et al.,NatProtocols2019</a>.'
SUMMARY="In situ RNA–linker–DNA ligation, on-bead copying and strand release yield a single-stranded contact chimera. Circularization and BamHI relinearization reverse the readable RNA/DNA end order before PCR."
TOP='NNAGATCGGAAGAGCGTCGTGTAGGGAGGATCCGTTCAGACGTGTGCTCTTCCGATCT'
BOTTOM='GATCGGAAGAGCACACGTCTGAACGGATCCTCCCTACACGACGCTCT'
CUT='TCGTGTAGGGAGGATCCGTTCAGACGTGTGCTCT'
UNIVERSAL='AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
INDEX_PREFIX='CAAGCAGAAGACGGCATACGAGAT'
INDEX_SUFFIX='GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
SAMPLE=feature('imargi.sample','sample_index','unknown',note='source six-base index family; selected member and acquired cycles unassigned')
SEQ_PRIMERS=()
SHOW_READ_CONTENT=False
READ_SPECS=tuple(sp.read_spec(r,100 if r.startswith('Read') else 'unassigned','RNA-derived end' if r=='Read 1' else 'GenomicDNA end; sourceAluI CT prefix' if r=='Read 2' else 'Index acquisition unassigned; physicali7 family6nt', 'include' if r.startswith('Read') else 'undetermined') for r in ('Read 1','Read 2','Index 1 (i7)','Index 2 (i5)'))
SEQUENCING_ENDING='PNAS PE100 HiSeq4000. Index acquisition unassigned; two linker NN bases are ligation bias, not a UMI.'
def primer(name,seq):
 sc=Scene();sc.strand(name,[Segment(name,seq[:-1]),Segment('terminalT',seq[-1])],label='sourcePCR DNA oligo');sc.junction(name,name,'terminalT',text='source terminalC–T phosphorothioate bond');return sc

def top_oligo(activated=False):
 sc=Scene();sc.strand('TOP',[Segment('TOPthroughPS_C',TOP[:-5]),Segment('GAafterPS','GA'),Segment('physicalinternalbiotinT','T'),Segment('TOPterminalCT','CT')],mod5='App' if activated else 'p',label='sourceTOP internalbiotin attached to physicalT')
 sc.junction('TOP','TOPthroughPS_C','GAafterPS',text='source C–G phosphorothioate bond');return sc

def bottom_oligo():
 sc=Scene();sc.strand('BOTTOM',[Segment('sourceBOTTOM',BOTTOM)],mod5='p',label='source5′PBOTTOM; no assignedPS');return sc

def cut_binding():
 start=BOTTOM.index(revcomp(CUT))
 sc=Scene();sc.strand('circlelocalBOTTOM',[Segment('BOTTOMprefix',BOTTOM[:start]),Segment('cutbinding',BOTTOM[start:start+len(CUT)]),Segment('BOTTOMsuffix',BOTTOM[start+len(CUT):])],label='localpiece of samecircle; full fold not asserted')
 sc.anneal('blockedCut',[Segment('sourceCUT',CUT),Segment('additionalinvertedT','T')],to='circlelocalBOTTOM',pair=('sourceCUT','cutbinding'),unpaired=('additionalinvertedT',),label='sourceadditional3′InvdT; nonextendable')
 return sc

def index_primer():
 sc=Scene();sc.strand('index primer',[Segment('P7 sourceprefix',INDEX_PREFIX),Segment('selected6ntsampleindex','N'*6,feature=SAMPLE),Segment('sourceannealingarm',INDEX_SUFFIX)],label='source6ntfamily; memberunassigned, terminalC–T PS');return sc

def cut_primer():
 sc=Scene();sc.strand('Cut_oligo',[Segment('sourcecutoligo',CUT),Segment('additional physical inverted terminaldT','T')],label='3′inverted dT blocks extension; not copied');return sc

def linker_scene():
 sc=Scene();sc.strand('top',[Segment('ligationbiasNN','NN'),Segment('unpairedRNA-side tail',TOP[2:10]),Segment('pairedcore',TOP[10:-1]),Segment('3primeT','T')],mod5='App',label='activated linker: original PS after ...TCC; internalbiotin on physicalT in ...GATCT')
 sc.anneal('bottom',[Segment('pairedBOTTOM',BOTTOM)],to='top',pair=('pairedBOTTOM','pairedcore'),mod5='p',label='source5′Pbottom strand')
 sc.ligation_end('top','3primeT',kind='cohesive',polarity='3-prime',owner='adapter',text='source3′T DNA-ligation end');sc.ligation_end('top','ligationbiasNN',kind='single-strand',polarity='5-prime',owner='adapter',text='sourceactivated5′App RNA-ligation end');return sc

def released_parts():
 return [role('genomicDNA original strand'),Segment('original sourcegenomicDNA3′dAtail','A'),Segment('BOTTOM',BOTTOM),Segment('copied RNA-side linker tail',revcomp(TOP[2:10])),role('copiedNN ligationbias; notUMI',length_bp=2),role('RNA-derived cDNA')]

def native_state(stage):
 if stage in ('released','linear','circle','annealed'):
  parts=released_parts()
  if stage=='linear':
   cut=BOTTOM.index('GGATCC')+1
   parts=[Segment('BamHIright halfBOTTOM',BOTTOM[cut:]),*parts[3:],*parts[:2],Segment('BamHIleft halfBOTTOM',BOTTOM[:cut])]
  if stage in ('circle','annealed'):
   rows=list(circle_rows(Construct(parts,name='source5′P single strand circle'),closure_label='ssDNA covalent closure; sourceCircLigase'))
   if stage=='annealed':rows.extend(cut_binding().rows())
   residual=Scene();residual.strand('residuals',[role('unreactedlinearSS, possibleRNAhydrolysates/oligofragments; no physicalseparation atcircleclosure')],label='independent residualpopulations');rows.extend(residual.rows())
   return MolecularState('single-strand contact circle'+(' plus localCutoligo; BamHIsite paired, globalfoldunassigned' if stage=='annealed' else ''),tuple(rows))
  sc=Scene();sc.strand('SS contact',parts,mod5='p',label='actual single-strand '+stage+' precursor');residues='possible originalRNA/oligo fragments; hydrolysis/recovery bounded by physical purification' if stage=='released' else 'possible uncircularizedSS, residualcircle/Cutfragments, RNA/oligofragments; physicalrecovery unassigned';sc.strand('possible residues',[role(residues)],label='independent residualclass');return sc
 sc=Scene()
 if stage=='library':
  sc.strand('amplified native product',[Segment('sourceuniversalprimer',UNIVERSAL),role('copiedNN ligationbias',length_bp=2),role('RNAderivedcDNA'),role('genomicDNA'),Segment('copiedindexannealingarm',revcomp(INDEX_SUFFIX)),Segment('copied6ntindex','N'*6,feature=SAMPLE),Segment('copiedP7',revcomp(INDEX_PREFIX))],label='conditional PCRdaughter footprint; full native assembly/firstcycle bounded')
  sc.strand('remaining templates',[role('possible originalSS/circle/offtarget/primer lineages; recovery unassigned')],label='independent residualorigins');return sc
 dna=role('originalAluIgenomicDNA' if stage!='raw' else 'nativegenomicDNA')
 if stage in ('tailed','rna_joined'):
  sc.strand('DNA',[dna,Segment('3primeA','A')],label='originalgenomic3′dAend')
  comp=complement_segments([dna]);sc.anneal('DNAopposite',comp,to='DNA',pair=(comp[0].name,dna.name),unpaired=(comp[0].name,));sc.overhang('DNA','3primeA',polarity='3-prime',owner='substrate',text='source3′dA')
 elif stage not in ('joined','rt','phosphorylated'):sc.strand('DNA',[dna],label='independent genomicDNA')
 if stage in ('rna_joined','joined','rt','phosphorylated'):
  sc.strand('RNA-top',[role('originalfragmentedRNA; physicalU'),Segment('TOPoriginal PS/internalbiotin retained',TOP),*([role('originalgenomicDNA joined throughdA/Tend')] if stage in ('joined','rt','phosphorylated') else [])],label='original RNA–linker'+('–DNA' if stage in ('joined','rt','phosphorylated') else ''))
  if stage in ('rt','phosphorylated'):sc.strand('copiedstrand',released_parts(),mod5='p' if stage=='phosphorylated' else '',label='oppositeoriginalgenomicDNA–BOTTOM–copiedtail–cDNA; '+('source5′P installed' if stage=='phosphorylated' else '5′P installed later'))
  else:sc.strand('BOTTOM',[Segment('BOTTOM',BOTTOM)],mod5='p',label='independent annealedbottomstrand; global extentunassigned')
 else:sc.strand('RNA',[role('nativeRNA withphysicalU; '+('3′OH' if stage in ('rna_oh','tailed') else 'sourcefragmentendunassigned'))],label='independentRNA')
 sc.strand('otherclasses',[role('incompleteRNA-linkers/freeDNA/RNA/oligos; sourcewash-specific recovery unresolved')],label='independent offproducts, not automatically discarded');return sc

def rna_ligation_geometry():return AdapterLigationGeometry.known('single-strand',linker_scene(),statement='Source activated5′App linker joins RNA3′OH viaRnl2trKQ; original linker3′T remains available.')
def dna_ligation_geometry():return AdapterLigationGeometry.known('cohesive',native_state('tailed'),linker_scene(),statement='Source genomic3′A and linker3′T cohesive ends; exactsealed strandextent bounded.')
FINAL_LIBRARY=Construct([role('nativeP5PCRend'),role('RNAderivedcDNA'),role('contactgenomicDNA'),role('sourcei7index',length_bp=6,molecular_feature=SAMPLE),role('nativeP7PCRend')],name='iMARGI native PCR reference')
FINAL_CAPTION='INFERRED — source-specific component/identifier reference; actual PCR input is BamHI-relinearized single-strand DNA, not this final readout.'

CAVEAT = 'The author protocol prints the modified bivalent linker, cut oligo and PCR primers. The workflow follows Yan et al. 2019 Methods and the Wu et al. 2019 Procedure (Steps 73–114): exonuclease biotin removal, streptavidin capture, on-bead RT, 5′ phosphorylation, NaOH release, CircLigase circularization, Cut_oligo/BamHI relinearization, PCR and gel size selection. Native copied-linker tail and BamHI strand order are derived from source oligos. Global fold, first-cycle landing and complete native endpoint assembly remain bounded; original PS/biotin and blocked Cut_oligo are not silently copied. Published citation discrepancy: Yan et al. 2019 reference 10 cites <a href="https://doi.org/10.1016/j.cub.2017.01.068">DOI 10.1016/j.cub.2017.01.068</a> for Sridhar et al. MARGI, but that DOI is the erratum; the <a href="https://doi.org/10.1016/j.cub.2017.01.011">original article</a> has DOI 10.1016/j.cub.2017.01.011. The model uses the iMARGI detailed-protocol oligos. Original Sridhar 2017 MARGI Document S1 (Supplemental Experimental Procedures, lines 204–215) prints a top linker ending ...TCCGA/BiodT/CT without a C* mark; later Wu 2019 iMARGI Nature Protocols (JATS lines 275–285) prints ...TCC*GA/iBiodT/CT and defines C* as phosphorothioated C. The bases agree, but the printed modifications differ. This model uses the later iMARGI-specific C* oligo; the 2019 PNAS paper\'s statement that it reused the MARGI linker sequence does not establish which modification was in its experimental oligo. Apparent reagent typo in Wu et al. 2019 detailed protocol: stock-preparation text (extracted JATS lines 395–397) says ATP for a reagent labeled 10 mM dATP, while its dA-tailing reaction table (lines 673–685) and Yan et al. 2019 Methods (line 203) specify dATP. This model uses a 3′-dA tail from the reaction descriptions; the reagent actually used in the experiment remains unverified. Wu et al. 2019 also names the PCR master mix inconsistently: Step 104 (pilot PCR, JATS line 1159) says “2× NEBNext HC PCR Master Mix”, whereas Step 106 (line 1173 table) and the reagent list (line 259) give NEBNext High-Fidelity 2× PCR Master Mix (M0541). The model uses the High-Fidelity mix; “HC” appears to be a typo.'
def oligos():
    return ['<b>bivalent linker top</b> — 5′-/Phos/NNAGATCGGAAGAGCGTCGTGTAGGGAGGATCCGTTCAGACGTGTGCTCTTCC*GA/iBiodT/CT-3′ (C* is phosphorothioated)', '<b>bivalent linker bottom</b> — 5′-/Phos/GATCGGAAGAGCACACGTCTGAACGGATCCTCCCTACACGACGCTCT-3′', '<b>BamHI cut oligo</b> — 5′-TCGTGTAGGGAGGATCCGTTCAGACGTGTGCTCT/3InvdT/-3′', '<b>NEBNext universal PCR</b> — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC*T-3′', '<b>NEBNext index PCR</b> — 5′-CAAGCAGAAGACGGCATACGAGAT[6 bp barcode]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-S-T-3′']
REVIEWED=('chromatin', 'restriction_enzyme', 'alui_cut_chromatin', 'fragmentation', 'rnase_i', 'fragmented_contacts', 'rna_fragmentation', 't4_pnk_phosphatase', 'rna_3oh', 'rna_3_dephosphorylation', 'klenow_exo_minus', 'da_tailed_contacts', 'dna_a_tailing', 'bivalent_linker', 'rna_ligase', 'rna_joined', 'rna_ligation', 'washed_rna_joined', 'free_linker', 'free_linker_wash', 'dna_ligase', 'contact_chimera', 'dna_ligation', 'proteinase_k', 'decrosslinked', 'reverse_crosslinks', 'phenol_chloroform', 'purified_nucleic_acids', 'extract_nucleic_acids', 'exo_i_t4_dna_polymerase', 'biotin_removed', 'remove_unligated_biotin', 'streptavidin_c1', 'captured_chimera', 'unbound_nucleic_acids', 'streptavidin_pulldown', 'reverse_transcriptase', 'cdna_chimera', 'reverse_transcription', 't4_pnk_kinase', 'phosphorylated_cdna', 'phosphorylate_5_ends', 'denaturing_buffer', 'released_ssdna', 'bead_bound_strand', 'release_ssdna', 'purified_ssdna', 'hydrolysed_rna', 'hydrolyse_and_purify', 'circligase', 'circular_ssdna', 'circularize', 'cut_oligo', 'cut_oligo_annealed', 'anneal_cut_oligo', 'bamhi', 'relinearized_ssdna', 'bamhi_relinearize', 'p5_primer', 'p7_primer', 'polymerase', 'pilot_aliquot', 'production_aliquot', 'pilot_split', 'pilot_pcr_product', 'pilot_pcr', 'amplified_library', 'library_pcr', 'library', 'size_rejected', 'size_selection')
REVIEWED += ('recovered_linear','post_bamhi_waste','post_bamhi_silane','linker_5p_top','linker_bottom','adenylation_enzyme','activated_top','activate_top','anneal_linker','washed_rt','rt_wash_waste','rt_washes','washed_phos','phos_wash_waste','phos_washes')

def protocol_graph():
    source = 'Defining source: https://doi.org/10.1073/pnas.1819788116; detailed protocol: https://doi.org/10.1038/s41596-019-0229-4'
    graph = ProtocolGraph('iMARGI', authority=source)
    chromatin = graph.input('chromatin', 'formaldehyde-cross-linked, SDS-permeabilized nuclei', MaterialKind.BIOLOGICAL, primary=True, state=native_state('raw'))
    nuclease = graph.input('restriction_enzyme', 'AluI (with RNasinPlus, CutSmart)', MaterialKind.CATALYST)
    cut = graph.product('alui_cut_chromatin', 'AluI-digested chromatin in nuclei', MaterialKind.MOLECULAR, state=native_state('cut'))
    graph.react('fragmentation', 'Digest chromatin with AluI overnight at 37 °C', StepKind.DIGESTION, ((chromatin, InputRole.TEMPLATE), (nuclease, InputRole.CATALYST)), (cut,))
    rnase = graph.input('rnase_i', 'RNase I (1:10, 3 min, 37 °C)', MaterialKind.CATALYST)
    fragmented = graph.product('fragmented_contacts', 'fragmented RNA and AluI-cut chromatin', MaterialKind.MOLECULAR, state=native_state('fragmented'))
    graph.react('rna_fragmentation', 'Briefly fragment chromatin-associated RNA with RNase I', StepKind.FRAGMENTATION, ((cut, InputRole.TEMPLATE), (rnase, InputRole.CATALYST)), (fragmented,))
    pnk = graph.input('t4_pnk_phosphatase', 'T4 PNK in pH 6.5 phosphatase buffer', MaterialKind.CATALYST)
    rna_oh = graph.product('rna_3oh', 'RNA fragments with 3′-OH ends', MaterialKind.MOLECULAR, state=native_state('rna_oh'))
    graph.react('rna_3_dephosphorylation', 'Remove RNA 3′ phosphates in situ', StepKind.DEPHOSPHORYLATION, ((fragmented, InputRole.TEMPLATE), (pnk, InputRole.CATALYST)), (rna_oh,))
    klenow = graph.input('klenow_exo_minus', 'Klenow fragment (3′→5′ exo−) with 0.1 mM dATP', MaterialKind.CATALYST)
    tailed = graph.product('da_tailed_contacts', 'RNA 3′-OH ends and dA-tailed genomic DNA ends', MaterialKind.MOLECULAR, state=native_state('tailed'))
    graph.react('dna_a_tailing', 'Add a 3′-dA tail to AluI-cut genomic DNA in situ', StepKind.TAILING, ((rna_oh, InputRole.TEMPLATE), (klenow, InputRole.CATALYST)), (tailed,))
    raw_top=graph.input('linker_5p_top','printed5′P TOP oligo withoriginalPS/internalbiotin',MaterialKind.OLIGO,state=top_oligo())
    raw_bottom=graph.input('linker_bottom','printed5′P BOTTOM oligo',MaterialKind.OLIGO,state=bottom_oligo())
    adenylation=graph.input('adenylation_enzyme','MthRNA ligase/5′DNAadenylationkit; detailedconditionsnotreread',MaterialKind.CATALYST)
    activated=graph.product('activated_top','5′App activatedTOP retains originalPS/internalbiotin',MaterialKind.OLIGO,state=top_oligo(True))
    graph.react('activate_top','Adenylate source5′P linkerTOP; conditions bounded',StepKind.MAGIC,((raw_top,InputRole.TEMPLATE),(adenylation,InputRole.CATALYST)),(activated,),reason='Source enzyme and activation are positive; detailed adenylation conditions outside inspected sourcepass.')
    linker=graph.product('bivalent_linker','5′AppTOP annealed to5′PBOTTOM',MaterialKind.ADAPTER,state=linker_scene())
    graph.react('anneal_linker','Anneal source asymmetric bivalent linker',StepKind.HYBRIDIZATION,((activated,InputRole.TEMPLATE),(raw_bottom,InputRole.STRUCTURAL)),(linker,))
    rna_ligase = graph.input('rna_ligase', 'T4 RNA ligase 2, truncated KQ (20% PEG 8000)', MaterialKind.CATALYST)
    rna_joined = graph.product('rna_joined', 'RNA-end-ligated bivalent linker contacts', MaterialKind.MOLECULAR, state=native_state('rna_joined'))
    graph.react('rna_ligation', 'Ligate the RNA end of the bivalent linker in situ', StepKind.LIGATION, ((tailed, InputRole.TEMPLATE), (linker, InputRole.STRUCTURAL), (rna_ligase, InputRole.CATALYST)), (rna_joined,), geometry=rna_ligation_geometry())
    washed = graph.product('washed_rna_joined', 'washed nuclei with RNA-ligated linker', MaterialKind.MOLECULAR, state=native_state('rna_joined'))
    free_linker = graph.product('free_linker', 'excess unligated linker', MaterialKind.WASTE)
    graph.react('free_linker_wash', 'Wash nuclei five times to remove excess free linker', StepKind.CLEANUP, ((rna_joined, InputRole.TEMPLATE),), (washed, free_linker))
    graph.discard(free_linker, 'removed in PNK-buffer washes of intact nuclei')
    dna_ligase = graph.input('dna_ligase', 'T4 DNA ligase (2 mL proximity-ligation volume)', MaterialKind.CATALYST)
    chimera = graph.product('contact_chimera', 'RNA–linker–DNA chimeras plus incomplete RNA–linker products', MaterialKind.MOLECULAR, state=native_state('joined'))
    graph.react('dna_ligation', 'Ligate the bridge DNA end to nearby genomic DNA', StepKind.LIGATION, ((washed, InputRole.TEMPLATE), (linker, InputRole.STRUCTURAL), (dna_ligase, InputRole.CATALYST)), (chimera,), geometry=dna_ligation_geometry())
    protease = graph.input('proteinase_k', 'Proteinase K extraction buffer (1% SDS, 65 °C, 3 h)', MaterialKind.CATALYST)
    reversed_ = graph.product('decrosslinked', 'cross-link-reversed nucleic acids', MaterialKind.MOLECULAR, state=native_state('joined'))
    graph.react('reverse_crosslinks', 'Lyse nuclei and reverse cross-links', StepKind.DEPROTEINIZATION, ((chimera, InputRole.TEMPLATE), (protease, InputRole.CATALYST)), (reversed_,))
    phenol = graph.input('phenol_chloroform', 'phenol:chloroform:isoamyl alcohol, ethanol precipitation', MaterialKind.CONDITION)
    nucleic = graph.product('purified_nucleic_acids', 'purified DNA and RNA', MaterialKind.MOLECULAR, state=native_state('joined'))
    graph.react('extract_nucleic_acids', 'Extract and precipitate DNA and RNA', StepKind.EXTRACTION, ((reversed_, InputRole.TEMPLATE), (phenol, InputRole.CONDITION)), (nucleic,))
    exo = graph.input('exo_i_t4_dna_polymerase', 'Exonuclease I, then T4 DNA polymerase with dATP/dGTP at 12 °C', MaterialKind.CATALYST)
    debiotinylated = graph.product('biotin_removed', 'chimeras keep biotin; biotin cleaved from linkers not ligated to DNA', MaterialKind.MOLECULAR, state=native_state('joined'))
    graph.react('remove_unligated_biotin', 'Remove biotin from linkers not ligated to DNA', StepKind.DIGESTION, ((nucleic, InputRole.TEMPLATE), (exo, InputRole.CATALYST)), (debiotinylated,))
    streptavidin = graph.input('streptavidin_c1', 'Dynabeads MyOne Streptavidin C1, high-salt biotin washes', MaterialKind.SELECTION)
    captured = graph.product('captured_chimera', 'bead-bound biotinylated RNA–linker–DNA chimeras', MaterialKind.MOLECULAR, state=native_state('joined'))
    unbound = graph.product('unbound_nucleic_acids', 'de-biotinylated RNA–linker, unligated RNA/DNA and washed-off nucleic acids', MaterialKind.WASTE)
    graph.react('streptavidin_pulldown', 'Pull down biotinylated chimeras', StepKind.AFFINITY_CAPTURE, ((debiotinylated, InputRole.TEMPLATE), (streptavidin, InputRole.SELECTION)), (captured, unbound))
    graph.discard(unbound, 'supernatant and high-salt washes of the streptavidin capture')
    rt = graph.input('reverse_transcriptase', 'SuperScript III, 50 °C, on beads', MaterialKind.CATALYST)
    cdna = graph.product('cdna_chimera', 'bead-bound duplex with cDNA extended from the linker bottom strand', MaterialKind.MOLECULAR, state=native_state('rt'))
    graph.react('reverse_transcription', 'Reverse-transcribe the RNA from the linker bottom strand', StepKind.REVERSE_TRANSCRIPTION, ((captured, InputRole.TEMPLATE), (linker, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (cdna,))
    rt_wash=graph.product('washed_rt','magnetically recoveredonbead RTproducts;2BW+PNKwash',MaterialKind.MOLECULAR,state=cdna.state)
    rt_waste=graph.product('rt_wash_waste','RTsupernatant andwash fractions',MaterialKind.WASTE)
    graph.react('rt_washes','Magnetically recoverbeads; removeRTsupernatant;2BW+PNKwash',StepKind.CLEANUP,((cdna,InputRole.TEMPLATE),),(rt_wash,rt_waste));graph.discard(rt_waste,'Physicalsupernatant/washes; molecularclasspurity not presumed.')
    kinase = graph.input('t4_pnk_kinase', 'T4 PNK with 1 mM ATP (pH 7.6)', MaterialKind.CATALYST)
    phos = graph.product('phosphorylated_cdna', 'bead-bound duplex with 5′-phosphorylated DNA ends', MaterialKind.MOLECULAR, state=native_state('phosphorylated'))
    graph.react('phosphorylate_5_ends', 'Phosphorylate DNA 5′ ends on beads', StepKind.PHOSPHORYLATION, ((rt_wash, InputRole.TEMPLATE), (kinase, InputRole.CATALYST)), (phos,))
    phos_wash=graph.product('washed_phos','onbead phosphorylatedproducts after2BWwashes',MaterialKind.MOLECULAR,state=phos.state)
    phos_waste=graph.product('phos_wash_waste','twoBWwash fractions',MaterialKind.WASTE)
    graph.react('phos_washes','Wash phosphorylatedbeads twice inBW',StepKind.CLEANUP,((phos,InputRole.TEMPLATE),),(phos_wash,phos_waste));graph.discard(phos_waste,'PhysicalBWwashfractions; no presumedpurity.')
    naoh = graph.input('denaturing_buffer', '0.1 M NaOH, 0.1 mM EDTA', MaterialKind.CONDITION)
    released = graph.product('released_ssdna', 'released ssDNA: genomic DNA–linker bottom strand–cDNA', MaterialKind.MOLECULAR, state=native_state('released'))
    bead_strand = graph.product('bead_bound_strand', 'biotinylated top-strand/RNA side retained on beads', MaterialKind.WASTE)
    graph.react('release_ssdna', 'Release the cDNA-containing strand with NaOH', StepKind.STRAND_SEPARATION, ((phos_wash, InputRole.TEMPLATE), (naoh, InputRole.CONDITION)), (released, bead_strand))
    graph.discard(bead_strand, 'beads magnetically removed; supernatant kept')
    ssdna = graph.product('purified_ssdna', 'purified single-stranded chimera DNA', MaterialKind.MOLECULAR, state=native_state('released'))
    hydrolysate = graph.product('hydrolysed_rna', 'alkali-hydrolysed RNA and SILANE flow-through', MaterialKind.WASTE)
    graph.react('hydrolyse_and_purify', 'Hydrolyse residual RNA at 98 °C, neutralize and purify on SILANE beads', StepKind.CLEANUP, ((released, InputRole.TEMPLATE),), (ssdna, hydrolysate))
    graph.discard(hydrolysate, 'not retained by SILANE bead purification')
    circligase = graph.input('circligase', 'CircLigase ssDNA ligase, 60 °C 4 h', MaterialKind.CATALYST)
    circle = graph.product('circular_ssdna', 'circularized single-stranded chimera', MaterialKind.MOLECULAR, state=native_state('circle'))
    graph.react('circularize', 'Circularize the released ssDNA', StepKind.CIRCULARIZATION, ((ssdna, InputRole.TEMPLATE), (circligase, InputRole.CATALYST)), (circle,), geometry=CircularizationGeometry.known('single-strand intramolecular ligation', statement='Wu et al. 2019 Steps 89–98: the released ssDNA is 5′-phosphorylated on beads and circularized with CircLigase ssDNA ligase.'))
    cut_oligo = graph.input('cut_oligo', 'BamHI Cut_oligo with 3′ inverted dT', MaterialKind.OLIGO, state=cut_primer())
    annealed = graph.product('cut_oligo_annealed', 'circle with a local Cut_oligo duplex over the linker BamHI site', MaterialKind.MOLECULAR, state=native_state('annealed'))
    graph.react('anneal_cut_oligo', 'Anneal the Cut_oligo across the linker BamHI site', StepKind.HYBRIDIZATION, ((circle, InputRole.TEMPLATE), (cut_oligo, InputRole.STRUCTURAL)), (annealed,))
    bamhi = graph.input('bamhi', 'BamHI, 37 °C 1 h; SILANE bead purification', MaterialKind.CATALYST)
    linear = graph.product('relinearized_ssdna', 'BamHI-relinearized ssDNA: half bottom linker–copied linker tail–cDNA–genomicDNA–other half', MaterialKind.MOLECULAR, state=native_state('linear'))
    graph.react('bamhi_relinearize', 'Relinearize the circle with BamHI, splitting the linker', StepKind.DIGESTION, ((annealed, InputRole.TEMPLATE), (bamhi, InputRole.CATALYST)), (linear,))
    recovered_linear=graph.product('recovered_linear','SILANE-recovered relinearizedSS plus possibleco-recoveredclasses',MaterialKind.MOLECULAR,state=linear.state)
    digest_waste=graph.product('post_bamhi_waste','SILANEflowthrough/washfractions',MaterialKind.WASTE)
    graph.react('post_bamhi_silane','Purify BamHIdigest onSILANE beforepilot/finalPCR',StepKind.CLEANUP,((linear,InputRole.TEMPLATE),),(recovered_linear,digest_waste));graph.discard(digest_waste,'Physicalnonretained fractions; no classpurity assumed.')
    p5 = graph.input('p5_primer', 'NEBNext Universal PCR primer', MaterialKind.OLIGO, state=primer('universal', UNIVERSAL))
    p7 = graph.input('p7_primer', 'NEBNext Index PCR primer (6-nt index unassigned)', MaterialKind.OLIGO, state=index_primer())
    polymerase = graph.input('polymerase', 'NEBNext High-Fidelity 2× PCR Master Mix; finalcycles selected by pilot minus one, not fixed12', MaterialKind.CATALYST)
    pilot = graph.product('pilot_aliquot', '5µl relinearized product for pilotPCR', MaterialKind.MOLECULAR, state=linear.state)
    production = graph.product('production_aliquot', 'remaining actual relinearized ssDNA for finalPCR', MaterialKind.MOLECULAR, state=linear.state)
    graph.react('pilot_split', 'Allocate 5µl to pilot PCR; retain production template', StepKind.PARTITION, ((recovered_linear, InputRole.TEMPLATE),), (pilot, production), changes_molecule=False)
    pilot_product = graph.product('pilot_pcr_product', 'pilot amplified DNA for cycle titration', MaterialKind.QC, state=native_state('library'))
    graph.react('pilot_pcr', 'Pilot amplification: sample at10/13/15/18/22cycles', StepKind.PCR, ((pilot, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL), (p7, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (pilot_product,), geometry=PCRGeometry.unavailable(pilot.state, reason='Actual SS pilot precursor; same bounded first-cycle geometry.'))
    graph.mark_terminal(Terminal(pilot_product.id, Disposition.QC, 'Pilot agarose QC determines finalcycle choice minus one; not production pooled material.'))
    amplified = graph.product('amplified_library', FINAL_LIBRARY.name + ' (before size selection)', MaterialKind.LIBRARY, state=native_state('library'))
    graph.react('library_pcr', 'PCR-amplify the relinearized ssDNA from its half-linker ends', StepKind.PCR, ((production, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL), (p7, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (amplified,), geometry=PCRGeometry.unavailable(production.state, reason='Actual relinearized single-strand precursor. First-cycle completion/native full primer landing geometry bounded; final library is not the template.'))
    library = graph.product('library', FINAL_LIBRARY.name, MaterialKind.LIBRARY, state=native_state('library'))
    out_of_range = graph.product('size_rejected', 'primer dimers and fragments outside 250–1000 bp', MaterialKind.WASTE)
    graph.react('size_selection', 'MinElute purify and cut 250–1000 bp from a 2% agarose E-gel', StepKind.CLEANUP, ((amplified, InputRole.TEMPLATE),), (library, out_of_range))
    graph.discard(out_of_range, 'left in the gel outside the excised 250–1000 bp range')
    graph.sequence(library, SequencingPlan(library=FINAL_LIBRARY, primers=(), read_lengths={'Read 1': 100, 'Read 2': 100}, acquired_roles=('Read 1', 'Read 2'), read_uses={'Read 1': 'include', 'Read 2': 'include'},limitation='Native full endpoint/loaded sequencing-primer geometry unassigned; source PE100 explicit.'))
    graph.declare_reviewed_scope(REVIEWED, source)
    return graph.require_publishable()

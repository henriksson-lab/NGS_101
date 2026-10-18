"""PANDORA source-owned native pretreatment and historical-kit boundaries."""
from chemdraw import Construct,Scene,feature
from protocol_model import molecule,role
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import ProtocolGraph,MaterialKind,InputRole,StepKind,SequencingPlan,Terminal,Disposition
import seqprimers as sp
TITLE='PANDORA-seq — pan-small-RNA sequencing'
NOTES='01_pandora_seq.html'
SOURCE='Shi et al., Nature Cell Biology(2021), <a href="https://doi.org/10.1038/s41556-021-00652-7">doi:10.1038/s41556-021-00652-7</a>; published PMC8236090 Methods RNA isolation/treatments/library preparation.'
SUMMARY='Pre-size-selected small RNA undergoes terminal repair and internal demethylation, then the source NEBNext E7330S adapter/RT/PCR route.'
CAVEAT='Selected PNK-first/AlkB-second route; the source also compares reversed treatment order. Initial15–50nt selection precedes treatment; library preparation reselects15–45nt. The exact oligos listed as2025 comparators belong to a later E7330 manual, not verified2021 molecular states. The Takara ligation-validation assay is a separate experiment. Historical kit geometry remains unassigned.'
SEQ_PRIMERS=()
READ_LENGTHS=None
SEQUENCING_INTRO='Source SE50 generally; sperm-head libraries useSE75. Sample-specific index acquisition/width and exact historical sequencing-primer sites remain unassigned.'
SEQUENCING_ENDING=SEQUENCING_INTRO
SEQUENCING_UNAVAILABLE='Native E7330S kit end/junction and primer-binding geometry have not been reconstructed for the2021 version.'
READ_SPECS=(sp.read_spec('Read 1','50;75 for sperm heads','small-RNA sequence after kit adapter; exact site unassigned','include'),sp.read_spec('Index 1 (i7)',None,'Acquisition and selected index species/width unassigned','undetermined'),sp.read_spec('Index 2 (i5)',None,'Presence/acquisition unassigned for historical kit branch','undetermined'))
FINAL_LIBRARY=Construct([role('E7330S flow-cell/read end; historical bases unassigned'),role('15–45nt small-RNA-derived DNA'),role('opposite E7330S end; sample index placement unassigned',molecular_feature=feature('pandora_sample','sample_index','unknown'))],name='PANDORA source-limited historical E7330S endpoint')
FINAL_CAPTION='Expected kit endpoint only; actual first-strand cDNA is the PCR input. No2025 index sequence or width is assigned.'
oligos = lambda: ('2025 comparator only — NEBNext E7330 3′ SR adapter — 5′-rAppAGATCGGAAGAGCACACGTCT-NH2-3′',
 '2025 comparator only — NEBNext E7330 5′ SR adapter — '
 '5′-rGrUrUrCrArGrArGrUrUrCrUrArCrArGrUrCrCrGrArCrGrArUrC-3′',
 '2025 comparator only — NEBNext E7330 SR RT primer — 5′-AGACGTGTGCTCTTCCGATCT-3′',
 '2025 comparator only — NEBNext E7330 SR PCR primer — '
 '5′-AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCG*A-3′',
 "2025 comparator only — NEBNext E7330 index 1 PCR primer (kit example; experiment's index choice unstated) — "
 '5′-CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC*T-3′')

def rna_scene(label,*,repaired=False,demethylated=False,adapters=0):
 sc=Scene();parts=[]
 if adapters==2:parts.append(role('2021 E7330S5′SR adapter; bases unassigned'))
 parts.append(role('small RNA; AlkB-treated, residual modification distribution unassigned' if demethylated else 'small RNA; internal modifications heterogeneous'))
 if adapters:parts.append(role('2021 E7330S3′SR adapter; bases/modifications unassigned'))
 sc.strand('RNA',parts,label=label,mod5='P' if repaired and adapters<2 else '',mod3='OH' if repaired and adapters==0 else '')
 if adapters:sc.junction('RNA',parts[-2].name,parts[-1].name,'RNA ligation')
 if adapters==2:sc.junction('RNA',parts[0].name,parts[1].name,'RNA ligation')
 return sc

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 total=g.input('total_rna','isolated total RNA',MaterialKind.BIOLOGICAL,state=rna_scene('native total RNA'),primary=True)
 size=g.input('initial_page','15%ureaPAGE;15–50nt excision,diffusion4°C overnight,ethanol precipitation',MaterialKind.SELECTION)
 rna=g.product('small_rna','15–50nt RNA fraction before treatment',MaterialKind.MOLECULAR,state=rna_scene('preselected15–50nt RNA'))
 excluded=g.product('initial_offsize','RNA outside initial selected gel window',MaterialKind.WASTE)
 g.react('initial_size_select','Select15–50nt RNA before enzymatic treatment',StepKind.SELECTION,((total,InputRole.TEMPLATE),(size,InputRole.SELECTION)),(rna,excluded));g.discard(excluded,'Not included in excised source gel fraction; reduces intact tRNA/rRNA substrates for degradation artifacts.')
 pnk=g.input('pnk','T4PNK10U,1mMATP,50µL;37°C20min',MaterialKind.CATALYST)
 repaired=g.product('repaired_rna','PNK-treated RNA targets',MaterialKind.MOLECULAR,state=rna_scene('PNK-treated RNA',repaired=True))
 g.react('end_repair','Repair3′P/cyclic-P and add5′P',StepKind.END_REPAIR,((rna,InputRole.TEMPLATE),(pnk,InputRole.CATALYST)),(repaired,))
 recovery=g.input('pnk_recovery','Add500µLTRIzol;source RNA-isolation procedure',MaterialKind.SELECTION)
 clean=g.product('pnk_clean','RNA recovered afterPNK',MaterialKind.MOLECULAR,state=repaired.state)
 loss=g.product('pnk_waste','TRIzol phases/precipitation washes',MaterialKind.WASTE)
 g.react('recover_pnk','Isolate RNA between the selected treatments',StepKind.EXTRACTION,((repaired,InputRole.TEMPLATE),(recovery,InputRole.SELECTION)),(clean,loss));g.discard(loss,'Source RNA extraction phase/supernatant transfers; molecular partition unspecified.')
 alkb=g.input('alkb','AlkB4µg/ml,Fe/αKG/ascorbate,200ngRNA/50µL;37°C30min',MaterialKind.CATALYST)
 demeth=g.product('demethylated_rna','PNK/AlkB-treated RNA',MaterialKind.MOLECULAR,state=rna_scene('PNK/AlkB-treated RNA',repaired=True,demethylated=True))
 g.react('demethylate','Reduce internal methylation barriers withAlkB',StepKind.CHEMICAL_CONVERSION,((clean,InputRole.TEMPLATE),(alkb,InputRole.CATALYST)),(demeth,))
 ar=g.input('alkb_recovery','Add500µLTRIzol;RNA-isolation procedure',MaterialKind.SELECTION)
 recovered=g.product('treated_clean','recovered doubly treated RNA',MaterialKind.MOLECULAR,state=demeth.state)
 aw=g.product('alkb_waste','postAlkB extraction/precipitation fractions',MaterialKind.WASTE)
 g.react('recover_alkb','Isolate RNA afterAlkB treatment',StepKind.EXTRACTION,((demeth,InputRole.TEMPLATE),(ar,InputRole.SELECTION)),(recovered,aw));g.discard(aw,'Physically excluded source extraction/precipitation fractions.')
 second_size=g.input('library_page','Library Methods:PAGE15–45nt stripe,recoverRNA',MaterialKind.SELECTION)
 selected=g.product('library_rna','15–45nt treated RNA for library preparation',MaterialKind.MOLECULAR,state=demeth.state)
 off=g.product('library_offsize','RNA outside15–45nt excision',MaterialKind.WASTE)
 g.react('library_size_select','Reselect15–45nt RNA for source library preparation',StepKind.SELECTION,((recovered,InputRole.TEMPLATE),(second_size,InputRole.SELECTION)),(selected,off));g.discard(off,'Not present in excised library-preparation gel stripe.')
 a3=g.input('three_adapter','2021 NEBNextE7330S3′SR adapter; historical physical oligo unassigned',MaterialKind.ADAPTER)
 l3=g.input('three_ligation_system','Source3′adapter system;70°C2min,25°C1h or sperm16°C18h',MaterialKind.CATALYST)
 lig3=g.product('three_ligated','3′adapter-ligated RNA',MaterialKind.MOLECULAR,state=rna_scene('3′adapted RNA',repaired=True,demethylated=True,adapters=1))
 g.react('three_ligation','Ligate source3′adapter before RT-primer annealing',StepKind.LIGATION,((selected,InputRole.TEMPLATE),(a3,InputRole.STRUCTURAL),(l3,InputRole.CATALYST)),(lig3,),geometry=AdapterLigationGeometry.unavailable('Historical2021 E7330S adapter end/modification geometry is not verified; later2025 manual and Takara validation assay are not substituted.'))
 primer=g.input('rt_primer','HistoricalE7330S RT primer;physical sequence unassigned',MaterialKind.OLIGO)
 anneal=g.input('anneal_conditions','75°C5min,37°C15min,15°C15min',MaterialKind.CONDITION)
 annealed=rna_scene('RNA/RT-primer mixture',repaired=True,demethylated=True,adapters=1);annealed.strand('RTprimer',[role('separate source RT primer; binding geometry unassigned')],label='RNA-bound/free primer population unresolved')
 primed=g.product('primer_annealed','3′adapted RNA plus RT primer before5′ligation',MaterialKind.MOLECULAR,state=annealed)
 g.react('anneal_rt_primer','Add and anneal RT primer before5′adapter mix',StepKind.HYBRIDIZATION,((lig3,InputRole.TEMPLATE),(primer,InputRole.STRUCTURAL),(anneal,InputRole.CONDITION)),(primed,))
 a5=g.input('five_adapter','2021 E7330S5′SR adapter mix;physical sequence unassigned',MaterialKind.ADAPTER)
 l5=g.input('five_ligation_system','Source5′adapter mix;70°C2min,25°C1h',MaterialKind.CATALYST)
 adapted=rna_scene('two-adapter RNA plus RT primer',repaired=True,demethylated=True,adapters=2);adapted.strand('RTprimer',[role('separate source RT primer; sequence unassigned')],label='independent RT primer')
 dual=g.product('adapted_rna','two-adapter RNA/RT-primer population',MaterialKind.MOLECULAR,state=adapted)
 g.react('five_ligation','Ligate source5′adapter after RT-primer addition',StepKind.LIGATION,((primed,InputRole.TEMPLATE),(a5,InputRole.STRUCTURAL),(l5,InputRole.CATALYST)),(dual,),geometry=AdapterLigationGeometry.unavailable('Historical kit5′adapter physical sequence/terminal chemistry and exact substrate pairing remain unassigned.'))
 rt=g.input('reverse_transcriptase','Source first-strand system;70°C2min,50°C1h',MaterialKind.CATALYST)
 ss=molecule('actual ss first-strand target',(role('source RT primer retained in cDNA'),role('copied3′adapter span beyond primer, extent unassigned'),role('small-RNA-derived cDNA'),role('copied5′adapter')))
 cdna=g.product('first_strand','actual source first-strand cDNA target',MaterialKind.MOLECULAR,state=ss)
 residual=g.product('rt_rna','RNA templates after RT',MaterialKind.MOLECULAR,state=rna_scene('retained RNA templates',repaired=True,demethylated=True,adapters=2))
 g.react('reverse_transcribe','Reverse-transcribe adapted small RNA',StepKind.REVERSE_TRANSCRIPTION,((dual,InputRole.TEMPLATE),(primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,residual));g.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,'RT does not consume RNA; the defining Methods does not assign a template-removal reaction here. Graph follows cDNA targets into PCR; RNA recovery is unresolved.'))
 f=g.input('pcr_forward','2021 PCR Primer Cocktail forward species unassigned',MaterialKind.OLIGO)
 r=g.input('pcr_reverse','2021 cocktail indexed primer species/width unassigned',MaterialKind.OLIGO)
 pol=g.input('pcr_mix','PCR MasterMix;11–22cycles94/62/70°C',MaterialKind.CATALYST)
 amp=g.product('amplified','amplified historical kit endpoint',MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
 g.react('library_pcr','Amplify actual first-strand cDNA target population',StepKind.PCR,((cdna,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(cdna.state,reason='Actual ss first-strand target, not completed library. Historical2021 cocktail/primer sites and native first-copy geometry remain unassigned.'))
 page=g.input('final_page','Source postPCR PAGE purification;exact selected window unassigned',MaterialKind.SELECTION)
 lib=g.product('library','PAGE-purified PANDORA library',MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 waste=g.product('final_offsize','postPCR fractions outside selected gel band',MaterialKind.WASTE)
 g.react('purify_library','PAGE-purify source amplified library',StepKind.SELECTION,((amp,InputRole.TEMPLATE),(page,InputRole.SELECTION)),(lib,waste));g.discard(waste,'Outside physically excised postPCR gel band.')
 g.sequence(lib,SequencingPlan(library=FINAL_LIBRARY,acquired_roles=('Read 1',),limitation=SEQUENCING_ENDING))
 g.declare_reviewed_scope(('total_rna','initial_page','initial_size_select','pnk','end_repair','pnk_recovery','recover_pnk','alkb','demethylate','alkb_recovery','recover_alkb','library_page','library_size_select','three_adapter','three_ligation_system','three_ligation','rt_primer','anneal_conditions','anneal_rt_primer','five_adapter','five_ligation_system','five_ligation','reverse_transcriptase','reverse_transcribe','pcr_forward','pcr_reverse','pcr_mix','library_pcr','final_page','purify_library','library'),SOURCE)
 return g.require_publishable()

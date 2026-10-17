"""CoBATCH molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct,MolecularState,Row,Scene,Workflow,duplex_rows,feature,revcomp
from chromatin_epigenetics import antibody_enzyme_rows
from tagmentation import loaded_tn5_adapter
from pcr import PCRGeometry
from protocol_graph import Disposition, Terminal
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind)

PUBLISHED=("Wang et al. 2019 Mol Cell STAR Methods (mmc7.pdf pp. e4-e5: 'Fixed single-cell CoBATCH', "
           "'Nextera library preparation and sequencing strategy') and Methods S1 bench protocol (mmc6.pdf steps 11-27)")
METHODS_S1="Wang et al. 2019 Methods S1 (mmc6.pdf) steps 23-27"
ASSEMBLY="Wang et al. 2019 STAR Methods 'Quality control for the assembled PAT transposome' (mmc7.pdf p. e4)"
TITLE="CoBATCH — combinatorial barcoding by targeted tagmentation"
NOTES="01_cobatch.html"
SOURCE='Defining source: <a href="https://doi.org/10.1016/j.molcel.2019.07.015">Wang et al., <i>Molecular Cell</i> (2019)</a>.'
SUMMARY="Antibody-tethered Protein-A–Tn5 transfers a pair of well barcodes at the chromatin target; pooled cells are redistributed and receive a second i5/i7 well identity during PCR."
CAVEAT="CoBATCH cell identity has four indexed parts: paired T5/T7 transposome barcodes and paired PCR-well indexes. The original published Table S1 confirms the displayed adapter/PCR families. The published STAR Methods and Methods S1 (step 27) specify paired-end 150-base reads on HiSeq X Ten or NovaSeq 6000; the authors' GEO sample record GSM3705255 states that the first 21 bases of Read 1 and the first 27 bases of Read 2 are cut for single-cell splitting, so the transposome barcodes lie in those cycles. Neither the published article nor Methods S1 prints the loaded sequencing primers or index-read lengths, so read entry points are not established in this model. Lysis, KAPA gap-fill/PCR, well pooling with column purification and AMPure 200-1000 bp selection follow the published STAR Methods and Methods S1."
def cb(identifier,part): return feature(identifier,"cell_barcode","combinatorial",group="cobatch_cell",part=part,whitelist="published CoBATCH barcode plate")
T5_OUTER="TCCACGC"
T5_INNER="GCGATCGAGGACGGC"
T7_OUTER="CTGTCCCTGTCC"
T7_INNER="CACCGTCTCCGCCTC"
T5_BARCODES=("TATAGCCT","ATAGAGGC","CCTATCCT","GGCTCTGA")
T7_BARCODES=("CGAGTAAT","TCTCCGGA","AATGAGCG","GGAATCTC")
TN5_ME_REV=nx.ME_RC

def t5_adapter(barcode=T5_BARCODES[0]):
 if barcode not in T5_BARCODES: raise ValueError("representative T5 barcode must be published")
 return [seg("S5",nx.S5,"s5"),seg("T5 outer connector",T5_OUTER),
         seg("T5 well barcode",barcode,"cbc",feature=cb("cobatch_t5","PAT-T5 well")),
         seg("T5 inner connector",T5_INNER),seg("T5 mosaic end",nx.ME,"me")]

def t7_adapter(barcode=T7_BARCODES[1]):
 if barcode not in T7_BARCODES: raise ValueError("representative T7 barcode must be published")
 return [seg("S7",nx.S7,"s7"),seg("T7 outer connector",T7_OUTER),
         seg("T7 well barcode",barcode,"cbc",feature=cb("cobatch_t7","PAT-T7 well")),
         seg("T7 inner connector",T7_INNER),seg("T7 mosaic end",nx.ME,"me")]

def loaded_scene(parts,label):
 sc=loaded_tn5_adapter(parts,mosaic_end=parts[-1].name,label=label,
                       bottom_name="common annealing primer",
                       bottom_label="5′-phosphorylated common ME bottom",
                       bottom_sequence=TN5_ME_REV,mod5="p")
 bc=next(s.name for s in parts if s.feature)
 sc.mark("Tn5 top",bc,"first-plate well identity")
 return sc

TABLE_S1="Wang et al. 2019 Table S1 (mmc2.xlsx)"

def transposome_strands():
 """Ordering strands of the representative loaded complexes, 5'->3' (sheet PA-Tn5_barcodes)."""
 return ("".join(s.top for s in t5_adapter()),"".join(s.top for s in t7_adapter()),TN5_ME_REV)

def transferred_fragment():
 right=t7_adapter()
 right_rc=[seg(s.name+" reverse complement",revcomp(s.top),s.tag,feature=s.feature)
           for s in reversed(right)]
 return Construct([*t5_adapter(),seg("targeted chromatin fragment","X"*42,placeholder=True),
                   *right_rc],name="gap-filled target-coded fragment")

TRANSFERRED=transferred_fragment()
FINAL_LIBRARY=Construct([
 seg("P5",il.P5,"p5"),seg("PCR-well i5","J"*8,"cbc",placeholder=True,feature=cb("cobatch_i5","PCR well i5")),seg("S5",nx.S5,"s5"),
 seg("T5 outer connector",T5_OUTER),seg("T5 well barcode",T5_BARCODES[0],"cbc",feature=cb("cobatch_t5","PAT-T5 well")),seg("T5 inner connector",T5_INNER),seg("left mosaic end",nx.ME,"me"),
 seg("targeted chromatin fragment","X"*42,placeholder=True),seg("right mosaic end reverse complement",nx.ME_RC,"me"),
 seg("T7 inner connector reverse complement",revcomp(T7_INNER)),seg("T7 well barcode reverse complement",revcomp(T7_BARCODES[1]),"cbc",feature=cb("cobatch_t7","PAT-T7 well")),seg("T7 outer connector reverse complement",revcomp(T7_OUTER)),seg("S7 reverse complement",nx.S7_RC,"s7"),
 seg("PCR-well i7 reverse complement","I"*8,"cbc",placeholder=True,feature=cb("cobatch_i7","PCR well i7")),seg("P7 reverse complement",il.P7_RC,"p7")],name="CoBATCH library")
# Published Table S1 prints adapters and PCR oligos. Adapter subsequences
# alone do not establish which sequencing primers were physically loaded.
SEQ_PRIMERS=()
FINAL_CAPTION="The first plate is encoded beside both mosaic ends; the second plate is read in i5 and i7. Their combination identifies a cell."
SEQUENCING_INTRO="Paired genomic reads map targeted chromatin fragments; index and inline cycles recover the two-round cell identity."
SEQUENCING_UNAVAILABLE="Published Table S1 establishes the adapter and i5/i7 PCR primer families. The published STAR Methods and Methods S1 do not print the loaded Read 1/Read 2/index sequencing oligos or index-read lengths, so their entry positions are not modelled. This is a bounded-pass observation, not an accepted source-absence claim."
SEQUENCING_ENDING="The published Methods S1 step 27 and STAR Methods specify paired-end 150-base reads on HiSeq X Ten or NovaSeq 6000. The completed duplex below is assembled from the original published adapter/PCR table; sequencing-entry geometry remains unresolved."
# Author-deposited GEO sample GSM3705255 (EC_CoBATCH sub-series GSE129327),
# Sample_data_processing: "For spliting single cell, we cut the first 21 bp of
# the read1 and the first 27 bp of the read2." Cycle counts only; read entry
# points are not stated there.
GEO_SPLIT_CYCLES={"Read 1":21,"Read 2":27}
def split_overrun(role):
 """Bases of the inner connector covered by the GEO split window IF the read began
 immediately 3' of S5 (Read 1) or S7 (Read 2); an inference, not a declared primer."""
 outer,barcode=(T5_OUTER,T5_BARCODES[0]) if role=="Read 1" else (T7_OUTER,T7_BARCODES[0])
 return GEO_SPLIT_CYCLES[role]-len(outer)-len(barcode)
READ_SPECS=(sp.read_spec("Read 1",150,"chromatin read; first 21 cycles used for cell splitting (GEO GSM3705255); exact entry unresolved","extract"),sp.read_spec("Read 2",150,"paired chromatin read; first 27 cycles used for cell splitting (GEO GSM3705255); exact entry unresolved","extract"))
READ_LENGTHS={"Read 1":150,"Read 2":150}
SHOW_READ_CONTENT=True

def oligos(): return (
 f"<b>T5-1</b> — 5′-{''.join(s.top for s in t5_adapter())}-3′",
 f"<b>T7-2</b> — 5′-{''.join(s.top for s in t7_adapter())}-3′",
 f"<b>common annealing primer</b> — 5′-/5Phos/{TN5_ME_REV}-3′",
 f"<b>representative i5 PCR primer</b> — 5′-{il.P5}[8-nt PCR-well i5]{nx.S5}-3′",
 f"<b>representative i7 PCR primer</b> — 5′-{il.P7}[8-nt PCR-well i7]{nx.S7}-3′",
)

def workflow():
 start=Scene.duplex([seg("native chromatin","X"*42,placeholder=True)],label="antibody-bound chromatin")
 wf=Workflow(MolecularState("Antibody-bound chromatin in first-plate wells",tuple(start.rows())))
 wf.react("Anneal each T5 or T7 oligo with the MErev bottom (95 °C 5 min, 0.1 °C/s to 75 °C, cool to RT) and load pA–Tn5 (25 °C, 60 min)",
          [*loaded_scene(t5_adapter(),"representative PAT-T5").rows(),
           *loaded_scene(t7_adapter(),"representative PAT-T7").rows()],
          name="First-plate barcoded pA–Tn5",
          note="The representative T5-1 and T7-2 sequences are drawn; the source table supplies the full whitelist.")
 wf.react("Tether through antibody and activate tagmentation (staggered gaps are filled later, at 72 °C in the PCR well)",
          duplex_rows(TRANSFERRED,label="T5/T7 target-coded fragment"),
          name="First-plate-indexed chromatin fragments")
 wf.react("Pool cells, sort 20-25 per second-plate well, lyse with proteinase K (65 °C 5-6 h, 85 °C 15 min), gap-fill and PCR with second-well i5/i7 primers, then pool wells for column and AMPure cleanup",
          duplex_rows(FINAL_LIBRARY,label="four-part CoBATCH library"),
          name="CoBATCH sequencing library",
          note="The two transposome barcodes plus the two PCR-well indexes jointly identify one cell.")
 return wf

def protocol_graph():
 source="Defining source: https://doi.org/10.1016/j.molcel.2019.07.015"
 graph=ProtocolGraph("CoBATCH",authority=source)
 chromatin=graph.input("chromatin","antibody-bound chromatin",MaterialKind.BIOLOGICAL,
                       primary=True,state=Scene.duplex([seg("native chromatin","X"*42,placeholder=True)]))
 transposomes=graph.input("barcoded_transposomes","T5/T7 first-plate pA-Tn5 complexes "
                          "(T5-1 and T7-2 top strands with the common annealing primer drawn; annealed 95 °C 5 min, "
                          "0.1 °C/s to 75 °C, cooled to RT; 37.5 µM adaptor + 37.5 µM PAT, 25 °C 60 min)",
                          MaterialKind.ADAPTER,state=transposome_strands(),source=TABLE_S1+"; "+ASSEMBLY)
 tn5=graph.input("pa_tn5","protein-A–Tn5 transposase",MaterialKind.CATALYST)
 tagged=graph.product("tagged_fragments","gapped first-plate-indexed chromatin fragments",
                      MaterialKind.MOLECULAR,state=TRANSFERRED)
 tagmentation_waste=graph.product("tagmentation_waste","candidate non-library tagmentation products",
                                  MaterialKind.MOLECULAR)
 graph.react("targeted_tagmentation","Tether pA-Tn5 and transfer the T5/T7 well adapters",
             StepKind.TAGMENTATION,
             ((chromatin,InputRole.TEMPLATE),(transposomes,InputRole.STRUCTURAL),
              (tn5,InputRole.CATALYST)),(tagged,tagmentation_waste))
 graph.mark_terminal(Terminal(tagmentation_waste.id, Disposition.UNRESOLVED, 'Candidate non-library end classes are carried into PCR. Their end arrangement does not establish whether they amplify, and no class-specific physical separation is represented.'))
 allocation=graph.input("second_plate_allocation","pool, DAPI-stain, strain and sort 20-25 cells per second-plate well",
                        MaterialKind.SELECTION,source=PUBLISHED)
 allocated=graph.product("allocated_fragments","tagged chromatin assigned to second-plate wells",
                         MaterialKind.MOLECULAR,state=TRANSFERRED)
 excluded=graph.product("allocation_loss","cells not sorted into an assayed second-plate well",
                        MaterialKind.WASTE)
 graph.react("redistribute","EDTA-stop, pool cells and sort them into second-plate wells",
             StepKind.SAMPLING,
             ((tagged,InputRole.TEMPLATE),(allocation,InputRole.SELECTION)),(allocated,excluded),
             changes_molecule=False,source=PUBLISHED)
 graph.discard(excluded,"material not sorted into an assayed second-plate well is not sequenced")
 lysis=graph.input("proteinase_k_lysis","4 µl SDS/proteinase K lysis, 65 °C 5-6 h reverse cross-linking, 85 °C 15 min inactivation, 1 µl 1.8% Triton X-100 quench at 37 °C 1 h",
                   MaterialKind.CONDITION,source=PUBLISHED)
 released=graph.product("released_fragments","released gapped first-plate-indexed fragments",
                        MaterialKind.MOLECULAR,state=TRANSFERRED)
 graph.react("release","Lyse, reverse cross-links and release tagmented DNA in the PCR well",
             StepKind.DEPROTEINIZATION,
             ((allocated,InputRole.TEMPLATE),(lysis,InputRole.CONDITION)),(released,),
             changes_molecule=False,source=PUBLISHED)
 polymerase=graph.input("polymerase","2X KAPA master mix (72 °C gap fill, then PCR)",MaterialKind.CATALYST,
                        source=PUBLISHED)
 filled=graph.product("gap_filled_fragments","gap-filled first-plate-indexed fragments",
                      MaterialKind.MOLECULAR,state=TRANSFERRED)
 graph.react("gap_fill","Fill the staggered transposase gaps (72 °C, 5 min in the PCR mix)",StepKind.EXTENSION,
             ((released,InputRole.TEMPLATE),(polymerase,InputRole.CATALYST)),(filled,),source=PUBLISHED)
 i5=graph.input("i5_primer","PCR-well i5 primer (Nextera i5-N5xx family)",MaterialKind.OLIGO,
                state=Construct([seg("P5",il.P5),seg("PCR-well i5","N"*8,feature=cb("cobatch_i5","PCR well i5")),
                                 seg("S5",nx.S5)]),source=TABLE_S1+", sheet Nextera_index")
 i7=graph.input("i7_primer","PCR-well i7 primer (Nextera i7-N7xx family)",MaterialKind.OLIGO,
                state=Construct([seg("P7",il.P7),seg("PCR-well i7","N"*8,feature=cb("cobatch_i7","PCR well i7")),
                                 seg("S7",nx.S7)]),source=TABLE_S1+", sheet Nextera_index")
 amplified=graph.product("amplified_library",FINAL_LIBRARY.name+" (unpurified PCR)",MaterialKind.MOLECULAR,
                         state=FINAL_LIBRARY)
 graph.react("second_well_pcr","Amplify 98 °C 45 s, then 10-18 cycles of 98 °C 15 s / 63 °C 30 s / 72 °C 1 min, final 72 °C 5 min, adding the second-well i5/i7 identities",StepKind.PCR,
             ((filled,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),
              (i7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(amplified,),
             geometry=PCRGeometry.unavailable(
                 FINAL_LIBRARY,reason="The source establishes the indexed PCR ends, but the page does not model exact primer annealing geometry."),
             source=PUBLISHED)
 column=graph.input("tiangen_column","combine all wells and purify on a TIANGEN DNA purification column (50 µl elution)",
                    MaterialKind.SELECTION,source=METHODS_S1)
 pooled=graph.product("pooled_library",FINAL_LIBRARY.name+" (pooled, column-purified)",MaterialKind.MOLECULAR,
                      state=FINAL_LIBRARY)
 column_waste=graph.product("column_waste","PCR-reaction components not retained by the column",MaterialKind.WASTE)
 graph.react("pool_purify","Combine the second-plate PCR wells and column-purify",StepKind.CLEANUP,
             ((amplified,InputRole.TEMPLATE),(column,InputRole.SELECTION)),(pooled,column_waste),
             changes_molecule=False,source=METHODS_S1)
 graph.discard(column_waste,"discarded in the column flow-through and washes")
 beads=graph.input("ampure","1X AMPure XP cleanup, then 0.5X/0.5X double-sided selection",
                   MaterialKind.SELECTION,source=PUBLISHED)
 library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 rejected=graph.product("rejected_sizes","PCR-reaction components and fragments outside 200-1000 bp",
                        MaterialKind.WASTE)
 graph.react("size_select","AMPure XP cleanup and 200-1000 bp size selection",StepKind.SELECTION,
             ((pooled,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(library,rejected),
             source=PUBLISHED)
 graph.discard(rejected,"removed by bead cleanup and double-sided size selection")
 graph.sequence(library,SequencingPlan(
     library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
     acquired_roles=("Read 1","Read 2"),read_lengths=READ_LENGTHS,
     read_uses={"Read 1":"include","Read 2":"include"},
     limitation=SEQUENCING_UNAVAILABLE))
 graph.annotate(SourceDiscrepancy(
     "barcoded_transposomes","Published Table S1 equivalence notes are offset by one row",
     (SourceClaim("Table S1 (mmc2.xlsx) sheet Mosaic-A&B, cells B3-B5",
                  "Mosaic-A TCGTCGGCAGCGTC+ME, Mosaic-B GTCTCGTGGGCTCGG+ME, common annealing primer 5Phos/CTGTCTCTTATACACATCT"),
      SourceClaim("Same sheet, cells C4-C6",
                  "'equally to Tn5ME-A_FC-121-1030' beside Mosaic-B, 'equally to Tn5ME-B_FC-121-1031' beside the common annealing primer, 'equally to TnExtMErev_primer' on an otherwise empty row")),
     "The model uses the printed sequences: the common annealing primer is the phosphorylated ME reverse strand (the T5/T7 rows say they anneal with TnExtMErev_primer).",
     "Read from the original publisher workbook cells; by sequence, Mosaic-A carries the Tn5ME-A/S5 arm and Mosaic-B the Tn5ME-B/S7 arm.",
     "The offset is in the published sheet itself; the STAR Methods only refer to Table S1 and do not restate the equivalences. The Mosaic-A/B oligos are not drawn on this single-cell page."))
 graph.declare_reviewed_scope((
     "chromatin","barcoded_transposomes","pa_tn5","targeted_tagmentation",
     "second_plate_allocation","redistribute","proteinase_k_lysis","release",
     "polymerase","gap_fill","i5_primer","i7_primer","second_well_pcr",
     "amplified_library","tiangen_column","pool_purify","pooled_library","ampure","size_select","library"),source)
 return graph.require_publishable()

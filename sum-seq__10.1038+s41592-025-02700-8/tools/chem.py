"""SUM-seq sample pre-indexing and overloaded droplet RNA/ATAC profiling."""
from multibranch import illumina_branch
from protocol_extensions import analyte_split, indexed_molecule
from chemdraw import Row
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)
TITLE="SUM-seq — sample-indexed, overloaded single-cell RNA and ATAC"
NOTES="01_sum-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41592-025-02700-8">Llorens-Rico et al., <i>Nature Methods</i> (2025)</a>.'
SUMMARY="ATAC fragments and RNA are separately sample-indexed in fixed nuclei, then a second 10x droplet barcode resolves overloaded nuclei; the sample-index–droplet-index pair identifies each cell across both modalities."
CELL_PARTS=(("sample index",8),("droplet barcode",16))
RNA,RP=illumina_branch("SUM-seq RNA","RNA-derived cDNA",cell_parts=CELL_PARTS,umi=10)
ATAC,AP=illumina_branch("SUM-seq ATAC","accessible genomic DNA",architecture="nextera",cell_parts=CELL_PARTS)
FINAL_LIBRARIES=(("RNA library",RNA,RP,"RT sample index plus droplet barcode identifies RNA molecules.","Illumina primers on the RNA branch."),("ATAC library",ATAC,AP,"Barcoded-Tn5 sample index plus the same droplet barcode identifies ATAC fragments.","Nextera primers on the ATAC branch."))
def sections():
    return [("Install modality-specific sample indexes",[Row(chunks=[("RNA: oligo(dT) RT ** sample index", "cbc", False)]),Row(chunks=[("ATAC: accessible DNA ** barcoded Tn5 sample index", "cbc", False)])],"Each nucleus receives matching sample provenance on its two analytes."),("Tagment RNA/cDNA hybrids",[Row(chunks=[("RNA–cDNA hybrid → Tn5 primer site", "me", False)])],"This makes the RNA branch compatible with the later microfluidic barcode transfer."),("Overload nuclei into 10x droplets",indexed_molecule(payload="sample-indexed RNA or ATAC molecule",barcode_parts=(("sample index",8),("droplet barcode",16),),umi=10).rows(),"The second barcode resolves multiple nuclei sharing a droplet when combined with sample index."),("Break droplets and split the library",analyte_split("RNA-specific amplification","ATAC-specific amplification"),"Branch-specific amplification retains the shared sample × droplet identity."),]

def protocol_graph():
    source="Defining source: https://doi.org/10.1038/s41592-025-02700-8"
    graph=ProtocolGraph("SUM-seq",authority=source)
    nuclei=graph.input("nuclei","fixed nuclei containing RNA and chromatin",
                       MaterialKind.BIOLOGICAL,primary=True)
    rna_analyte=graph.product("rna_analyte","nuclear RNA branch",MaterialKind.MOLECULAR)
    atac_analyte=graph.product("atac_analyte","nuclear chromatin branch",MaterialKind.MOLECULAR)
    graph.react("modality_split","Expose RNA and chromatin analytes",StepKind.SPLIT,
                ((nuclei,InputRole.TEMPLATE),),(rna_analyte,atac_analyte),changes_molecule=False)
    rt_primer=graph.input("indexed_rt_primer","sample-indexed oligo(dT) RT primer",
                          MaterialKind.OLIGO)
    rt=graph.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
    indexed_rna=graph.product("indexed_rna_cdna","sample-indexed RNA–cDNA hybrids",
                              MaterialKind.MOLECULAR)
    graph.react("indexed_rt","Reverse-transcribe RNA with the sample-indexed oligo(dT)",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna_analyte,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(indexed_rna,))
    atac_adapter=graph.input("atac_transposome","sample-indexed ATAC transposome",
                             MaterialKind.ADAPTER)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    indexed_atac=graph.product("indexed_atac","sample-indexed accessible-DNA fragments",
                               MaterialKind.MOLECULAR)
    atac_waste=graph.product("atac_waste","non-library ATAC tagmentation products",
                             MaterialKind.WASTE)
    graph.react("atac_tagmentation","Tagment chromatin with sample-indexed Tn5",
                StepKind.TAGMENTATION,
                ((atac_analyte,InputRole.TEMPLATE),(atac_adapter,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(indexed_atac,atac_waste))
    graph.discard(atac_waste,"molecules without usable paired transposase insertions do not amplify")
    rna_adapter=graph.input("rna_transposome","RNA-branch Tn5 primer adapter",
                            MaterialKind.ADAPTER)
    tagged_rna=graph.product("tagged_rna_cdna","tagmented sample-indexed RNA/cDNA hybrids",
                             MaterialKind.MOLECULAR)
    rna_waste=graph.product("rna_tagmentation_waste","non-library RNA-branch products",
                            MaterialKind.WASTE)
    graph.react("rna_tagmentation","Tagment the RNA/cDNA hybrids",StepKind.TAGMENTATION,
                ((indexed_rna,InputRole.TEMPLATE),(rna_adapter,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged_rna,rna_waste))
    graph.discard(rna_waste,"RNA/cDNA molecules without usable primer sites do not amplify")
    combined=graph.product("combined_modalities","sample-indexed RNA and ATAC molecules",
                           MaterialKind.MOLECULAR)
    graph.react("pool_modalities","Combine the indexed analytes for droplet loading",
                StepKind.POOL,
                ((tagged_rna,InputRole.TEMPLATE),(indexed_atac,InputRole.TEMPLATE)),(combined,),
                changes_molecule=False)
    bead=graph.input("droplet_bead_oligo","10x droplet barcode oligo",MaterialKind.ADAPTER)
    droplet_products=graph.product("droplet_products","sample × droplet indexed molecules",
                                   MaterialKind.MOLECULAR)
    graph.react("droplet_barcoding","Overload nuclei and transfer the droplet barcode",
                StepKind.ADAPTER_ATTACHMENT,
                ((combined,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL)),(droplet_products,))
    rna_precursor=graph.product("rna_precursor","RNA-specific amplification branch",
                                MaterialKind.MOLECULAR)
    atac_precursor=graph.product("atac_precursor","ATAC-specific amplification branch",
                                 MaterialKind.MOLECULAR)
    graph.react("library_split","Break droplets and split by modality",StepKind.SPLIT,
                ((droplet_products,InputRole.TEMPLATE),),(rna_precursor,atac_precursor),
                changes_molecule=False)
    pcr=graph.input("pcr_polymerase","branch-specific PCR polymerase",MaterialKind.CATALYST)
    rna_primers=graph.input("rna_pcr_primers","RNA branch PCR-primer pair",MaterialKind.OLIGO)
    atac_primers=graph.input("atac_pcr_primers","ATAC branch PCR-primer pair",MaterialKind.OLIGO)
    rna_library=graph.product("rna_library",RNA.name,MaterialKind.LIBRARY,state=RNA)
    atac_library=graph.product("atac_library",ATAC.name,MaterialKind.LIBRARY,state=ATAC)
    graph.react("rna_pcr","Amplify the RNA library",StepKind.PCR,
                ((rna_precursor,InputRole.TEMPLATE),(rna_primers,InputRole.STRUCTURAL),
                 (pcr,InputRole.CATALYST)),(rna_library,),
                geometry=PCRGeometry.unavailable(RNA,reason="The branch helper represents a complete library rather than separate PCR landing objects.",symmetric_species=True))
    graph.react("atac_pcr","Amplify the ATAC library",StepKind.PCR,
                ((atac_precursor,InputRole.TEMPLATE),(atac_primers,InputRole.STRUCTURAL),
                 (pcr,InputRole.CATALYST)),(atac_library,),
                geometry=PCRGeometry.unavailable(ATAC,reason="The branch helper represents a complete library rather than separate PCR landing objects.",symmetric_species=True))
    for material,library,primers in ((rna_library,RNA,RP),(atac_library,ATAC,AP)):
        graph.sequence(material,SequencingPlan(
            library=library,primers=primers,acquired_roles=tuple(p.role for p in primers),
            read_uses={p.role:("demultiplex" if "Index" in p.role else "include") for p in primers}))
    graph.declare_reviewed_scope((
        "nuclei","modality_split","indexed_rt_primer","reverse_transcriptase","indexed_rt",
        "atac_transposome","tn5","atac_tagmentation","rna_transposome","rna_tagmentation",
        "pool_modalities","droplet_bead_oligo","droplet_barcoding","library_split",
        "pcr_polymerase","rna_pcr_primers","atac_pcr_primers","rna_pcr","atac_pcr",
        "rna_library","atac_library"),source)
    return graph.require_publishable()

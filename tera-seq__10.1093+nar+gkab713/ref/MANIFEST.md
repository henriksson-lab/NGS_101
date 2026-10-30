# Source retrieval

Run `python3 tools/get_sources.py --slugs tera-seq__10.1093+nar+gkab713`. Originals and extracted twins live under `/husky/henriksson/NGS101_papers/sources/tera-seq__10.1093+nar+gkab713/`, outside the checkout.

- nar_gkab713_PMC8599856.xml — Published final article PMC8599856; inspected Methods59–73 and Results147/163: RNApreparation, sixend-selectionroutes, kit/acquisition and actual retainedfractions.
- nar_gkab713_supp_gkab713_supplemental_files_Supplementary_Table_1-Updated3.xlsx — Published SupplementaryTable1; inspected Whole adapter table; raw XLSX underlined rich-text runs examined for RTA bottom binding arm.
- nar_gkab713_supp_gkab713_supplemental_files_Supplement_TERA-Seq.pdf — Published supplementary information; inspected Identity/contents and SupplementaryNote firsttwo pages only; adapter-detection thresholds are computational, not chemistry.

Source checklist and cap are recorded in `catalogue/source_reviews/tera-seq__10.1093+nar+gkab713.json`. Filenames alone are not provenance.

# Source manifest

Originals live outside the checkout at `/husky/henriksson/NGS101_papers/sources/cyphyrna-seq__10.1080+15476286.2021.1999105`. Re-obtain using `python3 tools/get_sources.py --slugs cyphyrna-seq__10.1080+15476286.2021.1999105`, then run `tools/scrape_primers.py` on that directory.

Inspected originals:

- [15476286.2021.1999105_PMC8782182.xml](https://pmc-oa-opendata.s3.amazonaws.com/PMC8782182.1/PMC8782182.1.xml): Published Version of Record; Fig1/overview87–119, Methods185–247, selected branch preparative/library operations219–247.
- [15476286.2021.1999105_supp_KRNB_A_1999105_SM4973_20211015_Supplementary_Data_RNA_Biology_final.docx](https://pmc-oa-opendata.s3.amazonaws.com/PMC8782182.1/KRNB_A_1999105_SM4973.zip#Supplementary information/20211015_Supplementary Data_RNA Biology_final.docx): Published SupplementaryData TableS5 active oligos; original WordXML cells inspected, table5 row46/globalrow73 reversePCR primer and adjacent Index containing six bases. Original TableS5 chemistry modifications retained.

Deferred/blocked targets:

- Published Supplementary File1.xlsx and other data/validation sections: Available but uninspected within cap, not an unavailable source.
- Cited tRNA-blocking predecessor reference52 and historical named enzyme/vendor originals: Directly linked follow-up deferred; no full recipe inheritance or modern kit defaults.

The source-review JSON records actual manifest URLs and availability. This is a capped selected-source reading, not full publication review or confirmed absence.

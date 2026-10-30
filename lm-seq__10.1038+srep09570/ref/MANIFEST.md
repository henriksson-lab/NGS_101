# LM-Seq source retrieval

Published defining article: https://doi.org/10.1038/srep09570 (PMC4381617). Published supplement: https://www.nature.com/articles/srep09570#Sec10.

Use `python3 tools/get_sources.py --doi 10.1038/srep09570`. Originals and text twins belong in the external `sources/lm-seq__10.1038+srep09570/` directory, never in this checkout. The selected source pass used `srep09570_PMC4381617.xml.txt` and `srep09570_supp_srep09570-s1.pdf.txt`, checking actual Figure1 and SI original tables. Published main and its supplement are two originals; no vendor/predecessor originals were opened in this pass.

The local HZG885 source-guard exception handles terminal `ddC` notation only: the shared parser strips that appended nucleotide, while the physical model retains its C. Direct local whole modified-ordering matching uses published main Methods91; no alternate source sequence or shared parser change is supplied.

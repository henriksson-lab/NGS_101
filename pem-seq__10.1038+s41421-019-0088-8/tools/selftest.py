#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
from checks import Check
from chem import (PEM_BA_LOWER, PEM_BA_UP, PEM_BIO_RAG1AB,
                              PEM_I5_R1_RAG1AB, PEM_I7_B01, PEM_P5_I5, PEM_P7)

c = Check()
c.section("PEM-seq Supplementary Table S11 transcription")
c("bridge strands", (PEM_BA_UP, PEM_BA_LOWER), (
    "CCACGCGTGCTCTACANNNNTNNNANNNTNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGT",
    "TGTAGAGCACGCGTGGNNNNNN"))
c("tagging and I7-B01 primers", (PEM_P5_I5, PEM_P7, PEM_I7_B01), (
    "AATGATACGGCGACCACCGAGATCTACACACACTCTTTCCCTACACGACGC",
    "CAAGCAGAAGACGGCATACGAGAT",
    "CAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGC"))
c("RAG1A/B example primers", (PEM_BIO_RAG1AB, PEM_I5_R1_RAG1AB), (
    "GGACTGCTGGAGATTGCTCCAGAGAGGGTTTCCCCTCAAAGGAATCCTTC",
    "ACTCTTTCCCTACACGACGCTCTTCCGATCTAGCGAGGATCTCACCCGGAACAGC"))
c.report()

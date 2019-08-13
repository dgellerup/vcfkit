#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#TEST FOR PULL REQUEST
"""
vcfkit was created and is maintained by Dane Gellerup (https://github.com/dgellerup).

This module is intended to make reading, analyzing, manipulating, and writing VCF
files in Python easier, faster, and accessible to more users.

"""
import os
import gzip
import pandas as pd
import csv
from Bio.Seq import Seq
from Bio.SeqUtils import GC
import matplotlib.pyplot as plt

class VcfFile:

    def __init__(self, vcf):

        with open(vcf, 'r') as vcffile:
            lines = vcffile.readlines()

        metaList = [line[2:].strip() for line in lines if line.startswith("##")]

        self.info = {}
        self.filter = {}
        self.format= {}

        for item in metaList:
            if item.startswith("fileformat"):
                self.fileformat = item.split("=")[-1]
            elif item.startswith("fileDate"):
                self.fileDate = item.split("=")[-1]
            elif item.startswith("source"):
                self.source = item.split("=")[-1]
            elif item.startswith("contig"):
                self.contig = item.split("contig=")[-1]
                if "species" in item:
                    self.contigSpecies = item.strip(">").split("species=")[-1].split(",")[0]
                if "length" in item:
                    self.contigLength = int(item.strip(">").split("length=")[-1].split(",")[0])
                if "assembly" in item:
                    self.contigAssembly = item.strip(">").split("assembly=")[-1].split(",")[0]
                if "taxonomy" in item:
                    self.contigTaxonomy = item.strip(">").split("taxonomy=")[-1].split(",")[0]
            elif item.startswith("reference"):
                self.reference = item.split("=")[-1]
            elif item.startswith("phasing"):
                self.phasing = item.split("=")[-1]
            elif item.startswith("INFO"):
                iden = item.split("ID=")[-1].split(",")[0]
                self.info[iden] = item.strip(">").split("<")[-1].split(",")[1:]
            elif item.startswith("FILTER"):
                iden = item.split("ID=")[-1].split(",")[0]
                self.filter[iden] = item.strip(">").split("<")[-1].split(",")[1:]
            elif item.startswith("FORMAT"):
                iden = item.split("ID=")[-1].split(",")[0]
                self.format[iden] = item.strip(">").split("<")[-1].split(",")[1:]

        vcfdf = pd.DataFrame([line.split("\t") for line in lines if not line.startswith("##")])
        columns = list(vcfdf.iloc[0])
        last = columns[-1]
        last = last.strip()
        columns.pop()
        columns.append(last)
        vcfdf.columns = columns
        vcfdf.drop(vcfdf.index[0], inplace=True)

        self.vcfdf = vcfdf


    def listChromosomes(self):

        return list(set(self.vcfdf['#CHROM']))


    def getChromosomes(self, chromList):

        return self.vcfdf[self.vcfdf['#CHROM'].isin(chromList)]


    def getSNPs(self):

        return self.vcfdf[self.vcfdf['REF'].str.len() < 2]


    def getMicrosatellites(self):

        return self.vcfdf[self.vcfdf['ID'].str.startswith('micro')]


    def qFilter(self, passed='PASS'):

        if passed == 'PASS':
            return self.vcfdf[self.vcfdf['FILTER'] == 'PASS']
        elif passed == 'FAIL':
            return self.vcfdf[self.vcfdf['FILTER'] != 'PASS']
        else:
            return "Please pass string argument 'PASS' or 'FAIL'"


    def getPosition(self, chrom, pos):

        if type(chrom) == int:

            chromdf = self.vcfdf[pd.to_numeric(self.vcfdf['#CHROM']) == chrom]

        else:

            return "Please enter an int for chromosome as the first argument."

        if type(pos) == int:

            return chromdf[pd.to_numeric(chromdf['POS']) == pos]

        elif len(pos) == 2:

            return chromdf[(pd.to_numeric(chromdf['POS']) >= min(pos)) & (pd.to_numeric(chromdf['POS']) <= max(pos))]

        else:

            return "Please enter either one location or a list containing a range [start, stop] of locations."

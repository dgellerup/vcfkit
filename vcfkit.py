#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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


    common_Info_keys = {"AA" : "ancestral allele",

                        "AC" : "allele count in genotypes, for each ALT allele, in the same order as listed",
                        "AD" : "Total read depth for each allele",
                        "ADF" : "Read depth for each allele on the forward strand",
                        "ADR" : "Read depth for each allele on the reverse strand",
                        "AF" : "allele frequency for each ALT allele in the same order as listed: use this when estimated from primary data, not called genotypes",
                        "AN" : "total number of alleles in called genotypes",
                        "BQ" : "RMS base quality at this position",
                        "CIGAR" : "cigar string describing how to align an alternate allele to the reference allele",
                        "DB" : "dbSNP membership",
                        "DP" : "combined depth across samples, e.g. DP=154",
                        "END" : "end position of the variant described in this record (esp. for CNVs)",
                        "H2" : "membership in hapmap2",
                        "H3" : "HapMap3 membership",
                        "MQ" : "RMS mapping quality, e.g. MQ=52",
                        "MQ0" : "Number of MAPQ == 0 reads covering this record",
                        "NS" : "Number of samples with data",
                        "SB" : "strand bias at this position",
                        "SOMATIC" : "indicates that the record is a somatic mutation, for cancer genomics",
                        "VALIDATED" : "validated by follow-up experiment",
                        "1000G" : "1000 Genomes membership"}

    
    
    common_Format_keys = {"AD" : "Read depth for each allele",
                          "ADF" : "Read depth for each allele on the forward strand", 
                          "ADR" : "Read depth for each allele on the reverse strand",
                          "DP" : "Read depth",
                          "EC" : "Expected alternate allele counts",
                          "FT" : "Filter indicating if this genotype was 'called'",
                          "GL" : "Genotype likelihoods",
                          "GP" : "Genotype posterior probabilities",
                          "GQ" : "Conditional genotype quality",
                          "GT" : "Genotype",
                          "HQ" : "Haplotype quality",
                          "MQ" : "RMS mapping quality",
                          "PL" : "Phred-scaled genotype likelihoods rounded to the closest integer", 
                          "PQ" : "Phasing quality",
                          "PS" : "Phase set"}



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


    def getKey(self, key):
        
        if key in common_Info_key:
            
            print(f"{key} : {VcfFile.common_Info_keys.get(key)}")
            
        elif key in common_Format_key:
            
            print(f"{key} : {VcfFile.common_Format_keys.get(key)}")
            
        else:
            
            return "Please enter a valid common key."
        
    
    def QualifyFilter(self, Q_score):
        
        if Q_score == int:
            
            return df.loc[df['QUAL'] >= Q_score]
        
        else:
            
            return "Please enter a Q score as an int."


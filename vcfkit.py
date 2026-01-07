#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
vcfkit was created and is maintained by Dane Gellerup (https://github.com/dgellerup).

This module is intended to make reading, analyzing, manipulating, and writing VCF
files in Python easier, faster, and accessible to more users.

"""
import pandas as pd

class VcfFile:

    common_keys = {"AA" : "ancestral allele",
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

    def __init__(self, vcf):

        with open(vcf, 'r') as vcffile:
            lines = vcffile.readlines()

        metaList = [line[2:].strip() for line in lines if line.startswith("##")]

        self.contigs = {}
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
                contig_dict = {}
                contig_id = int(item.split("ID=")[-1].split(",")[0])
                contig_dict = {}
                contig_dict["contig"] = item.split("contig=")[-1]
                if "species" in item:
                    contig_dict["species"] = item.strip(">").split("species=")[-1].split(",")[0]
                if "length" in item:
                    contig_dict["length"] = int(item.strip(">").split("length=")[-1].split(",")[0])
                if "assembly" in item:
                    contig_dict["assembly"] = item.strip(">").split("assembly=")[-1].split(",")[0]
                if "taxonomy" in item:
                    contig_dict["taxonomy"] = item.strip(">").split("taxonomy=")[-1].split(",")[0]
                self.contigs[contig_id] = contig_dict
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

        vcfdf["#CHROM"] = pd.to_numeric(vcfdf["#CHROM"], errors="coerce").astype("Int64")
        vcfdf["POS"] = pd.to_numeric(vcfdf["POS"], errors="coerce").astype("Int64")

        self.vcfdf = vcfdf


    def list_chromosomes(self):

        return list(set(self.vcfdf['#CHROM']))


    def get_chromosomes(self, chromList):

        return self.vcfdf[self.vcfdf['#CHROM'].isin(chromList)]


    def get_snps(self):

        return self.vcfdf[self.vcfdf['REF'].str.len() < 2]


    def get_microsatellites(self):

        return self.vcfdf[self.vcfdf['ID'].str.startswith('micro')]


    def q_filter(self, passed='PASS'):

        if passed == 'PASS':
            return self.vcfdf[self.vcfdf['FILTER'] == 'PASS']
        elif passed == 'FAIL':
            return self.vcfdf[self.vcfdf['FILTER'] != 'PASS']
        else:
            return "Please pass string argument 'PASS' or 'FAIL'"


    def get_position(self, chrom: int, pos: int):

        # Handle case where user passes strings instead of integers, i.e. "20" vs 20
        chrom, pos = int(chrom), int(pos)

        try:
            chromdf = self.vcfdf[pd.to_numeric(self.vcfdf['#CHROM']) == chrom]
            return chromdf[pd.to_numeric(chromdf['POS']) == pos]
        except Exception as e:
            print(e)
        

    def get_position_slice(self, chrom: int, start: int, end: int):

        # Handle case where user passes strings instead of integers, i.e. "20" vs 20
        chrom, start, end = int(chrom), int(start), int(end)

        try:
            chromdf = self.vcfdf[pd.to_numeric(self.vcfdf['#CHROM']) == chrom]
            return chromdf[(pd.to_numeric(chromdf['POS']) >= start) & (pd.to_numeric(chromdf['POS']) <= end)]
        except TypeError as e:
            print(e)


    def get_contig_info(self, chrom: int):

        return self.contigs.get(chrom)


    def view_common_keys(cls):
        [print(f"{k}:\t{v}") for k, v in cls.common_keys.items()]


    def get_common_key(self, key):

        print(f"{key} : {VcfFile.common_keys.get(key)}")

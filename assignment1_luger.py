#!/usr/bin/python3

import mysql.connector
import mysql.connector
import numpy as np
import pysam
##
## Requirements:
##   chr21.bam and chr21.bam.bai have to be in the same
##   directory as the program.
##
##  specification:
##   Gene "DYRK1A" was used.

__author__ = "Martina Marina Luger"

class Assignment1:

    def __init__(self, geneName, fileName, bamFile):
        ## Your gene of interest
        self.gene = geneName
        self.file_name = fileName
        self.bam_file = bamFile
        self.alignfile = pysam.AlignmentFile(self.bam_file, "rb")
        self.genome_reference = 'hg38'

    def download_gene_coordinates(self):
        print("Connecting to UCSC to fetch data")

        ## Open connection
        cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep', passwd='password',
                                      db=self.genome_reference)

        ## Get cursor
        cursor = cnx.cursor()

        ## Build query fields
        query_fields = ["refGene.name2",
                        "refGene.name",
                        "refGene.chrom",
                        "refGene.txStart",
                        "refGene.txEnd",
                        "refGene.strand",
                        "refGene.exonCount",
                        "refGene.exonStarts",
                        "refGene.exonEnds"]

        ## Build query
        query = "SELECT DISTINCT %s \
                 FROM refGene \
                 WHERE refGene.name2 = 'DYRK1A'" \
                % ",".join(query_fields)

        ## Execute query
        cursor.execute(query)

        ## Write to file
        with open(self.file_name, "w") as fh:
            for row in cursor:
                fh.write(str(row) + "\n")

        ## Close cursor & connection
        cursor.close()
        cnx.close()

        self.read = self.get_read()

        self.gene_symbol = self.read[0]
        self.gene_name = self.read[1]
        self.chromosome = self.read[2]
        self.gene_start = int(self.read[3])
        self.gene_end = int(self.read[4])
        self.strand = self.read[5]
        self.number_of_exons = self.read[6]
        self.exon_starts = self.read[7]
        self.exon_ends = self.read[8]
        print("Done fetching data")

    def get_read(self):
        ## Use UCSC file
        read = ""

        with open(self.file_name, "r") as fh:
            for row in fh:
                read = row

        ## Remove unwanted chars
        unwanted_chars = ["(", ")", "'"]
        for char in unwanted_chars:
            read = read.replace(char, '')

        ## Split string into tuple
        read = read.split(" ")

        for i in range(0, 7):
            read[i] = read[i].replace(",", "")

        read[7] = read[7].split(",")
        read[8] = read[8].split(",")
        read[7][0] = read[7][0].replace("b", "")
        read[8][0] = read[8][0].replace("b", "")

        # print(read)

        return read

    def get_coordinates_of_gene(self, value):

        ## Return corresponding value depending on val
        if value == "start":
            return int(self.read[3])
        elif value == "end":
            return int(self.read[4])
        else:
            print("Error in ", self, ".get_coordinates_of_gene!")
            print("Wrong value supplied during call of function. (Parameter: val)")
            return 100

    def get_gene_chromosome(self):
        return self.chromosome

    def get_gene_symbol(self):
        return self.gene_symbol

    def get_sam_header(self):
        header = self.alignfile.header["HD"]

        head = ""

        for key in header:
            head += key + ": " + header[key] + "\t"

        print(head)
        return head

    def get_properly_paired_reads_of_gene(self):
        self.reads = list(self.alignfile.fetch(self.chromosome, self.gene_start, self.gene_end))
        self.proper_reads = len([i for i in self.reads if i.is_proper_pair])

    def get_gene_reads_with_indels(self):
        rd_indel = []
        for i in self.reads:
            if not i.is_unmapped:
                cig = i.cigartuples
                for (operation, length) in cig:
                    if (operation == 1) or (operation == 2):
                        rd_indel.append(i)
        self.reads_indel = len(rd_indel)

    def calculate_total_average_coverage(self):
        chr_len = [i["LN"] for i in self.alignfile.header["SQ"] if i["SN"] == self.chromosome][0]
        self.chr_avg_cov = round\
        (np.mean(self.alignfile.count_coverage(self.chromosome, start=0, stop=chr_len)), 2)

    def calculate_gene_average_coverage(self):
        self.gene_avg_cov = round\
        (np.mean(self.alignfile.count_coverage(self.chromosome, start=self.gene_start, stop=self.gene_end)), 2)

    def get_number_mapped_reads(self):
        self.reads_mapped = 0
        for i in self.reads:
            if not i.is_unmapped:
                self.reads_mapped += 1

    def get_region_of_gene(self):
        return "\n\tChromosome: " + self.chromosome + \
               "\n\tStart:      " + str(self.gene_start) + \
               "\n\tEnd:        " + str(self.gene_end)

    def get_number_of_exons(self):
        return self.number_of_exons

    def print_summary(self):
        region = self.get_region_of_gene()
        self.get_properly_paired_reads_of_gene()
        self.get_gene_reads_with_indels()
        self.calculate_gene_average_coverage()
        self.calculate_total_average_coverage()
        self.get_number_mapped_reads()
        header = self.get_sam_header()

        print("\n---------- RESULTS ----------")
        print("Gene Symbol:       ", self.gene_symbol)
        print("Start of gene:     ", self.gene_start)
        print("End of gene:       ", self.gene_end)
        print("Number of exons:   ", self.number_of_exons)
        print("Region of gene:    ", region)
        print("Paired reads:      ", self.proper_reads)
        print("Read indels:       ", self.reads_indel)
        print("Total avg coverage:", self.chr_avg_cov)
        print("Gene avg coverage: ", self.gene_avg_cov)
        print("Mapped reads:      ", self.reads_mapped)
        print("SAM Header:        ")
        print("  ", header)
        print("-----------------------------", end="\n\n")


#---------- RESULTS ----------
#Gene Symbol:        DYRK1A
#Start of gene:      37366933
#End of gene:        37515376
#Number of exons:    11
#Region of gene:
#        Chromosome: chr21
#       Start:      37366933
#        End:        37515376
#Paired reads:       50077
#Read indels:        1916
#Total avg coverage: 9.56
#Gene avg coverage:  12.39
#Mapped reads:       50849
#SAM Header:
#   VN: 1.0      SO: unsorted
#-----------------------------


def main():
    print("Assignment 1")
    file_name = "Output1.txt"
    gene_name = "DYRK1A"
    bam_file = "chr21.bam"
    assignment1 = Assignment1(gene_name, file_name, bam_file)

    assignment1.download_gene_coordinates()
    assignment1.print_summary()

    print("Done with assignment 1 :)")


if __name__ == '__main__':
    main()

    

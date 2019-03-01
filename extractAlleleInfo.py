import vcf

def main():

    def calcFreqs(gt, alt_depth, read_depth):
        freq_list = []
        if isinstance(alt_depth, list):  # [8,4] + 0/1 means take first digit from list i.e. 8
            genotype = str(gt.data[0]).split("=")
            list_positions = str(genotype[0]).split("/")
            # if one allele is the reference, calculate the ratio from the alternate allele in the genotype
            if list_positions[0] == '0':
                p = int(list_positions[1]) - 1
                freq = alt_depth[p] / float(read_depth)
                freq_list.append(freq)
                freq_list.append(1 - freq)
            # If both alleles are different from the reference
            else:
                total = 0
                for i in alt_depth:
                    total = total + i  # calculate read depth of all alternate alleles and compare this to the total read depth
                if total != float(read_depth):
                    # add another freq which is the ref reads -> 3 frequencies
                    allele1_pos = int(list_positions[0]) - 1
                    allele2_pos = int(list_positions[1]) - 1
                    freq1 = alt_depth[allele1_pos] / float(read_depth)
                    freq2 = alt_depth[allele2_pos] / float(read_depth)
                    freq3 = (read_depth - total) / float(read_depth)
                    freq_list.append(freq1)
                    freq_list.append(freq2)
                    freq_list.append(freq3)
                else:
                    # need to get positions of read depths from the list
                    allele1_pos = int(list_positions[0]) - 1
                    allele2_pos = int (list_positions[1]) -1
                    freq1 = alt_depth[allele1_pos] / float(read_depth)
                    freq2 = alt_depth[allele2_pos] / float(read_depth)
                    freq_list.append(freq1)
                    freq_list.append(freq2)
        else:
            freq = alt_depth / float(read_depth)
            freq_list.append(freq)
            freq_list.append(1 - freq)

        #print(freq_list)
        return freq_list

# open vcf file and create vcf reader
# hardcoded list of samples and output file
    vcf_file = "/Users/samanthacampbell/PycharmProjects/vcfToCircos/all_variants_parallel.filtered.vcf"
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))

    samples = ['ms018LMex1sp5_S23_L001_R.bam', 'ms018LMex2sp5_S19_L001_R.bam', 'ms018LMex3sp5_S15_L001_R.bam', '01-SP20-18059045_R.bam', '02-SP20-18060043_R.bam', '03-SP20-18058049_R.bam']
    fileOut = "sp5vsp20_alleleFreq-PyCharm_NEW.bed"

    with open(fileOut, "w") as out:
        for record in vcf_reader:
            chr = record.CHROM
            if chr != "LmxM.00":
                position = record.POS
                for sample in samples:
                    gt = record.genotype(sample)
                    # if the genotype for this site in this sample is heterozygous we can calculate a frequency
                    if gt.is_het:
                        read_depth = gt.data[1]
                        alt_depth = gt.data[4]

                        frequencies = calcFreqs(gt, alt_depth, read_depth)
                        for i in frequencies:
                            outLine = chr + "_" + str(position) + "\t" + chr + "\t" + sample + "\t" + str(i) + "\n"
                            out.write(outLine)

if __name__ == '__main__':
        main()
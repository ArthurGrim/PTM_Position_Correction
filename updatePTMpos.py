#!/usr/bin/python
import pip



from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd
import sys
import re


def main(fasta_file, proteo_form_quant_output_file, output_file):

    #load fasta file as a list of tuples, (tuple[0]: sequence ID second; tuple[1]: sequence)
    with open(fasta_file, "rU") as handle:
        fasta = list(SimpleFastaParser(handle))

    #load proteoFormQuant output as pandas dataframe:
    PFQ_df = pd.read_excel(proteo_form_quant_output_file)

    # Update PTMs positions
    # /!\ doesn't verify if a peptide matches to several histones /!\

    n_pos_updated = 0
    peptide_position = []

    for i in range(len(PFQ_df)):

        peptides_sequence = PFQ_df.loc[i, "PeptideSequence"]

        for entry in fasta:
            if entry[1][0] == "M":
                histone_seq = entry[1][1:]
            else:
                histone_seq = entry[1]

            if peptides_sequence in histone_seq:
                peptide_position.append(histone_seq.find(peptides_sequence)) #find the starting position of the peptide in histone's sequence
                PTM_positions = re.findall(r'[0-9]+', str(PFQ_df.loc[i, "PTM_code"])) #extract position(s) of PTMs

                for pos in PTM_positions:
                    n_pos_updated += 1
                    updated_pos = int(pos) + peptide_position[-1]
                    PFQ_df.loc[i, "PTM_code"] = PFQ_df.loc[i, "PTM_code"].replace(pos, str(updated_pos))

    PFQ_df.insert(4, "PositionOfPeptideInHistoneSequence", peptide_position) #Add collumn with peptide start position in respective histone sequence
    PFQ_df.to_excel(output_file,index=False)
    print(n_pos_updated, "PTMs positions have been updated, modified file has been saved as: \' ",output_file," \'   ")

if __name__ == "__main__":
        main(sys.argv[1], sys.argv[2], sys.argv[3])

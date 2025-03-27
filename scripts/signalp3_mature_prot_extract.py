# extracts mature protein sequence using signalp3 results and the input protein fasta (same as the one you used for signalp3).
# only extracts sequences for those that have D-status (D ?) == Y
# and splits at Y-pos which is the combined cleavage site position.

# see explanation copied from signalp3 manual:
##################################################################################################################
# The C-score is the ``cleavage site'' score. For each position in the submitted sequence, a C-score is reported, which should only be significantly high at the cleavage site. Confusion is often seen with the position numbering of the cleavage site. When a cleavage site position is referred to by a single number, the number indicates the first residue in the mature protein, meaning that a reported cleavage site between amino acid 26-27 corresponds to that the mature protein starts at (and include) position 27.

# Y-max is a derivative of the C-score combined with the S-score resulting in a better cleavage site prediction than the raw C-score alone. This is due to the fact that multiple high-peaking C-scores can be found in one sequence, where only one is the true cleavage site. The cleavage site is assigned from the Y-score where the slope of the S-score is steep and a significant C-score is found.
##################################################################################################################

# to-do: add user-defined D-score threshold instead of relying on Y/N status  


from Bio import SeqIO
import argparse
import re

def extract_mature_protein(signalp3_out, input_fasta, output_fasta):
    # parse signalp3 output to get cleavage site position (Y-pos)
    y_positions = {}
    with open(signalp3_out, "r") as f:
        next(f)
        for line in f:
            cols = re.split(r'\s+', line.strip())
            name = cols[0]
            y_pos = int(cols[5])  # column for y-max combined cleavage site position
            d_score = float(cols[12])
            d_status = cols[13]
            if d_status == "Y":
                y_positions[name] = y_pos

    # extract sequences from fasta
    sequences = {record.id.split(" ")[0]: str(record.seq) for record in SeqIO.parse(input_fasta, "fasta")}

    # write out fasta with mature protein sequences i.e. split at Y-pos cleavage site
    with open(output_fasta, "w") as out:
        for name, pos in y_positions.items():
            if name in sequences:
                mature_seq = sequences[name][pos-1:] # signalp3 uses 1-based coords
                out.write(f">{name}\n{mature_seq}\n")


def main():
    parser = argparse.ArgumentParser(description="Extract mature protein sequences from SignalP3 output and FASTA file.")
    parser.add_argument("signalp3_out", help="Path to the SignalP3 output file.")
    parser.add_argument("input_fasta", help="Path to the input FASTA file.")
    parser.add_argument("output_fasta", help="Path to the output FASTA file.")
    args = parser.parse_args()
    
    extract_mature_protein(args.signalp3_out, args.input_fasta, args.output_fasta)

if __name__ == "__main__":
    main()

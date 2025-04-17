import os
import argparse
import pandas as pd



def filter_gff3_by_transcripts_and_genes(gff3_df, sec_ids):
    """
    use a list of IDs of infection transcripts with transdecoder-predicted ORFs encoding putative SPs (sec_ids) to filter the transdecoder genome-based gff3.
    also handles escaping of regular expression as stringtie transcript IDs contain '.' which is problematic.

    note that not everything in the transcript ID list could be found in the genome gff3 as some infection transcripts may not be mappable against the Pst104E genome.
    """
    escaped_sec_ids = [sec_id.replace(".", "\\.") for sec_id in sec_ids]
    transcript_filter = gff3_df[gff3_df["attr"].str.contains('|'.join(escaped_sec_ids)) & (gff3_df["feature"] != "gene")]
    
    sec_id_stripped = ["\\.".join(sec_id.split("\\.")[0:2])+"\\^" for sec_id in escaped_sec_ids]
    gene_filter = gff3_df[gff3_df["attr"].str.contains('|'.join(sec_id_stripped)) & (gff3_df["feature"] == "gene")]
    print("# of SP transcripts found in trandecoder genome gff3: ", len(gene_filter))
    return pd.concat([transcript_filter, gene_filter]).sort_index()


def main(transdecoder_genome_gff, sec_id_file, output_gff):
    with open(sec_id_file, "r") as f:
        lines = f.readlines()
        sec_id = [l.rstrip() for l in lines]

    df = pd.read_csv(transdecoder_genome_gff, sep="\t", header=None, comment="#")
    df.columns = ["chr", "source", "feature", "start", "stop", "score", "strand", "phase", "attr"]

    print("# of infection transcripts with transdecoder-predicted ORFs encoding an SP: ", len(sec_id))
    print("Filtering GFF3 file by secretome IDs...")
    results = filter_gff3_by_transcripts_and_genes(df, sec_id)
    results.to_csv(output_gff, sep="\t", header=None, index=None)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract secretome transcripts/genes using a text file containing transcript IDs.")
    parser.add_argument("-i", "--in_gff3", required=True, help="Input GFF3 file")
    parser.add_argument("--in_secretome_ids", required=True, help="Input secretome transcript id listed line by line in a text file")
    parser.add_argument("--out_gff3", required=True, help="Output secretome gff3 filename")
    args = parser.parse_args()

    main(args.in_gff3, args.in_secretome_ids, args.out_gff3)
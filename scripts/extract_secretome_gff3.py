import os
import argparse
import pandas as pd



def filter_gff3_by_transcripts_and_genes(gff3_df, sec_ids):
    """
    use transcript IDs from secretome CDS IDs to filter the original gff3 for other features, e.g. genes, UTRs, exons.
    """
    transcript_filter = gff3_df[gff3_df["attr"].str.contains('|'.join(sec_ids))]
    gene_ids = set(transcript_filter["attr"].str.extract(r"Parent=([^;]+)")[0])
    gene_filter = gff3_df[gff3_df["attr"].str.contains('|'.join(gene_ids)) & (gff3_df["feature"] == "gene")]
    return pd.concat([transcript_filter, gene_filter]).sort_index()



def main(gff_file, sec_id_file, output_gff):
    with open(sec_id_file, "r") as f:
        lines = f.readlines()
        sec_ids = [l.rstrip() for l in lines]
    df = pd.read_csv(gff_file, sep="\t", header=None, comment="#")
    df.columns = ["chr", "source", "feature", "start", "stop", "score", "strand", "phase", "attr"]
    results = filter_gff3_by_transcripts_and_genes(df, sec_ids)
    results.to_csv(output_gff, sep="\t", header=None, index=None)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract secretome transcripts/genes using a text file containing secretome IDs.")
    parser.add_argument("-i", "--in_gff3", required=True, help="Input GFF3 file")
    parser.add_argument("--in_secretome_ids", required=True, help="Input secretome gene id listed line by line in a text file")
    parser.add_argument("--out_gff3", required=True, help="Output secretome gff3 filename")
    args = parser.parse_args()

    main(args.in_gff3, args.in_secretome_ids, args.out_gff3)
    
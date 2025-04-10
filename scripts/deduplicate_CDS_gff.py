import os
import argparse
import pandas as pd


def extract_cds_coordinates(gene_id_group):
    """
    for every gene id, get all unique transcript IDs and their CDS coordinates as a list of tuples,
    and turn it into a dictionary.
    e.g. {'Pst198E_000202-T1': [(809529, 810372), (810511, 810590), (810659, 810771)],
          'Pst198E_000202-T2': [(809529, 810372), (810511, 810590), (812266, 812414)]}
    """
    coords_dict = {}
    transcript_ids = sorted(gene_id_group["transcript_ID"].unique())
    for t in transcript_ids:
        t_rows = gene_id_group.loc[gene_id_group["transcript_ID"] == t, ["start", "stop"]]
        coords = list(zip(t_rows["start"].values, t_rows["stop"].values))
        coords_dict[t] = coords
    return coords_dict


def filter_gff3_by_transcripts_and_genes(gff3_df, dedup_ids):
    """
    use transcript IDs from deduplicated CDS IDs to filter the original gff3 for other features, e.g. genes, UTRs, exons.
    also keeps all entries for tRNAs as these don't have a CDS hence not subjected to filtering.
    """
    transcript_filter = gff3_df[gff3_df["attr"].str.contains('|'.join(dedup_ids))]
    
    gene_ids = set(transcript_filter["attr"].str.extract(r"Parent=([^;]+)")[0])
    gene_filter = gff3_df[gff3_df["attr"].str.contains('|'.join(gene_ids)) & (gff3_df["feature"] == "gene")]

    tRNA_filter = gff3_df[gff3_df["attr"].str.contains(";product=tRNA-")]
    tRNA_ids = set(tRNA_filter["attr"].apply(lambda x: x.split(";Parent=")[-1].split(";product=tRNA")[0]))
    tRNA = gff3_df[gff3_df["attr"].str.contains('|'.join(tRNA_ids))]
    return pd.concat([transcript_filter, gene_filter, tRNA]).sort_index()




def main(input_gff, dedup_cds_output, dedup_gff3_output):
    df = pd.read_csv(input_gff, sep="\t", header=None, comment="#")
    df.columns = ["chr", "source", "feature", "start", "stop", "score", "strand", "phase", "attr"]
    df["attr"] = df["attr"].apply(lambda x: x.split(";Alias=")[0]) # remove alias in attributes as they have the old IDs

    # get transcript IDs from mRNA entries, sanity check
    transcript_idlist = df[df["feature"] == "mRNA"]["attr"].apply(lambda x: x.split(";Parent=")[0].split("ID=")[-1])
    assert len(list(transcript_idlist)) == len(set(transcript_idlist))

    # deduplicate CDS
    cds = df[df["feature"] == "CDS"].copy()
    cds["transcript_ID"] = cds["attr"].apply(lambda x: x.split(".cds;Parent=")[0].split("ID=")[-1])
    cds["gene_ID"] = cds["transcript_ID"].apply(lambda x: x.split("-T")[0])
    cds_coords_dict = cds.groupby("gene_ID").apply(extract_cds_coordinates).to_dict()

    updated_transcript_ids = {}
    for gene_id, tx_coords in cds_coords_dict.items():
        unique_coords = {}
        for tx_id, coords in tx_coords.items():
            coords_tuple = tuple(coords)
            if coords_tuple not in unique_coords:
                unique_coords[coords_tuple] = tx_id
        for tx_id, coords in tx_coords.items():
            coords_tuple = tuple(coords)
            updated_transcript_ids[tx_id] = unique_coords[coords_tuple]

    cds2 = cds.copy()
    cds2["dedup_transcript_ID"] = cds2["transcript_ID"].map(updated_transcript_ids)

    deduplicated_cds = cds2.drop_duplicates(subset=["gene_ID", "dedup_transcript_ID", "start", "stop"])

    deduplicated_cds.iloc[:,:9].to_csv(dedup_cds_output, index=None, sep="\t", header=None)

    # process the full GFF3 to keep only entries from the deduplicated CDS IDs (plus tRNA as the filtering here does not apply to them)
    df2 = df.copy()
    dedup_ids = set(deduplicated_cds["dedup_transcript_ID"])
    full_gff3_filtered = filter_gff3_by_transcripts_and_genes(df2, dedup_ids)
    full_gff3_filtered.iloc[:,:9].to_csv(dedup_gff3_output, sep="\t", header=None, index=None)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge CDS entries with identical coordinates across transcripts.")
    parser.add_argument("-i", "--input", required=True, help="Input GFF3 file")
    parser.add_argument("--dedup_cds_gff3_out", required=True, help="Output deduplicated CDS GFF3 file")
    parser.add_argument("--dedup_full_gff3_out", required=True, help="Output deduplicated full GFF3 file")
    args = parser.parse_args()

    main(args.input, args.dedup_cds_gff3_out, args.dedup_full_gff3_out)
from collections import Counter
from Bio import SeqIO
import argparse


def combine_unions_faa(union_idlist, signalp3_mature_faa, signalp4_mature_faa, signalp6_mature_faa, output_union_faa, report_file="report.txt"):
    # fetch complete sequences with start and stop codon
    signalp3_seq = {record.id.split(" ")[0]: str(record.seq) for record in SeqIO.parse(signalp3_mature_faa, "fasta") if str(record.seq).endswith("*")} 
    signalp4_seq = {record.id.split(" ")[0]: str(record.seq) + "*" for record in SeqIO.parse(signalp4_mature_faa, "fasta")}   # signalp4 lacks stop codon
    signalp6_seq = {record.id.split(" ")[0]: str(record.seq) for record in SeqIO.parse(signalp6_mature_faa, "fasta") if str(record.seq).endswith("*")}

    # prep report entries
    report = {"total_processed": 0, "one_version": 0, "multiple_versions": 0,
        "identical_seqs": 0, "different_seqs": 0, "majority_seq": 0}

    with open(output_union_faa, "w") as out, open(report_file, "w") as report_out:
        with open(union_idlist, "r") as f:
            for line in f:
                line = line.strip()

                # for each union id, check if detected in each signalP version
                present_in_signalp3 = line in signalp3_seq
                present_in_signalp4 = line in signalp4_seq
                present_in_signalp6 = line in signalp6_seq

                num_versions = sum([present_in_signalp3, present_in_signalp4, present_in_signalp6])

                # case: found in only one version (e.g. TFF) -> write out directly.
                if num_versions == 1:
                    if present_in_signalp3:
                        out.write(f">{line}\n{signalp3_seq[line]}\n")
                    elif present_in_signalp4:
                        out.write(f">{line}\n{signalp4_seq[line]}\n")
                    elif present_in_signalp6:
                        out.write(f">{line}\n{signalp6_seq[line]}\n")
                    report["one_version"] += 1

                # case: found in multiple versions (e.g. TFT, TTF, TTT) -> process further.
                
                elif num_versions > 1:
                    sequences = []
                    if present_in_signalp3:
                        sequences.append(("signalp3", signalp3_seq[line]))
                    if present_in_signalp4:
                        sequences.append(("signalp4", signalp4_seq[line]))
                    if present_in_signalp6:
                        sequences.append(("signalp6", signalp6_seq[line]))

                    # extract unique sequences and count their occurrences
                    seqs = [seq for _, seq in sequences]
                    seq_counts = Counter(seqs)

                    # if all sequences are identical -> write once.
                    if len(seq_counts) == 1:
                        out.write(f">{line}\n{seqs[0]}\n")
                        report["identical_seqs"] += 1

                    # if two out of three sequences are identical -> write majority without version tag.
                    # this only applies to cases where all three versions detect an SP (i.e. TTT)
                    elif num_versions == 3 and max(seq_counts.values()) == 2:
                        majority_seq = max(seq_counts, key=seq_counts.get)
                        out.write(f">{line}\n{majority_seq}\n")
                        report["majority_seq"] += 1

                    # if all sequences different -> write everything with version tags.
                    else:
                        for tag, seq in sequences:
                            out.write(f">{line}_{tag}\n{seq}\n")
                        report["different_seqs"] += 1
                    report["multiple_versions"] += 1
                report["total_processed"] += 1

        # write report summary to file
        report_out.write("Summary Report\n")
        report_out.write(f"total union IDs processed: {report['total_processed']}\n")
        report_out.write(f"found in one SignalP version: {report['one_version']}\n")
        report_out.write(f"found in multiple versions: {report['multiple_versions']}\n")
        report_out.write(f"\t- all sequences identical: {report['identical_seqs']}\n")
        report_out.write(f"\t- different sequences: {report['different_seqs']}\n")
        report_out.write(f"\t- majority used: {report['majority_seq']}\n")


def main():
    parser = argparse.ArgumentParser(description="combine SignalP mature protein sequences from versions 3, 4 and 6.")
    parser.add_argument("union_idlist", help="file containing union IDs (i.e. detected by any signalP version).")
    parser.add_argument("signalp3_mature_faa", help="fasta file of SignalP3 mature proteins.")
    parser.add_argument("signalp4_mature_faa", help="fasta file of SignalP4 mature proteins.")
    parser.add_argument("signalp6_mature_faa", help="fasta file of SignalP6 mature proteins.")
    parser.add_argument("output_union_faa", help="Output fasta file with combined sequences.")
    parser.add_argument("--report_file", default="report.txt", help="report file (default: report.txt).")

    args = parser.parse_args()
    combine_unions_faa(args.union_idlist, args.signalp3_mature_faa, args.signalp4_mature_faa, args.signalp6_mature_faa, args.output_union_faa, args.report_file)       


if __name__ == "__main__":
    main()
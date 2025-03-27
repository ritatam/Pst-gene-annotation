import os
import pandas as pd
import argparse

def filter_tm_proteins(phobius, tmhmm, signalp_union_id, output): 

    phobius_df = pd.read_csv(phobius, delim_whitespace=True, header=None, names=["sequence_id", "TM", "SP", "prediction"])
    tmhmm_df = pd.read_csv(tmhmm, sep="	", header=None, names=["sequence_id", "len", "expected_num_TM_aa", "expected_num_TM_aa_first60", "num_predicted_TM_helices_by_Nbest", "topology"])

     # make sure both tables have the same length
    assert len(phobius_df) == len(tmhmm_df)
    
    phobius_tm_list = list(phobius_df[phobius_df["TM"] != 0]["sequence_id"])
    tmhmm_tm_list = list(tmhmm_df[tmhmm_df["num_predicted_TM_helices_by_Nbest"] != "PredHel=0"]["sequence_id"])
    
    # find entries that have a TM helice detected by both phobius and tmhmm (aka intersection).
    intersection_list = [] 
    for n in set(phobius_tm_list).intersection(tmhmm_tm_list):
        if "_signalp" in n:
            intersection_list.append(n.split("_signalp")[0])
        else:
            intersection_list.append(n)
    exclusion_list = list(set(intersection_list))  
    
     # filter mature protein ids to retain those without TM as our final SP list
    with open(output, "w") as out:
        with open(signalp_union_id, "r") as f:
            for line in f:
                line = line.rstrip()
                if line not in exclusion_list: 
                    print(line, file=out)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--phobius", required=True, help="Path to the Phobius output file")
    parser.add_argument("--tmhmm", required=True, help="Path to the TMHMM output file")
    parser.add_argument("--signalp_union_id", required=True, help="Path to the SignalP mature protein ID list")
    parser.add_argument("--output", required=True, help="Path to the output filtered ID list")
    args = parser.parse_args()
    
    filter_tm_proteins(args.phobius, args.tmhmm, args.signalp_union_id, args.output)

if __name__ == "__main__":
    main()

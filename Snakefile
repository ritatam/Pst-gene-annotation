import os
from glob import glob

configfile: "config.yaml"

REF = config["ref"]
PREASSEMBLED_TRANSCRIPTS = config["preassembled_transcripts_fasta"]
PREASSEMBLED_TRANSCRIPTS_GTF = config["preassembled_transcripts_gtf"]
SPORES_GTF = config["spores_gtf"]
INFECTIONS_GTF = config["infections_gtf"]
ONT_CDNA_FASTQ = config["ONT_cDNA_fastq"]
ONT_CDNA_FILT_Q = config["ONT_cDNA_filt_quality"]
OUTDIR = config["outdir"].rstrip("/")
FUNANNOTATE_DB_DIR = config["funannotate_db_dir"].rstrip("/")
GENEMARK_DIR = config["genemark_dir"].rstrip("/")
EGGNOG_DATA_DIR = config["eggnog_data_dir"].rstrip("/")
CODINGQUARRY_DIR = config["codingquarry_dir"]
FUNANNOTATE_SPECIES_NAME = config["funannotate_species_name"]
FUNANNOTATE_ISOLATE_NAME = config["funannotate_isolate_name"]

conda_init_cmd = "set +eu && . $(conda info --base)/etc/profile.d/conda.sh && conda activate "

# decide whether to do data filtering by user
if ONT_CDNA_FILT_Q == False:
    pass
else:
    assert int(ONT_CDNA_FILT_Q) > 0 , "Please provide a integer for the read filtering Q score. If no filtering is required, put False."
    tmp = f"{OUTDIR}/{os.path.basename(ONT_CDNA_FASTQ).split('.fastq.gz')[0]}.q{ONT_CDNA_FILT_Q}.fastq.gz"
    ONT_CDNA_FASTQ = tmp


SORTED_REF = f"{OUTDIR}/{os.path.splitext(os.path.basename(REF))[0]}.sorted.fasta"
SPORES_GFF3 = f"{OUTDIR}/codingquarryPM/{os.path.splitext(os.path.basename(SPORES_GTF))[0]}.gff3"
INFECTIONS_GFF3 = f"{OUTDIR}/codingquarryPM/{os.path.splitext(os.path.basename(INFECTIONS_GTF))[0]}.gff3"
ALL_TRANSCRIPTS_GFF3 = f"{OUTDIR}/{os.path.splitext(os.path.basename(PREASSEMBLED_TRANSCRIPTS_GTF))[0]}.gff3"
FUNANNOTATE_OUT_PREFIX = f"{'_'.join(config['funannotate_species_name'].split(' '))}_{config['funannotate_isolate_name']}"
infections_prefix = os.path.basename(INFECTIONS_GTF).split('.gtf')[0] 

#######################
##### RULES START #####
#######################


rule all:
    input:
        ONT_CDNA_FASTQ,
        f"{OUTDIR}/signalp/signalp.mature_prot.union.TM-filtered.IDlist"


rule chopper_filtering:
    input:
        ONT_CDNA_FASTQ
    params:
        ONT_CDNA_FILT_Q
    output:
        f"{OUTDIR}/{os.path.splitext(os.path.basename(ONT_CDNA_FASTQ))[0]}.q{ONT_CDNA_FILT_Q}.fastq.gz"
    threads:
        config["threads"]
    shell:
        conda_init_cmd+"{config[chopper_conda_env]}; "
        """
        chopper -q {params} -i {input} --threads {threads} | gzip > {output}
        """


rule funannotate_sort:
    input: REF
    output: SORTED_REF
    shell:
        conda_init_cmd+"{config[funannotate_conda_env]}; "
        """
        funannotate sort --simplify -i {input} -o {output}
        """


rule funannotate_train:
    input:
        ref = SORTED_REF,
        fastq = ONT_CDNA_FASTQ,
        preassembled_transcripts = PREASSEMBLED_TRANSCRIPTS
    output:
        rna_bam = f"{OUTDIR}/funannotate/training/funannotate_train.coordSorted.bam",
        pasa_gff = f"{OUTDIR}/funannotate/training/funannotate_train.pasa.gff3"
    params:
        outdir = directory(f"{OUTDIR}/funannotate")
    threads:
        config["threads"]
    shell:
        conda_init_cmd+"{config[funannotate_conda_env]}; "
        """
        funannotate train 
            -i {input.ref} 
            --nanopore_cdna {input.fastq} 
            --trinity {input.preassembled_transcripts} 
            --jaccard_clip 
            --species '{config[funannotate_species_name]}' 
            --isolate {config[funannotate_isolate_name]}
            --no_trimmomatic
            --cpus {threads}
            --memory 100G
            --pasa_db mysql
            -o {params.outdir}
        """


rule convert_gtf_to_gff3:
    input: 
        spores_gtf = SPORES_GTF, 
        infections_gtf = INFECTIONS_GTF,
        all_gtf = PREASSEMBLED_TRANSCRIPTS_GTF
    output: 
        spores_gff3 = SPORES_GFF3,
        infections_gff3 = INFECTIONS_GFF3,
        all_gff3 = ALL_TRANSCRIPTS_GFF3
    shell:
        conda_init_cmd+"{config[gffread_conda_env]}; "
        "gffread {input.spores_gtf} -o {output.spores_gff3}; "
        "gffread {input.infections_gtf} -o {output.infections_gff3}; "
        "gffread {input.all_gtf} -o {output.all_gff3}"



rule codingquarryPM_spores:
    input: 
        ref = SORTED_REF,
        spores_gff3 = SPORES_GFF3
    params:
        spores_outdir = f"{OUTDIR}/codingquarryPM/spores",
        cq_gff3_list = expand(f"{OUTDIR}/codingquarryPM/spores/out/{{types}}.gff3", types=["DubiousSet", "PGN_predictedPass", "PredictedPass"])
    output:
        f"{OUTDIR}/codingquarryPM/{FUNANNOTATE_ISOLATE_NAME}.CodingQuarry-PM.spores.gff3"
    shell:
        conda_init_cmd+"{config[codingquarry_conda_env]}; "
        """
        export PATH=$PATH:{config[signalp4_dir]}; 
        export PATH=$PATH:{config[codingquarry_dir]}; 
        mkdir -p {params.spores_outdir}; 
        cd {params.spores_outdir}; 
        cp -r {config[codingquarry_dir]}/QuarryFiles .; 
        run_CQ-PM_unstranded.sh {input.ref} {input.spores_gff3} && cat {params.cq_gff3_list} > {output}
        """


rule codingquarryPM_infections:
    input: 
        ref = SORTED_REF,
        infections_gff3 = INFECTIONS_GFF3
    params:
        infections_outdir = f"{OUTDIR}/codingquarryPM/infections",
        cq_gff3_list = expand(f"{OUTDIR}/codingquarryPM/infections/out/{{types}}.gff3", types=["DubiousSet", "PGN_predictedPass", "PredictedPass"])
    output:
        f"{OUTDIR}/codingquarryPM/{FUNANNOTATE_ISOLATE_NAME}.CodingQuarry-PM.infections.gff3"
    shell:
        conda_init_cmd+"{config[codingquarry_conda_env]}; "
        """
        export PATH=$PATH:{config[signalp4_dir]}; 
        export PATH=$PATH:{config[codingquarry_dir]}; 
        mkdir -p {params.infections_outdir}; 
        cd {params.infections_outdir}; 
        cp -r {config[codingquarry_dir]}/QuarryFiles .; 
        run_CQ-PM_unstranded.sh {input.ref} {input.infections_gff3} && cat {params.cq_gff3_list} > {output}
        """


rule funannotate_predict:
    input:
        ref = SORTED_REF,
        fastq = ONT_CDNA_FASTQ,
        preassembled_transcripts = PREASSEMBLED_TRANSCRIPTS,
        rna_bam = f"{OUTDIR}/funannotate/training/funannotate_train.coordSorted.bam",
        pasa_gff = f"{OUTDIR}/funannotate/training/funannotate_train.pasa.gff3",
        infections_transcripts = f"{OUTDIR}/codingquarryPM/{FUNANNOTATE_ISOLATE_NAME}.CodingQuarry-PM.infections.gff3",
        spores_transcripts = f"{OUTDIR}/codingquarryPM/{FUNANNOTATE_ISOLATE_NAME}.CodingQuarry-PM.spores.gff3",
        all_transcripts = ALL_TRANSCRIPTS_GFF3
    params:
        protein_evidence = expand(config["protein_evidence"]),
        weights = "augustus:4 hiq:6 genemark:1 pasa:10 codingquarry:0 snap:1 glimmerhmm:1 proteins:1 transcripts:2",
        outdir = f"{OUTDIR}/funannotate"
    output:
        predicted_gff3 = f"{OUTDIR}/funannotate/predict_results/{FUNANNOTATE_OUT_PREFIX}.gff3"
    threads:
        config["threads"]
    shell:
        conda_init_cmd+"{config[funannotate_conda_env]}; "
        """
        funannotate predict 
            -i {input.ref}
            -o {params.outdir}
            --species '{config[funannotate_species_name]}' 
            --isolate {config[funannotate_isolate_name]}
            --transcript_evidence {input.preassembled_transcripts}
            --rna_bam {input.rna_bam}
            --pasa_gff {input.pasa_gff}
            --other_gff {input.infections_transcripts}:20 {input.spores_transcripts}:10
            --transcript_alignments {input.all_transcripts}
            --protein_evidence {params.protein_evidence}
            --weights {params.weights}
            --optimize_augustus
            --repeats2evm
            --ploidy 1
            --cpus {threads}
        """


rule funannotate_update:
    input:
        ref = SORTED_REF,
        predicted_gff3 = f"{OUTDIR}/funannotate/predict_results/{FUNANNOTATE_OUT_PREFIX}.gff3"
    params:
        outdir = f"{OUTDIR}/funannotate"
    output:
        updated_gff3 = f"{OUTDIR}/funannotate/update_results/{FUNANNOTATE_OUT_PREFIX}.gff3"
    threads:
        config["threads"]
    shell:
        conda_init_cmd+"{config[funannotate_conda_env]}; "
        """
        funannotate update 
            -f {input.ref}
            -o {params.outdir}
            -g {input.predicted_gff3}
            --nanopore_cdna
            --trinity
            --jaccard_clip
            --no_trimmomatic
            --cpus {threads}
            --memory 100G
            --pasa_db mysql
            --species '{config[funannotate_species_name]}' 
            --isolate {config[funannotate_isolate_name]}
            --min_protlen 50
            --alt_transcripts 1
        """



rule transdecoder_preprocess_infection_transcripts:
    input:
        ref = SORTED_REF,
        infection_transcripts_gtf = INFECTIONS_GTF
    output:
        infection_transcripts_fasta = f"{OUTDIR}/SP_rescue/transdecoder/{infections_prefix}.fasta",
        infection_transcripts_gff3 = f"{OUTDIR}/SP_rescue/transdecoder/{infections_prefix}.gff3"
    params:
        transdecoder_util_dir = f"{config['funannotate_conda_env']}/opt/transdecoder/util"
    shell:
        conda_init_cmd+"{config[funannotate_conda_env]}; "
        """
        export PATH=$PATH:{params}; 
        gtf_genome_to_cdna_fasta.pl {input.infection_transcripts_gtf} {input.ref} > {output.infection_transcripts_fasta}; 
        gtf_to_alignment_gff3.pl {input.infection_transcripts_gtf} > {output.infection_transcripts_gff3}
        """


rule transdecoder:
    input:
        f"{OUTDIR}/SP_rescue/transdecoder/{infections_prefix}.fasta"
    output:
        directory(f"{OUTDIR}/SP_rescue/transdecoder/{infections_prefix}.fasta.transdecoder_dir"),
        f"{OUTDIR}/SP_rescue/transdecoder/{infections_prefix}.fasta.transdecoder.gff3",
        f"{OUTDIR}/SP_rescue/transdecoder/{infections_prefix}.fasta.transdecoder.pep"
    params:
        f"{OUTDIR}/SP_rescue/transdecoder"
    shell:
        conda_init_cmd+"{config[funannotate_conda_env]}; "
        """
        TransDecoder.LongOrfs -m 50 -t {input} -O {params}; 
        TransDecoder.Predict -t {input} --single_best_only -O {params}
        """

# generate genome-based transcript gff3 to merge detected SP annotations with funannotate gff3 later 
rule transdecoder_generate_genome_based_gff3:
    input:
        infections_transcripts_transdecoder_gff3 = f"{OUTDIR}/SP_rescue/transdecoder/{infections_prefix}.fasta.transdecoder.gff3",
        infection_transcripts_gff3 = f"{OUTDIR}/SP_rescue/transdecoder/{infections_prefix}.gff3",
        infection_transcripts_fasta = f"{OUTDIR}/SP_rescue/transdecoder/{infections_prefix}.fasta"
    output:
        f"{OUTDIR}/SP_rescue/transdecoder/{infections_prefix}.fasta.transdecoder.genome.gff3"
    params:
        transdecoder_util_dir = f"{config['funannotate_conda_env']}/opt/transdecoder/util"
    shell:
        conda_init_cmd+"{config[funannotate_conda_env]}; "
        """
        export PATH=$PATH:{params}; 
        cdna_alignment_orf_to_genome_orf.pl {input.infections_transcripts_transdecoder_gff3} {input.infection_transcripts_gff3} {input.infection_transcripts_fasta} > {output}
        """


rule filter_transdecoder_complete_orfs:
    input:
        transdecoder_pep = f"{OUTDIR}/SP_rescue/transdecoder/{infections_prefix}.fasta.transdecoder.pep",
        transdecoder_genome_gff3 = f"{OUTDIR}/SP_rescue/transdecoder/{infections_prefix}.fasta.transdecoder.genome.gff3"
    output:
        transdecoder_pep_complete = f"{OUTDIR}/SP_rescue/transdecoder/{infections_prefix}.fasta.transdecoder.pep.complete.fasta",
        transdecoder_pep_complete_ids_list = f"{OUTDIR}/SP_rescue/transdecoder/{infections_prefix}.fasta.transdecoder.pep.complete.IDlist"
    shell:
        conda_init_cmd+"{config[agat_conda_env]}; "
        # get complete transdecoder protein sequences
        "bash scripts/filter_transdecoder_complete_orfs.sh {input.transdecoder_pep} {output.transdecoder_pep_complete}; "
        # extract all headers from the complete protein sequences fasta
        "awk '/^>/ {{print substr($1,2)}}' {output.transdecoder_pep_complete} > {output.transdecoder_pep_complete_ids_list} "


rule signalp3_transdecoder:
    input:
        f"{OUTDIR}/SP_rescue/transdecoder/{infections_prefix}.fasta.transdecoder.pep.complete.fasta"
    params:
        signalp3 = config["signalp3_bin"]
    output:
        flatten = temp(directory(f"{OUTDIR}/SP_rescue/signalp/flatten")),
        signalp3_out = f"{OUTDIR}/SP_rescue/signalp/signalp3.out",
        mature_prot = f"{OUTDIR}/SP_rescue/signalp/signalp3.mature_prot.faa"
    threads:
        config["threads"]-2
    shell:
        conda_init_cmd+"{config[signalp6_conda_env]}; "
        """
        seqkit split --by-id {input} --by-id-prefix '' -O {output.flatten}; 
        echo running signalp3-nn euk; 
        echo 'name	Cmax	pos	?	Ymax	pos	?	Smax	pos	?	Smean	?	D	?' > {output.signalp3_out}; 
        find {output.flatten} -name "*.fasta" | parallel -j {threads} "{params.signalp3} -t euk -short -m nn {{}} | sed '1,2d' >> {output.signalp3_out}"; 
        conda deactivate; 
        """+
        conda_init_cmd+"{config[gffread_conda_env]}; "
        "python scripts/signalp3_mature_prot_extract.py {output.signalp3_out} {input} {output.mature_prot}"


rule signalp4_transdecoder:
    input:
        f"{OUTDIR}/SP_rescue/transdecoder/{infections_prefix}.fasta.transdecoder.pep.complete.fasta"
    params:
        signalp4 = config["signalp4_bin"]
    output:
        mature_prot_tmp1 = temp(f"{OUTDIR}/SP_rescue/signalp/signalp4.mature_prot.tmp1"),
        mature_prot_tmp2 = temp(f"{OUTDIR}/SP_rescue/signalp/signalp4.mature_prot.tmp2"),
        mature_prot = f"{OUTDIR}/SP_rescue/signalp/signalp4.mature_prot.faa",
        signalp4_out = f"{OUTDIR}/SP_rescue/signalp/signalp4.out"
    shell:
        "{params.signalp4} -f short -t euk -s best -m {output.mature_prot_tmp1} {input} > {output.signalp4_out}; "
        "awk '/^>/ {{print $1}} !/^>/ {{print}}' {output.mature_prot_tmp1} > {output.mature_prot_tmp2}; " # name cleanup
        """
        awk '!/^>/ {{printf "%s",$0; n="\n"}} /^>/ {{print n$0; n=""}} END {{printf "%s",n}}' {output.mature_prot_tmp2} > {output.mature_prot}
        """ # unwraps fasta


# slow-sequential, cpu mode 
rule signalp6_transdecoder:
    input:
        f"{OUTDIR}/SP_rescue/transdecoder/{infections_prefix}.fasta.transdecoder.pep.complete.fasta"
    params:
        outdir = f"{OUTDIR}/SP_rescue/signalp/signalp6-slow-sequential",
        results = f"{OUTDIR}/SP_rescue/signalp/signalp6-slow-sequential/prediction_results.txt",
        mature_prot = f"{OUTDIR}/SP_rescue/signalp/signalp6-slow-sequential/processed_entries.fasta"
    output:
        mature_prot = f"{OUTDIR}/SP_rescue/signalp/signalp6.mature_prot.faa",
        signalp6_out = f"{OUTDIR}/SP_rescue/signalp/signalp6.out"
    shell:
        conda_init_cmd+"{config[signalp6_conda_env]}; "
        "signalp6 --fastafile {input} --organism eukarya --format txt --mode slow-sequential --output_dir {params.outdir}; "
        # names cleanup
        """
        awk '/^>/ {{print $1}} !/^>/ {{print}}' {params.mature_prot} > {output.mature_prot}; 
        awk -F'\t' '!/^#/ {{split($1, a, " "); $1=a[1]}}1' OFS='\t' {params.results} > {output.signalp6_out}; 
        """ 


rule take_signalp_union_transdecoder:
    input:
        f"{OUTDIR}/SP_rescue/signalp/signalp3.mature_prot.faa",
        f"{OUTDIR}/SP_rescue/signalp/signalp4.mature_prot.faa",
        f"{OUTDIR}/SP_rescue/signalp/signalp6.mature_prot.faa"
    params:
        f"{OUTDIR}/SP_rescue/signalp/signalp*.mature_prot.faa"
    output:
        f"{OUTDIR}/SP_rescue/signalp/signalp.mature_prot.union.IDlist"
    run:
        """
        cat {params} | awk 'sub(/^>/, "")' | sort | uniq > {output}
        """



rule combine_signalp_unions_mature_proteins_fasta_transdecoder:
    input:
        f"{OUTDIR}/SP_rescue/signalp/signalp3.mature_prot.faa",
        f"{OUTDIR}/SP_rescue/signalp/signalp4.mature_prot.faa",
        f"{OUTDIR}/SP_rescue/signalp/signalp6.mature_prot.faa",
        union_idlist = f"{OUTDIR}/SP_rescue/signalp/signalp.mature_prot.union.IDlist"
    params:
        mature_prot_file_list = expand(f"{OUTDIR}/SP_rescue/signalp/signalp{{versions}}.mature_prot.faa", versions=["3", "4", "6"]),
        report_file = f"{OUTDIR}/SP_rescue/signalp/union_report.txt"
    output:
        f"{OUTDIR}/SP_rescue/signalp/signalp.mature_prot.union.faa"
    shell:
        conda_init_cmd+"{config[gffread_conda_env]}; "
        "python scripts/combine_signalp_unions_mature_prot_faa.py {input.union_idlist} {params.mature_prot_file_list} {output} --report_file {params.report_file}"


rule tmmhm_transdecoder:     #v2.0c
    input:
        f"{OUTDIR}/SP_rescue/signalp/signalp.mature_prot.union.faa"
    output:
        f"{OUTDIR}/SP_rescue/signalp/transmembrane_domain/tmhmm.txt"
    params:
        wdir = f"{OUTDIR}/SP_rescue/signalp/transmembrane_domain"
    shell:
        "export PATH=$PATH:{config[tmhmm_dir]}; "+
        conda_init_cmd+"{config[funannotate_conda_env]}; "
        "tmhmm --workdir {params.wdir} --short {input} > {output}"


rule phobius_transdecoder:       #phobius.pl v1.01
    input:
        f"{OUTDIR}/SP_rescue/signalp/signalp.mature_prot.union.faa"
    output:
        f"{OUTDIR}/SP_rescue/signalp/transmembrane_domain/phobius.txt"
    shell:
        "export PATH=$PATH:{config[phobius_dir]}; "+
        # awk to remove stop codon * in all sequences because phobius doesn't like it! why even, i want to bark. 
        """
        awk '/^>/ {{print}} !/^>/ {{sub(/\*$/, ""); print}}' {input} | phobius.pl -short | tail -n +2 > {output}
        """


rule filter_transdecoder_SP_for_absence_of_TM_helice:
    input:
        tmhmm = f"{OUTDIR}/SP_rescue/signalp/transmembrane_domain/tmhmm.txt",
        phobius = f"{OUTDIR}/SP_rescue/signalp/transmembrane_domain/phobius.txt",
        union_idlist = f"{OUTDIR}/SP_rescue/signalp/signalp.mature_prot.union.IDlist"
    output:
        f"{OUTDIR}/SP_rescue/signalp/signalp.mature_prot.union.TM-filtered.IDlist",
    shell:
        conda_init_cmd+"{config[gffread_conda_env]}; "
        "python scripts/transmembrane_filter.py  --phobius {input.phobius} --tmhmm {input.tmhmm} --signalp_union_id {input.union_idlist} --output {output}"


rule extract_transdecoder_SP_noTM_genes_gff3:    # plus some attribute and name clean up
    input:
        transdecoder_genome_gff3 = f"{OUTDIR}/SP_rescue/transdecoder/{infections_prefix}.fasta.transdecoder.genome.gff3",
        signalp_union_idlist = f"{OUTDIR}/SP_rescue/signalp/signalp.mature_prot.union.TM-filtered.IDlist"
    output:
        tmp1 = temp(f"{OUTDIR}/SP_rescue/transdecoder/{infections_prefix}.fasta.transdecoder.genome.SP_TM-filtered.gff3.tmp1"),
        tmp2 = temp(f"{OUTDIR}/SP_rescue/transdecoder/{infections_prefix}.fasta.transdecoder.genome.SP_TM-filtered.gff3.tmp2"),
        gff3 = f"{OUTDIR}/SP_rescue/transdecoder/{infections_prefix}.fasta.transdecoder.genome.SP_TM-filtered.gff3"
    shell: 
        conda_init_cmd+"{config[gffread_conda_env]}; "
        "agat_sp_filter_feature_from_keep_list.pl --gff {input.transdecoder_genome_gff3} --keep_list {input.signalp_union_idlist} -o {output.tmp1}; "
        "agat_sp_manage_attributes.pl -f {output.tmp1} --att Name -o {output.tmp2}; "
        "sed -e 's/\^.*\^-//g' -e 's/\^.*\^+//g' {output.tmp2} > {output.gff3}"


rule merge_transdecoderSP_with_funannotate:
    input:
        funannotate_update_gff3 = f"{OUTDIR}/funannotate/update_results/{FUNANNOTATE_OUT_PREFIX}.gff3",
        transdecoder_genome_SP_gff3 = f"{OUTDIR}/SP_rescue/transdecoder/{infections_prefix}.fasta.transdecoder.genome.SP_TM-filtered.gff3"
    output:
        f"{OUTDIR}/SP_rescue/funannotate_update.transdecoder_SP_TM-filtered.merged.gff3"
    shell:
        conda_init_cmd+"{config[gffread_conda_env]}; "
        "agat_sp_merge_annotations.pl -f {input.funannotate_update_gff3} -f {input.transdecoder_genome_SP_gff3} -o {output}"


rule renumber_gene_models_after_merge:
    input:
        ref = SORTED_REF,
        merged_gff3 = f"{OUTDIR}/SP_rescue/funannotate_update.transdecoder_SP_TM-filtered.merged.gff3"
    output:
        f"{OUTDIR}/SP_rescue/funannotate_update.transdecoder_SP_TM-filtered.merged.renumbered.gff3"
    shell:
        conda_init_cmd+"{config[funannotate_conda_env]}; "
        "funannotate gff-rename -g {input.merged_gff3} -f {input.ref} -o {output} --locus_tag {config[funannotate_locus_tag]}"


rule get_protein_sequences_from_merged_gff3:
    input:
        ref = SORTED_REF,
        merged_gff3 = f"{OUTDIR}/SP_rescue/funannotate_update.transdecoder_SP_TM-filtered.merged.renumbered.gff3"
    output:
        faa = f"{OUTDIR}/SP_rescue/funannotate_update.transdecoder_SP_TM-filtered.merged.renumbered.protein.faa",
        symlink = f"{OUTDIR}/signalp/funannotate.transdecoder-SP.combined.protein.faa"
    params:
        f"{OUTDIR}/signalp"
    shell:
        conda_init_cmd+"{config[funannotate_conda_env]}; "
        "funannotate util gff2prot -g {input.merged_gff3} -f {input.ref} > {output.faa}; "
        "mkdir -p {params}; "
        "ln -sr {output.faa} {output.symlink}"


rule signalp3:
    input:
        f"{OUTDIR}/signalp/funannotate.transdecoder-SP.combined.protein.faa"
    params:
        signalp3 = config["signalp3_bin"]
    output:
        flatten = temp(directory(f"{OUTDIR}/signalp/flatten")),
        signalp3_out = f"{OUTDIR}/signalp/signalp3.out",
        mature_prot = f"{OUTDIR}/signalp/signalp3.mature_prot.faa"
    threads:
        config["threads"]
    shell:
        conda_init_cmd+"{config[signalp6_conda_env]}; "
        """
        seqkit split --by-id {input} --by-id-prefix '' -O {output.flatten}; 
        echo running signalp3-nn euk; 
        echo 'name	Cmax	pos	?	Ymax	pos	?	Smax	pos	?	Smean	?	D	?' > {output.signalp3_out}; 
        find {output.flatten} -name "*.faa" | parallel -j {threads} "{params.signalp3} -t euk -short -m nn {{}} | sed '1,2d' >> {output.signalp3_out}"; 
        conda deactivate; 
        """+
        conda_init_cmd+"{config[gffread_conda_env]}; "
        "python scripts/signalp3_mature_prot_extract.py {output.signalp3_out} {input} {output.mature_prot}"


rule signalp4:
    input:
        f"{OUTDIR}/signalp/funannotate.transdecoder-SP.combined.protein.faa"
    params:
        signalp4 = config["signalp4_bin"]
    output:
        mature_prot_tmp1 = temp(f"{OUTDIR}/signalp/signalp4.mature_prot.tmp1"),
        mature_prot_tmp2 = temp(f"{OUTDIR}/signalp/signalp4.mature_prot.tmp2"),
        mature_prot = f"{OUTDIR}/signalp/signalp4.mature_prot.faa",
        signalp4_out = f"{OUTDIR}/signalp/signalp4.out"
    shell:
        """
        {params.signalp4} -f short -t euk -s best -m {output.mature_prot_tmp1} {input} > {output.signalp4_out}; 
        awk '/^>/ {{print $1}} !/^>/ {{print}}' {output.mature_prot_tmp1} > {output.mature_prot_tmp2}; 
        awk '!/^>/ {{printf "%s", $0; n="\\n"}} /^>/ {{print n$0; n=""}} END {{printf "%s",n}}' {output.mature_prot_tmp2} > {output.mature_prot}
        """ # unwraps fasta


# slow-sequential, cpu mode 
rule signalp6:
    input:
        f"{OUTDIR}/signalp/funannotate.transdecoder-SP.combined.protein.faa"
    params:
        outdir = f"{OUTDIR}/signalp/signalp6-slow-sequential",
        results = f"{OUTDIR}/signalp/signalp6-slow-sequential/prediction_results.txt",
        mature_prot = f"{OUTDIR}/signalp/signalp6-slow-sequential/processed_entries.fasta"
    output:
        mature_prot = f"{OUTDIR}/signalp/signalp6.mature_prot.faa",
        signalp6_out = f"{OUTDIR}/signalp/signalp6.out"
    shell:
        conda_init_cmd+"{config[signalp6_conda_env]}; "
        "signalp6 --fastafile {input} --organism eukarya --format txt --mode slow-sequential --output_dir {params.outdir}; "
        # names cleanup
        """
        awk '/^>/ {{print $1}} !/^>/ {{print}}' {params.mature_prot} > {output.mature_prot}; 
        awk -F'\t' '!/^#/ {{split($1, a, " "); $1=a[1]}}1' OFS='\t' {params.results} > {output.signalp6_out}; 
        """ 


rule take_signalp_union:
    input:
        expand(f"{OUTDIR}/signalp/signalp{{versions}}.mature_prot.faa", versions=["3", "4", "6"])
    params:
        f"{OUTDIR}/signalp/signalp*.mature_prot.faa"
    output:
        f"{OUTDIR}/signalp/signalp.mature_prot.union.IDlist"
    shell:
        """
        cat {params} | awk 'sub(/^>/, "")' | sort | uniq > {output}
        """


rule combine_signalp_unions_mature_proteins_fasta:
    input:
        union_idlist = f"{OUTDIR}/signalp/signalp.mature_prot.union.IDlist",
        mature_prot_file_list = expand(f"{OUTDIR}/signalp/signalp{{versions}}.mature_prot.faa", versions=["3", "4", "6"])
    params:
        report_file = f"{OUTDIR}/signalp/union_report.txt"
    output:
        f"{OUTDIR}/signalp/signalp.mature_prot.union.faa"
    shell:
        conda_init_cmd+"{config[gffread_conda_env]}; "
        "python scripts/combine_signalp_unions_mature_prot_faa.py {input.union_idlist} {input.mature_prot_file_list} {output} --report_file {params.report_file}"


rule tmmhm:     #v2.0c
    input:
        f"{OUTDIR}/signalp/signalp.mature_prot.union.faa"
    output:
        f"{OUTDIR}/signalp/transmembrane_domain/tmhmm.txt"
    params:
        wdir = f"{OUTDIR}/signalp/transmembrane_domain"
    shell:
        "export PATH=$PATH:{config[tmhmm_dir]}; "+
        conda_init_cmd+"{config[funannotate_conda_env]}; "
        "tmhmm --workdir {params.wdir} --short {input} > {output}"


rule phobius:       #phobius.pl v1.01
    input:
        f"{OUTDIR}/signalp/signalp.mature_prot.union.faa"
    output:
        f"{OUTDIR}/signalp/transmembrane_domain/phobius.txt"
    shell:
        "export PATH=$PATH:{config[phobius_dir]}; "+
        # awk to remove stop codon * in all sequences because phobius doesn't like it! why even, i want to bark. 
        """
        awk '/^>/ {{print}} !/^>/ {{sub(/\*$/, ""); print}}' {input} | phobius.pl -short | tail -n +2 > {output}
        """


rule filter_SP_for_absence_of_TM_helice:
    input:
        tmhmm = f"{OUTDIR}/signalp/transmembrane_domain/tmhmm.txt",
        phobius = f"{OUTDIR}/signalp/transmembrane_domain/phobius.txt",
        union_idlist = f"{OUTDIR}/signalp/signalp.mature_prot.union.IDlist"
    output:
        f"{OUTDIR}/signalp/signalp.mature_prot.union.TM-filtered.IDlist",
    shell:
        conda_init_cmd+"{config[gffread_conda_env]}; "
        "python scripts/transmembrane_filter.py  --phobius {input.phobius} --tmhmm {input.tmhmm} --signalp_union_id {input.union_idlist} --output {output}"


# 



# extract complete genes parent and all child features
    # "agat_sp_filter_feature_from_keep_list.pl --gff {input.transdecoder_genome_gff3} --keep_list {output.transdecoder_pep_complete_ids_list} -o {output.transdecoder_genome_gff3_complete}"




# rule bambu:

# rule DEseq2:
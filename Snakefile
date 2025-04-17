import os
from glob import glob

configfile: "config.yaml"

REF = config["ref"]
PREASSEMBLED_TRANSCRIPTS = config["preassembled_transcripts_fasta"]
PREASSEMBLED_TRANSCRIPTS_GTF = config["preassembled_transcripts_gtf"]
SPORES_GTF = config["spores_gtf"]
INFECTIONS_GTF = config["infections_gtf"]
ONT_CDNA_FILT_Q = config["ONT_cDNA_filt_quality"]
OUTDIR = config["outdir"].rstrip("/")
FUNANNOTATE_DB_DIR = config["funannotate_db_dir"].rstrip("/")
GENEMARK_DIR = config["genemark_dir"].rstrip("/")
EGGNOG_DATA_DIR = config["eggnog_data_dir"].rstrip("/")
CODINGQUARRY_DIR = config["codingquarry_dir"]
FUNANNOTATE_SPECIES_NAME = config["funannotate_species_name"]
FUNANNOTATE_ISOLATE_NAME = config["funannotate_isolate_name"]

conda_init_cmd = "set +eu && . $(conda info --base)/etc/profile.d/conda.sh && conda activate "

DATA_OUTDIR = f"{OUTDIR}/0_data"
FUNANNOTATE_OUTDIR = f"{OUTDIR}/1_funannotate"
CODINGQUARRY_OUTDIR = f"{OUTDIR}/2_codingquarryPM"
TRANSDECODER_SP_RESCUE_OUTDIR = f"{OUTDIR}/3_transdecoder_SP_rescue"
SIGNALP_OUTDIR = f"{OUTDIR}/4_signalp"
TM_OUTDIR = f"{OUTDIR}/5_transmembrane_domain"

SORTED_REF = f"{DATA_OUTDIR}/{os.path.splitext(os.path.basename(REF))[0]}.sorted.fasta"
SPORES_GFF3 = f"{CODINGQUARRY_OUTDIR}/{os.path.splitext(os.path.basename(SPORES_GTF))[0]}.gff3"
INFECTIONS_GFF3 = f"{CODINGQUARRY_OUTDIR}/{os.path.splitext(os.path.basename(INFECTIONS_GTF))[0]}.gff3"
ALL_TRANSCRIPTS_GFF3 = f"{DATA_OUTDIR}/{os.path.splitext(os.path.basename(PREASSEMBLED_TRANSCRIPTS_GTF))[0]}.gff3"
FUNANNOTATE_OUT_PREFIX = f"{'_'.join(config['funannotate_species_name'].split(' '))}_{config['funannotate_isolate_name']}"
infections_prefix = os.path.basename(INFECTIONS_GTF).split('.gtf')[0] 

#######################
##### RULES START #####
#######################


rule all:
    input:
        #f"{SIGNALP_OUTDIR}/funannotate.transdecoder-SP.combined.genes.gff3"
        # f"{TRANSDECODER_SP_RESCUE_OUTDIR}/funannotate_update.transdecoder_SP_TM-filtered.merged.gff3"
        f"{OUTDIR}/{FUNANNOTATE_OUT_PREFIX}.genes.secretome.gff3",
        f"{OUTDIR}/{FUNANNOTATE_OUT_PREFIX}.genes.gff3"


if type(ONT_CDNA_FILT_Q) == int:
    assert int(ONT_CDNA_FILT_Q) > 0
    ONT_CDNA_FASTQ = f"{DATA_OUTDIR}/{os.path.basename(config['ONT_cDNA_fastq']).split('.fastq.gz')[0]}.q{ONT_CDNA_FILT_Q}.fastq.gz"

    rule chopper_filtering:
        input:
            config["ONT_cDNA_fastq"]
        params:
            ONT_CDNA_FILT_Q
        output:
            ONT_CDNA_FASTQ
        threads:
            config["threads"]
        shell:
            conda_init_cmd+"{config[chopper_conda_env]}; "
            """
            chopper -q {params} -i {input} --threads {threads} | gzip > {output}
            """
else:
    ONT_CDNA_FASTQ = config["ONT_cDNA_fastq"]




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
        rna_bam = f"{FUNANNOTATE_OUTDIR}/training/funannotate_train.coordSorted.bam",
        pasa_gff = f"{FUNANNOTATE_OUTDIR}/training/funannotate_train.pasa.gff3"
    params:
        outdir = directory(FUNANNOTATE_OUTDIR)
    threads:
        config["threads"]
    shell:
        conda_init_cmd+"{config[funannotate_conda_env]}; "
        "funannotate train -i {input.ref}  --nanopore_cdna {input.fastq} --trinity {input.preassembled_transcripts}  --jaccard_clip  --species '{config[funannotate_species_name]}'  --isolate {config[funannotate_isolate_name]} --no_trimmomatic --cpus {threads} --memory 100G --pasa_db mysql -o {params.outdir}"


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
        spores_outdir = f"{CODINGQUARRY_OUTDIR}/spores",
        cq_gff3_list = expand(f"{CODINGQUARRY_OUTDIR}/spores/out/{{types}}.gff3", types=["DubiousSet", "PGN_predictedPass", "PredictedPass"])
    output:
        f"{CODINGQUARRY_OUTDIR}/{FUNANNOTATE_ISOLATE_NAME}.CodingQuarry-PM.spores.gff3"
    shell:
        conda_init_cmd+"{config[codingquarry_conda_env]}; "
        """
        export PATH=$PATH:{config[signalp4_dir]}; 
        export PATH=$PATH:{config[codingquarry_dir]}; 
        mkdir -p {params.spores_outdir}; 
        cp -r {config[codingquarry_dir]}/QuarryFiles {params.spores_outdir}/.; 
        export QUARRY_PATH={params.spores_outdir}/QuarryFiles; 
        cd {params.spores_outdir}; 
        run_CQ-PM_unstranded.sh {input.ref} {input.spores_gff3} && cat {params.cq_gff3_list} > {output}
        """


rule codingquarryPM_infections:
    input: 
        ref = SORTED_REF,
        infections_gff3 = INFECTIONS_GFF3
    params:
        infections_outdir = f"{CODINGQUARRY_OUTDIR}/infections",
        cq_gff3_list = expand(f"{CODINGQUARRY_OUTDIR}/infections/out/{{types}}.gff3", types=["DubiousSet", "PGN_predictedPass", "PredictedPass"])
    output:
        f"{CODINGQUARRY_OUTDIR}/{FUNANNOTATE_ISOLATE_NAME}.CodingQuarry-PM.infections.gff3"
    shell:
        conda_init_cmd+"{config[codingquarry_conda_env]}; "
        """
        export PATH=$PATH:{config[signalp4_dir]}; 
        export PATH=$PATH:{config[codingquarry_dir]}; 
        mkdir -p {params.infections_outdir}; 
        cp -r {config[codingquarry_dir]}/QuarryFiles {params.infections_outdir}/.; 
        export QUARRY_PATH={params.infections_outdir}/QuarryFiles; 
        cd {params.infections_outdir}; 
        run_CQ-PM_unstranded.sh {input.ref} {input.infections_gff3} && cat {params.cq_gff3_list} > {output}
        """


rule funannotate_predict:
    input:
        ref = SORTED_REF,
        fastq = ONT_CDNA_FASTQ,
        preassembled_transcripts = PREASSEMBLED_TRANSCRIPTS,
        rna_bam = f"{FUNANNOTATE_OUTDIR}/training/funannotate_train.coordSorted.bam",
        pasa_gff = f"{FUNANNOTATE_OUTDIR}/training/funannotate_train.pasa.gff3",
        infections_transcripts = f"{CODINGQUARRY_OUTDIR}/{FUNANNOTATE_ISOLATE_NAME}.CodingQuarry-PM.infections.gff3",
        spores_transcripts = f"{CODINGQUARRY_OUTDIR}/{FUNANNOTATE_ISOLATE_NAME}.CodingQuarry-PM.spores.gff3",
        all_transcripts = ALL_TRANSCRIPTS_GFF3
    params:
        protein_evidence = expand(config["protein_evidence"]),
        weights = "augustus:4 hiq:6 genemark:1 pasa:10 codingquarry:0 snap:1 glimmerhmm:1 proteins:1 transcripts:2",
        outdir = f"{FUNANNOTATE_OUTDIR}",
        funannotate_db_dir = FUNANNOTATE_DB_DIR,
        genemark_dir = GENEMARK_DIR

    output:
        predicted_gff3 = f"{FUNANNOTATE_OUTDIR}/predict_results/{FUNANNOTATE_OUT_PREFIX}.gff3"
    threads:
        config["threads"]
    shell:
        conda_init_cmd+"{config[funannotate_conda_env]}; "
        "export FUNANNOTATE_DB={params.funannotate_db_dir}; "
        "export GENEMARK_PATH={params.genemark_dir}; "
        "export PATH=$PATH:{params.genemark_dir}; "
        "funannotate predict -i {input.ref} -o {params.outdir} --species '{config[funannotate_species_name]}' --isolate {config[funannotate_isolate_name]} --transcript_evidence {input.preassembled_transcripts} --rna_bam {input.rna_bam} --pasa_gff {input.pasa_gff} --other_gff {input.infections_transcripts}:20 {input.spores_transcripts}:10 --transcript_alignments {input.all_transcripts} --protein_evidence {params.protein_evidence} --weights {params.weights} --optimize_augustus --repeats2evm --ploidy 1 --cpus {threads} --name {config[funannotate_locus_tag]}"


rule funannotate_update:
    input:
        ref = SORTED_REF,
        fastq = ONT_CDNA_FASTQ,
        predicted_gff3 = f"{FUNANNOTATE_OUTDIR}/predict_results/{FUNANNOTATE_OUT_PREFIX}.gff3",
        preassembled_transcripts = PREASSEMBLED_TRANSCRIPTS,
    params:
        outdir = f"{FUNANNOTATE_OUTDIR}"
    output:
        updated_gff3 = f"{FUNANNOTATE_OUTDIR}/update_results/{FUNANNOTATE_OUT_PREFIX}.gff3"
    threads:
        config["threads"]
    shell:
        conda_init_cmd+"{config[funannotate_conda_env]}; "
        "funannotate update -f {input.ref} -o {params.outdir} -g {input.predicted_gff3} --nanopore_cdna {input.fastq} --trinity {input.preassembled_transcripts} --jaccard_clip --no_trimmomatic --cpus {threads} --memory 100G --pasa_db mysql --species '{config[funannotate_species_name]}'  --isolate {config[funannotate_isolate_name]} --min_protlen 50 --alt_transcripts 1 --name {config[funannotate_locus_tag]}"


rule transdecoder_preprocess_infection_transcripts:
    input:
        ref = SORTED_REF,
        infection_transcripts_gtf = INFECTIONS_GTF
    output:
        infection_transcripts_fasta = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transdecoder/{infections_prefix}.fasta",
        infection_transcripts_gff3 = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transdecoder/{infections_prefix}.gff3"
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
        f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transdecoder/{infections_prefix}.fasta"
    output:
        directory(f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transdecoder/{infections_prefix}.fasta.transdecoder_dir"),
        f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transdecoder/{infections_prefix}.fasta.transdecoder.gff3",
        f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transdecoder/{infections_prefix}.fasta.transdecoder.pep"
    params:
        f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transdecoder"
    shell:
        conda_init_cmd+"{config[funannotate_conda_env]}; "
        """
        TransDecoder.LongOrfs -m 50 -t {input} -O {params}; 
        TransDecoder.Predict -t {input} --single_best_only -O {params}
        """

# generate genome-based transcript gff3 to merge detected SP annotations with funannotate gff3 later 
rule transdecoder_generate_genome_based_gff3:
    input:
        infections_transcripts_transdecoder_gff3 = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transdecoder/{infections_prefix}.fasta.transdecoder.gff3",
        infection_transcripts_gff3 = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transdecoder/{infections_prefix}.gff3",
        infection_transcripts_fasta = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transdecoder/{infections_prefix}.fasta"
    output:
        f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transdecoder/{infections_prefix}.fasta.transdecoder.genome.gff3"
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
        transdecoder_pep = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transdecoder/{infections_prefix}.fasta.transdecoder.pep",
        transdecoder_genome_gff3 = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transdecoder/{infections_prefix}.fasta.transdecoder.genome.gff3"
    output:
        transdecoder_pep_complete = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transdecoder/{infections_prefix}.fasta.transdecoder.pep.complete.fasta",
        transdecoder_pep_complete_ids_list = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transdecoder/{infections_prefix}.fasta.transdecoder.pep.complete.IDlist"
    shell:
        conda_init_cmd+"{config[agat_conda_env]}; "
        # get complete transdecoder protein sequences
        "bash scripts/filter_transdecoder_complete_orfs.sh {input.transdecoder_pep} {output.transdecoder_pep_complete}; "
        # extract all headers from the complete protein sequences fasta
        "awk '/^>/ {{print substr($1,2)}}' {output.transdecoder_pep_complete} > {output.transdecoder_pep_complete_ids_list} "


rule signalp3_transdecoder:
    input:
        f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transdecoder/{infections_prefix}.fasta.transdecoder.pep.complete.fasta"
    params:
        signalp3 = config["signalp3_bin"]
    output:
        flatten = temp(directory(f"{TRANSDECODER_SP_RESCUE_OUTDIR}/signalp/flatten")),
        signalp3_out = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/signalp/signalp3.out",
        mature_prot = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/signalp/signalp3.mature_prot.faa"
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
        f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transdecoder/{infections_prefix}.fasta.transdecoder.pep.complete.fasta"
    params:
        signalp4 = config["signalp4_bin"]
    output:
        mature_prot_tmp1 = temp(f"{TRANSDECODER_SP_RESCUE_OUTDIR}/signalp/signalp4.mature_prot.tmp1"),
        mature_prot_tmp2 = temp(f"{TRANSDECODER_SP_RESCUE_OUTDIR}/signalp/signalp4.mature_prot.tmp2"),
        mature_prot = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/signalp/signalp4.mature_prot.faa",
        signalp4_out = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/signalp/signalp4.out"
    shell:
        """
        {params.signalp4} -f short -t euk -s best -m {output.mature_prot_tmp1} {input} > {output.signalp4_out}; 
        awk '/^>/ {{print $1}} !/^>/ {{print}}' {output.mature_prot_tmp1} > {output.mature_prot_tmp2}; 
        awk '!/^>/ {{printf "%s", $0; n="\\n"}} /^>/ {{print n$0; n=""}} END {{printf "%s",n}}' {output.mature_prot_tmp2} > {output.mature_prot}
        """ # unwraps fasta



# slow-sequential, cpu mode 
rule signalp6_transdecoder:
    input:
        f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transdecoder/{infections_prefix}.fasta.transdecoder.pep.complete.fasta"
    params:
        outdir = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/signalp/signalp6-slow-sequential",
        results = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/signalp/signalp6-slow-sequential/prediction_results.txt",
        mature_prot = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/signalp/signalp6-slow-sequential/processed_entries.fasta"
    output:
        mature_prot = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/signalp/signalp6.mature_prot.faa",
        signalp6_out = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/signalp/signalp6.out"
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
        f"{TRANSDECODER_SP_RESCUE_OUTDIR}/signalp/signalp3.mature_prot.faa",
        f"{TRANSDECODER_SP_RESCUE_OUTDIR}/signalp/signalp4.mature_prot.faa",
        f"{TRANSDECODER_SP_RESCUE_OUTDIR}/signalp/signalp6.mature_prot.faa"
    params:
        f"{TRANSDECODER_SP_RESCUE_OUTDIR}/signalp/signalp*.mature_prot.faa"
    output:
        f"{TRANSDECODER_SP_RESCUE_OUTDIR}/signalp/signalp.mature_prot.union.IDlist"
    shell:
        """
        cat {params} | awk 'sub(/^>/, "")' | sort | uniq > {output}
        """



rule combine_signalp_unions_mature_proteins_fasta_transdecoder:
    input:
        f"{TRANSDECODER_SP_RESCUE_OUTDIR}/signalp/signalp3.mature_prot.faa",
        f"{TRANSDECODER_SP_RESCUE_OUTDIR}/signalp/signalp4.mature_prot.faa",
        f"{TRANSDECODER_SP_RESCUE_OUTDIR}/signalp/signalp6.mature_prot.faa",
        union_idlist = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/signalp/signalp.mature_prot.union.IDlist"
    params:
        mature_prot_file_list = expand(f"{TRANSDECODER_SP_RESCUE_OUTDIR}/signalp/signalp{{versions}}.mature_prot.faa", versions=["3", "4", "6"]),
        report_file = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/signalp/union_report.txt"
    output:
        f"{TRANSDECODER_SP_RESCUE_OUTDIR}/signalp/signalp.mature_prot.union.faa"
    shell:
        conda_init_cmd+"{config[gffread_conda_env]}; "
        "python scripts/combine_signalp_unions_mature_prot_faa.py {input.union_idlist} {params.mature_prot_file_list} {output} --report_file {params.report_file}"


rule tmmhm_transdecoder:     #v2.0c
    input:
        f"{TRANSDECODER_SP_RESCUE_OUTDIR}/signalp/signalp.mature_prot.union.faa"
    output:
        f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transmembrane_domain/tmhmm.txt"
    params:
        wdir = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transmembrane_domain"
    shell:
        "export PATH=$PATH:{config[tmhmm_dir]}; "+
        conda_init_cmd+"{config[funannotate_conda_env]}; "
        "tmhmm --workdir {params.wdir} --short {input} > {output}"


rule phobius_transdecoder:       #phobius.pl v1.01
    input:
        f"{TRANSDECODER_SP_RESCUE_OUTDIR}/signalp/signalp.mature_prot.union.faa"
    output:
        f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transmembrane_domain/phobius.txt"
    shell:
        "export PATH=$PATH:{config[phobius_dir]}; "+
        # awk to remove stop codon * in all sequences because phobius doesn't like it! why even, i want to bark. 
        """
        awk '/^>/ {{print}} !/^>/ {{sub(/\*$/, ""); print}}' {input} | phobius.pl -short | tail -n +2 > {output}
        """


rule filter_transdecoder_SP_for_absence_of_TM_helice:
    input:
        tmhmm = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transmembrane_domain/tmhmm.txt",
        phobius = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transmembrane_domain/phobius.txt",
        union_idlist = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/signalp/signalp.mature_prot.union.IDlist"
    output:
        f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transdecoder.signalp.mature_prot.union.TM-filtered.IDlist",
    shell:
        conda_init_cmd+"{config[gffread_conda_env]}; "
        "python scripts/transmembrane_filter.py  --phobius {input.phobius} --tmhmm {input.tmhmm} --signalp_union_id {input.union_idlist} --output {output}"


rule extract_transdecoder_SP_noTM_genes_gff3:    # plus some attribute and name clean up
    input:
        transdecoder_genome_gff3 = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transdecoder/{infections_prefix}.fasta.transdecoder.genome.gff3",
        signalp_union_idlist = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transdecoder.signalp.mature_prot.union.TM-filtered.IDlist"
    output:
        tmp1 = temp(f"{TRANSDECODER_SP_RESCUE_OUTDIR}/{infections_prefix}.fasta.transdecoder.genome.SP_TM-filtered.gff3.tmp1"),
        tmp2 = temp(f"{TRANSDECODER_SP_RESCUE_OUTDIR}/{infections_prefix}.fasta.transdecoder.genome.SP_TM-filtered.gff3.tmp2"),
        gff3 = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transdecoder_SP_TM-filtered.gff3"
    params:
        TRANSDECODER_SP_RESCUE_OUTDIR
    shell: 
        "python scripts/extract_transdecoder_SP_transcripts.py --in_gff3 {input.transdecoder_genome_gff3} --in_secretome_ids {input.signalp_union_idlist} --out_gff3 {output.tmp1}; "+
        conda_init_cmd+"{config[gffread_conda_env]}; "
        "cd {params}; "
        "agat_sp_manage_attributes.pl -f {output.tmp1} --att Name -o {output.tmp2}; "
        "sed -e 's/\^.*\^-//g' -e 's/\^.*\^+//g' {output.tmp2} > {output.gff3}"


# this step merges funannotate_update gff3 and transdecoder SP gff3 together. 
# note that the agat merge step only merges completely identical isoforms including the non-coding exons structures
rule merge_transdecoder_SP_gff3_with_funannotate_gff3:
    input:
        ref = SORTED_REF,
        funannotate_update_gff3 = f"{FUNANNOTATE_OUTDIR}/update_results/{FUNANNOTATE_OUT_PREFIX}.gff3",
        transdecoder_genome_SP_gff3 = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/transdecoder_SP_TM-filtered.gff3"
    output:
        tmp = temp(f"{TRANSDECODER_SP_RESCUE_OUTDIR}/funannotate_update.transdecoder_SP_TM-filtered.merged.tmp"),
        gff3 = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/funannotate_update.transdecoder_SP_TM-filtered.merged.gff3"
    params:
        TRANSDECODER_SP_RESCUE_OUTDIR
    shell:
        conda_init_cmd+"{config[gffread_conda_env]}; "
        "cd {params}; "
        "agat_sp_merge_annotations.pl -f {input.funannotate_update_gff3} -f {input.transdecoder_genome_SP_gff3} -o {output.tmp}; "
        "conda deactivate; "+
        conda_init_cmd+"{config[funannotate_conda_env]}; "
        "funannotate gff-rename -g {output.tmp} -f {input.ref} -o {output.gff3} --locus_tag {config[funannotate_locus_tag]} --numbering {config[funannotate_numbering_start]}"


# deduplicate isoforms that have identical CDS. keeps the first transcript ID (i.e. smallest number after "-T") that it encounters in every unique CDS coordinate set.
# outputs two files:
# (1) deduplicated CDS gff3, useful for expression analysis; 
# (2) full gff3 deduplicated using the CDS ids from (1), thus including other features eg genes, exons, UTRs and tRNA entries.
rule deduplicate_cds_gff3:
    input:
        f"{TRANSDECODER_SP_RESCUE_OUTDIR}/funannotate_update.transdecoder_SP_TM-filtered.merged.gff3"
    output:
        cds = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/dedup/funannotate_update.transdecoder_SP_TM-filtered.merged.dedup.cds.gff3",
        full = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/dedup/funannotate_update.transdecoder_SP_TM-filtered.merged.dedup.full.gff3"
    shell:
        "python scripts/deduplicate_CDS_gff.py -i {input} --dedup_cds_gff3_out {output.cds} --dedup_full_gff3_out {output.full}"


rule get_protein_and_cds_sequences_from_merged_gff3:
    input:
        ref = SORTED_REF,
        dedup_cds_gff3 = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/dedup/funannotate_update.transdecoder_SP_TM-filtered.merged.dedup.cds.gff3",
        dedup_full_gff3 = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/dedup/funannotate_update.transdecoder_SP_TM-filtered.merged.dedup.full.gff3"
    output:
        folded_ref = f"{DATA_OUTDIR}/{os.path.splitext(os.path.basename(REF))[0]}.sorted.folded.fasta",
        prot_tmp = temp(f"{TRANSDECODER_SP_RESCUE_OUTDIR}/dedup/funannotate_update.transdecoder_SP_TM-filtered.merged.dedup.protein.tmp"),
        cds_tmp = temp(f"{TRANSDECODER_SP_RESCUE_OUTDIR}/dedup/funannotate_update.transdecoder_SP_TM-filtered.merged.dedup.cds.tmp"),
        prot = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/dedup/funannotate_update.transdecoder_SP_TM-filtered.merged.dedup.protein.faa",
        cds = f"{TRANSDECODER_SP_RESCUE_OUTDIR}/dedup/funannotate_update.transdecoder_SP_TM-filtered.merged.dedup.cds.fna",
        symlink_prot = f"{SIGNALP_OUTDIR}/funannotate.transdecoder-SP.combined.protein.faa",
        symlink_cds = f"{SIGNALP_OUTDIR}/funannotate.transdecoder-SP.combined.cds.fna",
        symlink_cds_gff3 = f"{SIGNALP_OUTDIR}/funannotate.transdecoder-SP.combined.cds.gff3",
        symlink_full_gff3 = f"{SIGNALP_OUTDIR}/funannotate.transdecoder-SP.combined.genes.gff3"
    params:
        f"{SIGNALP_OUTDIR}"
    shell:
        conda_init_cmd+"{config[agat_conda_env]}; "+
        """
        fold {input.ref} > {output.folded_ref}; 
        cd {params}; 
        agat_sp_extract_sequences.pl --gff {input.dedup_cds_gff3} -f {output.folded_ref} -o {output.prot_tmp} -t cds -p ; 
        agat_sp_extract_sequences.pl --gff {input.dedup_cds_gff3} -f {output.folded_ref} -o {output.cds_tmp} -t cds ;
        awk '!/^>/ {{printf "%s", $0; n="\\n"}} /^>/ {{print n$0; n=""}} END {{printf "%s",n}}' {output.prot_tmp} > {output.prot}; 
        awk '!/^>/ {{printf "%s", $0; n="\\n"}} /^>/ {{print n$0; n=""}} END {{printf "%s",n}}' {output.cds_tmp} > {output.cds}; 
        """ #unwrap fasta
        "mkdir -p {params}; "
        "ln -sr {output.prot} {output.symlink_prot}; "
        "ln -sr {output.cds} {output.symlink_cds}; "
        "ln -sr {input.dedup_cds_gff3} {output.symlink_cds_gff3}; "
        "ln -sr {input.dedup_full_gff3} {output.symlink_full_gff3}"


rule signalp3:
    input:
        f"{SIGNALP_OUTDIR}/funannotate.transdecoder-SP.combined.protein.faa"
    params:
        signalp3 = config["signalp3_bin"]
    output:
        flatten = temp(directory(f"{SIGNALP_OUTDIR}/flatten")),
        signalp3_out = f"{SIGNALP_OUTDIR}/signalp3.out",
        mature_prot = f"{SIGNALP_OUTDIR}/signalp3.mature_prot.faa"
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
        f"{SIGNALP_OUTDIR}/funannotate.transdecoder-SP.combined.protein.faa"
    params:
        signalp4 = config["signalp4_bin"]
    output:
        mature_prot_tmp1 = temp(f"{SIGNALP_OUTDIR}/signalp4.mature_prot.tmp1"),
        mature_prot_tmp2 = temp(f"{SIGNALP_OUTDIR}/signalp4.mature_prot.tmp2"),
        mature_prot = f"{SIGNALP_OUTDIR}/signalp4.mature_prot.faa",
        signalp4_out = f"{SIGNALP_OUTDIR}/signalp4.out"
    shell:
        """
        {params.signalp4} -f short -t euk -s best -m {output.mature_prot_tmp1} {input} > {output.signalp4_out}; 
        awk '/^>/ {{print $1}} !/^>/ {{print}}' {output.mature_prot_tmp1} > {output.mature_prot_tmp2}; 
        awk '!/^>/ {{printf "%s", $0; n="\\n"}} /^>/ {{print n$0; n=""}} END {{printf "%s",n}}' {output.mature_prot_tmp2} > {output.mature_prot}
        """ # unwraps fasta


# slow-sequential, cpu mode 
rule signalp6:
    input:
        f"{SIGNALP_OUTDIR}/funannotate.transdecoder-SP.combined.protein.faa"
    params:
        outdir = f"{SIGNALP_OUTDIR}/signalp6-slow-sequential",
        results = f"{SIGNALP_OUTDIR}/signalp6-slow-sequential/prediction_results.txt",
        mature_prot = f"{SIGNALP_OUTDIR}/signalp6-slow-sequential/processed_entries.fasta"
    output:
        mature_prot = f"{SIGNALP_OUTDIR}/signalp6.mature_prot.faa",
        signalp6_out = f"{SIGNALP_OUTDIR}/signalp6.out"
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
        expand(f"{SIGNALP_OUTDIR}/signalp{{versions}}.mature_prot.faa", versions=["3", "4", "6"])
    params:
        f"{SIGNALP_OUTDIR}/signalp*.mature_prot.faa"
    output:
        f"{SIGNALP_OUTDIR}/signalp.mature_prot.union.IDlist"
    shell:
        """
        cat {params} | awk 'sub(/^>/, "")' | sort | uniq > {output}
        """


rule combine_signalp_unions_mature_proteins_fasta:
    input:
        union_idlist = f"{SIGNALP_OUTDIR}/signalp.mature_prot.union.IDlist",
        mature_prot_file_list = expand(f"{SIGNALP_OUTDIR}/signalp{{versions}}.mature_prot.faa", versions=["3", "4", "6"])
    params:
        report_file = f"{SIGNALP_OUTDIR}/union_report.txt"
    output:
        f"{SIGNALP_OUTDIR}/signalp.mature_prot.union.faa"
    shell:
        conda_init_cmd+"{config[gffread_conda_env]}; "
        "python scripts/combine_signalp_unions_mature_prot_faa.py {input.union_idlist} {input.mature_prot_file_list} {output} --report_file {params.report_file}"


rule tmmhm:     #v2.0c
    input:
        f"{SIGNALP_OUTDIR}/signalp.mature_prot.union.faa"
    output:
        f"{TM_OUTDIR}/tmhmm.txt"
    params:
        wdir = f"{TM_OUTDIR}"
    shell:
        "export PATH=$PATH:{config[tmhmm_dir]}; "+
        conda_init_cmd+"{config[funannotate_conda_env]}; "
        "tmhmm --workdir {params.wdir} --short {input} > {output}"


rule phobius:       #phobius.pl v1.01
    input:
        f"{SIGNALP_OUTDIR}/signalp.mature_prot.union.faa"
    output:
        f"{TM_OUTDIR}/phobius.txt"
    shell:
        "export PATH=$PATH:{config[phobius_dir]}; "+
        """
        awk '/^>/ {{print}} !/^>/ {{sub(/\*$/, ""); print}}' {input} | phobius.pl -short | tail -n +2 > {output}
        """


rule filter_SP_for_absence_of_TM_helice:
    input:
        tmhmm = f"{TM_OUTDIR}/tmhmm.txt",
        phobius = f"{TM_OUTDIR}/phobius.txt",
        union_idlist = f"{SIGNALP_OUTDIR}/signalp.mature_prot.union.IDlist"
    output:
        f"{SIGNALP_OUTDIR}/signalp.mature_prot.union.TM-filtered.IDlist",
    shell:
        conda_init_cmd+"{config[gffread_conda_env]}; "
        "python scripts/transmembrane_filter.py  --phobius {input.phobius} --tmhmm {input.tmhmm} --signalp_union_id {input.union_idlist} --output {output}"


rule extract_secretome_genes_gff3:
    input:
        gff3 = f"{SIGNALP_OUTDIR}/funannotate.transdecoder-SP.combined.genes.gff3",
        signalp_union_idlist = f"{SIGNALP_OUTDIR}/signalp.mature_prot.union.TM-filtered.IDlist"
    output:
        gff3 = f"{OUTDIR}/{FUNANNOTATE_OUT_PREFIX}.genes.secretome.gff3",
        cds_gff3 = f"{OUTDIR}/{FUNANNOTATE_OUT_PREFIX}.genes.secretome.cds.gff3"
    params:
        SIGNALP_OUTDIR
    shell:
        "python scripts/extract_secretome_gff3.py --in_gff3 {input.gff3} --in_secretome_ids {input.signalp_union_idlist} --out_gff3 {output.gff3}; "
        "grep -P '\tCDS\t' {output.gff3} > {output.cds_gff3}"
    

rule copy_all_genes_gff3:
    input:
        f"{SIGNALP_OUTDIR}/funannotate.transdecoder-SP.combined.genes.gff3"
    output:
        gff3 = f"{OUTDIR}/{FUNANNOTATE_OUT_PREFIX}.genes.gff3",
        cds_gff3 = f"{OUTDIR}/{FUNANNOTATE_OUT_PREFIX}.genes.cds.gff3"
    shell:
        "cp {input} {output.gff3}; "
        "grep -P '\tCDS\t' {output.gff3} > {output.cds_gff3}"



# rule bambu:

# rule DEseq2:
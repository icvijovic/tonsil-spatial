rule reads_to_fasta:
    input:
        "{base}/vdj/UMI_collapsed_on_tissue_read_list.csv.gz",
    output:
        "{base}/vdj/all_consensus_sequences.fasta",
    conda:
        "pacbio"#"../envs/pacbio.yaml",
    log:
        "{base}/logs/reads_to_fasta.log",
    resources:
        mem_mb="16000",
        threads=1,
    threads: 1,
    shell:
        "python scripts/collapsed_reads_to_fasta.py -i {input} -o {output} 2> {log}"

rule igblast:
    input:
        rules.reads_to_fasta.output,
    output:
        "{base}/vdj/igblast.tsv",
    conda:
        "pacbio"#"../envs/pacbio.yaml",
    resources:
        threads=16,
    params:
        organism="human",
        IGDBDIR=config["IGDBDIR"],
    shell:
        """
        ml load system libuv

        export IGDATA={params.IGDBDIR}
        
        igblastn \
                -germline_db_V {params.IGDBDIR}/database/imgt_{params.organism}_ig_v \
                -germline_db_D {params.IGDBDIR}/database/imgt_{params.organism}_ig_d \
                -germline_db_J {params.IGDBDIR}/database/imgt_{params.organism}_ig_j \
                -auxiliary_data {params.IGDBDIR}/optional_file/{params.organism}_gl.aux \
                -domain_system imgt \
                -ig_seqtype Ig \
                -organism {params.organism} \
                -outfmt 19 \
                -query {input} \
                -out {output} \
                -num_threads {resources.threads}
        """


rule filter_and_annotate:
    input:
        read_info=rules.reads_to_fasta.input,
        igblast=rules.igblast.output,
    output:
        "{base}/vdj/igblast_filtered_annotated.tsv.gz",
    conda:
        "pacbio"#"../envs/pacbio.yaml"
    log:
        "{base}/logs/filter_vdj.log",
    resources:
        mem_mb="131000",
    params:
        scripts=config["vdj_scripts"],
    shell:
        "python {params.scripts}/filter_igblast_and_annotate.py "
        "{input.igblast} "
        "-read_info {input.read_info} "
        "-outdir {wildcards.base}/vdj  "
        "--verbose "
        "> {log} "
        "2> {log}"



rule call_germlines:
    input:
        rules.filter_and_annotate.output,
    output:
        germlines="{base}/grmlin/igblast_filtered_annotated_v_germlines.tsv",
        preprocessed="{base}/grmlin/igblast_filtered_annotated_preprocessed.tsv.gz",
    log:
        "{base}/logs/grmlin.log",
    conda:
        "germline_inference"#../envs/grmlin.yaml"
    params:
        organism="human", 
        grmlin=config["grmlin"],
        IGDBDIR=config["IGDBDIR"],
    resources:
        mem_mb="131000",
        time="12:00:00",
    shell:
        "{params.grmlin}/grmlin "
        "{input} "
        "-annotate {params.IGDBDIR}/database/imgt_{params.organism}_ig_v "
        "-outdir {wildcards.base}/grmlin "
        "--verbose "
        "-max_sequences 500000 "
        "> {log}"


# the checkpoint that shall trigger re-evaluation of the DAG
checkpoint prepare_cdr3_groups_for_distance_evaluation:
    input:
        rules.call_germlines.output.preprocessed,
    output:
        groups=directory("{base}/lineage_clustering/cdr3"),
    wildcard_constraints:
        donor="TBd[0-9]",
    log:
        "{base}/logs/cluster_lineages/prepare_for_clustering.log",
    conda:
        "pacbio"#../envs/pacbio.yaml"
    params:
        scripts=config["vdj_scripts"],
    resources:
        mem_mb="64000",
        time="24:00:00",
    shell:
        """
        mkdir -p {wildcards.base}/lineage_clustering/cdr3 
        
        python {params.scripts}/prepare_distance_matrices.py \
        {input} \
        -outdir {wildcards.base}/lineage_clustering/cdr3 \
        -samplename tonsil_vdjs \
        2> {log}
        
        """


rule fasta_to_hamming_distance_matrix:
    input:
        "{base}/lineage_clustering/cdr3/{group}.fasta",
    output:
        "{base}/lineage_clustering/cdr3/{group}.npy",
    log:
        "{base}/logs/cluster_lineages/matrix_calc/cdr3/{group}.log",
    conda:
        "pacbio" #"../envs/pacbio.yaml"
    params:
        scripts=config["vdj_scripts"],
    resources:
        mem_mb="128000",
        time="24:00:00",
    shell:
        "python {params.scripts}/calc_matrix.py "
        "{input} "
        "{output} "
        "-hamming"


rule fasta_to_levenshtein_distance_matrix:
    input:
        "{base}/lineage_clustering/templated/{group}_{seq_element}.fasta",
    output:
        "{base}/lineage_clustering/templated/{group}_{seq_element}.npy",
    log:
        "{base}/logs/cluster_lineages/matrix_calc/templated/{group}_{seq_element}.log",
    conda:
        "pacbio"#../envs/pacbio.yaml"
    params:
        scripts=config["scripts"],
    resources:
        mem_mb="128000",
        time="48:00:00",
    shell:
        "python {params.scripts}/calc_matrix.py "
        "{input} "
        "{output} "


def aggregate_cdr3_npy_files(wildcards):
    checkpoint_output = checkpoints.prepare_cdr3_groups_for_distance_evaluation.get(
        **wildcards
    ).output.groups

    # small_groups = expand("{base}/lineage_clustering/cdr3/{group}_cdr3.npy",
    #        base=wildcards.base,

    #       group=glob_wildcards(os.path.join(checkpoint_output, "{group}_cdr3.npy")).group)

    large_groups = expand(
        "{base}/lineage_clustering/cdr3/{group}_cdr3.npy",
        base=wildcards.base,
        group=glob_wildcards(
            os.path.join(checkpoint_output, "{group}_cdr3.fasta")
        ).group,
    )
    return large_groups


def aggregate_templated_npy_files(wildcards):
    checkpoint_output = checkpoints.cluster_cdr3s.get(**wildcards).output.clusters
    v_files = expand(
        "{base}/lineage_clustering/templated/{group}_templated_v.npy",
        base=wildcards.base,
        group=glob_wildcards(
            os.path.join(checkpoint_output, "{group}_templated_v.fasta")
        ).group,
    )
    j_files = expand(
        "{base}/lineage_clustering/templated/{group}_templated_j.npy",
        base=wildcards.base,
        group=glob_wildcards(
            os.path.join(checkpoint_output, "{group}_templated_j.fasta")
        ).group,
    )
    return v_files + j_files


# an aggregation over all produced clusters
checkpoint cluster_cdr3s:
    input:
        airr=rules.call_germlines.output.preprocessed,
        npy=aggregate_cdr3_npy_files,
    output:
        airr="{base}/lineage_clustering/templated/tonsil_vdjs_unique_vdjs_cdr3_clusters.tsv.gz",
        clusters=directory("{base}/lineage_clustering/templated/"),

    log:
        "{base}/logs/cluster_cdr3s/cluster_cdr3s.log",
    conda:
        "pacbio"#"../envs/pacbio.yaml"
    params:
        scripts=config["scripts"],
    resources:
        mem_mb="131000",
        time="24:00:00",
    shell:
        """
        mkdir -p {output.clusters}
        rm {output.clusters}/*

        python {params.scripts}/cluster_sequences_based_on_cdr3_identity_uint8.py \
        {input.airr} \
        -outdir {wildcards.base}/lineage_clustering/templated \
        -matrixdir {wildcards.base}/lineage_clustering/cdr3 \
        -samplename tonsil_vdjs \
        2> {log}
        
        """


rule cluster_templated_regions:
    input:
        airr=rules.call_germlines.output.preprocessed,
        cdr3_cluster_table="{base}/lineage_clustering/templated/tonsil_vdjs_unique_vdjs_cdr3_clusters.tsv.gz",
        npy=aggregate_templated_npy_files,
    output:
        airr="{base}/lineage_clustering/final_lineage_ids/tonsil_vdjs.tsv.gz",
    log:
        "{base}/logs/cluster_templated/tonsil_vdjs_cluster_templated.log",
    conda:
        "pacbio"#"../envs/pacbio.yaml"
    params:
        scripts=config["scripts"],
    resources:
        mem_mb="131000",
        time="12:00:00",
    shell:
        "python {params.scripts}/cluster_templated_regions_uint8.py "
        "-airr {input.airr} "
        "-cdr3clusters {input.cdr3_cluster_table} "
        "-outdir {wildcards.base}/lineage_clustering/final_lineage_ids/ "
        "-matrixdir {wildcards.base}/lineage_clustering/templated "
        "-samplename tonsil_vdjs "
        "2> {log}"


rule create_germline_db:
    input:
        rules.call_germlines.output.germlines,
    output:
        "{base}/germline_databases/initial/tonsil_vdjs_v_germlines.fasta",
    log:
        "{base}/logs/makedb/tonsil_vdjs_create_germ_db.log",
    conda:
        "pacbio" #"../envs/pacbio.yaml"
    params:
        scripts=config["scripts"],
    shell:
        "python {params.scripts}/create_germline_database.py "
        "{input} "
        "-outdir {wildcards.base}/germline_databases/initial "
        "-samplename tonsil_vdjs "
        "> {log}"


rule blast_to_germline:
    input:
        seqs=rules.cluster_templated_regions.output,
        db=rules.create_germline_db.output,
    output:
        "{base}/germline_db_vcall/initial/tonsil_vdjs.tsv.gz",
    params:
        scripts=config["scripts"],
    log:
        "{base}/logs/blast_to_germline/tonsil_vdjs_blast_to_germline.log",
    resources:
        mem_mb="65000",
        time="12:00:00",
    conda:
        "pacbio"#"../envs/pacbio.yaml"
    shell:
        "python {params.scripts}/blast_to_germline.py "
        "{input.seqs} "
        "-germline_db {input.db} "
        "-outdir {wildcards.base}/germline_db_vcall/initial "
        "-samplename tonsil_vdjs "
        "2> {log}"


rule polish_germlines:
    input:
        seqs=rules.blast_to_germline.output,
        db=rules.create_germline_db.output,
    output:
        db="{base}/germline_databases/final/tonsil_vdjs_v_germlines_polished.fasta",
    params:
        scripts=config["vdj_scripts"],
        organism='human',
        IGDBDIR=config["IGDBDIR"],
    log:
        "{base}/logs/polish_germlines/tonsil_vdjs_polish_germlines.log",
    conda:
        "pacbio"#"../envs/pacbio.yaml"
    shell:
        "python {params.scripts}/polish_germline_db.py "
        "{input.seqs} "
        "-germline_db {input.db} "
        "-imgt_db {params.IGDBDIR}/fasta/imgt_{params.organism}_ig_v.fasta "
        "-imgt_allele_info {params.IGDBDIR}/internal_data/{params.organism}/{params.organism}.ndm.imgt "
        "-outdir {wildcards.base}/germline_databases/final "
        "2> {log}"


rule realign_to_polished_germline:
    input:
        seqs=rules.cluster_templated_regions.output,
        db=rules.polish_germlines.output.db,
    output:
        "{base}/germline_db_vcall/final/tonsil_vdjs.tsv.gz",
    params:
        scripts=config["scripts"],
    resources:
        mem_mb="65000",
        time="12:00:00",
    log:
        "{base}/logs/blast_to_germline/tonsil_vdjs_blast_to_germline_polished.log",
    conda:
        "pacbio"#"../envs/pacbio.yaml"
    shell:
        "python {params.scripts}/blast_to_germline.py "
        "{input.seqs} "
        "-germline_db {input.db} "
        "-outdir {wildcards.base}/germline_db_vcall/final "
        "-samplename tonsil_vdjs "
        "2> {log}"


rule reads_to_fasta:
    input:
        "{base}/vdj/UMI_collapsed_on_tissue_read_list.csv.gz",
    output:
        "{base}/vdj/all_consensus_sequences.fasta",
    conda:
        "pacbio",#"../envs/pacbio.yaml",
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
        "pacbio",#"../envs/pacbio.yaml",
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
        "pacbio",#"../envs/pacbio.yaml"
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
        germlines="{base}/grmlin/combined_v_germlines.tsv",
        preprocessed="{base}/grmlin/vdjs_combined_preprocessed.tsv.gz",
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
        time="24:00:00",
    shell:
        "{params.grmlin}/grmlin "
        "{input} "
        "-annotate {params.IGDBDIR}/database/imgt_{params.organism}_ig_v "
        "-outdir {wildcards.base}/grmlin "
        "--verbose "
        "-max_sequences 500000 "
        "> {log}"



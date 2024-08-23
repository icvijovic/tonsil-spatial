import os
import pandas as pd
import glob

shell.prefix("set +euo pipefail;")
configfile: "config/path_config.yaml"
base = config["base"]


rule all:
    input:
        expand("{base}/vdj/igblast_filtered_annotated.tsv.gz", base=base)
    params:
        name="all",
        partition="owners",
    resources:
        threads=1,


include: "rules/vdjc.smk"

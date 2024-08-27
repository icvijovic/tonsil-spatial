import os
import pandas as pd
import glob

shell.prefix("set +euo pipefail;")
configfile: "config/path_config.yaml"
base = config["base"]


rule all:
    input:
        expand("{base}/germline_db_vcall/final/tonsil_vdjs.tsv.gz", base=base)
    params:
        name="all",
        partition="owners",
    resources:
        threads=1,


include: "rules/vdjc.smk"

import pandas as pd
import json

# Get  metadata and config
meta = pd.read_csv(config['meta'], sep="\t")

with open("config/map.config") as mapconfig:
    map_config = json.load(mapconfig)

# Get sequencing info
_cols = ["BC","F1","F2"]
meta_seq = meta.set_index("FLI")[_cols]

FLI = meta_seq.index.to_list()


def get_barcode(wildcards):
    return meta_seq.BC.loc[wildcards.FLI]

def get_fli(wildcards):
    return wildcards.FLI




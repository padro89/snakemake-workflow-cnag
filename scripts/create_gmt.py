import pandas as pd

reactome = pd.read_csv(snakemake.input["reactome"], sep = "\t")
reactome_filtered = reactome.iloc[:,[0,3,5]]
reactome_filtered.columns = ["ID","Pathway","Species"]

# with open(snakemake.output["gmt"], "w") as file:
#     print()
import pandas as pd

print(f"Creating GMT file.")
reactome = pd.read_csv(snakemake.input["reactome"], sep = "\t")
filtered = reactome.iloc[:,[0,3,5]]
filtered.columns = ["ID","Pathway","Species"]
species = filtered[filtered.Species.eq(snakemake.config["species"])]
if species.shape[0] == 0:
    print(f"No matches for {snakemake.config['species']}. Creating human GMT and trying to use g:Profiler2 for orthologous human genes")
    species = filtered[filtered.Species.eq("Homo sapiens")]
species = species.iloc[:,[0,1]]
grouped = species.groupby("Pathway")["ID"].apply(' '.join).reset_index()
genes = pd.DataFrame(grouped["ID"].str.split(" ", expand=True))
gmt = pd.concat([grouped["Pathway"], genes], axis = 1)
print(f"GMT created. Saving it to {snakemake.output['gmt']}")
gmt.to_csv(snakemake.output["gmt"], sep ="\t", header = False, index = False)
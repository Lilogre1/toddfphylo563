from pathlib import Path
import pandas as pd
from collections import defaultdict


ORTHOGROUPS_FILE = Path("orthofinder_out/Results_Apr27/Orthogroups/Orthogroups.tsv") 
ALIGNED_DIR = Path("aligned")
OUT_DIR = Path("aligned_renamed")
OUT_DIR.mkdir(exist_ok=True)


df = pd.read_csv(ORTHOGROUPS_FILE, sep="\t")

species_cols = list(df.columns[1:])


gene_to_species = {}

for _, row in df.iterrows():
    for sp in species_cols:
        if pd.isna(row[sp]):
            continue

        genes = str(row[sp]).split(", ")

        for g in genes:
            g = g.strip()
            if g:
                gene_to_species[g] = sp

print(f"Mapped {len(gene_to_species)} genes to species")

def parse_fasta(path):
    seqs = {}
    with open(path) as f:
        name = None
        seq = []

        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(seq)
                name = line[1:]
                seq = []
            else:
                seq.append(line)

        if name is not None:
            seqs[name] = "".join(seq)

    return seqs


for f in ALIGNED_DIR.glob("*_aligned.fa"):
    seqs = parse_fasta(f)

    out_file = OUT_DIR / f.name

    with open(out_file, "w") as out:
        for gene_id, seq in seqs.items():

            gene_id = gene_id.split()[0]

            if gene_id not in gene_to_species:
                raise ValueError(f"Gene not found in Orthogroups.tsv: {gene_id}")

            species = gene_to_species[gene_id]

            out.write(f">{species}\n{seq}\n")

    print(f"Wrote: {out_file}")

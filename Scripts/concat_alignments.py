from pathlib import Path
from collections import defaultdict

aligned_dir = Path("aligned_renamed")
files = sorted(aligned_dir.glob("*_aligned.fa"))

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


def get_species(header):

    gene = header.split()[0]

    return gene.split(".")[0]


def pad(seq, length):
    return seq + "-" * (length - len(seq))


seqs = defaultdict(str)
species_order = None

for f in files:
    current = {}

    fasta = parse_fasta(f)

    for header, seq in fasta.items():
        sp = get_species(header)
        current[sp] = seq

    if species_order is None:
        species_order = sorted(current.keys())
        for sp in species_order:
            seqs[sp] = ""

    # enforce consistent length
    lengths = [len(s) for s in current.values()]
    if len(set(lengths)) != 1:
        raise ValueError(f"Inconsistent sequence lengths in {f}")

    block_len = lengths[0]

    for sp in species_order:
        if sp not in current:
            current[sp] = "-" * block_len

    for sp in species_order:
        seqs[sp] += current[sp]

# sanity check
lengths = {len(v) for v in seqs.values()}
if len(lengths) != 1:
    raise ValueError(f"Final concatenation length mismatch: {lengths}")

with open(r"results/concatenated.fa", "w") as out:
    for sp in species_order:
        out.write(f">{sp}\n{seqs[sp]}\n")

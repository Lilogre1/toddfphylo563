# 2/10/2026
# Dataset Description

I am conducting a phylogenetic analysis of nematodes within the genus *Caenorhabditis* and several closely related taxa, including *Heterorhabditis*, *Mesorhabditis*, and *Rhabditophanes*. The dataset includes 26 species, each represented by a predicted proteome downloaded from WormBase ParaSite. There were a few duplicate species available on Wormbase, but I selected the most recent one available. The redundancy/duplicates were due to multiple different papers/journals contributing to Wormbase. These fasta files originate from published BioProjects and contain annotated protein coding sequences for each genome.

The goal is to identify shared genes across species and use them to reconstruct evolutionary relationships within this clade. Individual loci will be aligned and analyzed to infer both gene trees and a species tree.

---

# Data Extraction and Organization

All proteome files were downloaded as `.fa.gz` archives. To keep the workflow reproducible and organized, I created a dedicated subdirectory inside the project folder to store the extracted FASTA files.

```powershell
#create subfolder for extracted FASTA files
$dest = Join-Path $PWD "extracted"
if (!(Test-Path $dest)) {
    New-Item -ItemType Directory -Path $dest | Out-Null
}

#extract all .gz files into the subfolder
Get-ChildItem -Filter *.gz | ForEach-Object {
    $outfile = Join-Path $dest ($_.BaseName)
    $in = [System.IO.File]::OpenRead($_.FullName)
    $out = [System.IO.File]::Create($outfile)
    $gzip = New-Object System.IO.Compression.GzipStream(
        $in,
        [System.IO.Compression.CompressionMode]::Decompress
    )
    $gzip.CopyTo($out)
    $gzip.Dispose()
    $in.Dispose()
    $out.Dispose()
}
```

This produced a clean directory:

```text
Data/
    extracted/
        caenorhabditis_angaria.PRJNA51225.WBPS19.protein.fa
        caenorhabditis_auriculariae.PRJEB40642.WBPS19.protein.fa
        ...
```

---

# Quality Control (QC) Procedures

As these are already-annotated protein fasta files, I did some rudimentary QC involving sequence count checks and checking for non-protein or blank lines in the .fa files.

---

## 1. Sequence Count Check

This verifies that each proteome contains a reasonable number of protein‑coding genes.

```powershell
Get-ChildItem *.fa | ForEach-Object {
    $count = (Select-String -Pattern "^>" -Path $_.FullName).Count
    Write-Output "$($_.Name) : $count sequences"
}
```

This step helps identify incomplete proteomes

---

## 2. FASTA Formatting Check

The point of this scan is to detect:
- blank lines  
- missing headers  
- malformed entries  

```powershell
Get-ChildItem *.fa | ForEach-Object {
    Write-Output "Checking $_"
    if (Select-String -Pattern "^\s*$" -Path $_.FullName) {
        Write-Output "  Contains blank lines"
    }
    if (-not (Select-String -Pattern "^>" -Path $_.FullName)) {
        Write-Output "  No FASTA headers found"
    }
}
```

This ensures each file is structurally valid before downstream analysis. I did not receive any negative feedback from these checks, with there being well over 10000 sequences for each fasta file. In the future, I would like to apply a consistent renaming standard to FASTA headers to improve downstream compatibility.

### 

3/10/26

I downloaded and used the MAFFT dockerfile as it was the most easily accessible for this task. I was not able to use conda to install and WSL has issues on my local machine, so I used the docker container at https://hub.docker.com/r/staphb/mafft/ to do the following:
```powershell
docker pull staphb/mafft:7.450

```


## Maximum Likelihood Step Update (IQ-TREE)

### Overview

In this step, I attempted to perform maximum likelihood phylogenetic inference using IQ-TREE on a previously generated multiple sequence alignment (`OG0000000_aligned.fa`) produced by MAFFT.

---

### Command Used

First I pulled the docker container via:

```Powershell
docker pull staphb/iqtree2 
```

I ran IQ-TREE via Docker from the directory containing the alignment file:

```Powershell
docker run --rm -v C:/Users/toddf/toddfphylo563/Data/Results_Mar04_1/Orthogroup_Sequences:/data staphb/iqtree2 iqtree2 -s /data/OG0000000_aligned.fa -m MFP -bb 1000 -nt AUTO
```

This command mounts the working directory into the container and runs IQ-TREE with:

* bootstrap support (`-bb 1000`)
* automatic thread detection (`-nt AUTO`)

---

### Issues encountered

The run failed with the following error:

```Powershell
ERROR: Unknown sequence format, please use PHYLIP, FASTA, CLUSTAL, MSF, or NEXUS format
```

# 4/25: Linux and Everything Else
From here on out, I decided to download and use WSL subsystem for Windows to directly run all of the required programs, as opposed to relying on Docker. I decided to do this as it is easier to save progress between runs of complex programs (see: IQtree) when it can take hours to process all data.

## OrthoFinder
The Orthofinder workflow remained the same (though I neglected to include it before):

Image pulled from:
https://hub.docker.com/r/davidemms/orthofinder

Using:
```Powershell
docker pull davidemms/orthofinder
```

and ran with
```Powershell
docker run --rm -v "C:\Users\toddf\toddfphylo563\Data\extracted:/input" -v "C:\Users\toddf\toddfphylo563\Data\orthofinder_out:/output" davidemms/orthofinder orthofinder -f /input -o /output -t 8
```
## Linux
The workflow is ultimately very similar, though, with me beginning with MAFFT downloaded via pixi in a Linux Subsystem window:
```linux
wget -qO- https://pixi.sh/install.sh | sh
```
## MAFFT
```linux
pixi global install mafft
```
I then (temporarily) moved my project files into a mounted wsl folder, in order to provide it with write access.

Accordingly, I then ran MAFFT on the 20 output single-copy orthologues from the OrthoFinder output
```linux
for f in Single_Copy_Orthologue_Sequences/*.fa; do
    base=$(basename "$f" .fa)
    mafft --auto "$f" > "aligned/${base}_aligned.fa"
done
```

# 5/1
## Remapping
After writing aligning, I realized my aligned genes were using internal species ID codes and non-readable stuff. Accordingly, I decided to write a python script to rename the species properly within the single-copy orthologues. I ran it within wsl, using a pyenv installed and enabled by:

```linux
sudo apt update
sudo apt install python3-venv

python3 -m venv phylo_env
source phylo_env/bin/activate
```

speciesmap.py:
```python
from pathlib import Path
import pandas as pd
from collections import defaultdict


ORTHOGROUPS_FILE = "/mnt/c/Users/toddf/toddfphylo563/Data/orthofinder_out/Results_Apr27/Orthogroups/Orthogroups.tsv"
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
```
This remapped the orthologous genes' species tags, so when I eventually concatenate the ortholog fasta files (which I tried to do a little bit before this step-- but had some issue with it), it shows the correct species in the final output trees.

# 5/2
## Concatenation
After renaming the species tags in the orthologues, I wrote the following script to concatenate the orthologues into one.
```python
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

with open("concatenated.fa", "w") as out:
    for sp in species_order:
        out.write(f">{sp}\n{seqs[sp]}\n")

```

## IQ Tree and NJ
After concatenating, I decided to use IQTree and Neighbor Joining because the former is a complex model and the latter is much simpler. I wanted to compare the results of the two.

Following the concatenation, I used IQTree on the Concatenated fasta file, and it selected LG+R4 as the model.
```Linux
iqtree2 -s concatenated.fa -m MFP -bb 1000 -safe
```

```Linux output
Akaike Information Criterion:           LG+R5
Corrected Akaike Information Criterion: LG+R5
Bayesian Information Criterion:         LG+R4
Best-fit model: LG+R4 chosen according to BIC

All model information printed to concatenated.fa.model.gz
CPU time for ModelFinder: 3175.433 seconds (0h:52m:55s)
Wall-clock time for ModelFinder: 3197.517 seconds (0h:53m:17s)
Generating 1000 samples for ultrafast bootstrap (seed: 453534)...

NOTE: 124 MB RAM (0 GB) is required!
Estimate model parameters (epsilon = 0.100)
1. Initial log-likelihood: -226505.812
2. Current log-likelihood: -160742.584
3. Current log-likelihood: -160521.536
4. Current log-likelihood: -160469.290
5. Current log-likelihood: -160424.476
6. Current log-likelihood: -160386.313
7. Current log-likelihood: -160354.500
8. Current log-likelihood: -160328.354
9. Current log-likelihood: -160307.354
10. Current log-likelihood: -160292.043
11. Current log-likelihood: -160283.001
12. Current log-likelihood: -160279.452
13. Current log-likelihood: -160277.044
14. Current log-likelihood: -160276.262
15. Current log-likelihood: -160275.885
16. Current log-likelihood: -160275.171
Optimal log-likelihood: -160275.060
Site proportion and rates:  (0.417,0.070) (0.260,0.496) (0.206,1.515) (0.117,4.526)
Parameters optimization took 16 rounds (43.422 sec)
Computing ML distances based on estimated model parameters... 0.039 sec
Computing BIONJ tree...
0.000 seconds
Log-likelihood of BIONJ tree: -141978.308
--------------------------------------------------------------------
|             INITIALIZING CANDIDATE TREE SET                      |
--------------------------------------------------------------------
Generating 98 parsimony trees... 1.635 second
Computing log-likelihood of 98 initial trees ... 24.428 seconds
Current best score: -141978.308

Do NNI search on 20 best initial trees
Estimate model parameters (epsilon = 0.100)
BETTER TREE FOUND at iteration 1: -141787.632
Iteration 10 / LogL: -141789.683 / Time: 0h:2m:21s
Iteration 20 / LogL: -141788.029 / Time: 0h:3m:14s
Finish initializing candidate tree set (1)
Current best tree score: -141787.632 / CPU time: 147.223
Number of iterations: 20
--------------------------------------------------------------------
|               OPTIMIZING CANDIDATE TREE SET                      |
--------------------------------------------------------------------
Iteration 30 / LogL: -141788.798 / Time: 0h:3m:53s (0h:9m:31s left)
Iteration 40 / LogL: -141787.883 / Time: 0h:4m:35s (0h:7m:10s left)
Iteration 50 / LogL: -141787.686 / Time: 0h:5m:19s (0h:5m:32s left)
Log-likelihood cutoff on original alignment: -141863.757
UPDATE BEST LOG-LIKELIHOOD: -141787.628
Iteration 60 / LogL: -141788.284 / Time: 0h:5m:57s (0h:4m:8s left)
Iteration 70 / LogL: -141788.841 / Time: 0h:6m:40s (0h:3m:0s left)
Iteration 80 / LogL: -141788.637 / Time: 0h:7m:22s (0h:1m:57s left)
UPDATE BEST LOG-LIKELIHOOD: -141787.627
Iteration 90 / LogL: -141788.754 / Time: 0h:8m:3s (0h:0m:59s left)
Iteration 100 / LogL: -141787.946 / Time: 0h:8m:44s (0h:0m:5s left)
Log-likelihood cutoff on original alignment: -141854.805
NOTE: Bootstrap correlation coefficient of split occurrence frequencies: 1.000
TREE SEARCH COMPLETED AFTER 102 ITERATIONS / Time: 0h:8m:52s

--------------------------------------------------------------------
|                    FINALIZING TREE SEARCH                        |
--------------------------------------------------------------------
Performs final model parameters optimization
Estimate model parameters (epsilon = 0.010)
1. Initial log-likelihood: -141787.627
2. Current log-likelihood: -141787.577
3. Current log-likelihood: -141787.549
4. Current log-likelihood: -141787.528
Optimal log-likelihood: -141787.515
Site proportion and rates:  (0.357,0.088) (0.348,0.609) (0.234,1.822) (0.062,5.359)
Parameters optimization took 4 rounds (7.402 sec)
BEST SCORE FOUND : -141787.515
Creating bootstrap support values...
Split supports printed to NEXUS file concatenated.fa.splits.nex
Total tree length: 7.143

Total number of iterations: 102
CPU time used for tree search: 481.184 sec (0h:8m:1s)
Wall-clock time used for tree search: 485.801 sec (0h:8m:5s)
Total CPU time used: 534.504 sec (0h:8m:54s)
Total wall-clock time used: 540.469 sec (0h:9m:0s)

Computing bootstrap consensus tree...
Reading input file concatenated.fa.splits.nex...
26 taxa and 65 splits.
Consensus tree written to concatenated.fa.contree
Reading input trees file concatenated.fa.contree
Log-likelihood of consensus tree: -141787.512

Analysis results written to:
  IQ-TREE report:                concatenated.fa.iqtree
  Maximum-likelihood tree:       concatenated.fa.treefile
  Likelihood distances:          concatenated.fa.mldist

Ultrafast bootstrap approximation results written to:
  Split support values:          concatenated.fa.splits.nex
  Consensus tree:                concatenated.fa.contree
  Screen log file:               concatenated.fa.log

Date and Time: Monday May  3 00:01:05 2026
```

This output the file 
concatenated.fa.treefile

which I re-rooted at the outgroup using the following script, rerooting_tree.py:
```python
from ete3 import Tree

t = Tree("concatenated.fa.treefile")

t.set_outgroup("rhabditophanes_kr3021")

t.write(outfile="rooted.treefile")
```
The new Newick tree is as follows:
```Newick IQTree Tree
(caenorhabditis_angaria:0.1578840406,(((caenorhabditis_auriculariae:0.3139830240,(caenorhabditis_parvicauda:0.4060425325,(heterorhabditis_bacteriophora:0.7473087454,((mesorhabditis_belari:0.2430178497,mesorhabditis_spiculigera:0.3358089989)100:0.3547388954,rhabditophanes_kr3021:1.5200132479)90:0.0785528294)100:0.1605080610)80:0.0466091351)100:0.1186013992,caenorhabditis_bovis:0.2410634839)88:0.0377340757,(((((caenorhabditis_becei:0.0675235884,(caenorhabditis_panamensis:0.0704412208,caenorhabditis_waitukubuli:0.0884971788)100:0.0242727941)100:0.0700039812,caenorhabditis_japonica:0.1843647796)83:0.0181272170,((((caenorhabditis_brenneri:0.0801094075,caenorhabditis_tropicalis:0.0998099043)100:0.0235878016,(((caenorhabditis_briggsae:0.0082673918,caenorhabditis_nigoni:0.0069717822)100:0.0788600085,((caenorhabditis_sinica:0.0703064580,caenorhabditis_tribulationis:0.0455868541)73:0.0083777419,caenorhabditis_zanzibari:0.0357364687)100:0.0539599183)100:0.0393019884,(caenorhabditis_latens:0.0143205989,caenorhabditis_remanei:0.0105675756)100:0.0717334021)100:0.0290901382)100:0.0326751146,caenorhabditis_elegans:0.1067970591)57:0.0119434688,caenorhabditis_inopinata:0.1223135069)100:0.0741032276)77:0.0244697481,caenorhabditis_sulstoni:0.1663692615)100:0.1634655699,caenorhabditis_uteleia:0.1978684135)100:0.0586912158)100:0.1576955918,caenorhabditis_quiockensis:0.0647318457);

```

After running this, I wrote an R script to do neighbor-joining
```r
library(ape)
library(phangorn)

#read concatenated protein alignment
alignment <- read.phyDat(
  "C:/Users/toddf/toddfphylo563/concatenated.fa",
  format = "fasta",
  type = "AA"
)

#ML distance matrix
dist_matrix <- dist.ml(alignment)

#neighbor joining tree
tree <- NJ(dist_matrix)

#root on outgroup
tree_rooted <- root(
  tree,
  outgroup = "rhabditophanes_kr3021",
  resolve.root = TRUE
)

# spread out tree visually
plot(
  tree_rooted,
  cex = 0.75,
  no.margin = TRUE,
  
)

#print newick text to console
write.tree(tree_rooted)

#newick file
write.tree(tree_rooted, file = "nj_rooted.treefile")

```
The Newick tree for the NJ tree was as follows:

```Newick Neighbor-Joining Tree
((((((((((caenorhabditis_inopinata:0.08021188621,caenorhabditis_elegans:0.07414184409):0.004431227529,(((((caenorhabditis_tribulationis:0.03578330731,caenorhabditis_zanzibari:0.03194864132):0.002017785265,caenorhabditis_sinica:0.05028905048):0.02679746658,(caenorhabditis_briggsae:0.004075015456,caenorhabditis_nigoni:0.008426233036):0.051658347):0.01175535392,(caenorhabditis_latens:0.01423703576,caenorhabditis_remanei:0.006926797894):0.05170809931):0.009549475991,(caenorhabditis_brenneri:0.05436340621,caenorhabditis_tropicalis:0.07113212588):0.01012281464):0.01027561732):0.02565343805,((((caenorhabditis_panamensis:0.05192246121,caenorhabditis_waitukubuli:0.05561648564):0.005184132109,caenorhabditis_becei:0.05393329461):0.03199446399,caenorhabditis_sulstoni:0.1037946519):0.004192856604,caenorhabditis_japonica:0.111671787):0.003639441717):0.05381312432,caenorhabditis_uteleia:0.1303892532):0.0117448394,(caenorhabditis_angaria:0.06808899306,caenorhabditis_quiockensis:0.06935503947):0.06407894919):0.007127270385,caenorhabditis_bovis:0.148988025):0.03069628836,caenorhabditis_auriculariae:0.1684583972):0.006944373942,caenorhabditis_parvicauda:0.2176447515):0.0450228074,((mesorhabditis_belari:0.149814455,mesorhabditis_spiculigera:0.1810145106):0.1146797136,heterorhabditis_bacteriophora:0.3365126412):0.0238874782):0,rhabditophanes_kr3021:0.488524261);
```

In R, I then compared the nj tree and the IQTree tree

```r
library(ape)

nj_tree <- read.tree(r"(C:\Users\toddf\toddfphylo563\nj_rooted.treefile)")
ml_tree <- read.tree(r"(\\wsl.localhost\Ubuntu\home\toddf\orthofinder_work\MAFFTIQTree\rooted.treefile)")

all.equal(nj_tree, ml_tree)
RF.dist(nj_tree, ml_tree)

```
Which output the following:
```R Output
[1] FALSE
[1] 10
```
And here are images of the respective plots, generated in R, with the image following it being the same plot in TreeViewer:
<img width="585" height="325" alt="image" src="https://github.com/user-attachments/assets/34e12cde-6ab5-4c73-ae52-4b3f7dee7b84" />
<img width="1713" height="971" alt="image" src="https://github.com/user-attachments/assets/57ffe562-97cc-401d-900e-b5683230b42d" />

<img width="567" height="356" alt="image" src="https://github.com/user-attachments/assets/d7575c3e-f48c-4f80-93d4-2ba52e181822" />
<img width="1237" height="690" alt="image" src="https://github.com/user-attachments/assets/d9ea6a7d-d925-4930-a26b-1b206c49910b" />


Software versions:
- IQ-TREE: v2.0.7 (Linux 64-bit build, Apr 1 2024)
- MAFFT: v7.526 (2024-04-26)
- OrthoFinder: v3.0.1b1
- R: v4.4.1 (2024-06-14, "Race for Your Life")
- ape: v5.8.1
- phangorn: v2.12.1

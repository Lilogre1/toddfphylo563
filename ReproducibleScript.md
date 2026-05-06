## Reproducible phylogenetics Workflow

---
### Software versions & Required Packages:


#### IQ-TREE (linux): v2.0.7 (Linux 64-bit build, Apr 1 2024)
#### MAFFT (linux): v7.526 (2024-04-26)
#### OrthoFinder (docker): v3.0.1b1
#### R (windows): v4.4.1 (2024-06-14, "Race for Your Life")
#### Python: v3.12.3
#### ape: v5.8.1
#### phangorn: v2.12.1
#### ete3: v3.1.3
#### pandas: v3.0.2
#### Docker (windows): v29.0.1, build eedd969
#### WSL (windows): v2.6.1.0
#### TreeViewer (windows): 2.2.0

---

### Input data

The species and their respective BioProject tags are as follows:

caenorhabditis_angaria: PRJNA51225

caenorhabditis_auriculariae: PRJEB40642

caenorhabditis_becei: PRJEB28243

caenorhabditis_bovis: PRJEB34497

caenorhabditis_brenneri: PRJNA20035

caenorhabditis_briggsae: PRJNA10731

caenorhabditis_elegans: PRJNA13758

caenorhabditis_inopinata: PRJDB5687

caenorhabditis_japonica: PRJNA12591

caenorhabditis_latens: PRJNA248912

caenorhabditis_nigoni: PRJNA384657

caenorhabditis_panamensis: PRJEB28259

caenorhabditis_parvicauda: PRJEB12595

caenorhabditis_quiockensis: PRJEB11354

caenorhabditis_remanei: PRJNA577507

caenorhabditis_sinica: PRJNA194557

caenorhabditis_sulstoni: PRJEB12601

caenorhabditis_tribulationis: PRJEB12608

caenorhabditis_tropicalis: PRJNA53597

caenorhabditis_uteleia: PRJEB12600

caenorhabditis_waitukubuli: PRJEB12602

caenorhabditis_zanzibari: PRJEB12596

heterorhabditis_bacteriophora: PRJNA13977

mesorhabditis_belari: PRJEB61636

mesorhabditis_spiculigera: PRJEB59059

rhabditophanes_kr3021: PRJEB1297

---

Place all downloaded WormBase ParaSite `.fa.gz` proteome files into:

```text
data/raw/
```

---

### Step 1 - extract proteomes

Run from project root:

```bash
mkdir -p data/extracted

for f in data/raw/*.gz; do
    out=$(basename "$f" .gz)
    gunzip -c "$f" > "data/extracted/$out"
done
```

---

### Step 2 - basic FASTA quality control

Count number of sequences per proteome:

```bash
for f in data/extracted/*.fa; do
    echo "$(basename "$f") : $(grep -c '^>' "$f") sequences"
done
```

Check for blank lines:

```bash
for f in data/extracted/*.fa; do
    if grep -q '^[[:space:]]*$' "$f"; then
        echo "Blank lines found in $f"
    fi
done
```

---

### Step 3 - infer orthologs with OrthoFinder

```bash
docker run --rm \
  -v "$PWD/data/extracted:/input" \
  -v "$PWD/orthofinder_out:/output" \
  davidemms/orthofinder \
  orthofinder -f /input -o /output -t 8
```

Single-copy orthologs will be produced in:

```text
orthofinder_out/Results_*/Single_Copy_Orthologue_Sequences/
```

---

### Step 4 - align single-copy orthologs with MAFFT

Adjust the Results_* in the following command to match the actual one you're intending to run. For me, this was Apr_27
```bash
mkdir -p aligned

for f in orthofinder_out/Results_*/Single_Copy_Orthologue_Sequences/*.fa; do 
    base=$(basename "$f" .fa)
    mafft --auto "$f" > "aligned/${base}_aligned.fa"
done
```

This creates:

```text
aligned_renamed/
```

Which should contain 20 aligned orthologues.

---

### Step 5 - rename ortholog headers

Again, speciesmap.py 's ORTHOGROUPS_FILE path to your orthofinder output's specific Results_* entry.

Run the following:

```bash
python3 scripts/speciesmap.py
```

This creates:

```text
aligned_renamed/
```

---

### Step 6 - concatenate ortholog alignments

Run:

```bash
python3 scripts/concatenate.py
```

Output:

```text
results/concatenated.fa
```

---

### Step 7 - maximum likelihood tree with IQ-TREE

```bash
cd results

iqtree2 -s concatenated.fa -m MFP -bb 1000 -safe
```

Best-fit model selected by BIC:

```text
LG+R4
```

---

### Step 8 - rooting the ML IQTree2 tree

Run:

```bash
python3 scripts/rerooting_tree.py
```

Outgroup:

```text
rhabditophanes_kr3021
```

Output:

```text
results/rooted.treefile
```

---

### Step 9 - neighbor-joining tree

# I struggled to working-directory-agnosticize the R script, so I simply wrote it in an adjustable .rmd file, found at toddfneighborjoining.rmd in Scripts

---

### Step 10 - compare trees

As seen in the above .rmd file, the Robinson–Foulds distance between the rooted Neighbor-Joining and IQ-TREE topologies was:

```text
RF distance = 10
```

### Step 11 - TreeViewer and Visual Comparisons

I then loaded my trees into TreeViewer and made images of them in my .rmd file.
Neighbor-Joining Tree (R):
<img width="585" height="325" alt="image" src="https://github.com/user-attachments/assets/34e12cde-6ab5-4c73-ae52-4b3f7dee7b84" />

Neighbor-Joining Tree (TreeViewer):
<img width="1713" height="971" alt="image" src="https://github.com/user-attachments/assets/57ffe562-97cc-401d-900e-b5683230b42d" />

IQTree2 ML Tree (R):
<img width="567" height="356" alt="image" src="https://github.com/user-attachments/assets/d7575c3e-f48c-4f80-93d4-2ba52e181822" />

IQTree2 ML Tree (TreeViewer):
<img width="1237" height="690" alt="image" src="https://github.com/user-attachments/assets/d9ea6a7d-d925-4930-a26b-1b206c49910b" />


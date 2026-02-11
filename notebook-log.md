# 2/10/2026
# Dataset Description

I am conducting a phylogenetic analysis of nematodes within the genus *Caenorhabditis* and several closely related taxa, including *Heterorhabditis*, *Mesorhabditis*, and *Rhabditophanes*. The dataset includes about 25 species, each represented by a predicted proteome downloaded from WormBase ParaSite. There are a few duplicate entries for species, due to multiple different papers/journals contributing to wormbase. These fasta  files originate from published BioProjects and contain annotated protein coding sequences for each genome.

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

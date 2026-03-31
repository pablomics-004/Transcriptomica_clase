from datetime import datetime
import subprocess as sb
import pandas as pd
import os
import re

# Ubicando la ejecución del script
if (pwd := os.getcwd()).split("/")[-1] == "src":
    os.chdir("..")
    print(f"Working directory changed to: {pwd}", flush=True)
else:
    print(f"Current working directory: {pwd}", flush=True)

# Descarga de datos SRA
file = "./data/GSE132040_MACA_Bulk_metadata.csv"
if not os.path.exists(file):
    print("Downloading metadata from GEO...", flush=True)
    sb.run(
        [ # Parámetros de ejecución de wget
            "wget", 
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132040/suppl/GSE132040_MACA_Bulk_metadata.csv", 
            "-O", 
            file
        ],
        check=True
    )

# Carga de metadatos
metadata = pd.read_csv(file)
del file

# Filtrado
patt_tissue = re.compile(r"^kidney(.?)+", re.IGNORECASE)
metadata = metadata[
    ((metadata["characteristics: age"] == "3") | (metadata["characteristics: age"] == "18")) &
    ((metadata["characteristics: sex"] == "m") & (metadata["source name"].astype(str).str.match(patt_tissue)))
]
SRR_files = list(metadata["raw file"]) # Extracción de los nombres de los archivos SRR
del pwd, patt_tissue, metadata

print(f"SRR files to download: {SRR_files}", flush=True)
print("-" * 50, flush=True)

sra_dir = "./data/sra/"
fastq_dir = "./data/fastq/"

os.makedirs(sra_dir, exist_ok=True)
os.makedirs(fastq_dir, exist_ok=True)
os.makedirs("./tmp/", exist_ok=True)

for SRR in SRR_files:

    sra_file = os.path.join(sra_dir, SRR, f"{SRR}.sra")
    fastq_1 = os.path.join(fastq_dir, f"{SRR}_1.fastq")
    fastq_2 = os.path.join(fastq_dir, f"{SRR}_2.fastq")
    fastq_single = os.path.join(fastq_dir, f"{SRR}.fastq")

    if not os.path.exists(sra_file):
        print(f"[{datetime.now().strftime('%d-%H:%M:%S')}] Downloading {SRR}...", flush=True)
        sb.run(["prefetch", SRR, "--output-directory", sra_dir], check=True)
    
    print(f"[{datetime.now().strftime('%d-%H:%M:%S')}] Processing {SRR}...", flush=True)

    if not (
        os.path.exists(fastq_single) or
        (os.path.exists(fastq_1) and os.path.exists(fastq_2))
    ):
        sb.run([
            "fasterq-dump", sra_file,
            "--split-files",
            # "--skip-technical",
            "--threads", "5",
            "--outdir", fastq_dir,
            "-t", "./tmp/"
        ], check=True)
        print(f"[{datetime.now().strftime('%d-%H:%M:%S')}] Finished processing {SRR}.", flush=True)
        print("-" * 50, flush=True)

print(f"[{datetime.now().strftime('%d-%H:%M:%S')}] All SRR files have been processed.", flush=True)
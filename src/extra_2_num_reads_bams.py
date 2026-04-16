#!/usr/bin/env python3

from datetime import datetime
import subprocess as sb
import pandas as pd
import numpy as np
import os

start_time = datetime.now()

def elapsed():
    return str(datetime.now() - start_time).split(".")[0]

# Ubicando la ejecución del script
if (pwd := os.getcwd()).split("/")[-1] == "src":
    os.chdir("..")
    print(f"[{elapsed()}] El directorio de trabajo cambió a: {os.getcwd()}", flush=True)
else:
    print(f"[{elapsed()}] Directorio de trabajo actual: {pwd}", flush=True)

def get_dirs_list(dir: str = "results") -> list[str]:
    if not os.path.isdir(dir):
        msg = f"[FILES] The given path don't correspond to a directory: {dir}"
        raise ValueError(msg)
    dirs = [path_dir for f in sorted(os.listdir(dir)) if os.path.isdir(path_dir := os.path.join(dir, f))]
    print(f"[{elapsed()}] Directorios encontrados en {dir}: {len(dirs)}", flush=True)
    return dirs

def get_files_list(dir: str = "results") -> list[str]:
    if not os.path.isdir(dir):
        msg = f"[FILES] The given path don't correspond to a directory: {dir}"
        raise ValueError(msg)
    files = [file for f in sorted(os.listdir(dir)) if os.path.isfile(file := (os.path.join(dir, f)))]
    print(f"[{elapsed()}] Archivos en {dir}: {len(files)}", flush=True)
    return files

def capture_num_reads(files_list, size):
    result = np.zeros(size, dtype=np.int64)

    print(f"[{elapsed()}] Contando reads en {len(files_list)} archivos...", flush=True)

    for i, file in enumerate(files_list):
        print(f"[{elapsed()}]   ({i+1}/{len(files_list)}) Procesando: {os.path.basename(file)}", flush=True)

        args = ["samtools", "view", "-c", "-f", "64", file] if "pe" in file.lower() else ["samtools", "view", "-c", file]
        
        output = sb.run(args, capture_output=True, text=True, check=True)
        result[i] = int(output.stdout.strip())

        print(f"[{elapsed()}]   → Reads: {result[i]}", flush=True)

    print(f"[{elapsed()}] Conteo terminado", flush=True)
    return result

def main():

    print(f"[{elapsed()}] Inicio del script", flush=True)

    # Alineadores con archivos BAM
    bam_types = ("BAM", "BAM_FILTERED_SORTED")
    aligners_dirs = [dir for dir in get_dirs_list("./results/") if "SALMON" not in dir]

    print(f"[{elapsed()}] Alineadores detectados: {len(aligners_dirs)}", flush=True)

    aligners_paths = {
        os.path.basename(aligner) : {
            bam_types[0] : get_files_list(os.path.join(aligner, bam_types[0])),
           bam_types[1] : get_files_list(os.path.join(aligner, bam_types[1]))
        }
        for aligner in aligners_dirs
    }

    # Dimensiones de la matriz
    l = len(aligners_paths[os.path.basename(aligners_dirs[0])][bam_types[1]])
    n = len(aligners_dirs)*l
    m = 2

    print(f"[{elapsed()}] Dimensiones matriz: ({n}, {m}) | Bloque: {l}", flush=True)

    M = np.zeros((n, m), dtype=np.int64)

    # Labels para la matriz
    columns = ["raw", "processed"]
    index = []

    # Captura y generación de valores de lecturas en cada BAM
    for i, aligner in enumerate(aligners_paths.keys()):

        print(f"\n[{elapsed()}] Procesando alineador {i+1}/{len(aligners_paths)}: {aligner}", flush=True)

        # Archivos de interés
        bam = aligners_paths[aligner][bam_types[0]]
        bam_filtered = aligners_paths[aligner][bam_types[1]]

        assert l == len(bam), f"El número de archivos BAM crudos del alineador {aligner} no coincide con el esperado ({l})"
        assert l == len(bam_filtered), f"El número de archivos BAM procesados del alineador {aligner} no coincide con el esperado ({l})"
        
        # Slicing
        start = i*l
        end = start+l

        print(f"[{elapsed()}] Filas asignadas: {start}:{end}", flush=True)

        M[start:end, 0] = capture_num_reads(bam, l)
        M[start:end, 1] = capture_num_reads(bam_filtered, l)

        # Obtención de nombres
        index.extend(map(lambda x: os.path.splitext(os.path.basename(x))[0], bam))

        print(f"[{elapsed()}] Alineador {aligner} terminado", flush=True)

    # Generación y guardado del CSV
    print(f"\n[{elapsed()}] Guardando resultados en CSV...", flush=True)
    pd.DataFrame(M, index=index, columns=columns).to_csv("./results/bam_reads.csv", sep=",")

    print(f"[{elapsed()}] Proceso finalizado", flush=True)

if __name__ == "__main__":
    main()
#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os

def get_files_list(dir: str = "results") -> list[str]:
    """
    Devuelve una lista con las rutas completas de los archivos contenidos en el directorio dado, ordenados alfabéticamente.
    """
    if not os.path.isdir(dir):
        msg = f"[FILES] The given path don't correspond to a directory: {dir}"
        raise ValueError(msg)
    return [os.path.join(dir, f) for f in sorted(os.listdir(dir))]

def get_srr_basename(str):
    """
    Devuelve el nombre base de un archivo sin la extensión ni el sufijo de alineador.
    """
    return os.path.basename(str).split("e")[0] + "e"

def standarize_time(time_str: str) -> float:
    """
    Convierte un tiempo en formato "h:m:s" a segundos.
    """
    h, m, s = map(float, time_str.strip().split(":"))
    return h * 3600 + m * 60 + s

def extract_hisat2_info(summary_file) -> np.ndarray:
    """
    Extrae el porcentaje de lecturas alineadas y el tiempo de ejecución a partir del archivo de resumen de HISAT2.
    """
    appeared_uniq, appeared_more = 0, 0
    info = np.empty(4, dtype=np.float64) # [aligned_uniq, >1, overall_rate, time_sec]

    with open(summary_file, "r") as f:
        for line in f.readlines():
            line = line.strip()
            # Reporte de valores
            if not appeared_uniq and "exactly 1 time" in line:
                info[0] = float(line.split()[1].replace("%", "").replace("(", "").replace(")", ""))
                appeared_uniq += 1

            if not appeared_more and ">1 times" in line:
                info[1] = float(line.split()[1].replace("%", "").replace("(", "").replace(")", ""))
                appeared_more += 1
            
            if "overall alignment rate" in line:
                info[2] = float(line.split()[0].replace("%", ""))
            
            if "Time searching" in line:
                info[3] = standarize_time(line.split()[-1])

    return info

def extract_star_info(log_file) -> np.ndarray:
    """
    Extrae el porcentaje de lecturas alineadas y el tiempo de ejecución a partir del archivo de log de STAR.
    """
    info = np.empty(4, dtype=np.float64) # [aligned_uniq, >1, overall, time_sec]

    with open(log_file, "r") as f:
        for line in f.readlines():
            line = line.strip()
            if "Started mapping on" in line:
                t0 = standarize_time(line.split()[-1])
            if "Finished on" in line:
                t1 = standarize_time(line.split()[-1])
                info[3] = t1 - t0

            if "Uniquely mapped reads %" in line:
                info[0] = float(line.split()[-1].replace("%", ""))

            if "% of reads mapped to multiple loci" in line:
                info[1] = float(line.split()[-1].replace("%", ""))
                info[2] = info[0] + info[1]

    return info

# Código principal
def main():
    # Obtener la lista de archivos de resumen de HISAT2 y log de STAR
    alignment_dir = {"hisat2": "time_align", "star": "Logs"}
    aligners = list(alignment_dir.keys())
    files = sum(
        (
            get_files_list(
                os.path.join("results", (aligner + f"_{sufix}"), alignment_dir[aligner])
            )
            for aligner in aligners
            for sufix in ("pe", "se")
        ),
        []
    )

    # Extracción de información para cada archivo y almacenarla en un DataFrame
    results = {
        get_srr_basename(f): extract_hisat2_info(f) if "hisat2" in f else extract_star_info(f)
        for f in files
    }

    # Guardado de los resultados en un archivo CSV
    pd.DataFrame(
        results,
        index=["aligned_uniq", "aligned_mult", "overall_rate", "time_sec"]
    ).T.to_csv("./results/aligners_performance_table.csv", sep=",")

if __name__ == "__main__":
    main()
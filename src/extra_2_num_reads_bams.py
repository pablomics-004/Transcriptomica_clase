#!/usr/bin/env python3

import subprocess as sb
import pandas as pd
import numpy as np
import os

# Ubicando la ejecución del script
if (pwd := os.getcwd()).split("/")[-1] == "src":
    os.chdir("..")
    print(f"El directorio de trabajo cambió a: {pwd}", flush=True)
else:
    print(f"Directorio de trabajo actual: {pwd}", flush=True)

def capture_num_reads(files_list):
    return np.array([
        int(sb.run([
            "samtools view", "-c", "-f", 64, file # PE
        ])) if "pe" in file else int(sb.run([
            "samtools view", "-c", file #SE
        ]))
        for file in files_list
    ], dtype=np.int32)

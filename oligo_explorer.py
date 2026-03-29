"""
CIENBIO Oligo Explorer v2.0
Diseño de oligos degenerados para PCR a partir de secuencias de proteínas.

Correcciones respecto a v1.0:
  - back_table reconstruida desde forward_table (todos los codones por AA)
  - Parámetros de iniciación SantaLucia 1998 completos (AT y GC)
  - Fórmula Tm corregida: usa CT/4 para oligos no autocomplementarios
  - MUSCLE resuelto dinámicamente (no ruta hardcodeada)
  - Límite de explosión combinatoria configurable
  - Archivo temporal eliminado correctamente
  - Parámetros configurables interactivamente
  - Estructura modular con main()

Autor: Augusto Manubens — CIENBIO SpA
Referencia Tm: SantaLucia (1998) PNAS 95:1460-1465
"""

import sys
import os
import math
import shutil
import tempfile
import itertools
import subprocess
from collections import defaultdict
from datetime import datetime

from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
import pandas as pd


# ─────────────────────────────────────────────
#  PARÁMETROS NEAREST-NEIGHBOR (SantaLucia 1998, Tabla 2)
#  ΔH en kcal/mol, ΔS en cal/mol·K
# ─────────────────────────────────────────────
NN_VALUES = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2), "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}

# ─────────────────────────────────────────────
#  BACK TABLE COMPLETA (todos los codones por AA)
#  FIX #1: back_table de BioPython devuelve un solo str por AA.
#  Se reconstruye desde forward_table para incluir todos los codones sinónimos.
# ─────────────────────────────────────────────
def construir_back_table():
    tabla = CodonTable.unambiguous_dna_by_name["Standard"]
    back = defaultdict(list)
    for codon, aa in tabla.forward_table.items():
        back[aa].append(codon)
    return dict(back)

BACK_TABLE = construir_back_table()


# ─────────────────────────────────────────────
#  MUSCLE
# ─────────────────────────────────────────────
def get_muscle_path():
    """Resuelve la ruta de MUSCLE dinámicamente. FIX #4."""
    # 1. En PATH del sistema
    path = shutil.which("muscle")
    if path:
        return path
    # 2. Junto al script
    local = os.path.join(os.path.dirname(os.path.abspath(__file__)), "muscle")
    if os.path.exists(local):
        return local
    local_exe = local + ".exe"
    if os.path.exists(local_exe):
        return local_exe
    raise FileNotFoundError(
        "MUSCLE no encontrado.\n"
        "  · Instálalo (https://drive5.com/muscle5/) y agrégalo al PATH, o\n"
        "  · Coloca el ejecutable junto a este script."
    )


def alinear_muscle(fasta_path, muscle_path):
    """Ejecuta MUSCLE y devuelve la ruta del alineamiento temporal. FIX #6."""
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".aln")
    tmp.close()
    subprocess.run(
        [muscle_path, "-align", fasta_path, "-output", tmp.name],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return tmp.name


# ─────────────────────────────────────────────
#  ANÁLISIS DE ALINEAMIENTO
# ─────────────────────────────────────────────
def calcular_identidad_global(alignment):
    total_id = 0
    total_positions = 0
    nseqs = len(alignment)
    for i in range(alignment.get_alignment_length()):
        columna = [record.seq[i] for record in alignment]
        if "-" in columna:
            continue
        max_id = max(columna.count(res) for res in set(columna))
        total_id += max_id / nseqs
        total_positions += 1
    return (total_id / total_positions) * 100 if total_positions else 0.0


def calcular_conservacion_columnas(alignment):
    n = len(alignment)
    conservacion = []
    for i in range(alignment.get_alignment_length()):
        columna = [record.seq[i] for record in alignment]
        if "-" in columna:
            conservacion.append(0.0)
        else:
            max_id = max(columna.count(aa) for aa in set(columna))
            conservacion.append(max_id / n)
    return conservacion


def calcular_consenso_con_criterio(alignment, umbral):
    consenso = ""
    for i in range(alignment.get_alignment_length()):
        columna = [record.seq[i] for record in alignment]
        aa_frecuente = max(set(columna), key=columna.count)
        if columna.count(aa_frecuente) / len(columna) >= umbral:
            consenso += aa_frecuente
        else:
            consenso += "-"
    return consenso


# ─────────────────────────────────────────────
#  DISEÑO DE OLIGOS
# ─────────────────────────────────────────────
def codones_por_aminoacido(aa_block):
    """Devuelve TODOS los codones sinónimos para cada AA. FIX #1."""
    return [BACK_TABLE.get(aa, ["NNN"]) for aa in aa_block]


def calcular_tm_dg_santalucia(seq, Na=0.05, CT=2.5e-7):
    """
    Calcula Tm y ΔH por el método nearest-neighbor (SantaLucia 1998).

    FIX #2: Parámetros de iniciación completos (AT y GC).
    FIX #3: Fórmula Tm usa CT/4 para oligos no autocomplementarios.

    Args:
        seq: secuencia de ADN (5'→3')
        Na:  concentración de Na+ en M (default: 50 mM)
        CT:  concentración total del oligo en M (default: 250 nM)

    Returns:
        (Tm en °C, ΔH en kcal/mol)
    """
    seq = seq.upper()
    total_dh = 0.0
    total_ds = 0.0

    # Suma de contribuciones nearest-neighbor
    for i in range(len(seq) - 1):
        pair = seq[i:i+2]
        if pair in NN_VALUES:
            dh, ds = NN_VALUES[pair]
            total_dh += dh
            total_ds += ds

    # Correcciones de iniciación (SantaLucia 1998, Tabla 2)
    # FIX #2: incluye terminales AT (antes ausentes)
    for terminal in (seq[0], seq[-1]):
        if terminal in "GC":
            total_dh += 0.1
            total_ds += -2.8
        elif terminal in "AT":
            total_dh += 2.3
            total_ds += 4.1

    # Tm para oligos no autocomplementarios: ΔH / (ΔS + R·ln(CT/4))
    # FIX #3: CT/4 en lugar de CT directamente
    R = 1.987  # cal/mol·K
    tm_kelvin = (1000 * total_dh) / (total_ds + R * math.log(CT / 4))

    # Corrección de sal (Marmur-Schildkraut-Doty simplificada)
    salt_correction = 16.6 * math.log10(Na)

    return round(tm_kelvin - 273.15 + salt_correction, 2), round(total_dh, 2)


# ─────────────────────────────────────────────
#  UTILIDADES DE INPUT
# ─────────────────────────────────────────────
def pedir_float(prompt, default, minval=None, maxval=None):
    while True:
        raw = input(f"{prompt} [{default}]: ").strip()
        val = float(raw) if raw else default
        if minval is not None and val < minval:
            print(f"  ⚠️  Debe ser ≥ {minval}")
            continue
        if maxval is not None and val > maxval:
            print(f"  ⚠️  Debe ser ≤ {maxval}")
            continue
        return val


def pedir_int(prompt, default, minval=None, maxval=None):
    while True:
        raw = input(f"{prompt} [{default}]: ").strip()
        val = int(raw) if raw else default
        if minval is not None and val < minval:
            print(f"  ⚠️  Debe ser ≥ {minval}")
            continue
        if maxval is not None and val > maxval:
            print(f"  ⚠️  Debe ser ≤ {maxval}")
            continue
        return val


# ─────────────────────────────────────────────
#  MAIN
# ─────────────────────────────────────────────
def main():
    print("=" * 60)
    print("  CIENBIO Oligo Explorer v2.0")
    print("  Diseño de oligos degenerados desde proteínas")
    print("=" * 60)
    print()

    # ── Inputs ──────────────────────────────
    fasta_input = input("📂 Ruta archivo FASTA: ").strip().strip('"')
    if not os.path.exists(fasta_input):
        print(f"❌ Archivo no encontrado: '{fasta_input}'")
        sys.exit(1)

    bloque_longitud = pedir_int("🔢 Longitud del bloque en AA (ej. 5)", default=5, minval=2, maxval=20)

    print()
    print("── Parámetros de conservación ──────────────────────────")
    umbral_identidad       = pedir_float("🎯 Umbral identidad global % (ej. 50)", default=50, minval=0, maxval=100)
    umbral_conservacion    = pedir_float("📊 Umbral conservación local 0–1 (ej. 0.85)", default=0.85, minval=0, maxval=1)

    print()
    print("── Parámetros termodinámicos ────────────────────────────")
    Na_mM  = pedir_float("🧪 [Na+] en mM (ej. 50)", default=50, minval=1, maxval=1000)
    CT_nM  = pedir_float("🧬 Concentración oligo CT en nM (ej. 250)", default=250, minval=1)
    Na_M   = Na_mM / 1000
    CT_M   = CT_nM / 1e9

    print()
    print("── Límites de seguridad ─────────────────────────────────")
    max_combinaciones = pedir_int(
        "⚠️  Máx. combinaciones por bloque (ej. 1000)", default=1000, minval=1
    )
    max_oligos_total = pedir_int(
        "📋 Máx. oligos en output total (ej. 10000)", default=10000, minval=1
    )

    # ── Output dir ──────────────────────────
    output_dir = os.path.dirname(os.path.abspath(fasta_input))

    # ── MUSCLE ──────────────────────────────
    print()
    print("🔧 Buscando MUSCLE...")
    try:
        muscle_path = get_muscle_path()
        print(f"   ✅ Encontrado: {muscle_path}")
    except FileNotFoundError as e:
        print(f"❌ {e}")
        sys.exit(1)

    print("🧬 Ejecutando alineamiento con MUSCLE...")
    aligned_file = alinear_muscle(fasta_input, muscle_path)

    try:
        alignment = AlignIO.read(aligned_file, "fasta")
    finally:
        os.unlink(aligned_file)   # FIX #6: limpieza garantizada

    print(f"   ✅ {len(alignment)} secuencias × {alignment.get_alignment_length()} posiciones\n")
    for record in alignment:
        print(f"   >{record.id}\n   {record.seq}\n")

    # ── Identidad global ─────────────────────
    identidad = calcular_identidad_global(alignment)
    print(f"🔬 Identidad global promedio: {identidad:.2f}%")
    if identidad < umbral_identidad:
        print(f"⚠️  Identidad insuficiente (< {umbral_identidad}%). Abortando.")
        sys.exit(0)

    # ── Conservación por columna ─────────────
    print("📊 Calculando conservación por columna...")
    conservacion = calcular_conservacion_columnas(alignment)
    conservacion_path = os.path.join(output_dir, "conservacion_columnas.xlsx")
    pd.DataFrame({
        "Posición": list(range(1, len(conservacion) + 1)),
        "Conservación (%)": [round(c * 100, 2) for c in conservacion],
    }).to_excel(conservacion_path, index=False)
    print(f"   ✅ Guardado: {conservacion_path}")

    # ── Consenso y diseño de oligos ──────────
    print("🔍 Generando oligos degenerados por bloque conservado...")
    consenso = calcular_consenso_con_criterio(alignment, umbral_conservacion)

    oligos = []
    bloques_omitidos = 0
    limite_alcanzado = False

    for i in range(len(consenso) - bloque_longitud + 1):
        if limite_alcanzado:
            break
        bloque = consenso[i:i + bloque_longitud]
        if "-" in bloque:
            continue

        codon_list = codones_por_aminoacido(bloque)

        # FIX #5: control de explosión combinatoria
        n_comb = math.prod(len(c) for c in codon_list)
        if n_comb > max_combinaciones:
            print(f"   ⚠️  Pos {i+1}-{i+bloque_longitud} ({bloque}): {n_comb} combinaciones → omitido")
            bloques_omitidos += 1
            continue

        for comb in itertools.product(*codon_list):
            if len(oligos) >= max_oligos_total:
                print(f"   ⚠️  Límite de {max_oligos_total} oligos alcanzado. Output truncado.")
                limite_alcanzado = True
                break
            oligo_fwd = "".join(comb)
            oligo_rev = str(Seq(oligo_fwd).reverse_complement())
            tm, dh = calcular_tm_dg_santalucia(oligo_fwd, Na=Na_M, CT=CT_M)
            oligos.append({
                "Posición (AA)": f"{i+1}-{i+bloque_longitud}",
                "Aminoácidos":   bloque,
                "Oligo Forward": oligo_fwd,
                "Oligo Reverso": oligo_rev,
                "Longitud (nt)": len(oligo_fwd),
                "Tm (°C)":       tm,
                "ΔH (kcal/mol)": dh,
                "N combinaciones bloque": n_comb,
            })

    # ── Output ──────────────────────────────
    print()
    if oligos:
        filename = f"oligos_degenerados_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx"
        out_path = os.path.join(output_dir, filename)
        df = pd.DataFrame(oligos)
        df.to_excel(out_path, index=False)
        print(f"✅ {len(oligos)} oligos generados → {out_path}")
        if bloques_omitidos:
            print(f"   ℹ️  {bloques_omitidos} bloques omitidos por alta degeneración (>{max_combinaciones} comb.)")
    else:
        print("⚠️  No se generaron oligos.")
        print("   Revisa los umbrales de conservación o el bloque_longitud.")


if __name__ == "__main__":
    main()
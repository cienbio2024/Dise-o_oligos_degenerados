# CIENBIO Oligo Explorer v2.0

Script para el diseño de oligos degenerados para PCR a partir de secuencias de proteínas múltiples alineadas.

Desarrollado por **Augusto Manubens — [CIENBIO SpA](https://cienbio.cl)**

---

## ¿Qué hace?

A partir de un archivo FASTA con dos o más secuencias de aminoácidos:

1. Alinea las secuencias (MUSCLE o BioPython automáticamente)
2. Calcula la conservación por columna
3. Identifica bloques conservados según umbrales configurables
4. Genera **todos los oligos degenerados posibles** (forward y reverso) para cada bloque
5. Calcula **Tm y ΔH** por el método nearest-neighbor (SantaLucia 1998)
6. Exporta resultados a Excel

---

## Requisitos

```bash
pip install biopython pandas openpyxl
```

**MUSCLE (opcional pero recomendado para datasets grandes)**
Descarga el ejecutable desde https://github.com/rcedgar/muscle/releases/latest
y colócalo junto al script como `muscle.exe` (Windows) o `muscle` (Linux/Mac).
Si MUSCLE no está disponible, el script usa automáticamente el alineador interno de BioPython.

---

## Uso

```bash
python oligo_explorer.py
```

El script es interactivo y solicita todos los parámetros:

| Parámetro | Descripción | Default |
|---|---|---|
| Ruta FASTA | Archivo con ≥2 secuencias de aminoácidos | — |
| Longitud bloque | Número de AA por bloque de oligo | 5 |
| Umbral identidad global | % mínimo de identidad para continuar | 50% |
| Umbral conservación local | Fracción mínima de conservación por columna | 0.85 |
| [Na⁺] | Concentración de sal en mM | 50 mM |
| CT oligo | Concentración total del oligo en nM | 250 nM |
| Máx. combinaciones/bloque | Límite para evitar explosión combinatoria | 1000 |
| Máx. oligos total | Límite total del archivo de salida | 10000 |

---

## Archivos de salida

| Archivo | Contenido |
|---|---|
| `conservacion_columnas.xlsx` | Conservación (%) por posición del alineamiento |
| `oligos_degenerados_YYYYMMDD_HHMMSS.xlsx` | Todos los oligos: forward, reverso, Tm, ΔH |

---

## Fundamento científico

### Diseño de oligos degenerados
Los oligos se diseñan desde la secuencia consenso de aminoácidos, expandiendo cada AA a **todos sus codones sinónimos** según el código genético estándar (tabla de BioPython `CodonTable.unambiguous_dna_by_name["Standard"]`). Las combinaciones se generan por producto cartesiano.

### Cálculo de Tm (SantaLucia 1998)
```
Tm = ΔH / (ΔS + R·ln(CT/4)) − 273.15 + 16.6·log₁₀([Na⁺])
```
- Parámetros nearest-neighbor de SantaLucia (1998) PNAS 95:1460–1465
- Correcciones de iniciación para terminales GC y AT
- Fórmula para oligos **no autocomplementarios** (CT/4)
- Corrección de sal por Marmur-Schildkraut-Doty

---

## Cambios en v2.0 respecto a v1.0

| # | Problema en v1.0 | Solución en v2.0 |
|---|---|---|
| 🔴 | `back_table` devolvía 1 codón por AA → sin degeneración real | Reconstruida desde `forward_table` (todos los codones sinónimos) |
| 🔴 | Fórmula Tm usaba CT directo en lugar de CT/4 | Corregido según SantaLucia 1998 |
| 🔴 | Faltaban parámetros de iniciación para terminales AT | Agregados (+2.3 kcal/mol, +4.1 cal/mol·K) |
| 🟠 | Ruta de MUSCLE hardcodeada (solo funcionaba en 1 PC) | Resolución dinámica: PATH → carpeta del script → fallback BioPython |
| 🟠 | Sin límite de combinaciones (posible cuelgue con bloques degenerados) | Límite configurable por bloque y total |
| 🟠 | Archivo `.aln` temporal nunca eliminado | Eliminación garantizada con `finally` |
| 🟡 | Parámetros hardcodeados en el código | Todos los parámetros son inputs interactivos con defaults |
| 🟡 | Código sin estructura modular | `main()` + `if __name__ == "__main__":` |

---

## Licencia

MIT © 2026 Augusto Manubens — CIENBIO SpA

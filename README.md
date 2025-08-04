# 🔬 ANNOMAF vs Funcotator Comparison: Test Framework

Este repositorio contiene un framework de comparación entre herramientas de anotación genómica, específicamente **ANNOMAF**, **Funcotator** y **Oncotator**, aplicado a muestras WES, WGS y paneles personalizados.

---

## 📁 Estructura General

```
.
├── data-test          # Muestras VCF para test
├── programs           # Pipelines y programas utilizados
├── test_WES           # Comparación con muestras WES
├── test_WGS           # Comparación con muestras WGS
└── test_panels        # Comparación con paneles clínicos
```

---

## 1️⃣ Descarga de Muestras

Desde el directorio `data-test`, ejecuta el script de descarga:

```bash
bash input_dwnld.sh
```

Esto generará carpetas con archivos de prueba como:

```
data-test/
├── Panels/
│   └── CC004BTUCH.mutect2.filtered_VEP.ann.vcf
├── WES/
│   ├── WES_IL_1.bwa.strelka.vcf.gz
│   └── WES_NV_1.bwa.muTect2.vcf.gz
├── WGS/
│   ├── WGS_FD_1.bwa.muTect2.vcf
│   └── WGS_NV_2.bwa.strelka.vcf.gz
```

> Descomprime los archivos `.vcf.gz` con el siguiente comando (conserva los originales en caso de):

```bash
gzip -dk *.gz
```

---

## 2️⃣ Preparación del Test

> ⚠️ *Nota: Ants de ejecutar cualquier archivo revisa las rutas que utiliza para configurarlas.

Entra al folder correspondiente (por ejemplo, `test_WES`, `test_WGS`, etc.) y ejecuta el script de configuración:

```bash
chmod +x setup_annomaf_test_casesWES.sh
./setup_annomaf_test_casesWES.sh
```

Esto creará la siguiente estructura de carpetas para realizar las pruebas:

```
test_WES/
├── T1_byCaller/
│   ├── annomaf/
│   │   ├── output/
│   │   ├── params.json
│   │   └── samplesheet.csv
│   ├── funcotator/
│   └── oncotator/
├── T2_multiCaller_IL1/
├── run_annomaf.sh
├── run_funcotator.sh
└── setup_annomaf_test_casesWES.sh
```

---

## 3️⃣ Ejecución

### 🔹 3.1 Ejecutar ANNOMAF

Ejecuta el siguiente script para realizar las pruebas con **ANNOMAF**:

```bash
bash run_annomaf.sh
```

Este script:

- Ejecuta **ANNOMAF** utilizando `nextflow`.
- Mide el uso de **CPU**, **RAM** y el **tiempo de ejecución**.
- Genera tres archivos con los resultados:
  - `annomaf_report.txt` con el reporte de ejecución.
  - `annomaf_stdout.log` con los pasos detallados.
  - `annomaf_time.log` con los tiempos de ejecución.

Ejemplo de salida en `annomaf_report.txt`:

```
Test: T1_byCaller
Execution time (s): 237
VCF file size: 1.9M
Total variants: 7702
Unique genes: 85
Actionable variants: 0
CPU: AMD Ryzen 5 5600 6-Core Processor
Total RAM (MB): 15957
Max RSS: 502868 KB
CPU usage: 9%
```

También se mostrará un resumen en la terminal:

```
🔁 Iniciando ejecución de casos ANNOMAF...
⚙️ Procesando: T1_byCaller
📊 [Resumen - T1_byCaller]
⏱️ Tiempo: 237s
💾 VCF size: 1.9M
📈 Variantes: 7702
🧬 Genes únicos: 85
🎯 Accionables: 0
🧠 RAM: 15957MB / Máximo RSS: 502868 KB
⚙️ CPU: AMD Ryzen 5 5600 6-Core Processor
```

---

### 🔹 3.2 Ejecutar Funcotator

Ejecuta el siguiente script para realizar las pruebas con **Funcotator**:

```bash
bash run_funcotator.sh
```

> ⚠️ *Nota: Actualmente este script puede presentar errores con algunas muestras WES y WGS.* 

Cuando se corrijan o se agreguen más detalles, los mismos serán documentados aquí.

---

## 📌 Pendientes

- [ ] Agregar ejemplos de salida de Funcotator.
- [ ] Documentar comparaciones numéricas cruzadas (drive).

---

## 💬 Créditos

Proyecto coordinado por [@conifva](https://github.com/conifva), con colaboración de [NOMBRE DE TU AMIGA].

---

#!/bin/bash

set -e  # Abortar si falla algo

# Ruta fija al archivo de parámetros
BASE_DIR="/home/conifva/MathBioLab/git_works/OneDrive_1_20-5-2025"
PARAM_FILE="$BASE_DIR/parameters_file.json"

# Casos de test definidos
TARGET_CASES=("T1_byCaller" "T2_multiCaller_IL1" "T3_byInstitution" "T4_replicates_FD" "T5_lightestFile" "T6_heaviestFile")

# Directorio donde se ejecuta
EXEC_DIR="$(pwd)"
echo "🔁 Creando estructura de test en: $EXEC_DIR"

for CASE in "${TARGET_CASES[@]}"; do
  echo "📁 Configurando $CASE..."

  # Crear carpetas si no existen
  [[ -d "$EXEC_DIR/$CASE/annomaf/output" ]] || mkdir -p "$EXEC_DIR/$CASE/annomaf/output"
  [[ -d "$EXEC_DIR/$CASE/funcotator" ]]     || mkdir -p "$EXEC_DIR/$CASE/funcotator"
  [[ -d "$EXEC_DIR/$CASE/oncotator" ]]      || mkdir -p "$EXEC_DIR/$CASE/oncotator"

  # Copiar param.json si no existe
  if [[ ! -f "$EXEC_DIR/$CASE/annomaf/params.json" ]]; then
    cp "$PARAM_FILE" "$EXEC_DIR/$CASE/annomaf/params.json"
  else
    echo "⚠️  params.json ya existe para $CASE, no se sobreescribe."
  fi
done

echo "🧪 Verificando creación de samplesheet.csv..."

create_samplesheet() {
  local path="$1"
  local content="$2"

  if [[ -f "$path" ]]; then
    echo "⚠️  samplesheet.csv ya existe en $(dirname "$path"), no se sobreescribe."
  else
    echo "$content" > "$path"
    echo "✅ Creado: $path"
  fi
}

# T1: Test por caller
create_samplesheet "$EXEC_DIR/T1_byCaller/annomaf/samplesheet.csv" \
"sample,variantcaller,path,lb_bam
FD1,mutect2,../../../data-test/WES/WES_FD_1.bwa.muTect2.vcf,no"

# T2: Comparación de callers
# T2: Test con múltiples callers para IL1, ver si pasa el mismo id pero sin considerar el caller (VERIFICAR)
create_samplesheet "$EXEC_DIR/T2_multiCaller_IL1/annomaf/samplesheet.csv" \
"sample,variantcaller,path,lb_bam
IL1,mutect2,../../../data-test/WES/WES_IL_1.bwa.muTect2.vcf,no
IL1,sniper,../../../data-test/WES/WES_IL_1.bwa.somaticSniper.vcf,no
IL1,strelka,../../../data-test/WES/WES_IL_1.bwa.strelka.vcf,no"

# T3: Comparación entre instituciones
create_samplesheet "$EXEC_DIR/T3_byInstitution/annomaf/samplesheet.csv" \
"sample,variantcaller,path,lb_bam
FD1,mutect2,../../../data-test/WES/WES_FD_1.bwa.muTect2.vcf,no
IL1,mutect2,../../../data-test/WES/WES_IL_1.bwa.muTect2.vcf,no
NV1,mutect2,../../../data-test/WES/WES_NV_1.bwa.muTect2.vcf,no"

# T4: Réplicas de misma institución
create_samplesheet "$EXEC_DIR/T4_replicates_FD/annomaf/samplesheet.csv" \
"sample,variantcaller,path,lb_bam
FD1,mutect2,../../../data-test/WES/WES_FD_1.bwa.muTect2.vcf,no
FD2,mutect2,../../../data-test/WES/WES_FD_2.bwa.muTect2.vcf,no
FD3,mutect2,../../../data-test/WES/WES_FD_3.bwa.muTect2.vcf,no"

# T5: Archivo más liviano
create_samplesheet "$EXEC_DIR/T5_lightestFile/annomaf/samplesheet.csv" \
"sample,variantcaller,path,lb_bam
LIGHT,sniper,../../../data-test/WES/WES_FD_2.bwa.somaticSniper.vcf,no"

# T6: Archivo más pesado
create_samplesheet "$EXEC_DIR/T6_heaviestFile/annomaf/samplesheet.csv" \
"sample,variantcaller,path,lb_bam
HEAVY,strelka,../../../data-test/WES/WES_FD_1.bwa.strelka.vcf,no"

echo ""
echo "✅ ¡Listo! Las estructuras de los tests se han generado sin sobrescribir archivos existentes."

echo ""
echo "📌 Para correr ANNOMAF en uno de los tests manualmente:"
echo "cd T1_byCaller/annomaf && nextflow run ../../../programs/nf-core-annomaf/ \\"
echo "  -params-file params.json --input samplesheet.csv --outdir output \\"
echo "  -profile docker --max_memory '16.GB' --max_cpus 8 --max_time '12.h'"

#!/bin/bash

SCRIPT_DIR="$(dirname "$(realpath "$0")")"
echo "🔧 Paso 1: Preparando archivos FASTA adicionales para Funcotator..."

FASTA_ORIGINAL="/home/conifva/MathBioLab/git_works/OneDrive_1_20-5-2025/references/fasta_hg38/Homo_sapiens_assembly38.fasta"
DEST_DIR="/home/conifva/MathBioLab/publi_ANNOM_comparison/data-test/input_adicional_funcotator"
FASTA_DEST="$DEST_DIR/Homo_sapiens_assembly38.fasta"
DICT_DEST="${FASTA_DEST%.fasta}.dict"
GATK_JAR="/home/conifva/MathBioLab/publi_ANNOM_comparison/programs/funcotator/gatk-4.6.2.0/gatk-package-4.6.2.0-local.jar"
DATASOURCES="/home/conifva/MathBioLab/publi_ANNOM_comparison/programs/funcotator/gatk-4.6.2.0/funcotator_dataSources.v1.8.hg38.20230908s"

mkdir -p "$DEST_DIR"

# Copiar y preparar FASTA si no existe
if [ ! -f "$FASTA_DEST" ]; then
    echo "📄 Copiando FASTA..."
    cp "$FASTA_ORIGINAL" "$FASTA_DEST"
fi

if [ ! -f "$FASTA_DEST.fai" ]; then
    echo "🛠️  Generando .fai con samtools..."
    if ! command -v samtools &> /dev/null; then
        echo "❌ Error: 'samtools' no está instalado. Instálalo con: sudo apt install samtools"
        exit 1
    fi
    samtools faidx "$FASTA_DEST"
fi

if [ ! -f "$DICT_DEST" ]; then
    echo "🛠️  Generando .dict con GATK..."
    java -jar "$GATK_JAR" CreateSequenceDictionary -R "$FASTA_DEST" -O "$DICT_DEST"
fi

echo "✅ FASTA listo en: $FASTA_DEST"
echo ""

ejecutar_funcotator() {
    local vcf=$1
    local sample=$2
    local caller=$3
    local output_dir=$4
    local output_maf="$output_dir/${sample}.${caller}.maf"
    local stdout_log="$output_dir/funcotator_stdout_${sample}_${caller}.log"
    local report_txt="$output_dir/funcotator_report_${sample}_${caller}.txt"
    local time_log="$output_dir/funcotator_time_${sample}_${caller}.log"

    if [ -f "$output_maf" ]; then
        echo "🟡 Ya existe: $output_maf — Omitiendo"
        return
    fi

    echo "🔄 Ejecutando Funcotator para $sample con $caller → $output_maf"
    mkdir -p "$output_dir"

    start_time=$(date +%s)

    # Ejecutar con monitoreo si /usr/bin/time existe
    if [ -x "/usr/bin/time" ]; then
        /usr/bin/time -v -o "$time_log" java -Dsamjdk.use_async_io_read_samtools=false \
             -Dsamjdk.use_async_io_write_samtools=true \
             -Dsamjdk.use_async_io_write_tribble=false \
             -Dsamjdk.compression_level=2 \
             -jar "$GATK_JAR" Funcotator \
             --variant "$vcf" \
             --reference "$FASTA_DEST" \
             --ref-version hg38 \
             --data-sources-path "$DATASOURCES" \
             --output "$output_maf" \
             --output-file-format MAF \
             &> "$stdout_log"
    else
        echo "⚠️  /usr/bin/time no disponible, ejecutando sin monitoreo detallado"
        java -Dsamjdk.use_async_io_read_samtools=false \
             -Dsamjdk.use_async_io_write_samtools=true \
             -Dsamjdk.use_async_io_write_tribble=false \
             -Dsamjdk.compression_level=2 \
             -jar "$GATK_JAR" Funcotator \
             --variant "$vcf" \
             --reference "$FASTA_DEST" \
             --ref-version hg38 \
             --data-sources-path "$DATASOURCES" \
             --output "$output_maf" \
             --output-file-format MAF \
             &> "$stdout_log"
    fi

    end_time=$(date +%s)
    duration=$((end_time - start_time))

    # Info máquina
    cpu_model=$(grep -m1 "model name" /proc/cpuinfo | cut -d: -f2 | xargs)
    mem_total=$(grep MemTotal /proc/meminfo | awk '{print $2}')
    mem_total_mb=$((mem_total / 1024))

    # Tamaño VCF
    vcf_size=$(du -h "$vcf" | cut -f1)
    total_variants=$(zgrep -v '^#' "$vcf" 2>/dev/null | wc -l || grep -v '^#' "$vcf" | wc -l)

    # Genes y accionables
    if [ -f "$output_maf" ]; then
        unique_genes=$(cut -f13 "$output_maf" | grep -v "Hugo_Symbol" | sort | uniq | wc -l)
        accionables=$(grep -ic "actionable" "$output_maf")
    else
        unique_genes="NA"
        accionables="NA"
    fi

    # Extraer métricas de time_log si existe
    if [ -f "$time_log" ]; then
        max_rss=$(grep "Maximum resident set size" "$time_log" | awk '{print $6 " KB"}')
        cpu_usage=$(grep "Percent of CPU this job got" "$time_log" | awk '{print $NF}')
    else
        max_rss="NA"
        cpu_usage="NA"
    fi

    # Mostrar por terminal
    echo "📊 [Resumen - $sample / $caller]"
    echo "⏱️ Tiempo de ejecución: ${duration}s"
    echo "💾 VCF size: $vcf_size"
    echo "📈 Variantes totales: $total_variants"
    echo "🧬 Genes únicos (MAF): $unique_genes"
    echo "🎯 Variantes accionables: $accionables"
    echo "🖥️ CPU: $cpu_model"
    echo "🧠 RAM total: ${mem_total_mb}MB"
    echo "📊 Máximo RSS: $max_rss"
    echo "⚙️ CPU usage: $cpu_usage"
    echo ""

    # Guardar en archivo
    {
      echo "Sample: $sample"
      echo "Caller: $caller"
      echo "Execution time (s): $duration"
      echo "VCF file size: $vcf_size"
      echo "Total variants: $total_variants"
      echo "Unique genes (MAF): $unique_genes"
      echo "Actionable variants: $accionables"
      echo "CPU: $cpu_model"
      echo "Total RAM (MB): $mem_total_mb"
      echo "Max RSS (KB): $max_rss"
      echo "CPU usage (%): $cpu_usage"
    } > "$report_txt"
}

# Procesar casos
for case_dir in T1_byCaller  T2_multiCaller_IL1  T3_byInstitution  T4_replicates_FD  T5_lightestFile  T6_heaviestFile; do
    sheet="$SCRIPT_DIR/$case_dir/annomaf/samplesheet.csv"
    output_dir="$SCRIPT_DIR/$case_dir/funcotator"

    echo "📁 Procesando $case_dir..."

    if [ ! -f "$sheet" ]; then
        echo "⚠️  No se encontró $sheet — Omitiendo caso"
        continue
    fi

    tail -n +2 "$sheet" | while IFS=',' read -r sample caller vcf_path lb_bam; do
        vcf_dir="$(dirname "$sheet")"
        original_vcf_path="$vcf_dir/$vcf_path"
        filtered_vcf_path="${original_vcf_path%.vcf}.filtered.vcf.gz"

        # Usar el .filtered.vcf.gz si existe, si no usar el original
        if [ -f "$filtered_vcf_path" ]; then
            full_vcf_path="$filtered_vcf_path"
        else
            full_vcf_path="$original_vcf_path"
        fi

        ejecutar_funcotator "$full_vcf_path" "$sample" "$caller" "$output_dir"
    done
    echo ""
done

echo "🎉 Todos los casos procesados."

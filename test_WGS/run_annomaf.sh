#!/bin/bash

SCRIPT_DIR="$(dirname "$(realpath "$0")")"
ANNOMAF_PIPELINE="/home/conifva/MathBioLab/publi_ANNOM_comparison/programs/nf-core-annomaf"

echo "🔁 Iniciando ejecución de casos ANNOMAF..."

# Función para ejecutar ANNOMAF y registrar métricas
ejecutar_annomaf() {
    local annomaf_dir=$1
    local test_dir="$(dirname "$annomaf_dir")"
    local test_name="$(basename "$test_dir")"
    local sheet="samplesheet.csv"
    local output_dir="output"
    local report_txt="annomaf_report.txt"
    local time_log="annomaf_time.log"
    local stdout_log="annomaf_stdout.log"

    echo "⚙️ Procesando: $test_name"

    if [ ! -f "$annomaf_dir/$sheet" ]; then
        echo "⚠️  No se encontró $sheet en $annomaf_dir — Omitiendo caso"
        return
    fi

    (
        cd "$annomaf_dir" || exit

        rm -rf "$output_dir"
        mkdir -p "$output_dir"

        start_time=$(date +%s)

        if [ -x "/usr/bin/time" ]; then
            /usr/bin/time -v -o "$time_log" \
            nextflow run "$ANNOMAF_PIPELINE" \
              -params-file "params.json" \
              --input "$sheet" \
              --outdir "$output_dir" \
              -profile docker \
              --max_memory '10.GB' --max_cpus 8 --max_time '20.h' \
              &> "$stdout_log"
        else
            echo "⚠️  /usr/bin/time no disponible, ejecutando sin monitoreo detallado"
            nextflow run "$ANNOMAF_PIPELINE" \
              -params-file "params.json" \
              --input "$sheet" \
              --outdir "$output_dir" \
              -profile docker \
              --max_memory '10.GB' --max_cpus 8 --max_time '20.h' \
              &> "$stdout_log"
        fi

        end_time=$(date +%s)
        duration=$((end_time - start_time))

        # Info sistema
        cpu_model=$(grep -m1 "model name" /proc/cpuinfo | cut -d: -f2 | xargs)
        mem_total=$(grep MemTotal /proc/meminfo | awk '{print $2}')
        mem_total_mb=$((mem_total / 1024))

        # VCF original
        sample_vcf=$(tail -n +2 "$sheet" | head -n 1 | cut -d',' -f3)
        vcf_size=$(du -h "$sample_vcf" | cut -f1)
        total_variants=$(grep -v '^#' "$sample_vcf" | wc -l)

        # Buscar output .tsv
        output_maf=$(find "$output_dir" -type f -iname '*.maf' -path '*/vcf2maf/*' | head -n 1)
        if [ -f "$output_maf" ]; then
            unique_genes=$(cut -f13 "$output_maf" | grep -v "Hugo_Symbol" | sort | uniq | wc -l)
            accionables=$(grep -ic "actionable" "$output_maf")
        else
            unique_genes="NA"
            accionables="NA"
        fi

        # Extraer métricas
        if [ -f "$time_log" ]; then
            max_rss=$(grep "Maximum resident set size" "$time_log" | awk '{print $6 " KB"}')
            cpu_usage=$(grep "Percent of CPU this job got" "$time_log" | awk '{print $NF}')
        else
            max_rss="NA"
            cpu_usage="NA"
        fi

        echo "📊 [Resumen - $test_name]"
        echo "⏱️ Tiempo: ${duration}s"
        echo "💾 VCF size: $vcf_size"
        echo "📈 Variantes: $total_variants"
        echo "🧬 Genes únicos: $unique_genes"
        echo "🎯 Accionables: $accionables"
        echo "🧠 RAM: ${mem_total_mb}MB / Máximo RSS: $max_rss"
        echo "⚙️ CPU: $cpu_model"
        echo ""

        {
          echo "Test: $test_name"
          echo "Execution time (s): $duration"
          echo "VCF file size: $vcf_size"
          echo "Total variants: $total_variants"
          echo "Unique genes: $unique_genes"
          echo "Actionable variants: $accionables"
          echo "CPU: $cpu_model"
          echo "Total RAM (MB): $mem_total_mb"
          echo "Max RSS: $max_rss"
          echo "CPU usage: $cpu_usage"
        } > "$report_txt"
    )
}

# Buscar carpetas del tipo T*/annomaf
for dir in "$SCRIPT_DIR"/T*/annomaf; do
    ejecutar_annomaf "$dir"
done

echo "✅ Todos los casos ANNOMAF ejecutados correctamente."

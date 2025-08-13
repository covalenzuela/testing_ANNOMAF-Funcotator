#!/usr/bin/env bash
set -Eeuo pipefail

################################
# Config por defecto (pensado para ejecutar 1 nivel arriba de Panels/)
################################
PANELS_DIR="Panels"
ORIG_DIR="$PANELS_DIR/panel_test_vcf_RAW"
DEST_DIR="$PANELS_DIR/panel_test_vcf_hg38"
CHAIN="$PANELS_DIR/hg19ToHg38.over.chain.gz"
FASTA="/home/conifva/MathBioLab/git_works/OneDrive_1_20-5-2025/references/fasta_hg38/Homo_sapiens_assembly38.fasta"  # ⚠️ AJUSTA A TU ENTORNO

FORCE=0
VERBOSE=0

################################
# Utilidades y mensajes
################################
GREEN="\033[0;32m"; YELLOW="\033[1;33m"; RED="\033[0;31m"; NC="\033[0m"
log()  { echo -e "${GREEN}[$(date +'%F %T')]${NC} $*"; }
warn() { echo -e "${YELLOW}[$(date +'%F %T')] ⚠ ${NC} $*" >&2; }
err()  { echo -e "${RED}[$(date +'%F %T')] ✖ ${NC} $*" >&2; }

trap 'err "Fallo en la línea $LINENO. Abortando."' ERR

show_help() {
cat <<EOF
Uso: $(basename "$0") [opciones]

Convierte VCFs de hg19 a hg38 con CrossMap y limpia variantes no mapeadas.
Pensado para ejecutarse 1 nivel arriba de 'Panels/', usando por defecto:
  - Input : \$PANELS_DIR/panel_test_vcf_RAW
  - Output: \$PANELS_DIR/panel_test_vcf_hg38
  - Chain : \$PANELS_DIR/hg19ToHg38.over.chain.gz
  - FASTA : $FASTA

Opciones:
  -h, --help                  Mostrar esta ayuda y salir.
  -i, --input-dir DIR         Directorio de entrada (VCF/VCF.gz). Por defecto: $ORIG_DIR
  -o, --output-dir DIR        Directorio de salida. Por defecto: $DEST_DIR
  --chain FILE                Ruta al archivo chain. Por defecto: $CHAIN
  --fasta FILE                Ruta al FASTA de hg38. Por defecto: $FASTA
  -f, --force                 Reprocesar aunque los outputs existan y estén actualizados.
  -v, --verbose               Verbose (muestra más detalle).

Notas:
  • El script es idempotente: si existe <base>_hg38.clean.vcf y es más nuevo que el input,
    no se vuelve a convertir a menos que uses --force.
  • Acepta .vcf y .vcf.gz como entrada; genera .vcf como salida.
  • Asegúrate de que CrossMap, grep y gzip estén disponibles.

Ejemplos:
  $(basename "$0")
  $(basename "$0") -i Panels/panel_test_vcf_RAW -o Panels/panel_test_vcf_hg38 \\
                   --chain Panels/hg19ToHg38.over.chain.gz \\
                   --fasta /ruta/a/Homo_sapiens_assembly38.fasta
EOF
}

################################
# Parseo de argumentos
################################
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) show_help; exit 0 ;;
    -i|--input-dir) ORIG_DIR="$2"; shift 2 ;;
    -o|--output-dir) DEST_DIR="$2"; shift 2 ;;
    --chain) CHAIN="$2"; shift 2 ;;
    --fasta) FASTA="$2"; shift 2 ;;
    -f|--force) FORCE=1; shift ;;
    -v|--verbose) VERBOSE=1; shift ;;
    *) err "Opción desconocida: $1"; exit 1 ;;
  esac
done

################################
# Validaciones
################################
command -v CrossMap >/dev/null || { err "No se encontró 'CrossMap' en PATH."; }
command -v grep >/dev/null || { err "No se encontró 'grep' en PATH."; }
command -v gzip >/dev/null || { warn "No se encontró 'gzip' en PATH (requerido para .vcf.gz)."; }

[[ -d "$ORIG_DIR" ]] || { err "No existe el directorio de entrada: $ORIG_DIR"; }
mkdir -p "$DEST_DIR"

[[ -f "$CHAIN" ]] || { err "No se encuentra chain file: $CHAIN"; }
[[ -f "$FASTA" ]] || { err "No se encuentra FASTA: $FASTA
⚠️ Ajusta la ruta con --fasta o edita el valor por defecto en el script."; }

if [[ $VERBOSE -eq 1 ]]; then
  set -x
fi

log "Entrada : $ORIG_DIR"
log "Salida  : $DEST_DIR"
log "Chain   : $CHAIN"
log "FASTA   : $FASTA"

################################
# Función: convertir un VCF
################################
convert_one() {
  local vcf="$1"
  local name base out_raw out_clean

  name="$(basename "$vcf")"
  base="${name%.vcf}"
  base="${base%.vcf.gz}"   # si venía .vcf.gz
  out_raw="$DEST_DIR/${base}_hg38.vcf"
  out_clean="$DEST_DIR/${base}_hg38.clean.vcf"

  # Idempotencia: si clean existe y es más nuevo que input, saltar (a menos que --force)
  if [[ $FORCE -eq 0 && -f "$out_clean" && "$out_clean" -nt "$vcf" ]]; then
    log "✓ Ya convertido/limpio (más nuevo que input): $(basename "$out_clean"). Omito."
    return 0
  fi

  log "🔄 Convirtiendo $name → $(basename "$out_raw") ..."
  # CrossMap puede leer .vcf o .vcf.gz directamente
  CrossMap vcf "$CHAIN" "$vcf" "$FASTA" "$out_raw"

  log "🧹 Limpiando variantes no mapeadas → $(basename "$out_clean") ..."
  # Eliminamos líneas con el aviso de no conversión manteniendo encabezados
  # Si no hay coincidencias, grep devuelve 1; por eso usamos '|| true'
  { grep -v "was not converted" "$out_raw" > "$out_clean"; } || true

  log "✅ Listo: $(basename "$out_clean")"
}

################################
# Loop de conversión
################################
shopt -s nullglob
inputs=( "$ORIG_DIR"/*.vcf "$ORIG_DIR"/*.vcf.gz )
shopt -u nullglob

if [[ ${#inputs[@]} -eq 0 ]]; then
  warn "No se encontraron .vcf ni .vcf.gz en $ORIG_DIR"
else
  for vcf in "${inputs[@]}"; do
    convert_one "$vcf"
  done
  log "🎉 Conversión completa. Archivos en: $DEST_DIR"
fi

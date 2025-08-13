#!/usr/bin/env bash
# parseLBP.sh - Script para procesar archivos LBP específicos (desde bam)
# chmod +x parseLBP.sh

set -Eeuo pipefail

#########################
# Configuración
#########################
THREADS=8
LBP_DIR="LBP"
LBP_FILES=(
  "LBP_IL_T_1ng_1.bwa.dedup.bam"
  "LBP_IL_T_10ng_1.bwa.dedup.bam"
  "LBP_IL_T_100ng_1.bwa.dedup.bam"
)
CHR="chr10"

#########################
# Utilidades
#########################
GREEN="\033[0;32m"; YELLOW="\033[1;33m"; RED="\033[0;31m"; NC="\033[0m"
log()  { echo -e "${GREEN}[$(date +'%F %T')]${NC} $*"; }
warn() { echo -e "${YELLOW}[$(date +'%F %T')] ⚠ ${NC} $*" >&2; }
err()  { echo -e "${RED}[$(date +'%F %T')] ✖ ${NC} $*" >&2; }

trap 'err "Error en línea $LINENO. Abortando."; exit 1' ERR

require() {
  command -v "$1" &>/dev/null || { err "No se encontró '$1' en PATH"; exit 127; }
}

show_help() {
cat <<EOF
Uso: parseLBP [opciones]

Opciones:
  -h, --help       Mostrar esta ayuda y salir.
  -c, --chr STR    Cromosoma o región a extraer (por defecto: $CHR).
  -t, --threads N  Número de hilos a usar (por defecto: $THREADS).

Descripción:
  Procesa únicamente los BAM listados en el array LBP_FILES dentro de la carpeta $LBP_DIR.
  Para cada archivo:
    1) Verifica integridad con 'samtools quickcheck'
    2) Crea índice si falta
    3) Extrae la región especificada (-c), ordena e indexa el resultado.
  El script es idempotente: no repite pasos si los outputs ya existen y están actualizados.

Ejemplos:
  parseLBP --chr chr10 --threads 4
  parseLBP -c chrX -t 16
EOF
}

#########################
# Parseo de argumentos
#########################
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      show_help
      exit 0
      ;;
    -c|--chr)
      CHR="$2"
      shift 2
      ;;
    -t|--threads)
      THREADS="$2"
      shift 2
      ;;
    *)
      err "Opción desconocida: $1"
      exit 1
      ;;
  esac
done

#########################
# Chequeos iniciales
#########################
require samtools
[[ -d "$LBP_DIR" ]] || { err "No existe el directorio $LBP_DIR"; exit 1; }

read -r -p "¿Procesar archivos LBP? (si/no): " confirmar
[[ "${confirmar:-no}" == "si" ]] || { warn "Operación cancelada por el usuario."; exit 0; }

log "Iniciando pipeline para ${#LBP_FILES[@]} archivo(s) en $LBP_DIR ..."
declare -gA SKIP

#########################
# 0) quickcheck
#########################
for bam in "${LBP_FILES[@]}"; do
  bam_path="$LBP_DIR/$bam"
  if [[ ! -f "$bam_path" ]]; then
    warn "No se encontró $bam en $LBP_DIR. Se omite."
    continue
  fi
  if ! samtools quickcheck -v "$bam_path"; then
    warn "$bam podría estar corrupto. Se omitirá."
    SKIP["$bam_path"]=1
  fi
done

#########################
# 1) Indexado
#########################
for bam in "${LBP_FILES[@]}"; do
  bam_path="$LBP_DIR/$bam"
  [[ -f "$bam_path" ]] || continue
  [[ ${SKIP["$bam_path"]+x} ]] && continue
  if [[ ! -e "${bam_path}.bai" && ! -e "${bam_path}.csi" ]]; then
    log "Indexando $bam ..."
    samtools index -@ "$THREADS" "$bam_path"
  else
    log "Índice presente para $bam."
  fi
done

#########################
# 2) Procesar región
#########################
for bam in "${LBP_FILES[@]}"; do
  bam_path="$LBP_DIR/$bam"
  [[ -f "$bam_path" ]] || continue
  [[ ${SKIP["$bam_path"]+x} ]] && continue

  base="${bam%.bam}"
  out="$LBP_DIR/${base}_${CHR}.sorted.bam"
  out_bai="${out}.bai"
  out_csi="${out}.csi"

  if [[ -f "$out" && ( -f "$out_bai" || -f "$out_csi" ) && "$out" -nt "$bam_path" ]]; then
    log "Output ya listo para $bam. Omitido."
    continue
  fi

  log "Procesando $bam ..."
  tmp_out="$(mktemp "${out}.tmp.XXXXXX")"
  samtools view -@ "$THREADS" -b "$bam_path" "$CHR" \
    | samtools sort -@ "$THREADS" -o "$tmp_out" -
  samtools index -@ "$THREADS" "$tmp_out"

  mv -f "$tmp_out" "$out"
  [[ -f "${tmp_out}.bai" ]] && mv -f "${tmp_out}.bai" "$out_bai"
  [[ -f "${tmp_out}.csi" ]] && mv -f "${tmp_out}.csi" "$out_csi"

  log "Listo: $(basename "$out")"
done

log "✅ Proceso completado."


#!/usr/bin/env bash
set -Eeuo pipefail

################################
# Colores y logs
################################
GREEN="\033[0;32m"; YELLOW="\033[1;33m"; RED="\033[0;31m"; NC="\033[0m"
log()  { echo -e "${GREEN}[$(date +'%F %T')]${NC} $*"; }
warn() { echo -e "${YELLOW}[$(date +'%F %T')] ⚠ ${NC} $*" >&2; }
err()  { echo -e "${RED}[$(date +'%F %T')] ✖ ${NC} $*" >&2; }

trap 'err "Error en línea $LINENO. Abortando."; exit 1' ERR

################################
# Rutas base robustas
################################
SCRIPT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PANELS_DIR="$SCRIPT_DIR/Panels"
RAW_ARCHIVE="$SCRIPT_DIR/panel_test_vcf_RAW.tar.gz"

################################
# Defaults (config)
################################
THREADS=8
CHR="chr10"
ASK_LBP_DOWNLOAD=1
ASK_LBP_PARSE=1

# Endpoints
WES_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/analysis/SNVs/vcfs/WES"
WGS_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/analysis/SNVs/vcfs/WGS"
LBP_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/LBP"

# Archivos
WES_FILES=(
  "WES_FD_1.bwa.muTect2.vcf.gz"
#  "WES_FD_1.bwa.somaticSniper.vcf.gz"
  "WES_FD_1.bwa.strelka.vcf.gz"
  "WES_FD_2.bwa.muTect2.vcf.gz"
  "WES_FD_2.bwa.somaticSniper.vcf.gz"
#  "WES_FD_2.bwa.strelka.vcf.gz"
  "WES_FD_3.bwa.muTect2.vcf.gz"
#  "WES_FD_3.bwa.somaticSniper.vcf.gz"
#  "WES_FD_3.bwa.strelka.vcf.gz"
  "WES_IL_1.bwa.muTect2.vcf.gz"
  "WES_IL_1.bwa.somaticSniper.vcf.gz"
  "WES_IL_1.bwa.strelka.vcf.gz"
#  "WES_IL_2.bwa.muTect2.vcf.gz"
#  "WES_IL_2.bwa.somaticSniper.vcf.gz"
#  "WES_IL_2.bwa.strelka.vcf.gz"
#  "WES_IL_3.bwa.muTect2.vcf.gz"
#  "WES_IL_3.bwa.somaticSniper.vcf.gz"
#  "WES_IL_3.bwa.strelka.vcf.gz"
  "WES_NV_1.bwa.muTect2.vcf.gz"
#  "WES_NV_1.bwa.somaticSniper.vcf.gz"
#  "WES_NV_1.bwa.strelka.vcf.gz"
#  "WES_NV_2.bwa.muTect2.vcf.gz"
#  "WES_NV_2.bwa.somaticSniper.vcf.gz"
#  "WES_NV_2.bwa.strelka.vcf.gz"
#  "WES_NV_3.bwa.muTect2.vcf.gz"
#  "WES_NV_3.bwa.somaticSniper.vcf.gz"
#  "WES_NV_3.bwa.strelka.vcf.gz"
)

WGS_FILES=(
  "WGS_FD_1.bwa.muTect2.vcf.gz"
  "WGS_FD_1.bwa.somaticSniper.vcf.gz"
#  "WGS_FD_1.bwa.strelka.vcf.gz"
  "WGS_FD_2.bwa.muTect2.vcf.gz"
#  "WGS_FD_2.bwa.somaticSniper.vcf.gz"
#  "WGS_FD_2.bwa.strelka.vcf.gz"
  "WGS_FD_3.bwa.muTect2.vcf.gz"
#  "WGS_FD_3.bwa.somaticSniper.vcf.gz"
#  "WGS_FD_3.bwa.strelka.vcf.gz"
  "WGS_IL_1.bwa.muTect2.vcf.gz"
  "WGS_IL_1.bwa.somaticSniper.vcf.gz"
  "WGS_IL_1.bwa.strelka.vcf.gz"
#  "WGS_IL_2.bwa.muTect2.vcf.gz"
#  "WGS_IL_2.bwa.somaticSniper.vcf.gz"
#  "WGS_IL_2.bwa.strelka.vcf.gz"
#  "WGS_IL_3.bwa.muTect2.vcf.gz"
#  "WGS_IL_3.bwa.somaticSniper.vcf.gz"
#  "WGS_IL_3.bwa.strelka.vcf.gz"
  "WGS_NV_1.bwa.muTect2.vcf.gz"
#  "WGS_NV_1.bwa.somaticSniper.vcf.gz"
#  "WGS_NV_1.bwa.strelka.vcf.gz"
#  "WGS_NV_2.bwa.muTect2.vcf.gz"
#  "WGS_NV_2.bwa.somaticSniper.vcf.gz"
  "WGS_NV_2.bwa.strelka.vcf.gz"
#  "WGS_NV_3.bwa.muTect2.vcf.gz"
#  "WGS_NV_3.bwa.somaticSniper.vcf.gz"
#  "WGS_NV_3.bwa.strelka.vcf.gz"
)

LBP_FILES=(
  "LBP_IL_T_1ng_1.bwa.dedup.bam"
  "LBP_IL_T_10ng_1.bwa.dedup.bam"
  "LBP_IL_T_100ng_1.bwa.dedup.bam"
)

################################
# Ayuda
################################
show_help() {
cat <<EOF
Uso: $(basename "$0") [opciones]

Descarga datasets WES/WGS y (opcionalmente) LBP. Además:
  • Si existe 'panel_test_vcf_RAW.tar.gz' junto a este script, lo descomprime en Panels/.
  • Si descargas LBP, puedes ejecutar al final el parseo (extract/sort/index) o hacerlo posteriormente
    con 'parseLBP' de forma independiente.

Opciones:
  -h, --help            Mostrar esta ayuda y salir.
  -t, --threads N       Hilos para parseo LBP (por defecto: $THREADS).
  -c, --chr STR         Región/cromosoma para parseo LBP (por defecto: $CHR).
  --yes-lbp             No preguntar por descarga LBP (asumir "sí").
  --no-parse            No preguntar por parseo LBP (omitir parseo).
  --auto-parse          No preguntar y ejecutar parseo LBP automáticamente.

Notas:
  • Si existe 'parseLBP' en PATH, se usa; si no, se intenta './parseLBP.sh'.
  • Si ninguno está, se usa una implementación interna idempotente con samtools.
  • El parseo opera sobre ./LBP y crea *_<chr>.sorted.bam + índice.

Ejemplos:
  $(basename "$0")
  $(basename "$0") --yes-lbp --auto-parse -t 16 -c chr10:1-50000000
EOF
}

################################
# Parseo de argumentos
################################
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) show_help; exit 0 ;;
    -t|--threads) THREADS="$2"; shift 2 ;;
    -c|--chr) CHR="$2"; shift 2 ;;
    --yes-lbp) ASK_LBP_DOWNLOAD=0; shift ;;
    --no-parse) ASK_LBP_PARSE=0; shift ;;
    --auto-parse) ASK_LBP_PARSE=0; AUTO_PARSE=1; shift ;;
    *) err "Opción desconocida: $1"; exit 1 ;;
  esac
done

################################
# Requisitos y carpetas
################################
command -v wget >/dev/null || { err "No se encontró 'wget' en PATH"; exit 127; }
mkdir -p "$PANELS_DIR" "$SCRIPT_DIR/WES" "$SCRIPT_DIR/WGS" "$SCRIPT_DIR/LBP"

################################
# Descompresión de paneles RAW (archivo al lado del script)
################################
if [[ -f "$RAW_ARCHIVE" ]]; then
  log "Descargando archivos Panels..."
  dest_dir="$PANELS_DIR/panel_test_vcf_RAW"
  if [[ -d "$dest_dir" && -n "$(ls -A "$dest_dir" 2>/dev/null)" ]]; then
    warn "Panels/panel_test_vcf_RAW ya existe con contenido. Omitiendo descompresión."
  else
    log "Detectado $(basename "$RAW_ARCHIVE"). Preparando descompresión en Panels/ ..."
    # Detectar primer entry del tar para saber si viene con carpeta o archivos sueltos
    first_entry="$(tar -tzf "$RAW_ARCHIVE" | head -1 || true)"
    if [[ -z "$first_entry" ]]; then
      warn "El tar.gz parece vacío. Omitiendo."
    elif [[ "$first_entry" == panel_test_vcf_RAW/* ]]; then
      # El tar ya trae la carpeta; extraer directo a Panels/
      tar -xzvf "$RAW_ARCHIVE" -C "$PANELS_DIR"
      log "✅ Descomprimido dentro de Panels/."
    else
      # El tar trae archivos sueltos; crear la carpeta destino y extraer ahí
      mkdir -p "$dest_dir"
      tar -xzvf "$RAW_ARCHIVE" -C "$dest_dir"
      log "✅ Descomprimido en $dest_dir."
    fi
  fi
else
  warn "No se halló $RAW_ARCHIVE; se omite descompresión de paneles RAW."
fi

################################
# Funciones comunes (descarga y parseo LBP)
################################
descargar_si_no_existe() {
  local url="$1" dest="$2" file
  file=$(basename "$url")
  if [[ ! -f "$dest/$file" ]]; then
    log "Descargando $file ..."
    wget -q --show-progress -P "$dest" "$url"
  else
    warn "Ya existe: $file, se omite descarga."
  fi
}

parse_lbp_internal() {
  command -v samtools >/dev/null || { err "No se encontró 'samtools' para parsear LBP"; return 127; }
  log "Iniciando parseo LBP interno con $THREADS hilos y región '$CHR' ..."
  declare -A SKIP=()

  # quickcheck
  for bam in "${LBP_FILES[@]}"; do
    local bam_path="$SCRIPT_DIR/LBP/$bam"
    [[ -f "$bam_path" ]] || { warn "No encontrado $bam_path. Omito."; continue; }
    if ! samtools quickcheck -v "$bam_path"; then
      warn "$bam podría estar corrupto. Se omitirá."
      SKIP["$bam_path"]=1
    fi
  done

  # index de entrada
  for bam in "${LBP_FILES[@]}"; do
    local bam_path="$SCRIPT_DIR/LBP/$bam"
    [[ -f "$bam_path" ]] || continue
    [[ ${SKIP["$bam_path"]+x} ]] && continue
    if [[ ! -e "${bam_path}.bai" && ! -e "${bam_path}.csi" ]]; then
      log "Indexando $(basename "$bam_path") ..."
      samtools index -@ "$THREADS" "$bam_path"
    else
      log "Índice presente para $(basename "$bam_path")."
    fi
  done

  # extraer región + sort + index (idempotente)
  for bam in "${LBP_FILES[@]}"; do
    local bam_path="$SCRIPT_DIR/LBP/$bam"
    [[ -f "$bam_path" ]] || continue
    [[ ${SKIP["$bam_path"]+x} ]] && continue

    local base="${bam%.bam}"
    local out="$SCRIPT_DIR/LBP/${base}_${CHR}.sorted.bam"
    local out_bai="${out}.bai"
    local out_csi="${out}.csi"

    if [[ -f "$out" && ( -f "$out_bai" || -f "$out_csi" ) && "$out" -nt "$bam_path" ]]; then
      log "Output ya listo: $(basename "$out"). Omito."
      continue
    fi

    log "Procesando $(basename "$bam_path") → $(basename "$out") ..."
    local tmp_out
    tmp_out="$(mktemp "${out}.tmp.XXXXXX")"

    samtools view -@ "$THREADS" -b "$bam_path" "$CHR" \
      | samtools sort -@ "$THREADS" -o "$tmp_out" -
    samtools index -@ "$THREADS" "$tmp_out"

    mv -f "$tmp_out" "$out"
    [[ -f "${tmp_out}.bai" ]] && mv -f "${tmp_out}.bai" "$out_bai"
    [[ -f "${tmp_out}.csi" ]] && mv -f "${tmp_out}.csi" "$out_csi"

    log "Listo: $(basename "$out") (+ índice)."
  done

  log "✅ Parseo LBP interno completado."
}

maybe_run_parse_lbp() {
  if command -v parseLBP >/dev/null 2>&1; then
    log "Ejecutando 'parseLBP' del sistema (-t $THREADS -c $CHR) ..."
    parseLBP --threads "$THREADS" --chr "$CHR"
  elif [[ -x "$SCRIPT_DIR/parseLBP.sh" ]]; then
    log "Ejecutando './parseLBP.sh' (-t $THREADS -c $CHR) ..."
    "$SCRIPT_DIR/parseLBP.sh" --threads "$THREADS" --chr "$CHR"
  else
    warn "No se encontró 'parseLBP' ni './parseLBP.sh'. Uso implementación interna."
    parse_lbp_internal
  fi
}

################################
# Descargas
################################
log "Descargando archivos WES..."
for f in "${WES_FILES[@]}"; do
  descargar_si_no_existe "$WES_URL/$f" "$SCRIPT_DIR/WES"
done

log "Descargando archivos WGS..."
for f in "${WGS_FILES[@]}"; do
  descargar_si_no_existe "$WGS_URL/$f" "$SCRIPT_DIR/WGS"
done

# LBP (gigantes)
if (( ASK_LBP_DOWNLOAD )); then
  echo -e "${RED}¿Descargar archivos LBP? Son >90 GB cada uno.${NC}"
  read -r -p "Escribe 'si' para continuar o Enter para omitir: " confirmar
  LBP_DO="$([[ "${confirmar:-}" == "si" ]] && echo 1 || echo 0)"
else
  LBP_DO=1
fi

if (( LBP_DO )); then
  log "Descargando archivos LBP..."
  for f in "${LBP_FILES[@]}"; do
    descargar_si_no_existe "$LBP_URL/$f" "$SCRIPT_DIR/LBP"
  done

  # ¿Parsear LBP?
  if [[ -n "${AUTO_PARSE:-}" ]]; then
    maybe_run_parse_lbp
  elif (( ASK_LBP_PARSE )); then
    read -r -p "¿Ejecutar parseo LBP ahora (extraer '$CHR', sort+index)? (si/no): " p
    [[ "${p:-no}" == "si" ]] && maybe_run_parse_lbp || log "Omitiendo parseo LBP."
  else
    log "Parseo LBP omitido (flag --no-parse)."
  fi
else
  warn "Omitiendo descarga de LBP."
fi

log "✅ Finalizado."
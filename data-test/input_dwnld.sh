#!/bin/bash

set -e

# Colores para mensajes
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # Sin color

# Crear carpetas
mkdir -p Panels WES WGS LBP

# URL base
WES_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/analysis/SNVs/vcfs/WES"
WGS_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/analysis/SNVs/vcfs/WGS"
LBP_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/LBP"

# Función para descargar si no existe
descargar_si_no_existe() {
    local url="$1"
    local dest="$2"
    local file=$(basename "$url")
    if [ ! -f "$dest/$file" ]; then
        echo -e "${GREEN}Descargando $file...${NC}"
        wget -q --show-progress -P "$dest" "$url"
    else
        echo -e "${RED}Ya existe: $file, se omite descarga.${NC}"
    fi
}

# Archivos WES
echo -e "${GREEN}Descargando archivos WES...${NC}"
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

for file in "${WES_FILES[@]}"; do
  descargar_si_no_existe "$WES_URL/$file" "WES"
done

# Archivos WGS
echo -e "${GREEN}Descargando archivos WGS...${NC}"
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

for file in "${WGS_FILES[@]}"; do
  descargar_si_no_existe "$WGS_URL/$file" "WGS"
done

# Archivos LBP (gigantes)
echo -e "${RED}¿Deseas descargar los archivos de biopsia líquida (LBP)? Son mayores a 90 GB cada uno.${NC}"
read -p "Escribe 'si' para continuar o presiona Enter para omitir: " confirmar

if [[ "$confirmar" == "si" ]]; then
  echo -e "${GREEN}Descargando archivos LBP...${NC}"
  LBP_FILES=(
    "LBP_IL_T_100ng_1.bwa.dedup.bam"
    "LBP_IL_T_10ng_1.bwa.dedup.bam"
    "LBP_IL_T_1ng_1.bwa.dedup.bam"
  )

  for file in "${LBP_FILES[@]}"; do
    descargar_si_no_existe "$LBP_URL/$file" "LBP"
  done
else
  echo -e "${RED}Omitiendo descarga de LBP.${NC}"
fi

echo -e "${GREEN}\nFinalizado.${NC}"
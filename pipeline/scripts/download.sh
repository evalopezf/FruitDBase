#!/bin/bash

INPUT=$1
OUTDIR=$2
THREADS=4   # puedes ajustar en sbatch

mkdir -p "$OUTDIR"

if [ -z "$INPUT" ] || [ -z "$OUTDIR" ]; then
  echo "Uso: $0 archivo.tsv output_dir"
  exit 1
fi

# Detectar columna URL
URL_COL=$(head -n 1 "$INPUT" | tr '\t' '\n' | grep -n -i "^url$" | cut -d: -f1)

if [ -z "$URL_COL" ]; then
  echo "❌ No se encontró columna 'url'"
  exit 1
fi

echo "✔ Columna URL en posición: $URL_COL"

# -----------------------------
# función de descarga
# -----------------------------
download_one() {
  realurl=$1
  outdir=$2

  realurl=$(echo "$realurl" | tr -d '\r' | xargs)

  if [ -z "$realurl" ]; then
    return
  fi

  filename=$(basename "$realurl")
  outpath="$outdir/$filename"

  if [ -f "$outpath" ]; then
    echo "[EXISTS] $filename"
    return
  fi

  echo "[DOWNLOADING] $realurl"

  curl -L --fail --retry 3 --retry-delay 2 -o "$outpath" "$realurl"

  if [ $? -ne 0 ]; then
    echo "[ERROR] $realurl"
    rm -f "$outpath"
  fi
}
export -f download_one

# -----------------------------
# procesar TSV
# -----------------------------
tail -n +2 "$INPUT" | while IFS=$'\t' read -r line; do

  url=$(echo "$line" | cut -f"$URL_COL")
  url=$(echo "$url" | tr -d '\r' | xargs)

  if [ -z "$url" ] || [ "$url" = "NA" ]; then
    echo "[SKIP] sin URL"
    continue
  fi

  # expandir ; (ENA) y ,
  echo "$url" | tr ';' '\n' | tr ',' '\n' | \
  xargs -n 1 -P "$THREADS" -I {} bash -c 'download_one "$@"' _ {} "$OUTDIR"

done

echo "✔ Descarga completada"
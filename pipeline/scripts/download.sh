#!/bin/bash

INPUT=$1
OUTDIR=$2
THREADS=4   

mkdir -p "$OUTDIR"

if [ -z "$INPUT" ] || [ -z "$OUTDIR" ]; then
  echo "Uso: $0 archivo.tsv output_dir"
  exit 1
fi

#Get URL
URL_COL=$(head -n 1 "$INPUT" | tr '\t' '\n' | grep -n -i "^url$" | cut -d: -f1)

if [ -z "$URL_COL" ]; then
  echo " Column 'url' not found."
  exit 1
fi

echo "✔ URL column in position: $URL_COL"

# -----------------------------
# Download Function
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
# Processing TSV
# -----------------------------
tail -n +2 "$INPUT" | while IFS=$'\t' read -r line; do

  url=$(echo "$line" | cut -f"$URL_COL")
  url=$(echo "$url" | tr -d '\r' | xargs)

  if [ -z "$url" ] || [ "$url" = "NA" ]; then
    echo "[SKIP] sin URL"
    continue
  fi


  echo "$url" | tr ';' '\n' | tr ',' '\n' | \
  xargs -n 1 -P "$THREADS" -I {} bash -c 'download_one "$@"' _ {} "$OUTDIR"

done

echo "✔ Download completed"
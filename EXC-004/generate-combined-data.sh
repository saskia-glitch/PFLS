#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# ============================================================
# EXC-004 — Combine RAW-DATA into COMBINED-DATA
# - Mac + Linux compatible (Bash 3.2+)
# - Output (stdout) at the end: ONLY filenames in COMBINED-DATA (sorted)
# - No warnings / no extra text
# ============================================================

die() { printf "ERROR: %s\n" "$*" >&2; exit 1; }

# ============================================================
# SECTION 1 — Paths
# ============================================================
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${1:-$SCRIPT_DIR}"

RAW_DIR="$BASE_DIR/RAW-DATA"
OUT_DIR="$BASE_DIR/COMBINED-DATA"
MAP_FILE="$BASE_DIR/sample-translation.txt"

[[ -d "$RAW_DIR" ]] || die "Cannot find directory: $RAW_DIR"
cd "$BASE_DIR"

# ============================================================
# SECTION 2 — Recreate COMBINED-DATA (delete if exists)
# ============================================================
case "$OUT_DIR" in
  */COMBINED-DATA) ;;
  *) die "Safety stop: OUT_DIR does not end with /COMBINED-DATA: $OUT_DIR" ;;
esac

if [[ -e "$OUT_DIR" ]]; then
  [[ -d "$OUT_DIR" ]] || die "Path exists but is not a directory: $OUT_DIR"
  rm -rf "$OUT_DIR"
fi
mkdir -p "$OUT_DIR"

# ============================================================
# SECTION 3 — sample-translation mapping (DNAxx -> COxx) via TSV
# - Bash 3.2 compatible (no associative arrays)
# - If mapping missing: infer CO## from DNA##
# ============================================================
TMPDIR_USE="${TMPDIR:-/tmp}"
map_tsv="$(mktemp "$TMPDIR_USE/sample_map.XXXXXX")"
cleanup() { rm -f "$map_tsv"; }
trap cleanup EXIT
: > "$map_tsv"

if [[ -f "$MAP_FILE" ]]; then
  while IFS= read -r line; do
    [[ -z "${line//[[:space:]]/}" ]] && continue
    [[ "$line" =~ ^[[:space:]]*# ]] && continue

    set -- $line
    [[ $# -lt 2 ]] && continue

    lib=""; cult=""
    for tok in "$@"; do
      [[ -z "$lib"  && "$tok" =~ ^DNA[0-9]+$ ]] && lib="$tok"
      [[ -z "$cult" && "$tok" =~ ^CO[0-9]+$  ]] && cult="$tok"
    done

    if [[ -n "$lib" && -n "$cult" ]]; then
      printf "%s\t%s\n" "$lib" "$cult" >> "$map_tsv"
    else
      printf "%s\t%s\n" "$1" "$2" >> "$map_tsv"
    fi
  done < "$MAP_FILE"
fi

culture_for_lib() {
  local lib="$1" cult=""
  if [[ -s "$map_tsv" ]]; then
    cult="$(awk -v L="$lib" 'BEGIN{FS="\t"} $1==L {print $2; exit}' "$map_tsv" 2>/dev/null || true)"
    [[ -n "$cult" ]] && { printf "%s\n" "$cult"; return; }
  fi
  if [[ "$lib" =~ ^DNA([0-9]+)$ ]]; then
    printf "CO%s\n" "${BASH_REMATCH[1]}"
    return
  fi
  printf "%s\n" "$lib"
}

# ============================================================
# SECTION 4 — FASTA bin number extraction (bin-03.fasta -> 3)
# ============================================================
bin_num_from_fasta() {
  local bn="$1"
  if [[ "$bn" =~ ^bin-0*([0-9]+)\.(fa|fasta|fna)$ ]]; then
    printf "%d\n" "$((10#${BASH_REMATCH[1]}))"
  else
    printf "%s\n" ""
  fi
}

# ============================================================
# SECTION 5 — Robust MAG-bin extraction from checkm.txt (fixed-width parsing)
# - Find header positions of "Completeness" and "Contamination"
# - Extract values via substr() using those positions
# - Determine MAG if completeness>=50 and contamination<=5
# - Output: unique sorted bin numbers (integers)
# ============================================================
build_mag_bins_from_checkm() {
  local checkm_file="$1"
  local out_file="$2"

  awk '
    function trim(s){ sub(/^[ \t]+/,"",s); sub(/[ \t]+$/,"",s); return s }

    BEGIN{
      cpos=0; tpos=0; spos=0;
    }

    {
      line=$0
      gsub(/\r/,"",line)

      # detect header line
      if (line ~ /Completeness/ && line ~ /Contamination/) {
        cpos = index(line, "Completeness")
        tpos = index(line, "Contamination")
        spos = index(line, "Strain heterogeneity")  # may be 0
        next
      }

      if (cpos==0 || tpos==0) next
      if (line ~ /^-+$/) next
      if (line ~ /^[[:space:]]*$/) next

      # first token = Bin Id
      if (!match(line, /^[[:space:]]*[^[:space:]]+/)) next
      id = trim(substr(line, RSTART, RLENGTH))
      if (id !~ /bin/) next

      # trailing digits = bin number
      if (!match(id, /[0-9]+$/)) next
      n = substr(id, RSTART, RLENGTH) + 0

      # fixed-width extract
      comp = trim(substr(line, cpos, tpos - cpos))
      if (spos > tpos) cont = trim(substr(line, tpos, spos - tpos))
      else            cont = trim(substr(line, tpos))

      if (!(comp ~ /^[0-9]+(\.[0-9]+)?$/ && cont ~ /^[0-9]+(\.[0-9]+)?$/)) next
      comp += 0; cont += 0

      if (comp >= 50 && cont <= 5) print n
    }
  ' "$checkm_file" | LC_ALL=C sort -n -u > "$out_file"
}

# ============================================================
# SECTION 6 — Reheader FASTA deflines (unique + culture-associated)
# Header format: >CO64_MAG_001|000001|<old header>
# ============================================================
reheader_fasta() {
  local in_fa="$1" out_fa="$2" prefix="$3"
  awk -v p="$prefix" '
    BEGIN{n=0}
    /^>/{
      n++
      h=$0; sub(/^>/,"",h)
      printf(">%s|%06d|%s\n", p, n, h)
      next
    }
    {print}
  ' "$in_fa" > "$out_fa"
}

# ============================================================
# SECTION 7 — Main loop over RAW-DATA/DNAxx
# - Copy checkm/gtdb files to COMBINED-DATA with culture prefix
# - Determine MAG bin numbers once per culture
# - Classify each bin FASTA by membership in MAG set
# - Write *_MAG_###.fa and *_BIN_###.fa with sequential numbering
# ============================================================
for libdir in "$RAW_DIR"/*; do
  [[ -d "$libdir" ]] || continue

  lib="$(basename "$libdir")"
  culture="$(culture_for_lib "$lib")"

  bins_dir="$libdir/bins"
  [[ -d "$bins_dir" ]] || die "Missing bins/ directory in: $libdir"

  checkm_path="$libdir/checkm.txt"
  gtdb_path="$libdir/gtdb.gtdbtk.tax"
  [[ -f "$checkm_path" ]] || die "Missing checkm.txt in: $libdir"
  [[ -f "$gtdb_path"   ]] || die "Missing gtdb.gtdbtk.tax in: $libdir"

  cp -f "$checkm_path" "$OUT_DIR/${culture}-CHECKM.txt"
  cp -f "$gtdb_path"   "$OUT_DIR/${culture}-GTDB-TAX.txt"

  fasta_files=( "$bins_dir"/*.fa "$bins_dir"/*.fasta "$bins_dir"/*.fna )

  # unbinned
  unbinned=""
  if [[ -f "$bins_dir/bin-unbinned.fasta" ]]; then
    unbinned="$bins_dir/bin-unbinned.fasta"
  else
    for f in "${fasta_files[@]}"; do
      b="$(basename "$f")"
      [[ "$b" =~ unbinned ]] && { unbinned="$f"; break; }
    done
  fi

  if [[ -n "$unbinned" ]]; then
    reheader_fasta "$unbinned" "$OUT_DIR/${culture}_UNBINNED.fa" "${culture}_UNBINNED"
  fi

  # Build MAG bin-number set
  mag_bins_file="$(mktemp "$TMPDIR_USE/mag_bins.XXXXXX")"
  build_mag_bins_from_checkm "$checkm_path" "$mag_bins_file"

  # Store "binNumber<TAB>filepath" for stable sorting/numbering
  mag_list_file="$(mktemp "$TMPDIR_USE/mag_list.XXXXXX")"
  bin_list_file="$(mktemp "$TMPDIR_USE/bin_list.XXXXXX")"
  : > "$mag_list_file"
  : > "$bin_list_file"

  for f in "${fasta_files[@]}"; do
    [[ -f "$f" ]] || continue
    [[ -n "$unbinned" && "$f" == "$unbinned" ]] && continue

    bn="$(basename "$f")"
    n="$(bin_num_from_fasta "$bn")"

    if [[ -z "$n" ]]; then
      printf "999999\t%s\n" "$f" >> "$bin_list_file"
      continue
    fi

    if grep -Fxq "$n" "$mag_bins_file"; then
      printf "%d\t%s\n" "$n" "$f" >> "$mag_list_file"
    else
      printf "%d\t%s\n" "$n" "$f" >> "$bin_list_file"
    fi
  done

  sort -n -k1,1 "$mag_list_file" > "$mag_list_file.sorted" || true
  sort -n -k1,1 "$bin_list_file" > "$bin_list_file.sorted" || true

  mag_i=0
  while IFS=$'\t' read -r _ ff; do
    [[ -z "$ff" ]] && continue
    mag_i=$((mag_i+1))
    num="$(printf "%03d" "$mag_i")"
    reheader_fasta "$ff" "$OUT_DIR/${culture}_MAG_${num}.fa" "${culture}_MAG_${num}"
  done < "$mag_list_file.sorted"

  bin_i=0
  while IFS=$'\t' read -r _ ff; do
    [[ -z "$ff" ]] && continue
    bin_i=$((bin_i+1))
    num="$(printf "%03d" "$bin_i")"
    reheader_fasta "$ff" "$OUT_DIR/${culture}_BIN_${num}.fa" "${culture}_BIN_${num}"
  done < "$bin_list_file.sorted"

  rm -f "$mag_bins_file" \
        "$mag_list_file" "$bin_list_file" \
        "$mag_list_file.sorted" "$bin_list_file.sorted"
done

# ============================================================
# SECTION 8 — Required output: ONLY filenames in COMBINED-DATA
# ============================================================
LC_ALL=C ls -1 "COMBINED-DATA" | LC_ALL=C sort

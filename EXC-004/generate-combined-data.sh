mkdir -p COMBINED-DATA

declare -A CULTURE_MAP

while IFS=$'\t' read -r library culture rest; do
    if [ -n "$library" ] && [ "$library" != "library" ]; then
        CULTURE_MAP[$library]=$culture
    fi
done < RAW-DATA/sample-translation.txt


for dir in RAW-DATA/*/; do

    library=$(basename "$dir")

    culture=${CULTURE_MAP[$library]}

    if [ -f "${dir}checkm.txt" ]; then
        cp "${dir}checkm.txt" "COMBINED-DATA/${culture}-CHECKM.txt"
    else
        echo "  WARNING: checkm.txt not found for $library"
    fi

    if [ -f "${dir}gtdb.gtdbtk.tax" ]; then
        cp "${dir}gtdb.gtdbtk.tax" "COMBINED-DATA/${culture}-GTDB-TAX.txt"

    else
        echo "  WARNING: gtdb.gtdbtk.tax not found for $library"
    fi

    mag_counter=1
    bin_counter=1

    for fasta in "${dir}bins/"*.fasta; do
        [ -e "$fasta" ] || continue

        bin_name=$(basename "$fasta" .fasta)
        if [[ "$bin_name" == *"unbinned"* ]]; then
            output_file="COMBINED-DATA/${culture}_UNBINNED.fa"

            awk -v culture="$culture" -v bin="UNBINNED" '{
                if ($0 ~ /^>/) {
                    print ">" culture "_" bin "_" substr($0, 2)
                } else {
                    print $0
                }
            }' "$fasta" > "$output_file"

            continue
        fi

        comp=""
        contam=""

        if [ -f "${dir}checkm.txt" ]; then
            read comp contam < <(awk -v bin="$bin_name" '$1 ~ bin && $1 != "Bin" {print $13, $14; exit}' "${dir}checkm.txt")
        fi

        is_mag=0
        if [ -n "$comp" ] && [ -n "$contam" ]; then
            if (( $(echo "$comp >= 50" | bc -l) )) && (( $(echo "$contam < 5" | bc -l) )); then
                is_mag=1
            fi
        fi

        if [ $is_mag -eq 1 ]; then
            seq_num=$(printf "%03d" $mag_counter)
            output_file="COMBINED-DATA/${culture}_MAG_${seq_num}.fa"
            mag_counter=$((mag_counter + 1))
        else
            seq_num=$(printf "%03d" $bin_counter)
            output_file="COMBINED-DATA/${culture}_BIN_${seq_num}.fa"
            bin_counter=$((bin_counter + 1))
        fi

        if [ $is_mag -eq 1 ]; then
            bin_identifier="MAG_${seq_num}"
        else
            bin_identifier="BIN_${seq_num}"
        fi

        awk -v culture="$culture" -v bin_id="$bin_identifier" '{
            if ($0 ~ /^>/) {
                print ">" culture "_" bin_id "_" substr($0, 2)
            } else {
                print $0
            }
        }' "$fasta" > "$output_file"

    done
done
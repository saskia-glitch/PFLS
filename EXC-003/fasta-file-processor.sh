# EXC-003


for f in $@
    do

        num_seq=$(grep '>' $f | wc -l)
        len_seq=$(grep '>' $f | wc -c)
        long_seq=0
        short_seq=0
        avg_seq=0
        GC_cont=0 #$(echo  | awk '{gc_count += gsub(/[GgCc]/, "", $1)} END {print gc_count}')


        echo "FASTA File Statistics:"
        echo "----------------------"
        echo "Number of sequences: $num_seq"
        echo "Total length of sequences: $len_seq"
        echo "Length of the longest sequence: $long_seq"
        echo "Length of the shortest sequence: $short_seq"
        echo "Average sequence length: $avg_seq"
        echo "GC Content (%): $GC_cont"
done


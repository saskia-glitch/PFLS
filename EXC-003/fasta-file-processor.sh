# EXC-003
# 06/02/26
# finding the Number of sequences, the total length of sequences, length of the longest and 
# shortest sequence, the avarage sequence length and the GC content

# iterating over f as the given filename by user
for f in $@
do 
    # get the Number of sequences
    num_seq=$(grep '>' $f | wc -l)

    # Get the total length of sequences
    len_seq=$(awk '!/>/{gc_count += gsub(/[AaTtGgCc]/, "")} END {print gc_count}' $f)

    # get the longest sequence
    long_seq=$(awk '
        /^>/ {
            if (c > 0 && c > max) 
            max = c
            c = 0
            next
        }
        {
            gsub(/[^AaTtCcGg]/, "", $0)
            c += length($0)
        }
        END {
            if (c > 0 && c > max) 
            max = c
            print max
        }
    ' "$f")
    # If a line starts with ">", it is a header line
    # Check if this length is larger than the current maximum (max):
    # max is not defined in first iteration!
    # if (c > 0 && c > max) max = c
    # new max found, so c = 0
    # new sequence starts

    # get the shortest sequence
    short_seq=$(awk '
        /^>/ {
            if (c > 0) {if (min == 0 || c < min) 
            min = c
            }
            c = 0
            next
        }
        {
            gsub(/[^AaTtCcGg]/, "", $0)
            c += length($0)
        }
        END {
            if (c > 0) {
            if (min == 0 || c < min) min = c
            }
            print min
        }
    ' "$f")

    # get the GC content
    GC_cont=$(awk '!/>/{gc_count += gsub(/[GgCc]/, "")} END {print gc_count}' $f)




    # formatting the output
    echo "FASTA File Statistics:"
    echo "----------------------"
    echo "Number of sequences: $num_seq"
    echo "Total length of sequences: $len_seq"
    echo "Length of the longest sequence: $long_seq"
    echo "Length of the shortest sequence: $short_seq"
    echo "Average sequence length:" $(echo "scale=3; $len_seq/$num_seq" | bc) # calculating the Avarage
    echo "GC Content (%):" $(echo "scale=3; $GC_cont *100 / $len_seq" | bc) # calculating the GC content in %
done



# How to GitHub:
# file into git directory
# terminal: git add .
# terminal: git commit -m "my changes"
# terminal: git push -u origin main
# check GitHub!
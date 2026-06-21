#!/bin/bash

#!/usr/bin/env bash

input="$1"

if [[ -z "$input" ]]; then
    echo "Usage: $0 assembly.fna"
    exit 1
fi

awk '
BEGIN {
    outfile = ""
}

/^>/ {
    if (match($0, /chromosome: *([^ ]+)/, arr)) {
        chr = arr[1]
        outfile = "chr_" chr ".fasta"
    } else {
        outfile = "unclassified.fasta"
    }

    print $0 > outfile
    next
}

{
    print $0 >> outfile
}
' "$input"

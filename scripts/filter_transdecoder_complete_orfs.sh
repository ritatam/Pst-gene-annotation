#!/bin/bash
awk 'BEGIN {RS=">";FS="\n";OFS=""} NR>1 {print ">"$1; $1=""; print}' "$1" | 
grep "complete" -A1 --no-group-separator > "$2"
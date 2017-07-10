#!/usr/bin/awk -f
BEGIN{
    OFS="\t"
}
{
    if (NR==FNR){
        size_arr[$1] = $2
    } else {
        gsub("ID=gene:", "", $9)
        end=$4+l-1
        end=end<size_arr[$1]?end:size_arr[$1];
        if ($4 != 1){
            print $1, $4-1, end, $9"_fwd", ".", "+"
        }
        if ($5 != size_arr[$1]){
            start=$5-l
            start=start<1?1:start
            print $1, start, $5, $9"_rev", ".", "-"
        }
    }
}

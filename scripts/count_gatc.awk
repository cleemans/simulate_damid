#!/usr/bin/awk -f
{
    gsub(/_[fr][we][dv]/, "", $4)
    if ($6=="+"){
        start[$4] = $1"\t"$2"\t"$3
        if ($2==$8){
            count_fwd[$4]++
        }
    } else if ($6=="-" && $3==$9){
        count_rev[$4]++
    }
}
END{
    OFS="\t"
    print "seqnames", "start", "end", "fwd", "rev", "nreads"
    for (key in start){
        fwd = (key in count_fwd)?count_fwd[key]:0
        rev = (key in count_rev)?count_rev[key]:0
        print start[key], fwd, rev, fwd + rev
    }
}

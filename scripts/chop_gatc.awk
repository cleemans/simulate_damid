#!/usr/bin/awk -f
{
    a = gensub(/(GATC[ACTGacgtNn]+)GATC[ACTGacgtNn]*/, "\\1", $0)
    print a
}

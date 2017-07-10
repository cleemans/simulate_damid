import itertools
import re
import subprocess


def run_shell(command):
    p = subprocess.Popen(['/bin/bash', '-c', command],
                         stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    try:
        outs, errs = p.communicate()
    except subprocess.TimeoutExpired:
        p.kill()
        outs, errs = p.communicate()
    return(outs)


fwd_command = ("<(bedtools bamtobed -i %s | "
               "awk '{"
               "  if (($4 ~ /fwd/ && $6==\"+\") || "
               "      ($4 ~ /rev/ && $6==\"-\")){"
               "      print $1\"\t\"$2\"\t\"$2+4\"\t\"$4\"\t\"$5\"\t\"$6"
               "  }}')")

rev_command = ("<(bedtools bamtobed -i %s | "
               "awk '{"
               "  if (($4 ~ /fwd/ && $6==\"+\") || "
               "      ($4 ~ /rev/ && $6==\"-\")){"
               "      print $1\"\t\"$3-4\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6"
               "  }}')")

command_fmt = "bedtools coverage -a {} -b %s" % snakemake.input.gatc[0]

input_list = snakemake.input.split + snakemake.input.mismatch
command_list = [command_fmt.format(fwd_command % f) for f in input_list]
command_list.extend([command_fmt.format(rev_command % f) for f in input_list])

count_fmt = '{};FWD_MATCH=%s/%s:{}/{};REV_MATCH=%s/%s:{}/{};MISMATCH={}' % tuple(snakemake.params[0] * 2)

out_list = [str(run_shell(cmd)).split('\\n') for cmd in command_list]

print(command_list[0])
with open(snakemake.output[0], 'w') as fout:
    for line_tuple in itertools.zip_longest(*out_list):
        line_split_list = [line.split('\\t') for line in line_tuple]
        print(line_split_list[0])
        gff_list = line_split_list[0][0:8]
        gff_attribute = line_split_list[0][8]
        count_list = [split[9] for split in line_split_list]
        mismatch = sum(int(count_list[i]) for i in c(2,3,6,7))
        count_str = count_fmt.format(gff_attribute, count_list[0],
                                     count_list[1], mismatch,
                                     count_list[4], count_list[5])
        gff_list.append(count_str)
        fout.write('\t'.join(gff_list))

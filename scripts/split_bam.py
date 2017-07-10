import pysam
import random
import itertools
import math
import re


sam_list = [pysam.AlignmentFile(b, 'rb', check_sq=False) for b in snakemake.input.bam]
sam_out = [pysam.AlignmentFile(out, "wb", template=sam_list[0]) for out in snakemake.output]
for line_tuple in itertools.zip_longest(*sam_list):
    match_list = [-math.inf if line.is_unmapped else line.get_tag('AS')
                  for line in line_tuple]
    if match_list[0] > match_list[1]:
        sam_out[0].write(line_tuple[0])
    elif match_list[1] > match_list[0]:
        sam_out[1].write(line_tuple[1])
    else:
        sam_out[2].write(line_tuple[round(random.random())])
[b.close() for b in sam_list]
[b.close() for b in sam_out]

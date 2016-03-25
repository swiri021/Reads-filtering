[Instruction]

This is a class for filtering out designated reads in BAM file.

In case of multi-pools in targeted sequencing or multiple reads-groups in sample region, you can filter out select filtering out specific read-group. 

If you want to remove useless reads in some amplicon region or specific region of reads, just import this code.

This code calculates a distance between a target and a background, and wipes out a target.

This code is only for Ion Torrent PGM, if you want to change it to use other platform, just let me know.



[USAGE]

from filtering_amplicon.py import FILTER_AMP
…
…
filter = FILTER_AMP(input_BAM, input_BED)
filter.write_proc()

[INPUT BAM]

Ion Torrent Suite BAM

[INPUT BED]

format :

chromosome[tab]start_position[tab]end_position[tab]description[tab]flag

flag = T : Filtering target
flag = B : Background read group to compare target (NOT filtering out reads)

example :

chr7	55242328	55242488	EGFR_EXON19	T
chr7	55242434	55242585	EGFR_EXON19	B
chr3	178951823	178951971	PIK3CA_EXON20	B
chr3	178951920	178952090	PIK3CA_EXON20	T
chr3	178952037	178952212	PIK3CA_EXON20	B



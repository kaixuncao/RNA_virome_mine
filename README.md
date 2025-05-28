This project is supplementary material to the article “Spiders are capsules of viral diversity and neglected vectors with zoonotic risks”. Here we provide, where possible, the script code used for assembly and RNA lookup and details of its implementation. In addition we provide the evolutionary tree from the Supplementary Material as supplementary material.

The specific steps are as follows:
1. the data were first quality controlled and assembled, and potential RdRP fragments were queried by homology;
2. after that, fragments ≥75% of the length were de-compared by the provided R language script (virus_selected.r) to ensure the integrity of the fragments;
3. The command to remove the corresponding fragments by bedtools is: 
bedtools getfasta -fi potential_virus_rmdup.fasta -bed $Folde/RdRp_grep_list.bed -fo$Folde/RdRp_list_75grep.fasta
4. Use the Ref sequences provided by ICTV as a reference for evolutionary tree construction and exclude branches on the outside of the evolutionary tree (branch positions are not trusted);
5. further refined evolutionary tree construction by IQ-Tree.

Here we express our special thanks to UriNeri, some of the code in this topic is borrowed from previous projects: https://github.com/UriNeri/RVMT

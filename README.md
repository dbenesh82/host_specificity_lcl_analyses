# Host specificity analysis

This repository analyzes patterns of host specificity in complex life cycle parasites. Complex life cycle parasites infect dissimilar hosts during their life, e.g. an invertebrate as first host and a vertebrate as second host. Do complex life cycle parasites infect more hosts than simple life cycle parasites? How does specificity vary within cycles (i.e. larval vs adult stage)? Is the ability to infect many different kinds of hosts costly? I use a [database](http://onlinelibrary.wiley.com/doi/10.1002/ecy.1680/suppinfo) of helminth (parasitic worm) life cycles to explore these questions. The analysis is broken in to two parts. First, there is a [script](get_taxonomy/get_clean_host_taxonomy.RMD) to download host taxonomy for calculating a host specificity index. Second, there is a [script](calc_specificity/calc_specificity_patterns.RMD) that analyzes patterns of host specificity.
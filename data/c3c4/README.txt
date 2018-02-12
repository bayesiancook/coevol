Alignment, tree and data about photosynthetic regimes are from the following original study: Kapralov MV, Smith JAC and Filatov DA (2012). Rubisco evolution in C4 eudicots: an analysis of Amaranthaceae sensu lato. PLoS One 7(12):e52974. They were kindly provided by Maxim Kapralov. If you use the files in this folder, please cite this original study.

To reproduce the analysis published by Parto and Lartillot, 2018:
- run the inddiffsel program, with the following command;
inddiffsel -d C4Amaranthaceae.ali -t C4Amaranthaceae.tree -x 1 1100 -ncond 3 c3c4dirmodel

- read the chain and output posterior probabilities for position- and amino-acid-specific differential selection, with the following command:
readinddiffsel -x 100 1 1100 c3c4dirmodel


Alignment, tree and data about photosynthetic regimes are from the following original study: Kapralov MV, Smith JAC and Filatov DA (2012). Rubisco evolution in C4 eudicots: an analysis of Amaranthaceae sensu lato. PLoS One 7(12):e52974. They were kindly provided by Maxim Kapralov. If you use the files in this folder, please cite this original study.

To reproduce the analysis published by Parto and Lartillot, 2018:

- run the inddiffsel program (differential selection model, such as described in Parto and Lartilot, 2018), with the following command:
inddiffsel -d C4Amaranthaceae.ali -t C4Amaranthaceae.tree -c C4Amaranthaceae.photoregime -x 1 1100 -ncond 3 c3c4dir31
inddiffsel -d C4Amaranthaceae.ali -t C4Amaranthaceae.tree -c C4Amaranthaceae.photoregime -x 1 1100 -ncond 2 c3c4dir21

the -ncond option specifies the number of conditions for the analysis; here, there are 2 states (c3=1, c4=2), such as specified in C4Amaranthaceae.photoregime.The 3-condition analysis attributes condition 1 and 2 for all maximal monophyletic c3 and c4 groups, and condition 0 for all interior branches. With -ncond 2, the analysis attributes condition 1 for all maximal monophyletic c4 groups and condition 0 for all other branches (thus all assumed c3).

- read the chain and output posterior probabilities for position- and amino-acid-specific differential selection, with the following command:
readinddiffsel -x 100 1 -1 -c 0.9 c3c4dir31

-x <burnin> <every> <until>
-c <cutoff>: writes out only those sites and amino-acids with posterior probability for a differential effect pp > 1-c or pp < c. If pp>1-c, then the differental effect is positive, if pp<c, the effect is negative.

- similarly you can run the diffom program (model with site- and condition-specific dN/dS), with either 2 or 3 conditions
diffom -d C4Amaranthaceae.ali -t C4Amaranthaceae.tree -c C4Amaranthaceae.photoregime -x 1 1100 -ncond 1 c3c4om11
diffom -d C4Amaranthaceae.ali -t C4Amaranthaceae.tree -c C4Amaranthaceae.photoregime -x 1 1100 -ncond 3 c3c4om31

- read the chains and output posterior mean omega values (and posterior probabilities for omega > 1)
readdiffom -x 100 1 -1 -c 0.9 c3c4om31


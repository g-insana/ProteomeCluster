
Script to create a presence/absence (and count) matrix where the columns
are proteomes, the rows are clusters and the values are either 0 if there are no proteins from a
column proteome in the row cluster or a positive integer corresponding to the number of proteins from
that proteome in the cluster. The matrix will be saved in tsv format.

Minimal required input is a filtered protein clusters file (-i)
  (like the output of the filter_clusters script) whose first two columns contains cluster_id and protein_ids (space separated entries in the form proteome_id:protein_id)
and a file containing a list of proteome_ids (-p)
  (assumed to be all the proteome_ids present in the protein clusters file).

Sample input file:
    cluster_id	protein_ids	proteins_count	proteomes_count	representative	seqlen_range
    0	35497:ENSSSCP00000055324|661 4698918:ENSSSCP00035021320|596 4698922:ENSSSCP00055011301|568	3	3	ENSSSCP00000055324|661	568-661

Sample output file:
    cluster	35497	4229143	4698918	4698922
          0	1	0	1	1

If the proteomes list file contains two columns, the second column will be assumed to contain proteome labels which will be used to replace the proteome_ids from the output file as column labels.
The order of proteomes in the proteomes list file also determines the orders of the columns in the output matrix.

An optional clusterslabel file can be specified (-l), which will be used to name the rows. This file should contain two columns: the first one containing cluster_ids and the second one containing the cluster labels.

An optional additional file (-a) can be specified. This file is assumed to contain singleton clusters (clusters with protein(s) from a single proteome). The second column of this file will be considered to contain cluster_ids, the third column protein_id and the fourth column containing proteome_ids.

If the -t option is passed, a column "proteomes_count" will be added (as 2nd column) with the number of proteomes having at least one protein for the cluster corresponding to each row (same as the 4th column of the clusters input file).

If the -c option is passed, a final row "clusters_per_proteome" will be added with the number of clusters present for each proteome (number of clusters in which the proteome is present with at least one protein).

An optional cluster members file can be created, specifying its filename (-m).

Example call:
    presence_matrix.py -i Protein_clusters_m2.tsv -o matrix.tsv -p proteomelabels.tsv -l clusterlabels.tsv -a add_singletons.acc -t

| Options | Values  | Help |
| ------- | ------- | ---- |
| <pre>-h --help</pre> | Flag. | show this help message and exit |
| <pre>-i --input\_file</pre> | INPUT_FILE<br/>Required.<br/>Type: file | Path to the input filefirst column cluster\_id, second column protein\_ids \(proteome\_id\:protein\_id space separated pairs\) |
| <pre>-p --proteomes</pre> | PROTEOMES<br/>Required.<br/>Type: file | Path to the list of proteomes or proteome labels mapping file\: either one \(proteome\_ids\) or two columns \(proteome\_id proteome\_label\) |
| <pre>-o --out\_file</pre> | OUT_FILE<br/>Required.<br/>Type: str | Path to the matrix output file |
| <pre>-m --members\_file</pre> | Optional.<br/>Type: str | Path to the optional cluster members file |
| <pre>-l --labels</pre> | Optional.<br/>Type: file | optional cluster labels mapping file\: two columns \(cluster\_id cluster\_label\)\; NOTE that if same label is applied to different cluster\_id, the corresponding rows will be merged |
| <pre>-a --additional</pre> | Optional.<br/>Type: file | optional .acc file for additional singleton clusters to add to the output matrix\: 2nd col \= cluster\_id, 4th col \= proteome\_id |
| <pre>-t --totals</pre> | Flag. | Add a column with total proteomes per cluster |
| <pre>-c --counts</pre> | Flag. | Add a final row with cluster counts per proteome |
| <pre>-s --strict</pre> | Flag. | Ensure totals match \(if -t used\) |

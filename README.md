# miR2Pathway
miR2Pathway is a PageRank-based method that can be used to rank disease risk of miRNA-mediated biological pathways. miR2Pathway can help explore how much miRNAs differentially influence the activity of biological pathways between two classes of phenotypes.

For information about the method, please read:

Chaoxing Li and Valentin Dinu. “miR2Pathway: A Novel Analytical Method to Discover MicroRNA-mediateding Dysregulated Pathways Involved in Disease”, 2017.

Original website: http://dinulab.org/tools/mir2pathway/

R code, sample data and sample scripts: 

miR2Pathway.R

The R script to perform miR2Pathway

miR2Pathway-example.R

Example usage of miR2Pathway.R

miR2Pathway-example-data.Rdata

Data used in miR2Pathway-example.R

# Usage

miR2Pathway <- function(mydata.gene,mydata.miR,genelist,name.genelist,miRlist,miRlist.full,N.miR,N.gene,N.path,Num.sample.normal,Num.sample.case,Pathway.database,cor.cutoff,N.parallel)

# Arguments
mydata.gene: 	a gene expression dataset (matrix). Rows represent genes, and columns represent samples (from control to case).
mydata.miR: 	a miRNA expression dataset (matrix). Rows represent miRNAs, and columns represent samples (from control to case)
genelist: 	a list of gene names (e.g., gene symobol).
name.genelist: 	a matrix of gene expression dataset (The order of gene names must be consistent with genelist). This matrix is 1*N, where N is the number of genes. The value of each element could be a ramdom value. See an example in Examples Section.
miRlist: 	miRNA in a standard format
miRlist.full: 	the full name of all miRNAs.

# Value

T.score: the total differential influence of all the miRNAs on the activity of a single pathway between control and case.

LengthOfPathway: gene counts of pathways.

# Examples

See miR2Pathway-example.R

An example for the matrix of name.genelist:

  100130426      100133144      100134869      10357   10431   136542

  Sample1       0            0               0            0         0        0

 

# Credits

Authors: Chaoxing Li and Valentin Dinu

If you use or modify the code, please cite:

“miR2Pathway: A Novel Analytical Method to Discover MicroRNA-mediateding Dysregulated Pathways Involved in Disease”, 2017.

Issues?

Please email Chaoxing Li <chaoxing@asu.edu> or Valentin Dinu <Valentin.Dinu@asu.edu> if you have any questions, concerns or feedback.





# miR2Pathway
miR2Pathway is a PageRank-based method that can be used to rank disease risk of miRNA-mediated biological pathways. miR2Pathway can help explore how much miRNAs differentially influence the activity of biological pathways between two classes of phenotypes.

<b> For information about the method, please read: </b>

Chaoxing Li and Valentin Dinu. “miR2Pathway: A Novel Analytical Method to Discover MicroRNA-mediateding Dysregulated Pathways Involved in Disease”, 2017.

Main website: http://dinulab.org/tools/mir2pathway/

<b> R code, sample data and sample scripts: </b>

The R script to perform miR2Pathway: <samp>   <br> miR2Pathway.R</samp>

Example usage of miR2Pathway.R: <samp> <br>   miR2Pathway-example.R</samp>

Data used in miR2Pathway_Example_Dataset.RData: <samp> <br>   miR2Pathway_Example_Dataset.RData</samp>



# Usage

  
<b> Example Usage (to be used with example data) </b>  

library(compiler) <br>
library(foreach) <br>
library(iterators) <br>
library(parallel) <br>
library(doParallel) <br>
library(AnnotationDbi) <br>
library(stats4) <br>
library(BiocGenerics) <br>
library(miRNAtap)  <br>
library(miRNAtap.db)  <br>
library(graphite) <br>
library(BiocGenerics) <br>
library(graph) <br>
library(igraph) <br>

<b>If you wish, you may download the example datasset files individually and run the code below to pre-process them <br> </b>

<code>genelist <- read.table("genelist.txt", header = TRUE, sep = "\t")</code><br>
<code>miRlist <- read.table("miRlist.txt", header = TRUE, sep = "\t")</code><br>
<code>miRlist.full <- read.table("miRlist.full", header = TRUE, sep = "\t")</code><br>

<code>mydata.miR <- read.table("mydata.miR", header = TRUE, sep = "\t")</code><br>
<code>names.genelist <- read.table("names.genelist", header = TRUE, sep = "\t")</code><br>
<code>mydata.gene <- read.table("mydata.gene", header = TRUE, sep = "\t")</code><br>

<code>mydata.gene2 <- mydata.gene[,-1]</code><br>
<code>rownames(mydata.gene2) <- mydata.gene[,1]</code><br>

<code>mydata.miR2 <- mydata.miR[,-1]</code><br>
<code>rownames(mydata.miR2) <- mydata.miR[,1]</code><br>

<code>names.genelist2 <- names.genelist[,-1]</code><br>
<code>rownames(names.genelist2) <- names.genelist[,1]</code><br>
 
<b>Or you may simply load the workspace file: miR2Pathway_Example_Dataset.RData <br></b>

<code>results<-miR2Pathway(mydata.gene=mydata.gene2,mydata.miR=mydata.miR2,genelist=genelist,name.genelist=names.genelist2,miRlist=miRlist,miRlist.full=miRlist.full,N.miR=4,N.gene=4,N.path=4,Num.sample.normal=4,Num.sample.case=4,Pathway.database=humanKEGG,cor.cutoff=-0.4,N.parallel=4)</code> 
   
   

# Arguments

<code> miR2Pathway <- function(mydata.gene,mydata.miR,genelist,name.genelist,miRlist,miRlist.full,N.miR,N.gene,N.path,Num.sample.normal,Num.sample.case,Pathway.database,cor.cutoff,N.parallel) </code>

<html>

<body>

<table>
    <tr>
    <th>Object</th>
     <th>Description</th>
   </tr>
   <tr>
    <td>mydata.gene</td>
    <td>a gene expression dataset (matrix). Rows represent genes, and columns represent samples (from control to case).</td>
  </tr>
  <tr>
   <tr>
    <td>mydata.gene</td>
    <td>a gene expression dataset (matrix). Rows represent genes, and columns represent samples (from control to case).</td>
  </tr>
  <tr>
    <td>mydata.miR</td>
    <td>a miRNA expression dataset (matrix). Rows represent miRNAs, and columns represent samples (from control to case).</td>
  </tr>
  <tr>
    <td>genelist</td>
    <td>a list of gene names (e.g., gene symobol).</td>
  </tr>
  <tr>
    <td>name.genelist</td>
    <td>a matrix of gene expression dataset (The order of gene names must be consistent with genelist). This matrix is 1*N, where N is the number of genes. The value of each element could be a ramdom value. See an example in Examples Section. </td>
  </tr>
  <tr>
    <td>miRlist</td>
    <td>miRNA in a standard format</td>
  </tr>
  <tr>
    <td>miRlist.full</td>
    <td>the full name of all miRNAs.</td>
  </tr>
</table>

</body>
</html>

  

# Value

<html>
<body>

<table>
    <tr>
    <th>Object</th>
    <th>Description</th>
   </tr>
   <td>LengthOfPathway</td>
   <td>gene counts of pathways</td>
</table>

</body>
</html>      
       

# Examples

See miR2Pathway-example.R <br>
An example for the matrix of name.genelist:
<html>
<body>

<table>
   <tr>
    <th>100130426</th>
    <th>100133144</th>
    <th>100134869</th>
    <th>10357</th>
    <th>10431</th>
    <th>136542</th>
    <th>145363</th>
  </tr>
   <tr>
    <td>Sample1 </td>
    <td>0 </td>
      <td>0 </td>
      <td>0 </td>
      <td>0 </td>
      <td>0 </td>
      <td>0 </td>
  </tr>
</table>

</body>
</html>

 
# Credits

<b> Authors: </b> Chaoxing Li and Valentin Dinu

<b> If you use or modify the code, please cite: </b>

“miR2Pathway: A Novel Analytical Method to Discover MicroRNA-mediateding Dysregulated Pathways Involved in Disease”, 2017.

<b> Issues? </b>

Please email Chaoxing Li <chaoxing@asu.edu> or Valentin Dinu <Valentin.Dinu@asu.edu> if you have any questions, concerns or feedback.






library(clusterProfiler)
library(org.Hs.eg.db)

geneOntologyEnrichment_updated <- function(geneList, geneType="ENTREZ_GENE_ID", annotations="ALL", pvalCutOff=1){
    if(!is(geneList)[1]=='character'){stop('The class of geneList must be character')}
    if(!is(geneType)[1]=='character'){stop('The class of geneType must be character')}
    if(!is(annotations)[1]=='character'){stop('The class of ontologies must be character')}
    if(!is.numeric(pvalCutOff)){stop("The class of pvalCutOff parameter must be numeric.")}
    
    if ( geneType == 'ENTREZ_GENE_ID')
    {gene.type <- 'entrezgene_id'}
    else if (geneType == 'GENE_SYMBOL')
    {gene.type <- 'external_gene_name'}
    else 
    {gene.type <- tolower(geneType)}
    
    cat('Getting gene symbols...')
    genes.annotations <- getGenesAnnotation(geneList,attributes=c("external_gene_name","entrezgene_id"),filter=gene.type)
    
    cat('Retrieving Gene Ontology terms related to the list of DEGs...\n')
    gos.data <- enrichGO(genes.annotations$entrezgene_id, 'org.Hs.eg.db', ont=annotations, pvalueCutoff=pvalCutOff)
    final.gos <- data.frame(gos.data)
    
    return(final.gos)
  }
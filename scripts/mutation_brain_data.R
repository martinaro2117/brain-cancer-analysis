

# CARGA DE LIBRERÍAS
library(KnowSeq)
library(data.table)

# LECTURA DE FICHEROS NECESARIOS
source("scripts/geneOntologyEnrichment_updated.R")
source("scripts/knowseqReport_updated.R")
source("scripts/convert_to_counts.R")


######################################
# PREPARACIÓN DE DATOS DE ENTRADA
######################################

# Preparamos los datos de entrada
path_rna <- "brain_data/"
path_muts <- "mutation_data/"
path_sample_sheet <- "Downloaded_data/gdc_manifest_sample_sheet_brain.tsv" 
my_outpath <- "cancer_analysis_data/_count_files"


######################################
#  PROCESAMIENTO RNA-SEQ
######################################

samples <- convert_to_counts(down_folder = path_rna, 
                             sample_sheet_file = path_sample_sheet, 
                             outpath = my_outpath,
                             column_name = "unstranded", 
                             filter_lines = 4)

# FILTRADO: SOLO TUMORES PRIMARIOS
# SE FILTRA POR Tumor.Descriptor == "Primary" (Elimina Recurrence y Normal, realmente solamente recurrencia pq normal no hay)
valid_index <- which(samples$Tumor.Descriptor == "Primary")

# SE RECORTA LA INFORMACIÓN DE LAS MUESTRAS
samples <- samples[valid_index, ]

print(paste("Pacientes tras filtro (Solo Primarios):", nrow(samples)))


########################################
# PROCESAMIENTO MUTACIONES (IDH1 / TP53)
########################################

# SE BUSCAN TODOS LOS ARCHIVOS .maf.gz DESCARGADOS
maf_files <- list.files(path = path_muts,
                        pattern = "\\.maf\\.gz$", 
                        full.names = TRUE, recursive = TRUE)

# SE LEEN TODOS Y SE UNEN EN UNA TABLA
list_muts <- lapply(maf_files, function(x) data.table::fread(x, skip = "Hugo_Symbol", fill = TRUE))
mutations_df <- as.data.frame(data.table::rbindlist(list_muts, fill = TRUE))

# CREACIÓN DE LAS CLASES
# FUNCIÓN PARA LIMPIAR IDs (TCGA-XX-XXXX...)
get_patient_id <- function(x) substr(x, 1, 12)

# PACIENTES DISPONIBLES EN RNA 
rna_patients <- samples$Case.ID

# PACIENTES CON MUTACIONES (SEGÚN EL ARCHIVO DESCARGADO)
mut_idh1 <- unique(get_patient_id(mutations_df$Tumor_Sample_Barcode[mutations_df$Hugo_Symbol == "IDH1"]))
mut_tp53 <- unique(get_patient_id(mutations_df$Tumor_Sample_Barcode[mutations_df$Hugo_Symbol == "TP53"]))

# SE ASIGNA LA CLASE A CADA PACIENTE DE RNA
etiquetas <- vector("character", length(rna_patients))

for (i in 1:length(rna_patients)) {
  p <- rna_patients[i]
  
  has_IDH1 <- p %in% mut_idh1
  has_TP53 <- p %in% mut_tp53
  
  if (has_IDH1 & has_TP53)  etiquetas[i] <- "IDH1_TP53"
  else if (has_IDH1)        etiquetas[i] <- "IDH1_Only"
  else if (has_TP53)        etiquetas[i] <- "TP53_Only"
  else                      etiquetas[i] <- "No_IDH1_TP53"
}


#########################
#   RESULTADO FINAL
#########################

# Conteo de Clases Generadas
print(table(etiquetas))         # hay q usar SMOTE?

# GENERACIÓN DEL CSV Y LA MATRIZ FINAL
print("Generando data_info.csv y countsMatrix...")

# SE CREA EL CSV APUNTANDO A LOS ARCHIVOS LIMPIOS
data_info <- data.frame(
  Run   = samples$File.Name,  # Nombre
  Path  = rep(my_outpath, nrow(samples)),  # Ruta limpia repetida
  Class = etiquetas  # Clases calculadas
)
write.csv(x = data_info, file = "data_info.csv", row.names = FALSE)

######################################
# ANÁLISIS DATOS DE EXPRESIÓN
######################################

# Cargamos y aunamos los ficheros count
# SE GENERA LA MATRIZ DE CONTEOS DESDE EL CSV
countsInfo <- countsToMatrix("data_info.csv", extension = "")

# SE GUARDAN LOS DATOS "CRUDOS" PERO YA ORGANIZADOS
save(countsMatrix, etiquetas, samples, file = "cancer_analysis_data/datos_info_counts.RData")
load("cancer_analysis_data/datos_info_counts.RData")

# Exportamos tanto la matriz de datos como los labels a nuevas variables
countsMatrix <- countsInfo$countsMatrix
labels <- countsInfo$labels
table(labels)

######################################
#    ANOTACIÓN Y NORMALIZACIÓN
######################################

# SE OBTIENE LA ANOTACIÓN (NOMBRES REALES Y %GC) PARA LOS GENES
# 'rownames(countsMatrix)' SON LOS IDENTIFICADORES (Ensembl IDs) DE LAS FILAS
myAnnotation <- getGenesAnnotation(rownames(countsMatrix))

# SE CALCULAN LOS VALORES DE EXPRESIÓN NORMALIZADOS
# SE USA LA MATRIZ DE CONTEOS CRUDOS Y LA INFORMACIÓN QUE SE ACABA DE DESCARGAR
geneExprMatrix <- calculateGeneExpressionValues(countsMatrix, 
                                                annotation = myAnnotation)

# LIMPIEZA FINAL
# SE ELIMINAN LAS FILAS (GENES) QUE NO SE PUDIERON IDENTIFICAR (NAs)
geneExprMatrix <- geneExprMatrix[!is.na(rownames(geneExprMatrix)), ]
save(geneExprMatrix, file="cancer_analysis_data/geneExpr_brain.RData")

######################################################
# ANÁLISIS DE CALIDAD (QA) Y ELIMINACIÓN DE OUTLIERS
######################################################

# Realizamos el analisis de calidad
# Esta función busca outliers mediante PCA y distancias entre muestras
# toRemoval = TRUE indica que queremos que detecte los que se deben borrar
# toPNG/toPDF = FALSE evita que guarde las gráficas en disco 
QAResults <- RNAseqQA(geneExprMatrix, toRemoval = TRUE, toPNG=FALSE, toPDF=FALSE)

# SE SELECCIONAN LAS MUESTRAS QUE PASARON EL FILTRO (LA MATRIZ LIMPIA)
# Seleccionar muestras que pasen el filtro (no outliers)
# QAResults$matrix ya contiene los datos sin las muestras malas
qualityMatrix <- QAResults$matrix
qualityLabels <- labels[-which(colnames(geneExprMatrix) %in% QAResults$outliers)]

table(qualityLabels)

save(qualityMatrix, qualityLabels, file="cancer_analysis_data/QAResults_brain.RData")
#load("cancer_analysis_data/QAResults_brain.RData")

######################################
# ELIMINACIÓN DEL EFECTO BATCH (SVA)
######################################

# Creamos el modelo SVA de variables surrogadas para tratar el efecto batch
# ELIMINACIÓN DEL EFECTO BATCH (RUIDO TÉCNICO)
# Se usa el método SVA (Surrogate Variable Analysis) para detectar y corregir
# variaciones ocultas que no son biológicas
batchMatrix <- batchEffectRemoval(qualityMatrix, qualityLabels, method = "sva")

rownames(batchMatrix) <- make.names(rownames(batchMatrix))

save(batchMatrix, file="cancer_analysis_data/batchMatrix_brain.RData")
#load("cancer_analysis_data/batchMatrix_brain.RData")

########################################
# VISUALIZACIÓN DE LOS DATOS CORREGIDOS
########################################

# Se definen colores para los 4 grupos
colores <- c("cornflowerblue", "indianred1", "gold", "forestgreen")

# Se asigna un color a cada paciente según su grupo (qualityLabels)
colores_pacientes <- colores[as.factor(qualityLabels)]

# Se dibuja el Boxplot 
boxplot(batchMatrix,
        col = colores_pacientes,  
        outline = FALSE,        
        xaxt = 'n',               
        main = "Distribución Normalizada (Sin Efecto Batch)",
        ylab = "Nivel de Expresión (TPM)",
        xlab = paste("Pacientes (Total:", ncol(batchMatrix), ")")) 

# Se añade la leyenda para saber qué es cada color
legend("topright", 
       legend = levels(as.factor(qualityLabels)), 
       fill = colores, 
       cex = 0.8,      
       bty = "n") 

########################################
# VISUALIZACIÓN DE LOS DATOS CORREGIDOS
########################################

# Separar train-test
require(CORElearn)
set.seed(1)
p=5 # 5 particiones
foldIDx <- cvGenStratified(qualityLabels,p)

indexTrn <- foldIDx !=1
indexTest <- foldIDx ==1

XTRN <- batchMatrix[,indexTrn] # aquí se trabajan todavía con los genes
YTRN <- qualityLabels[indexTrn]
XTEST <- batchMatrix[,indexTest] # aquí se trabajan todavía con los genes
YTEST <- qualityLabels[indexTest]

# Extraemos los genes diferencialmente destacados teniendo en cuenta las correccion mediante SVA. DEG=Diferential Expressed Genes
DEGsInfo <- DEGsExtraction(XTRN, YTRN, lfc = 2, cov=2, pvalue = 0.001)

# Extraemos la tabla de estadisticas de los genes diferencialmente expresados, asi como la matriz ya filtrada con dichos genes. 
topTable <- DEGsInfo$DEG_Results$DEGs_Table
#DEGsInfo$DEG_Results$MulticlassLFC te da el LFC por par de gen
DEGsMatrix <- DEGsInfo$DEG_Results$DEGs_Matrix

# Heatmap de expresion de 12 primeros DEGs
dataPlot(DEGsMatrix[1:12,], YTRN, mode = "heatmap", toPNG=FALSE, toPDF=FALSE)

######################################
# PREDICCIÓN DE BIOMARCADORES CON ML
######################################

# Se preparan tanto la matriz como los labels
XTRN <- t(DEGsMatrix)
YTRN <- YTRN

XTEST <- t(XTEST[rownames(DEGsMatrix),])
YTEST <- YTEST

# Llevamos a cabo un proceso de Seleccion de Caracteristicas
FSRanking <- featureSelection(XTRN, YTRN, mode = "mrmr", vars_selected = colnames(XTRN))
dataPlot(t(XTRN[,FSRanking[1:6]]), YTRN, mode = "genesBoxplot", toPNG=FALSE, toPDF=FALSE)

# Evaluamos los biomarcadores mediante un proceso de validacion cruzada
knn_trn_res <- knn_trn(XTRN, YTRN, vars_selected = names(FSRanking)) 
knn_results_tr <- rbind(knn_trn_res$accuracyInfo$meanAccuracy, knn_trn_res$sensitivityInfo$meanSensitivity,knn_trn_res$specificityInfo$meanSpecificity)
dataPlot(knn_results_tr, YTEST, legend = c("Mean Accuracy","Mean Sensitivity","Mean Specificity"), mode = "classResults", xlab="# Genes", ylab="Prediction Score")

knn_test_res <- knn_test(XTRN, YTRN, XTEST, YTEST, vars_selected = names(FSRanking), knn_trn_res$bestK) 
knn_results_ts <- rbind(knn_test_res$accVector, knn_test_res$sensVector,knn_test_res$specVector)
dataPlot(knn_results_ts, YTEST, legend = c("Mean Accuracy","Mean Sensitivity","Mean Specificity"), mode = "classResults", xlab="# Genes", ylab="Prediction Score")

dataPlot(t(XTRN[,names(FSRanking[1:10])]), YTRN, mode = "heatmap")
dataPlot(knn_trn_res$cfMats[[3]]$table, YTRN, mode = "confusionMatrix")
dataPlot(t(XTRN[,names(FSRanking[1:10])]), YTRN, mode = "genesBoxplot")

dataPlot(t(XTEST[,names(FSRanking[1:10])]), YTEST, mode = "heatmap")
dataPlot(knn_test_res$cfMats[[3]]$table, YTEST, mode = "confusionMatrix") #knn_test_res$cfMats
dataPlot(t(XTEST[,names(FSRanking[1:10])]), YTEST, mode = "genesBoxplot")

######################################
# ENRIQUECIMIENTO FUNCIONAL
######################################

# Cogemos los identificadores ENTREZ
entrezAnnotation <- getGenesAnnotation(names(FSRanking[1:10]), attributes = c("external_gene_name","entrezgene_id"), filter = "external_gene_name")
entrezGeneIds<- entrezAnnotation$entrezgene_id[!is.na(entrezAnnotation$entrezgene_id)]

# Se descarga informacion sobre los Gene Ontology
GOs <- geneOntologyEnrichment_updated(as.character(entrezGeneIds), geneType = "ENTREZ_GENE_ID")

# Se descarga informacion sobre los Pathways
pathways <- DEGsToPathways(entrezAnnotation$external_gene_name)

# Se descarga informacion sobre las enfermedades relacionadas
diseases <- DEGsToDiseases(entrezAnnotation$external_gene_name, getEvidences = TRUE)

# Lanzamos el report automatico
knowseqReport_updated(geneExprMatrix,labels,'knowSeq-report',clasifAlgs=c('knn'), qualityAnalysis = F,
                      getDiseases=TRUE, geneOntology=TRUE, maxGenes = 12)
# TCGA Brain Cancer – Mutational Stratification and RNA-seq Classification

Este repositorio contiene un **pipeline completo en R** para el análisis de datos de **RNA-seq y mutaciones** de cáncer cerebral obtenidos del **The Cancer Genome Atlas (TCGA)**.  
El objetivo principal es **estratificar pacientes según mutaciones clave (IDH1 y TP53)** y aplicar **selección de características y clasificación supervisada** para la identificación de biomarcadores.

El flujo de trabajo incluye:
- Procesamiento de RNA-seq
- Integración de mutaciones (MAF)
- Normalización y control de calidad
- Eliminación de efecto batch (SVA)
- Selección de genes diferencialmente expresados
- Clasificación con Machine Learning
- Enriquecimiento funcional y reporte automático

---

## 📂 Estructura del repositorio
├── brain_data/ # Datos RNA-seq descargados desde GDC
├── mutation_data/ # Archivos .maf.gz con mutaciones
├── Downloaded_data/
│ └── gdc_manifest_sample_sheet_brain.tsv
├── cancer_analysis_data/
│ ├── _count_files/ # Conteos procesados
│ ├── geneExpr_brain.RData
│ ├── QAResults_brain.RData
│ └── batchMatrix_brain.RData
├── scripts/
│ ├── convert_to_counts.R
│ ├── geneOntologyEnrichment_updated.R
│ └── knowseqReport_updated.R
├── main_analysis.R # Script principal
└── README.md

---

## 🧬 Datos

### RNA-seq
- Datos de expresión génica de cáncer cerebral (TCGA)
- Descargados desde el **Genomic Data Commons (GDC)**

### Mutaciones
- Archivos **`.maf.gz`**
- Se utilizan las mutaciones en **IDH1** y **TP53**
- Se combinan múltiples archivos MAF en una única tabla

---

## ⬇️ Descarga de datos desde GDC

Para descargar los datos usando los archivos *manifest* de TCGA es **imprescindible** usar el cliente oficial:

### 🔧 Requisito
- **`gdc-client.exe`** (Windows) o `gdc-client` (Linux/Mac)

Disponible en:  
https://gdc.cancer.gov/access-data/gdc-data-transfer-tool

### Ejemplo de descarga
```bash
gdc-client.exe download -m gdc_manifest.txt
```

El sample sheet (gdc_manifest_sample_sheet_brain.tsv) se usa posteriormente para enlazar muestras con pacientes.

---

## Librerías necesarias
```r
library(KnowSeq)
library(data.table)
library(CORElearn)
```

> Nota: KnowSeq es fundamental para normalización, QA, selección de características, clasificación y enriquecimiento funcional.

--- 

## Flujo de análisis
### Procesamiento RNA-seq

- Conversión de archivos descargados a conteos
- Filtrado de muestras solo Tumores Primarios
- Generación de matriz de conteos

### Procesamiento de mutaciones

- Lectura y unión de archivos .maf.gz
- Extracción de IDs de pacientes TCGA
- Creación de clases:
  - IDH1_TP53
  - IDH1_Only
  - TP53_Only
  - No_IDH1_TP53

### Normalización y anotación

- Anotación génica (Ensembl → nombre real, %GC)
- Normalización a valores de expresión (TPM)
- Eliminación de genes no anotados

### Control de calidad (QA)

- Detección de outliers mediante PCA y distancias
- Eliminación automática de muestras problemáticas

### Eliminación de efecto batch

- Corrección mediante SVA (Surrogate Variable Analysis)
- Reducción de ruido técnico no biológico

### Selección de genes y Machine Learning

- División Train / Test estratificada
- Extracción de genes diferencialmente expresados (DEGs)
- Selección de características con mRMR
- Clasificación usando k-NN
- Evaluación con:
  - Accuracy
  - Sensitivity
  - Specificity
  - Confusion Matrix

### Enriquecimiento funcional

- Gene Ontology (GO)
- Pathways
- Enfermedades asociadas
- Reporte automático con KnowSeq

### Licencia
Este proyecto se distribuye bajo licencia MIT.
Los datos utilizados pertenecen al The Cancer Genome Atlas (TCGA) y están sujetos a sus condiciones de uso.

convert_to_counts <- function(down_folder, sample_sheet_file, outpath,
                              column_name = "unstranded", 
                              filter_lines = 0) {
  
  # Leemos sample sheet para obtener todas las muestras
  samples <- read.table(sample_sheet_file, sep="\t", header=TRUE)
  
  # Crear path de salida si no existe
  if (!file.exists(outpath)){
    dir.create(outpath)
  }
  
  # Recorremos cada muestra
  nsamples <- nrow(samples)
  for (i in 1:nsamples){
      print(paste0("Converting ", samples$File.ID[i]," ..."))
      # Identificamos el fichero TSV de STAR
      file <- paste0(down_folder, "/", samples$File.ID[i], "/", samples$File.Name[i])
      # Leemos el archivo tsv en un data frame
      df <- read.table(file = file, sep = "\t", header = TRUE)
      # Seleccionamos las columnas 1 y 4 del data frame
      df <- df[, c(1, which(colnames(df) == column_name))]
      # Tomamos solo las filas desde la 5 hasta la última
      df <- df[(filter_lines+1):nrow(df), ]
      # Guardamos el data frame en formato csv, aunque le ponemos extensión counts. Da igual, sería cuestión de adaptar el código de después.
      counts_file <- paste0(outpath, "/",
                            sub(".tsv", ".counts", 
                                samples$File.Name[i], fixed = TRUE))
      write.table(df, file = counts_file, sep ="\t", row.names = FALSE,
                  col.names =FALSE,quote = FALSE)
      # Actualizar sample sheet con el nuevo fichero
      samples$File.Name[i] <- sub(".tsv", ".counts", samples$File.Name[i], fixed = TRUE)
  }
  
  return(samples)
}
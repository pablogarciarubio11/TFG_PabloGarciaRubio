
#Seurat Command List https://satijalab.org/seurat/articles/essential_commands.html 
#EXECUTE THIS ALWAYS AT THE START OF EACH R SESSION
#The directories used were both coding and non-coding
setwd("D:/TFG/Pytables_coding")

.libPaths()

#1.- PACKAGES & LIBRARIES
{
  #===================================================================
  #This is commented because the necessary packages were already installed, but if they are not installed, uncomment this part and execute it
  #Execute from here if not installed ...
  
  #install.packages('Seurat', repos = c('https://satijalab.r-universe.dev', 'https://cloud.r-project.org'))
  # SeuratDisk is not currently available on CRAN. You can install it from GitHub with:
  # NO funciona para obtener el seurat-disk en esta nueva version de Seurat 5.0
  
  # SeuratDisk is not currently available on CRAN. You can install it from GitHub with:
  # if (!requireNamespace("remotes", quietly = TRUE)) {
  #   install.packages("remotes")
  # }
  # remotes::install_github("mojaveazure/seurat-disk")
  # 
  # install.packages("SeuratObject")
  # 
  # # SingleCellExperiment is not currently available on CRAN. You can install it from :
  # if (!require("BiocManager", quietly = TRUE))
  #   install.packages("BiocManager")
  # library(BiocManager)
  #... until here. 
  #=======================================================================
  
  
  
  
  #BiocManager::install("zellkonverter")
  
  # En R:
  #install.packages("remotes")
  #No es necesaria esta parte
  
  #remotes::install_github("cellgeni/sceasy") MIRAR LUEGO ESTA PARTE SI HAY ALGÚN ERROR
  
  #BiocManager::install("SingleCellExperiment")
  
  #BiocManager::install("glmGamPoi")
  
  #BiocManager::install("SingleR")
  
  #BiocManager::install("plger/scDblFinder")
  
  #BiocManager::install("AnnotationHub")
  
  #BiocManager::install("ensembldb")
  
  #BiocManager::install("BPCells")
  
  #BiocManager::install("scCustomize")
  
  #BiocManager::install('glmGamPoi')
  
  #BiocManager::install("batchelor")
  
  #BiocManager::install("ensembldb")
  
  #BiocManager::install("SingleR")
  
  # SETUP ERRORs
  # SI EN ALGÚN MOMENTO OCURRE ESTE ERROR:
  # Error: vector memory exhausted (limit reached?)
  # The problem with Rstudio, and identifies a solution, shown below:
  # Step 1: Open terminal,
  # Step 2:
  #  cd ~
  #  touch .Renviron
  # open .Renviron
  # Step 3: Save the following as the first line of .Renviron:
  #  R_MAX_VSIZE=100Gb
  # Note: This limit includes both physical and virtual memory; so setting _MAX_VSIZE=16Gb on a machine with 16Gb of physical memory
  # may not prevent this error. You may have to play with this parameter, depending on the specs of your machine.
  
  
  
  #=========================================================================
  #EMPEZAR A EJECUTAR DESDE AQUÍ EN CASO DE QUE SEA NECESARIO. LO ANTERIOR HASTA EL PRIMER COMENTARIO 
  #NO HACERLO
  #AÑADIDO:
  #SOLO HACER ESTO EN CASO DE NO TENER LAS LIBRERÍAS INSTALADAS
  # install.packages("tidyverse")
  # install.packages("patchwork")
  # install.packages("scales")
  # 
  # BiocManager::install("SingleCellExperiment", ask = FALSE, update = FALSE)
  # BiocManager::install("bluster", ask = FALSE, update = FALSE)
  # BiocManager::install("scDblFinder", ask = FALSE, update = FALSE)
  # BiocManager::install("zellkonverter", ask = FALSE, update = FALSE)
  #==============================================================================
}

#2.START SETUP 
  {
  # LOAD R PACKAGES 
  # ESTO SÍ HACE FALTA HACERLO, PARA COMPROBAR QUE NO HAY NINGÚN PROBLEMA
  library(tidyverse)
  library(dplyr)
  library(Matrix)
  library(Seurat)
  library(SeuratObject)
  library(scDblFinder)
  library(sctransform)
  library(SingleCellExperiment)
  library(ggplot2)
  library(patchwork)
  library(sctransform)
  library(scales)
  library(zellkonverter)
}


#2.- CREATE SEURAT OBJECT 
 #   H5 from CellBender
  #  2.1 PYTABLES!!!
    #{
      # setwd('/Users/rodrigo/Desktop/scRNAseq/Mmu_MicroRef/Pytables')
      # Seurat 4.0.2 uses a dataloader Read10X_h5() which is not currently compatible with the CellBender output file format.
      # Hopefully Seurat will update its dataloader to ignore extra information in the future, but in the interim, we can use 
      # a super handy utility from PyTables to strip the extra CellBender information out of the output file so that Seurat 
      # can load it.From a python environment in which PyTables is installed, do the following at the command line
      # PYTABLES
      # BASH-colocate en la carpeta donde vas a guardar los h5
      # ptrepack --complevel 5 CellBender_Append_MiR_ALL_CTR_Brennan_filtered.h5:/matrix Global_Brennan_.h5:/matrix
    #TIENE QUE TENER LA ESCTURTURA DE 
    #ref = Global ... GLobalPrimiR... (sin espacios o guión bajos)
    #samp = Tabula o Brennan o Hou
    #cond = GSMXXXX_Ctr
     # }
    
    #2.2  INDIVIDUALMENTE (Sin tener que dividir en CellBender)
  #   {
  #   # 10X CellRanger .HDF5 format
  #     GlobalprimiR_Tabula_GSM7474475_1dpi <- Read10X_h5(
  #     filename = "/Users/rodrigo/Desktop/scRNAseq/Mmu_MicroRef/Pytables/GlobalprimiR_Tabula_GSM7474475_1dpi.h5",
  #     use.names = TRUE,
  #     unique.features = TRUE)
  #   
  #     GlobalprimiR_Tabula_GSM7474475_1dpi[1:10, 1:10]
  #     dim(GlobalprimiR_Tabula_GSM7474475_1dpi)
  #   # 78988 50343
  #     GlobalprimiR_Tabula_GSM7474475_1dpi <- CreateSeuratObject(counts = GlobalprimiR_Tabula_GSM7474475_1dpi, assay = "RNA", min.cells = 3, min.features = 200)
  #   dim(GlobalprimiR_Tabula_GSM7474475_1dpi)
  #   
  #   
  #   
  #   Global_Squair_GSM5961588_Ctr <- CreateSeuratObject(counts = Count_Global)
  #   dim(Global_Squair_GSM5961588_Ctr)
  #   saveRDS(GlobalPrimiRtronLink_Brennan_GSM5904825_Ctr, file = "GlobalPrimiRtronLink_Brennan_GSM5904825_Ctr.rds") # SAVE
  #   saveRDS(Global_Brennan, file = "Global_Brennan_GSM5904825_Ctr.rds") 
  # }
    
    
    # 2.4 COLECTIVAMENTE si los archivos son pequeños y deben dividirse para ejecutar CellBender entonces ir arriba y despues run
    # EJECUTAR ESTA PARTE DEL CÓDIGO
    {
      # Script integrado para listar y procesar archivos .h5
      # Directorio objetivo
      directorio <- "D:/TFG/Pytables_noncoding"
     
      # Verificar si el directorio existe
      if (!dir.exists(directorio)) {
        stop("El directorio no existe: ", directorio)
      }
      
      # Obtener todos los archivos del directorio
      todos_archivos <- list.files(directorio, full.names = FALSE)
      cat(todos_archivos)
      
      # Filtrar solo archivos .h5
      archivos_h5 <- todos_archivos[grepl("\\.h5$", todos_archivos, ignore.case = TRUE)]
      cat(archivos_h5)
      
      # Verificar si se encontraron archivos .h5
      if (length(archivos_h5) == 0) {
        stop("No se encontraron archivos .h5 en el directorio.")
      }
      
      # Mostrar los archivos encontrados
      cat("Archivos .h5 encontrados:\n")
      lista_c <- paste0("c(\"", paste(archivos_h5, collapse = "\", \""), "\")")
      cat(lista_c, "\n")
      cat("Total de archivos .h5 encontrados:", length(archivos_h5), "\n\n")
      
      # Establecer el directorio de trabajo para los archivos
      setwd(directorio)
      
      # Crear lista empleando el script automatizado
      # Lista de archivos a procesar (ahora generada automáticamente)
      files <- archivos_h5
      
      # Lista vacía para almacenar los objetos Seurat
      seurat_objects <- list()
      
      # Iterar sobre los archivos y procesarlos
      for (file in files) {
        cat("Procesando archivo:", file, "\n")
        
        # Leer el archivo HDF5 con Read10X_h5
        tryCatch({
          adata <- Read10X_h5(filename = file, use.names = TRUE, unique.features = TRUE)
          
          # Crear nombre de variable limpio (sin extensión .h5 y caracteres especiales)
          var_name <- gsub("\\.h5$", "", file)
          var_name <- gsub("[^A-Za-z0-9_]", "_", var_name)  # Reemplazar caracteres especiales
          var_name <- paste0("seurat_", var_name)  # Prefijo para identificar fácilmente
          
          # Crear objeto Seurat
          seurat_obj <- CreateSeuratObject(
            counts = adata, 
            project = gsub("\\.h5$", "", file),
            min.cells = 3,
            min.features = 200 #ESTo ES UNA RESTRICCIÓN DEL NÚMERO DE CÉLULAS/EVENTOS QUE VAMOS A EVALUAR
            #IMPORTANTE COMENTAR ESTO
          )
          
          # Asignar el objeto Seurat al entorno global con nombre dinámico
          assign(var_name, seurat_obj, envir = .GlobalEnv)
          
          # También guardarlo en la lista para acceso programático
          seurat_objects[[var_name]] <- seurat_obj
          
          cat("Archivo", file, "procesado exitosamente como objeto:", var_name, "\n")
          cat("Dimensiones:", dim(seurat_obj), "(genes x células)\n")
          
        }, error = function(e) {
          cat("Error procesando", file, ":", e$message, "\n")
        })
      }
      
      cat("\nProcesamiento completado para", length(files), "archivos.\n")
      cat("\nObjetos Seurat creados en el entorno:\n")
      for (obj_name in names(seurat_objects)) {
        cat("-", obj_name, "\n")
      }
      
      # Opcional: Mostrar resumen de todos los objetos
      cat("\n=== RESUMEN DE OBJETOS SEURAT ===\n")
      for (obj_name in names(seurat_objects)) {
        obj <- get(obj_name, envir = .GlobalEnv)
        cat("Objeto:", obj_name, "\n")
        cat("  Proyecto:", obj@project.name, "\n")
        cat("  Dimensiones:", dim(obj)[1], "genes x", dim(obj)[2], "células\n")
        cat("  Archivo origen:", paste0(gsub("^seurat_", "", obj_name), ".h5"), "\n\n")
        #cat("  Archivo origen:", paste0(gsub("seurat_", "", gsub("_", ".", obj_name)), ".h5"), "\n\n")
      }
    }
  
# 3.- DOUBLETS scDblFinder INCLUIR EN EN METADATA UNA COLUMNA CON DOUBLETS_SCDBLF
#     3.1 Individual.
#     {
#     # Extract count matrix, prepare a sce object to scDblFinder and execute
#     # Extraer matriz de cuentas
#     countsGlobal_Squair <- Global_Squair_GSM5961588_Ctr[["RNA"]]$counts
#     
#     # Crear objeto SingleCellExperiment
#     sceGlobal_Squair <- SingleCellExperiment(assays = list(counts = countsGlobal_Squair))
#     
#     # Ejecutar scDblFinder
#     sceGlobal_Squair <- scDblFinder(sceGlobal_Squair)
#     
#     # Ver tabla de clasificación
#     table(sceGlobal_Squair$scDblFinder.class)
#     
#     # Transferir los resultados de scDblFinder al objeto Seurat original
#     Global_Squair_GSM5961588_Ctr$DoubletClass <- sceGlobal_Squair$scDblFinder.class
#     Global_Squair_GSM5961588_Ctr$DoubletScore <- sceGlobal_Squair$scDblFinder.score
#     
#     # Opcional: ver los metadatos con los nuevos campos
#     View(Global_Squair_GSM5961588_Ctr@meta.data)
#     
#     # Guardar el objeto Seurat enriquecido con clasificación doublet/singlet
#     saveRDS(GlobalPrimiRtronLink_Brennan_GSM5904825_Ctr, file = "GlobalPrimiRtronLink_Brennan_GSM5904825_Ctr.rds")
#     
#     # Global_Tabula_GSM7474457_Con <- readRDS("Global_Tabula_GSM7474457_Con.rds") # LOAD
#     
#    
#   NOOOOOOOOO  # REMOVE DOUBLETS
#     # Global_Tabula_GSM7474457_Con_SinDbl <- subset(Global_Tabula_GSM7474457_Con, subset = scDblFinder.class == "singlet")
#     # View(Global_Tabula_GSM7474457_Con_SinDbl@meta.data)
#     # dim(Global_Tabula_GSM7474457_Con_SinDbl)
#     # 
#     # saveRDS(Global_Tabula_GSM7474457_Con_SinDbl, file = "Global_Tabula_GSM7474457_ConSinDbl.rds") # SAVE
#     # 
#     # Global_Tabula_GSM7474457_Con_SinDbl<- readRDS("Global_Tabula_GSM7474457_ConSinDbl.rds") # LOAD
#     }    
#     
    # .2 ITERACION.
    {
      # Script iterativo para detección de doublets con scDblFinder
      
      # Función para identificar objetos Seurat en el entorno
      get_seurat_objects <- function() {
        all_objects <- ls(envir = .GlobalEnv)
        seurat_objects <- c()
        for (obj_name in all_objects) {
          obj <- get(obj_name, envir = .GlobalEnv)
          if (inherits(obj, "Seurat") && grepl("^seurat_GSM", obj_name)) { #El ^ marca dónde inicia el texto, el resto es el patrón a filtrar
            #Se usa grepl porque grep solo devuelve la posición del vector obj_name en el que el elemento coincide con el patrón
            #y también puedes pedir el texto directamente con value = TRUE.
            #Por otro lado, grepl devuelve un TRUE o FALSE si coincide o no
            seurat_objects <- c(seurat_objects, obj_name)
          }
        }
        return(seurat_objects)
      }
      
      
      # Obtener lista de objetos Seurat
      seurat_object_names <- get_seurat_objects()
      
      cat("Objetos Seurat encontrados en el entorno:\n")
      for (i in seq_along(seurat_object_names)) {
        cat(paste0(i, ". ", seurat_object_names[i], "\n"))
      }
      cat("\nTotal de objetos Seurat:", length(seurat_object_names), "\n\n")
      
      processed_objects <- c() #En este caso es mejor usar un vector en vez de una lista, ya que no voy a guardar datos complejos
      failed_objects <- c()
      
      for (obj_name in seurat_object_names) {
        cat("==========================================\n")
        cat("Procesando objeto:", obj_name, "\n")
        cat("==========================================\n")
        
        tryCatch({
          seurat_obj <- get(obj_name, envir = .GlobalEnv)
          
          # Verificar si el assay RNA tiene múltiples layers y unirlos si es necesario
          if ("RNA" %in% names(seurat_obj@assays)) {
            assay_obj <- seurat_obj[["RNA"]]
            if ("Assay5" %in% class(assay_obj) && length(Layers(assay_obj)) > 1) {
              cat(" → RNA tiene múltiples layers, aplicando JoinLayers...\n")
              seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")
            }
          }
          
          # Verificar que el objeto tenga datos RNA
          if (!"RNA" %in% names(seurat_obj@assays)) {
            warning("El objeto ", obj_name, " no tiene ensayo RNA. Saltando...")
            failed_objects <<- c(failed_objects, obj_name)
            next
          }
          
          # Extraer matriz de cuentas
          cat("Extrayendo matriz de cuentas...\n")
          counts_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
          
          if (ncol(counts_matrix) == 0 || nrow(counts_matrix) == 0) {
            warning("Matriz de cuentas vacía para ", obj_name, ". Saltando...")
            failed_objects <<- c(failed_objects, obj_name)
            next
          }
          
          cat("Dimensiones de la matriz:", dim(counts_matrix), "(genes x células)\n")
          
          sce_obj <- SingleCellExperiment(assays = list(counts = counts_matrix))
          cat("Ejecutando scDblFinder...\n")
          sce_obj <- scDblFinder(sce_obj, verbose = TRUE)
          
          doublet_table <- table(sce_obj$scDblFinder.class)
          cat("Tabla de clasificación de doublets:\n")
          print(doublet_table)
          
          cat("Transfiriendo resultados al objeto Seurat...\n")
          seurat_obj$DoubletClass <- sce_obj$scDblFinder.class
          seurat_obj$DoubletScore <- sce_obj$scDblFinder.score
          
          assign(obj_name, seurat_obj, envir = .GlobalEnv)
          processed_objects <<- c(processed_objects, obj_name) #Se pone así por el mismo motivo que lo de failed_objects
          
          cat("✓ Objeto", obj_name, "procesado exitosamente!\n")
          cat("  - Doublets detectados:", doublet_table["doublet"], "\n")
          cat("  - Singlets detectados:", doublet_table["singlet"], "\n\n")
          
        }, error = function(e) {
          cat("✗ Error procesando", obj_name, ":", e$message, "\n\n")
          failed_objects <<- c(failed_objects, obj_name) #Se usa <<- porque el objeto failed_objects se encuentra fuera del bucle error, por lo que tiene que buscar en entornos externos 
        })
      }
      
      cat("=== RESUMEN FINAL ===\n")
      cat("Total de objetos procesados exitosamente:", length(processed_objects), "\n")
      cat("Total de objetos fallidos:", length(failed_objects), "\n\n")
      
      if (length(processed_objects) > 0) {
        cat("Objetos procesados exitosamente:\n")
        for (i in seq_along(processed_objects)) {
          cat(paste0(i, ". ", processed_objects[i], "\n"))
        }
      }
      
      if (length(failed_objects) > 0) {
        cat("\nObjetos que fallaron:\n")
        for (i in seq_along(failed_objects)) {
          cat(paste0(i, ". ", failed_objects[i], "\n"))
        }
      }
      
      cat("\n=== LISTA DE OBJETOS SEURAT PROCESADOS ===\n")
      if (length(processed_objects) > 0) {
        seurat_objects_processed <- as.list(processed_objects)
        names(seurat_objects_processed) <- processed_objects
        cat("seurat_objects <- c(\n")
        for (i in seq_along(processed_objects)) {
          if (i == length(processed_objects)) {
            cat(paste0('  "', processed_objects[i], '"\n'))
          } else {
            cat(paste0('  "', processed_objects[i], '",\n'))
          }
        }
        cat(")\n\n")
        seurat_objects <- processed_objects
        assign("seurat_objects", seurat_objects, envir = .GlobalEnv)
        cat("Variable 'seurat_objects' creada en el entorno global.\n")
      } else {
        cat("No se procesaron objetos exitosamente.\n")
      }
      
      cat("\n¡Procesamiento completado!\n")
      
}

# 4.- EXPORT to PYTHON
#     4.1 Individual. Export to Python to create an ANNDATA OBJECT
    # {
    # Global_Tabula_GSM7474457_Con
    # Global_Tabula_GSM7474487_Cru
    # Global_Tabula_GSM7474501_Ctr 
      
     #  counts_matrix <- GetAssayData(Global_Squair_GSM5961588_Ctr, assay='RNA', layer='counts')
     #  dim(counts_matrix)  # array([18044.,  2728.])
     #  writeMM(counts_matrix, file=paste0(file='/Users/rodrigo/Desktop/LAB/matrix_files/matrix_seurat_con.mtx')) 
     #  # El mensaje NULL sugiere que el archivo se ha escrito correctamente en el directorio
     #  
     #  # write gene names
     #  write.table(data.frame('gene'=rownames(counts_matrix)), file='/Users/rodrigo/Desktop/LAB/matrix_files/gene_names_con.csv', quote=F,row.names=F,col.names=F)
     #  
     #  # En el caso de que emplemos el método de integración CCA o RCPA para obtener embeding layers...write dimensional reduction matrix (PCA)
     #  # write.csv (Objeto_Seurat@reductions$pca@cell.embeddings, file='pca.csv', quote=F, row.names=F)
     #  
     #  View(Global_Squair_GSM5961588_Ctr@meta.data) 
     #  dim(Global_Squair_GSM5961588_Ctr@meta.data)  # array([2728,    5], dtype=int32)
     #  # REVISAR    Eliminar la columna 'column_to_remove' del metadata si es necesario y esta repetida. Tambien se puede en python
     #  # Objeto_Seurat@meta.data$X <- NULL
     #  
     #  # save metadata table:
     #  Global_Squair_GSM5961588_Ctr$barcode <- colnames(Global_Squair_GSM5961588_Ctr)
     #  dim(counts_matrix)
     #  # Global_Squair_GSM5961588_Ctr$UMAP_1 <- Global_Squair_GSM5961588_Ctr@reductions$umap@cell.embeddings[,1]
     #  # Global_Squair_GSM5961588_Ctr$UMAP_2 <- Global_Squair_GSM5961588_Ctr@reductions$umap@cell.embeddings[,2]
     #  write.csv(Global_Squair_GSM5961588_Ctr@meta.data, file='/Users/rodrigo/Desktop/LAB/matrix_files/Global_Squair_GSM5961588_Ctr.csv', quote=F, row.names=F)
     # }
     
      
       
    # 4.2 ITERACION. Export to Python to create an ANNDATA OBJECTs
    # IMPORTANTE GUARDAR  el file_list.csv generado. Importante para el próximo SCRIPT.
    # {
    #   # Cargar las librerías necesarias
    #   library(Seurat)
    #   library(Matrix)
    #   
    #   # Usar la lista de objetos Seurat generada por el script anterior
    #   # (En lugar de definir manualmente la lista)
    #   if (!exists("seurat_objects")) {
    #     stop("La lista 'seurat_objects' no existe. Ejecuta primero el script de scDblFinder.")
    #   }
    #   
    #   # Directorio de salida
    #   output_dir <- "D:/TFG/Anndatas_coding"
    #   
    #   # Crear directorio si no existe
    #   if (!dir.exists(output_dir)) {
    #     dir.create(output_dir, recursive = TRUE)
    #   }
    #   
    #   # Lista para almacenar los nombres de los archivos generados
    #   file_records <- data.frame(
    #     object_name = character(),
    #     metadata_file = character(),
    #     gene_names_file = character(),
    #     stringsAsFactors = FALSE
    #   )
    #   
    #   # Iterar sobre cada objeto usando la lista seurat_objects
    #   for (obj_name in seurat_objects) {
    #     # Verificar si el objeto Seurat existe en el entorno
    #     if (!exists(obj_name)) {
    #       print(paste("El objeto", obj_name, "no existe en el entorno. Saltando..."))
    #       next
    #     }
    #     
    #     seurat_obj <- get(obj_name)
    #     
    #     # Obtener la matriz de conteos
    #     counts_matrix <- GetAssayData(seurat_obj, assay='RNA', layer='counts')
    #     
    #     # Definir nombres de archivos de salida
    #     matrix_file <- paste0(output_dir, "matrix_seurat_", obj_name, ".mtx")
    #     gene_file <- paste0(output_dir, "gene_names_", obj_name, ".csv")
    #     metadata_file <- paste0(output_dir, obj_name, ".csv")
    #     
    #     # Guardar la matriz en formato .mtx
    #     writeMM(counts_matrix, file=matrix_file)
    #     
    #     # Guardar los nombres de los genes
    #     write.table(
    #       data.frame('gene' = rownames(counts_matrix)), 
    #       file=gene_file,
    #       quote=F, row.names=F, col.names=F
    #     )
    #     
    #     # Guardar metadata
    #     seurat_obj$barcode <- colnames(seurat_obj)
    #     
    #     # Opcional: agregar coordenadas UMAP si existen
    #     if ("umap" %in% names(seurat_obj@reductions)) {
    #       seurat_obj$UMAP_1 <- seurat_obj@reductions$umap@cell.embeddings[,1]
    #       seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[,2]
    #     }
    #     
    #     # Guardar la metadata
    #     write.csv(
    #       seurat_obj@meta.data, 
    #       file=metadata_file, 
    #       quote=F, row.names=F
    #     )
    #     
    #     # Agregar los archivos generados a la lista
    #     file_records <- rbind(file_records, data.frame(
    #       object_name = obj_name,
    #       metadata_file = basename(metadata_file),
    #       gene_names_file = basename(gene_file)
    #     ))
    #     
    #     print(paste("Procesado:", obj_name))
    #   }
    #   
    #   # Guardar la lista de archivos generados en files_list.csv
    #   files_list_path <- paste0(output_dir, "files_list.csv")
    #   write.csv(file_records, files_list_path, row.names=FALSE, quote=FALSE)
    #   
    #   print(paste("Archivo de lista generado:", files_list_path))
    # }
      
      
      
    #Voy a usar mejor este porque  es más robusto y porque puede tratar con Seurat v5 con lo de los layers
   # 4.3 ITERACION Archivos grandes (batch_0+batch_1). Export to Python to create an ANNDATA OBJECTs
    {
    # Cargar las librerías necesarias
    library(Seurat)
    library(Matrix)
    
    # Usar la lista de objetos Seurat generada por el script anterior
    # (En lugar de definir manualmente la lista)
    if (!exists("seurat_objects")) {
      stop("La lista 'seurat_objects' no existe. Ejecuta primero el script de scDblFinder.")
    }
    
    # Directorio de salida
    output_dir <- "D:/TFG/PRUEBA_Anndatas_coding"
    
    # Crear directorio si no existe
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Lista para almacenar los nombres de los archivos generados
    file_records <- data.frame(
      object_name = character(),
      object_folder = character(),#Añado la subcarpeta a las columnas del dataframe
      matrix_file = character(),#Añado la matrix a las columnas
      metadata_file = character(),
      gene_names_file = character(),
      stringsAsFactors = FALSE
    )
    
    # Iterar sobre cada objeto usando la lista seurat_objects
    for (obj_name in seurat_objects) {
      # Verificar si el objeto Seurat existe en el entorno
      if (!exists(obj_name)) {
        print(paste("El objeto", obj_name, "no existe en el entorno. Saltando..."))
        next
      }
      
      seurat_obj <- get(obj_name)
      
      # Para Seurat v5: Unir capas si hay múltiples capas
      if (length(Layers(seurat_obj, assay = "RNA")) > 1) {
        print(paste("Uniendo capas para objeto:", obj_name))
        seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")
      }
      
      # Obtener la matriz de conteos
      counts_matrix <- GetAssayData(seurat_obj, assay='RNA', layer='counts')
      
      # Crea la subcarpeta para este objeto Seurat
      #Guardo la subcarpeta como GSM55... en vez de seurat_GSM55... para que sea todo más claro para python
      sample_name <- sub("^seurat_", "", obj_name)
      object_dir <- file.path(output_dir, sample_name)
      if (!dir.exists(object_dir)) {
        dir.create(object_dir, recursive = TRUE)
      }
      
      # Definir nombres de archivos de salida
      # matrix_file <- paste0(output_dir, "matrix_seurat_", obj_name, ".mtx")
      # gene_file <- paste0(output_dir, "gene_names_", obj_name, ".csv")
      # metadata_file <- paste0(output_dir, obj_name, ".csv")
      #Lo he cambiado porque le falta un separador / entre el output_dir y la subcarpeta y también entre 
      #la subcarpeta y el directorio del output donde va a ir cada archivo
      #File.path me junta file.path("D:/TFG/Anndatas_coding", "GSM5514787", "matrix.mtx") para dar
      #[1] "D:/TFG/Anndatas_coding/GSM5514787/matrix.mtx"
      matrix_file <- file.path(object_dir, "matrix.mtx")
      gene_file <- file.path(object_dir, "gene_names.csv")
      metadata_file <- file.path(object_dir, "metadata.csv")
      
      # Guardar la matriz en formato .mtx
      writeMM(counts_matrix, file=matrix_file)
      
      # Guardar los nombres de los genes
      write.table(
        data.frame('gene' = rownames(counts_matrix)), 
        file=gene_file,
        quote=FALSE, row.names=FALSE, col.names=FALSE
      )
      
      # Guardar metadata
      seurat_obj$barcode <- colnames(seurat_obj)
      
      # Opcional: agregar coordenadas UMAP si existen
      if ("umap" %in% names(seurat_obj@reductions)) {
        seurat_obj$UMAP_1 <- seurat_obj@reductions$umap@cell.embeddings[,1]
        seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[,2]
      }
      
      # Guardar la metadata
      write.csv(
        seurat_obj@meta.data, 
        file=metadata_file, 
        quote=FALSE, row.names=FALSE
      )
      
      # Agregar los archivos generados a la lista
      # file_records <- rbind(file_records, data.frame(
      #   object_name = obj_name,
      #   metadata_file = basename(metadata_file),
      #   gene_names_file = basename(gene_file)
      # ))
      #Debido a que he añadido 2 columnas nuevas, una con la subcarpeta y otra con la matrix, hay que 
      #cambiar la forma en la que se agregan los archivos generados a la lista
      file_records <- rbind(file_records, data.frame(
        object_name = obj_name,
        object_folder = sample_name,
        matrix_file = file.path(sample_name, "matrix.mtx"),
        metadata_file = file.path(sample_name, "metadata.csv"),
        gene_names_file = file.path(sample_name, "gene_names.csv"),
        stringsAsFactors = FALSE
      ))
      
      print(paste("Procesado:", obj_name))
    }
    
    # Guardar la lista de archivos generados en files_list.csv
    #files_list_path <- paste0(output_dir, "files_list.csv") 
    #Esto anterior está mal, ya que no pone el separador / entre el nombre del directorio
    #y el nombre del archivo csv
    files_list_path <- file.path(output_dir, "files_list.csv")
    write.csv(file_records, files_list_path, row.names=FALSE, quote=FALSE)
    
    print(paste("Archivo de lista generado:", files_list_path))
    }  
    

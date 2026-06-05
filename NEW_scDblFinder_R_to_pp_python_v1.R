
#Seurat Command List https://satijalab.org/seurat/articles/essential_commands.html 
#EXECUTE THIS ALWAYS AT THE START OF EACH R SESSION
#The directories used were both coding and non-coding
setwd("D:/TFG/Pytables_coding")

.libPaths()

#1.- PACKAGES & LIBRARIES
{
  #===================================================================
  #This is commented because the necessary packages were already installed, but if they are not installed, uncomment this part and execute it
  #Execute if not installed 
  
  #install.packages('Seurat', repos = c('https://satijalab.r-universe.dev', 'https://cloud.r-project.org'))
  
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
  # To verify that everything works properly.
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
    
    # 2.4 Collectively
    {
      # Script integrated to be able to list and process .h5 files
      # Target directory
      directorio <- "D:/TFG/Pytables_noncoding"
     
      # Verify if the directory exists
      if (!dir.exists(directorio)) {
        stop("The directory does not exist: ", directorio)
      }
      
      # Obtain all the files from the directory
      todos_archivos <- list.files(directorio, full.names = FALSE)
      cat(todos_archivos)
      
      # Filter to find only the .h5 files
      archivos_h5 <- todos_archivos[grepl("\\.h5$", todos_archivos, ignore.case = TRUE)]
      cat(archivos_h5)
      
      # Verify that .h5 files were found
      if (length(archivos_h5) == 0) {
        stop(".h5 files were not found in the directory.")
      }
      
      # Show the found files
      cat(".h5 files found:\n")
      lista_c <- paste0("c(\"", paste(archivos_h5, collapse = "\", \""), "\")")
      cat(lista_c, "\n")
      cat("Total number of .h5 files found:", length(archivos_h5), "\n\n")
      
      # Establish the work directory for the files
      setwd(directorio)
      
      # Create a list 
      # File list to be processed (automatically generated)
      files <- archivos_h5
      
      # Empty list to store the Seurat objects
      seurat_objects <- list()
      
      # Iterate over the files and process them
      for (file in files) {
        cat("Processing file:", file, "\n")
        
        # Read the HDF5 file with Read10x_h5
        tryCatch({
          adata <- Read10X_h5(filename = file, use.names = TRUE, unique.features = TRUE)
          
          # Create the clean variable name (without the .h5 extension nor any special characters)
          var_name <- gsub("\\.h5$", "", file)
          var_name <- gsub("[^A-Za-z0-9_]", "_", var_name)  # Replace the special characteres
          var_name <- paste0("seurat_", var_name)  # Preffix to identify the objects easily
          
          # Create the Seurat object
          seurat_obj <- CreateSeuratObject(
            counts = adata, 
            project = gsub("\\.h5$", "", file),
            min.cells = 3,
            min.features = 200 
            #This restricts the number of cells or events that we are going to evaluate
          )
          
          # Assign the Seurat object to the global environment with a dynamic name
          assign(var_name, seurat_obj, envir = .GlobalEnv)
          
          # Store the object in a list to be able to access it
          seurat_objects[[var_name]] <- seurat_obj
          
          cat("File", file, "processed successfully as the object:", var_name, "\n")
          cat("Dimensions:", dim(seurat_obj), "(genes x cells)\n")
          
        }, error = function(e) {
          cat("Processing error", file, ":", e$message, "\n")
        })
      }
      
      cat("\nProcessing completed for", length(files), "files.\n")
      cat("\nSeurat objects created in the environment:\n")
      for (obj_name in names(seurat_objects)) {
        cat("-", obj_name, "\n")
      }
      
      cat("\n=== SUMMARY OF THE SEURAT OBJECTS ===\n")
      for (obj_name in names(seurat_objects)) {
        obj <- get(obj_name, envir = .GlobalEnv)
        cat("Object:", obj_name, "\n")
        cat("  Project:", obj@project.name, "\n")
        cat("  Dimensions:", dim(obj)[1], "genes x", dim(obj)[2], "células\n")
        cat("  Origin file:", paste0(gsub("^seurat_", "", obj_name), ".h5"), "\n\n")
      }
    }
  
# 3.- DOUBLETS scDblFinder INCLUDE IN THE METADATA A COLUMN WITH DOUBLETS_SCDBLF
    # .2 ITERATION
    {
      # Iterative script for the detection of doublets with scDblFinder
      
      # Function that identifies the Seurat objects in the environment
      get_seurat_objects <- function() {
        all_objects <- ls(envir = .GlobalEnv)
        seurat_objects <- c()
        for (obj_name in all_objects) {
          obj <- get(obj_name, envir = .GlobalEnv)
          if (inherits(obj, "Seurat") && grepl("^seurat_GSM", obj_name)) { 
            seurat_objects <- c(seurat_objects, obj_name)
          }
        }
        return(seurat_objects)
      }
      
      
      # Obtain the list of Seurat objects
      seurat_object_names <- get_seurat_objects()
      
      cat("Seurat objects found in the environment:\n")
      for (i in seq_along(seurat_object_names)) {
        cat(paste0(i, ". ", seurat_object_names[i], "\n"))
      }
      cat("\nTotal of Seurat objects:", length(seurat_object_names), "\n\n")
      
      processed_objects <- c() 
      failed_objects <- c()
      
      for (obj_name in seurat_object_names) {
        cat("==========================================\n")
        cat("Processing object:", obj_name, "\n")
        cat("==========================================\n")
        
        tryCatch({
          seurat_obj <- get(obj_name, envir = .GlobalEnv)
          
          # Verify if the assay RNA has multiple layers and unite them if necessary
          if ("RNA" %in% names(seurat_obj@assays)) {
            assay_obj <- seurat_obj[["RNA"]]
            if ("Assay5" %in% class(assay_obj) && length(Layers(assay_obj)) > 1) {
              cat(" → RNA contains multiple layers, applying JoinLayers...\n")
              seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")
            }
          }
          
          # Verify that the object contains RNA data
          if (!"RNA" %in% names(seurat_obj@assays)) {
            warning("El objeto ", obj_name, " no tiene ensayo RNA. Saltando...")
            failed_objects <<- c(failed_objects, obj_name)
            next
          }
          
          # Extract counts matrix
          cat("Extracting counts matrix...\n")
          counts_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
          
          if (ncol(counts_matrix) == 0 || nrow(counts_matrix) == 0) {
            warning("Counts matrix empty for ", obj_name, ". Skipping...")
            failed_objects <<- c(failed_objects, obj_name)
            next
          }
          
          cat("Matrix dimensions:", dim(counts_matrix), "(genes x cells)\n")
          
          sce_obj <- SingleCellExperiment(assays = list(counts = counts_matrix))
          cat("Executing scDblFinder...\n")
          sce_obj <- scDblFinder(sce_obj, verbose = TRUE)
          
          doublet_table <- table(sce_obj$scDblFinder.class)
          cat("Doublet classification table:\n")
          print(doublet_table)
          
          cat("Transfering the results to the Seurat object...\n")
          seurat_obj$DoubletClass <- sce_obj$scDblFinder.class
          seurat_obj$DoubletScore <- sce_obj$scDblFinder.score
          
          assign(obj_name, seurat_obj, envir = .GlobalEnv)
          processed_objects <<- c(processed_objects, obj_name) 
          
          cat("✓ Object", obj_name, "processed successfully!\n")
          cat("  - Doublets detected:", doublet_table["doublet"], "\n")
          cat("  - Singlets detected:", doublet_table["singlet"], "\n\n")
          
        }, error = function(e) {
          cat("✗ Error processing", obj_name, ":", e$message, "\n\n")
          failed_objects <<- c(failed_objects, obj_name) #Se usa <<- porque el objeto failed_objects se encuentra fuera del bucle error, por lo que tiene que buscar en entornos externos 
        })
      }
      
      cat("=== FINAL SUMMARY ===\n")
      cat("Total of successfully processed objects:", length(processed_objects), "\n")
      cat("Total of failed objects:", length(failed_objects), "\n\n")
      
      if (length(processed_objects) > 0) {
        cat("Objects successfully processed:\n")
        for (i in seq_along(processed_objects)) {
          cat(paste0(i, ". ", processed_objects[i], "\n"))
        }
      }
      
      if (length(failed_objects) > 0) {
        cat("\nObjects that failed:\n")
        for (i in seq_along(failed_objects)) {
          cat(paste0(i, ". ", failed_objects[i], "\n"))
        }
      }
      
      cat("\n=== LIST CONTAINING THE SEURAT OBJECTS PROCESSED ===\n")
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
        cat("Variable 'seurat_objects' created in the global environment.\n")
      } else {
        cat("No objects were processed successfully.\n")
      }
      
      cat("\nProcessing completed!\n")
      
}

# 4.- EXPORT to PYTHON
   # 4.3 ITERATION for big files (batch_0+batch_1). Export to Python to create an ANNDATA OBJECTs
    {
    # Load the necessary libraries
    library(Seurat)
    library(Matrix)
    
    # Uses the Seurat objects list generated previously
    if (!exists("seurat_objects")) {
      stop("La lista 'seurat_objects' no existe. Ejecuta primero el script de scDblFinder.")
    }
    
    # Output directory
    output_dir <- "D:/TFG/PRUEBA_Anndatas_coding"
    
    # Create the directory if it does not exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # List to store the names of the files generated
    file_records <- data.frame(
      object_name = character(),
      object_folder = character(),
      matrix_file = character(),
      metadata_file = character(),
      gene_names_file = character(),
      stringsAsFactors = FALSE
    )
    
    # Iterate over every object using the list seurat_objects
    for (obj_name in seurat_objects) {
      # Verify if the Seurat object exists in the environment
      if (!exists(obj_name)) {
        print(paste("The object", obj_name, "does not exist in the environment. Skipping..."))
        next
      }
      
      seurat_obj <- get(obj_name)
      
      # For Seurat v5: Unite the layers if there are multiple layers
      if (length(Layers(seurat_obj, assay = "RNA")) > 1) {
        print(paste("Uniendo capas para objeto:", obj_name))
        seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")
      }
      
      # Obtain the counts matrix
      counts_matrix <- GetAssayData(seurat_obj, assay='RNA', layer='counts')
      
      # Creates the subfolder for this Seurat Object
      # Saves the subfolder as GSM55... instead of seurat_GSM55... to simplify the Python analysis
      sample_name <- sub("^seurat_", "", obj_name)
      object_dir <- file.path(output_dir, sample_name)
      if (!dir.exists(object_dir)) {
        dir.create(object_dir, recursive = TRUE)
      }
      
      # Defining the output files names 
      # For example: "D:/TFG/Anndatas_coding/GSM5514787/matrix.mtx"
      matrix_file <- file.path(object_dir, "matrix.mtx")
      gene_file <- file.path(object_dir, "gene_names.csv")
      metadata_file <- file.path(object_dir, "metadata.csv")
      
      # Save the matrix in .mtx format
      writeMM(counts_matrix, file=matrix_file)
      
      # Save the gene names
      write.table(
        data.frame('gene' = rownames(counts_matrix)), 
        file=gene_file,
        quote=FALSE, row.names=FALSE, col.names=FALSE
      )
      
      # Save metadata
      seurat_obj$barcode <- colnames(seurat_obj)
      
      # Add UMAP coordinates if existing
      if ("umap" %in% names(seurat_obj@reductions)) {
        seurat_obj$UMAP_1 <- seurat_obj@reductions$umap@cell.embeddings[,1]
        seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[,2]
      }
      
      # Save metadata
      write.csv(
        seurat_obj@meta.data, 
        file=metadata_file, 
        quote=FALSE, row.names=FALSE
      )
      
      # Add the generated files to the list
      file_records <- rbind(file_records, data.frame(
        object_name = obj_name,
        object_folder = sample_name,
        matrix_file = file.path(sample_name, "matrix.mtx"),
        metadata_file = file.path(sample_name, "metadata.csv"),
        gene_names_file = file.path(sample_name, "gene_names.csv"),
        stringsAsFactors = FALSE
      ))
      
      print(paste("Processed:", obj_name))
    }
    
    # Save the list of generated files in files_list.csv
    files_list_path <- file.path(output_dir, "files_list.csv")
    write.csv(file_records, files_list_path, row.names=FALSE, quote=FALSE)
    
    print(paste("List file generated:", files_list_path))
    }  
    

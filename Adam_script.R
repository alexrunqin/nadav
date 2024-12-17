suppressPackageStartupMessages({
  library(CATALYST)
  library(methods)
  library(flowCore)
  library(ggplot2)
  library(SingleCellExperiment)
  
  library(DT)
  
  library(dplyr)
  library(tidyr)
  library(broom)
  library(stringr)
  library(janitor)
  library(readxl)
  
  library(ggiraph)
  library(viridis)
})


#data_folder <- "IVoCC PABLO, Blyth 2022"
data_folder <- "/dskh/nobackup/alexq/nadav"

create_dt <- function(x){
  DT::datatable(x,
                extensions = 'Buttons',
                options = list(dom = 'Blfrtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                               lengthMenu = list(c(10,25,50,-1),
                                                 c(10,25,50,"All"))))
}








library(CytoML)
library(flowWorkspace)

flowWsDataPath <- file.path(data_folder, "data")
wsp_files <- list.files(flowWsDataPath, pattern = ".wsp")
fcs_folders <- file.path(flowWsDataPath, 
                         c("Debarcoded fcs1", "Debarcoded FCS2", 
                           "Debarcoded FCS3", "Debarcoded FCS4"))

#' Fixing fcs with SampleID channel
#' ----------------
# library(cytoqc)
# rawfiles <- list.files(fcs_folders[4], ".fcs", recursive = TRUE, full.names = TRUE)
# cqc_data <- cqc_load_fcs(rawfiles)
# 
# check_res <- cqc_check(cqc_data, type = "channel")
# check_res
# 
# match_res <- cqc_match(check_res, ref = 2)
# match_res
# 
# fcs_file <- read.FCS(
#     # "data/Debarcoded FCS4/Patients/CT_023_TP3_CD45_104.fcs",
#     "data/Debarcoded FCS4/HC229/HC229_CD45_154_CT_023_TP3.fcs",
#     # "data/MHA_IVoCC_002/FCS/MHA_002_TP4_CD45_154.fcs",
#     column.pattern = "SampleID", # Remove sample ID as only some have it
#     invert.pattern = TRUE
# )
# 
# write.FCS(
#     fcs_file,
#     # "IVoCC PABLO, Blyth 2022//data/Debarcoded FCS4/Patients/CT_023_TP3_CD45_104.fcs"
#     # "IVoCC PABLO, Blyth 2022//data/Debarcoded FCS4/HC229/HC229_CD45_154_CT_023_TP3.fcs"
#     "IVoCC PABLO, Blyth 2022//data/MHA_IVoCC_002/FCS/MHA_002_TP4_CD45_154.fcs"
# )


#' ======================
#' Get the workspace files and label cells by gating
#' ----------------------

priority_cell_paths = c()

getGatedCells <- function(ws_file,
                          ws_sample_group_name,
                          fcs_files_loc) {
  # Parse workspace file
  ws <- CytoML::open_flowjo_xml(ws_file)
  # Get gating set
  gs <- CytoML::flowjo_to_gatingset(ws, 
                                    name = ws_sample_group_name, 
                                    path = fcs_files_loc)
  # Get list of files
  sample_list <- sampleNames(gs)
  # # Vector of fcs file sizes
  # sample_num <- as.numeric(unlist(lapply(strsplit(sample_list, "_"), 
  #                                        function(x) x[length(x)])))
  # names(sample_num) <- unlist(lapply(strsplit(sample_list, "_"), 
  #                                    function(x) paste(x[-length(x)], collapse = "_")))
  # Get annotations for each cell
  cellTypes_list_gs1 <- list()
  cellTypes_df_list_gs1 <- list()
  
  pb <- txtProgressBar(min = 0, max = length(sample_list), style = 3)
  for(s in 1:length(sample_list)) { # iterate through each sample
    leaf_nodes <- gs_get_pop_paths(gs[[s]], order="bfs")
    cellTypes_list_sample <- list()
    for (ln in leaf_nodes) { # iterate through each node
      cellTypes <- character()
      cellTypes<- gh_pop_get_indices(gs[[s]], ln)
      cellTypes_list_sample[[ln]] <- cellTypes
    }
    cellTypes_list_sample <- do.call(cbind, cellTypes_list_sample)
    # Unprioritise cell subset annotations
    idx <- which(colnames(cellTypes_list_sample) %in% c("root", "/singlets", "/first part", "/CD34+ singlets"))
    # Prioritise cell subset annotations
    idx_p <- which(colnames(cellTypes_list_sample) %in% priority_cell_paths)
    idx <- c(idx, seq_len(ncol(cellTypes_list_sample))[-c(idx,idx_p)], idx_p)
    # Pick the furthest classified annotation
    cellTypes <- apply(cellTypes_list_sample[, idx], 1, function(x) {
      res <- names(which(x))
      res[length(res)]
      # return(res)
    })
    cellTypes_list_gs1[[sample_list[s]]] <- cellTypes
    cellTypes_df_list_gs1[[sample_list[s]]] <- cellTypes_list_sample[, idx]
    
    setTxtProgressBar(pb, s)
  }
  close(pb)
  
  # A list of all the cell labels, sorted by fcs file
  names(cellTypes_list_gs1) <- unlist(lapply(strsplit(names(cellTypes_list_gs1), "_"), 
                                             function(x) paste(x[-length(x)], collapse = "_")))
  names(cellTypes_df_list_gs1) <- names(cellTypes_list_gs1)
  
  cell_type_df <- lapply(seq_along(cellTypes_list_gs1), function(i) {
    data.frame(mg_cell_path=cellTypes_list_gs1[[i]], 
               cell_index=1:length(cellTypes_list_gs1[[i]]),
               fcs_name=names(cellTypes_list_gs1)[i])
  }) %>% 
    bind_rows() %>%
    mutate(wsp_name=unlist(lapply(strsplit(ws_file, "/"), function(x) x[length(x)])))
  
  cell_type_df_full <- lapply(seq_along(cellTypes_list_gs1), function(i) {
    data.frame(mg_cell_path=cellTypes_list_gs1[[i]], 
               cell_index=1:length(cellTypes_list_gs1[[i]]),
               fcs_name=names(cellTypes_list_gs1)[i]) %>% 
      bind_cols(
        cellTypes_df_list_gs1[[i]]
      ) 
  })%>%
    bind_rows() %>%
    mutate(wsp_name=unlist(lapply(strsplit(ws_file, "/"), function(x) x[length(x)])))
  
  return(list(cell_type_df=cell_type_df,
              cell_type_df_full=cell_type_df_full,
              leaf_nodes=leaf_nodes))
}


# Get all leaf nodes and put into excel
mg_list_mha <- getGatedCells(
  ws_file = file.path("/dskh/nobackup/alexq/nadav/mha_data/MHA_IVoCC_002/BG_20210618_CB_NM.wsp"),
  fcs_files_loc = file.path("/dskh/nobackup/alexq/nadav/mha_data/MHA_IVoCC_002/FCS"),
  ws_sample_group_name = "All Samples"
)

mg_list <- lapply(1:4, 
                  function(x) getGatedCells(
                    ws_file = file.path(flowWsDataPath, wsp_files[x]),
                    fcs_files_loc = fcs_folders[x],
                    ws_sample_group_name = "All Samples"
                  )
)

saveRDS(
  mg_list,
  file=file.path(data_folder, "data/cell_type_mg_list_b4.RDS")
)

# saveRDS(
#   mg_list_mha,
#   file=file.path(data_folder, "data/cell_type_mg_list_mha.RDS")
# )




#' ==================================================
#' Read in data
#' --------------------------------------------------
fcsFiles <- list.files(file.path(data_folder, "data"),
                       pattern=".fcs",
                       recursive = TRUE,
                       full.names = TRUE) %>%
  .[!grepl("Don't USE",.)]

fcs_raw <- read.flowSet(fcsFiles, 
                        transformation = FALSE, 
                        truncate_max_range = FALSE)
# tmp <- read.flowSet(fcsFiles[161:167], 
#                         transformation = FALSE, 
#                         truncate_max_range = FALSE)
#,
# column.pattern = "SampleID", # Remove sample ID as only some have it
# invert.pattern = TRUE)

# 
# x <- colnames(read.flowSet(fcsFiles[160], 
#                         transformation = FALSE, 
#                         truncate_max_range = FALSE))
# 
# y <- colnames(read.flowSet(fcsFiles[161], 
#                         transformation = FALSE, 
#                         truncate_max_range = FALSE))


### Get marker information ----------------------------------
# pregating_channels <- c("Bead", "DNA1", "DNA2", "Dead", "Cell_length")
lineage_channels <- c("CD11c","CD56",
                      "CD8a",#"CD8A",
                      "CD57","CD27","CD19","CD45RA",
                      "CD4","KLRG1","CD45RO","CD31",
                      "CD194_CCR4",#"CD194 (CCR4)",
                      "CD197_CCR7",#"CD197 (CCR7)",
                      "CD14","APC","CD16","IgD","Ki67",
                      "CD25","CD3","CD38",
                      "IntB7",#"integrin beta7",
                      "CD127")

# Base the panel off the batch 1 .fcs
fcs_panel <- pData(parameters(fcs_raw[[which(fcsFiles == "/dskh/nobackup/alexq/nadav/data/Debarcoded fcs1/HC229/HC299_CHW_002_d180_CD45_104.fcs")]])) %>%
  select(name, desc) %>%
  mutate(marker = gsub("^[^_]+_", "", desc)) %>%
  mutate(marker_class = ifelse(marker %in% lineage_channels,
                               "type", ifelse(is.na(marker),
                                              "none",
                                              "state")))

### Get sample information ----------------------------------
### Read in flow data
samp_df <- data.frame(file_name_full=list.files(file.path(data_folder, "data"),
                                                pattern=".fcs",
                                                recursive = TRUE,
                                                full.names = FALSE)) %>%
  filter(!grepl("Don't USE",file_name_full)) %>%
  mutate(file_name_full = gsub("Decoded fcs\\/|Debarcoded FCS2\\/|Debarcoded FCS3\\/|Debarcoded FCS4\\/", 
                               "", file_name_full)) %>%
  mutate(sample_type = gsub( "\\/.*$", "", file_name_full ),
         file_name = gsub( ".*\\/", "", file_name_full )) %>%
  rowwise() %>% 
  mutate(pat_id = ifelse(sample_type == "HC229",
                         paste0(strsplit(file_name, split="_")[[1]][2:3], collapse="_"),
                         paste0(strsplit(file_name, split="_")[[1]][1:2], collapse="_"))) %>%
  ungroup()

saveRDS(
  lapply(mg_list, function(dat) {
    dat$cell_type_df_full %>%
      left_join(samp_df,
                by=c("fcs_name"="file_name"))
  }) %>%
    bind_rows(),
  file=file.path(data_folder, "data/cell_type_full_df_b4.RDS")
)

all_cells_df <- readRDS(file.path(data_folder, "data/cell_type_full_df_b4.RDS"))

# Convert to SingleCellExperiment
sce <- prepData(x=fcs_raw, 
                panel=fcs_panel,
                md=samp_df,
                cofactor = 5,
                panel_cols = list(channel="name",
                                  antigen="desc"),
                md_cols = list(file="file_name",
                               id="pat_id",
                               factors=c("sample_type", "file_name"))
)

# Edit colData
colData <- colData(sce)
colData$file_name <- as.character(colData$file_name)

### Join wsp cell types ----------------------------------
colData$mg_cell_path <- colData %>%
  as.data.frame() %>%
  group_by(file_name) %>%
  mutate(cell_index=1:n()) %>%
  left_join(all_cells_df %>%
              select(cell_index, fcs_name, mg_cell_path),
            by=c("cell_index"="cell_index",
                 "file_name"="fcs_name")) %>%
  dplyr::pull(mg_cell_path)

# Add batch information, taken from the workspace name
colData$batch <- colData %>%
  as.data.frame() %>%
  group_by(file_name) %>%
  mutate(cell_index=1:n()) %>%
  left_join(all_cells_df %>%
              select(cell_index, fcs_name, wsp_name),
            by=c("cell_index"="cell_index",
                 "file_name"="fcs_name")) %>%
  pull(wsp_name) %>%
  readr::parse_number()

# Add indexes for cells
colData$cell_index <- colData %>%
  as.data.frame() %>%
  group_by(file_name) %>%
  mutate(cell_index=1:n()) %>%
  pull(cell_index)

colData(sce) <- colData

# Switch the values for "156Gd_PE" "175Lu_CD279_PD1" for batch 2
sce_exprs_new = assay(sce, "exprs")
sce_exprs_new[c(which(rownames(sce)=="156Gd_CD279_PD1"), which(rownames(sce)=="175Lu_PE")), which(colData$batch==2)] <-
  sce_exprs_new[c(which(rownames(sce)=="175Lu_PE"), which(rownames(sce)=="156Gd_CD279_PD1")), which(colData$batch==2)]

assay(sce, "exprs") <- sce_exprs_new

saveRDS(sce, 
        file=file.path(data_folder, "data/sce_b4.rds"))

#' ==================================================
#' Read in data - MHA
#' --------------------------------------------------
fcsFiles <- list.files(file.path("/dskh/nobackup/alexq/nadav/mha_data/MHA_IVoCC_002"),
                       pattern=".fcs",
                       recursive = TRUE,
                       full.names = TRUE) %>%
  .[!grepl("Don't USE",.)]

fcs_raw <- read.flowSet(fcsFiles, 
                        transformation = FALSE, 
                        truncate_max_range = FALSE)


### Get marker information ----------------------------------
lineage_channels <- c("CD11c","CD56",
                      "CD8a",#"CD8A",
                      "CD57","CD27","CD19","CD45RA",
                      "CD4","KLRG1","CD45RO","CD31",
                      "CD194_CCR4",#"CD194 (CCR4)",
                      "CD197_CCR7",#"CD197 (CCR7)",
                      "CD14","APC","CD16","IgD","Ki67",
                      "CD25","CD3","CD38",
                      "IntB7",#"integrin beta7",
                      "CD127")

# Base the panel off the batch 1 .fcs
fcs_panel <- pData(parameters(fcs_raw[[which(fcsFiles == "/dskh/nobackup/alexq/nadav/mha_data/MHA_IVoCC_002/FCS/MHA_002_BMT_CD45_154.fcs")]])) %>%
  select(name, desc) %>%
  mutate(marker = gsub("^[^_]+_", "", desc)) %>%
  mutate(marker_class = ifelse(marker %in% lineage_channels,
                               "type", ifelse(is.na(marker),
                                              "none",
                                              "state")))

### Get sample information ----------------------------------
### Read in flow data
samp_df <- data.frame(file_name_full=list.files(file.path("/dskh/nobackup/alexq/nadav/mha_data/MHA_IVoCC_002"),
                                                pattern=".fcs",
                                                recursive = TRUE,
                                                full.names = FALSE)) %>%
  filter(!grepl("Don't USE",file_name_full)) %>%
  mutate(file_name_full = gsub("FCS\\/", 
                               "", file_name_full)) %>%
  mutate(sample_type = "Samples",
         file_name = gsub( ".*\\/", "", file_name_full )) %>%
  rowwise() %>% 
  mutate(pat_id = ifelse(sample_type == "HC229",
                         paste0(strsplit(file_name, split="_")[[1]][2:3], collapse="_"),
                         paste0(strsplit(file_name, split="_")[[1]][1:2], collapse="_"))) %>%
  ungroup()

saveRDS(
  mg_list_mha$cell_type_df_full %>%
    left_join(samp_df,
              by=c("fcs_name"="file_name")),
  file=file.path(data_folder, "data/cell_type_full_df_mha.RDS")
)

all_cells_df_mha <- readRDS(file.path(data_folder, "data/cell_type_full_df_mha.RDS"))

# Convert to SingleCellExperiment
sce_mha <- prepData(x=fcs_raw, 
                    panel=fcs_panel,
                    md=samp_df,
                    cofactor = 5,
                    panel_cols = list(channel="name",
                                      antigen="desc"),
                    md_cols = list(file="file_name",
                                   id="pat_id",
                                   factors=c("sample_type", "file_name"))
)

# Edit colData
colData <- colData(sce_mha)
colData$file_name <- as.character(colData$file_name)

### Join wsp cell types ----------------------------------
colData$mg_cell_path <- colData %>%
  as.data.frame() %>%
  group_by(file_name) %>%
  mutate(cell_index=1:n()) %>%
  left_join(all_cells_df_mha %>%
              select(cell_index, fcs_name, mg_cell_path),
            by=c("cell_index"="cell_index",
                 "file_name"="fcs_name")) %>%
  dplyr::pull(mg_cell_path)

# Add batch information, taken from the workspace name
colData$batch <- colData %>%
  as.data.frame() %>%
  group_by(file_name) %>%
  mutate(cell_index=1:n()) %>%
  left_join(all_cells_df_mha %>%
              select(cell_index, fcs_name, wsp_name),
            by=c("cell_index"="cell_index",
                 "file_name"="fcs_name")) %>%
  pull(wsp_name) %>%
  readr::parse_number()

# Add indexes for cells
colData$cell_index <- colData %>%
  as.data.frame() %>%
  group_by(file_name) %>%
  mutate(cell_index=1:n()) %>%
  pull(cell_index)

colData(sce_mha) <- colData

saveRDS(sce_mha, 
        file=file.path(data_folder, "data/sce_mha.rds"))


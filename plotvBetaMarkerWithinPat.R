plotvBetaMarkerWithinPat <- function(pat_id_val, all_pats=FALSE) {
  if (!all_pats) {
    vbeta_pops <- nadav_tcrvb_filt %>%
      # filter(sample_id %in% pat_id_val) %>%
      filter(fcs_name %in% pat_id_val) %>%
      pull(name)
    
  } else {
    vbeta_pops <- unique(nadav_tcrvb_filt$name)
  }
  
  vbeta_pop_paths <- colnames(all_cells_df)[
    grepl(paste0(paste0(vbeta_pops, "$"), collapse="|"), 
          colnames(all_cells_df)) &
      grepl("CD3\\+\\/non-NK", 
            colnames(all_cells_df))]
  
  # Subset 
  cd3_cell_ind <- all_cells_df %>%
    # filter patients with a v beta available
    {if(all_pats) . else filter(., fcs_name %in% pat_id_val)} %>%
    # filter CD3+ only
    filter(`/CD3+/non-NK T cells`) %>%
    # filter samples only %>%
    filter(sample_type %in% c("Samples", "Patients")) %>% 
    # get dpt from sample file 
    select(-dpt) %>%
    left_join(
      nadav_tcrvb_filt %>%
        {if(all_pats) . else filter(., fcs_name %in% pat_id_val)} %>%
        distinct(fcs_name, dpt, batch),
      by="fcs_name"
    ) %>%
    # select_if(function(col) sum(is.na(col)) == 0 ) %>%
    select(any_of(vbeta_pop_paths), mg_cell_path, cell_index, fcs_name,sample_type,dpt, pat_id) %>%
    left_join(
      colData(sce) %>% 
        as.data.frame() %>%
        mutate(sce_ind = 1:n()) %>%
        select(#sample_id, #dpt, 
          file_name, cell_index, sce_ind) ,
      by=c(#"pat_id"="sample_id",
        # "dpt"="dpt",
        "fcs_name"="file_name",
        "cell_index"="cell_index")
    ) %>%
    select_if(function(col) any(!is.na(col)))
  
  trbv_cols <- cd3_cell_ind %>%
    select(any_of(vbeta_pop_paths)) %>%
    colnames()
  
  pat_cell_df <- cd3_cell_ind %>%
    bind_cols(
      t(assay(sce, "exprs")) %>%
        as.data.frame() %>%
        select(all_of(nadav_panel_marker$`FCS name`)) %>%
        `colnames<-`(nadav_panel_marker$name) %>%
        dplyr::slice(cd3_cell_ind$sce_ind)
    ) %>%
    pivot_longer(cols=all_of(trbv_cols),
                 names_to="pop",
                 values_to="value") %>%
    filter(!is.na(value)) %>%
    mutate(value = ifelse(value, "vbeta", "rest")) %>%
    select(-c(mg_cell_path, sample_type)) 
  
  if (all_pats) {
    pop_comp_df <- pat_cell_df %>% 
      distinct(pop)
    
    pat_n_cells_df <- pat_cell_df %>% 
      group_by(pop) %>% 
      summarise(n_cells=n(), 
                n_vbeta = sum(value == "vbeta"),
                .groups='drop')
  } else {
    pop_comp_df <- pat_cell_df %>% 
      distinct(pop, fcs_name, dpt, pat_id)  %>%
      mutate(vbeta_name = sub(".*/", "", pop)) %>%
      # Filter samples with NA dpt, since they don't have a vbeta listed on Nadavs sheet
      filter(!is.na(dpt)) %>%
      # Need to filter only vbeta patient combinations
      filter(
        grepl(
          #legal vbeta name combinations
          nadav_tcrvb_filt %>% 
            mutate(match_str = paste(name,sample_id,dpt,sep="_")) %>% 
            pull(match_str) %>% 
            paste(., collapse="|"),
          paste(pop,pat_id,dpt,sep="_")
        )
      )
    
    pat_n_cells_df <- pat_cell_df %>% 
      group_by(fcs_name, dpt, pat_id, pop) %>% 
      summarise(n_cells=n(), 
                n_vbeta = sum(value == "vbeta"),
                .groups='drop')
  }
  
  num_markers = length(nadav_panel_marker$name)
  num_pops = nrow(pop_comp_df)
  
  MAGpop = matrix(nrow=num_pops,ncol=num_markers)
  MAGref = matrix(nrow=num_pops,ncol=num_markers)
  IQRpop = matrix(nrow=num_pops,ncol=num_markers)
  IQRref = matrix(nrow=num_pops,ncol=num_markers)
  
  for (i in seq_len(nrow(pop_comp_df))) {
    if (all_pats) {
      filt_dat <- pat_cell_df %>% 
        filter(pop == pop_comp_df$pop[i])
    } else {
      filt_dat <- pat_cell_df %>% 
        filter(fcs_name == pop_comp_df$fcs_name[i] & 
                 dpt  == pop_comp_df$dpt[i] &
                 pat_id == pop_comp_df$pat_id[i] & 
                 pop == pop_comp_df$pop[i])
    }
    MAGpop[i,] = filt_dat %>% 
      filter(value == "vbeta") %>% 
      select(all_of(nadav_panel_marker$name)) %>% 
      summarise_all(function(col) median(col, na.rm=TRUE)) %>%
      t() %>% t()
    
    IQRpop[i,] = filt_dat %>% 
      filter(value == "vbeta") %>% 
      select(all_of(nadav_panel_marker$name)) %>% 
      summarise_all(function(col) IQR(col, na.rm=TRUE)) %>%
      t() %>% t()
    
    MAGref[i,] = filt_dat %>% 
      filter(value == "rest") %>% 
      select(all_of(nadav_panel_marker$name)) %>% 
      summarise_all(function(col) median(col, na.rm=TRUE)) %>%
      t() %>% t()
    IQRref[i,] = filt_dat %>% 
      filter(value == "rest") %>% 
      select(all_of(nadav_panel_marker$name)) %>% 
      summarise_all(function(col) IQR(col, na.rm=TRUE)) %>%
      t() %>% t()
  }
  
  # Remove the line where there is NA, probably due to not filtering the 
  # vbetas for the particular days post
  remove_ind = which(is.na(rowSums(IQRpop)))
  if (!identical(remove_ind, integer(0))) {
    pop_comp_df <- pop_comp_df %>% slice(-remove_ind)
    MAGpop=MAGpop[-remove_ind,]
    IQRpop=IQRpop[-remove_ind,]
    MAGref=MAGref[-remove_ind,]
    IQRref=IQRref[-remove_ind,]
  }
  
  
  IQR_thresh <- IQR.thresh(MAGpop,MAGref,IQRpop,IQRref,num_markers) # nolint
  
  for (i in seq_len(num_markers)){
    IQRpop[,i] = pmax(IQRpop[,i],IQR_thresh)
    IQRref[,i] = pmax(IQRref[,i],IQR_thresh)
  }
  
  # Calculate MEM scores
  MAG_diff = MAGpop-MAGref
  MEM_matrix = abs(MAGpop-MAGref)+(IQRref/IQRpop)-1
  MEM_matrix[!(MAG_diff>=0)] <- (-MEM_matrix[!(MAG_diff>=0)])
  
  # Put MEM values on -10 to +10 scale
  scale_max = max(abs(MEM_matrix[,c(1:ncol(MEM_matrix)-1)]))
  if (nrow(MEM_matrix) == 1) {
    MEM_matrix = c((MEM_matrix[,c(1:ncol(MEM_matrix)-1)]/scale_max)*10,MEM_matrix[,ncol(MEM_matrix)]) %>%
      matrix(nrow=1)
  } else {
    MEM_matrix = cbind((MEM_matrix[,c(1:ncol(MEM_matrix)-1)]/scale_max)*10,MEM_matrix[,ncol(MEM_matrix)])
  }
  
  colnames(MEM_matrix) <- nadav_panel_marker$name
  colnames(MAGpop) <- nadav_panel_marker$name
  
  pop_mem_df <- pop_comp_df %>%
    bind_cols(MEM_matrix)%>%
    {if (all_pats) left_join(
      ., 
      pat_n_cells_df,
      by=c("pop"="pop")) else left_join(
        .,
        pat_n_cells_df,
        by=c("fcs_name"="fcs_name",
             "dpt"="dpt",
             "pat_id"="pat_id",
             "pop"="pop")
      )}
  
  pop_median_df <- pop_comp_df %>%
    bind_cols(MAGpop) %>%
    {if (all_pats) left_join(
      ., 
      pat_n_cells_df,
      by=c("pop"="pop")) else left_join(
        .,
        pat_n_cells_df,
        by=c("fcs_name"="fcs_name",
             "dpt"="dpt",
             "pat_id"="pat_id",
             "pop"="pop")
      )}
  
  return(list(pop_mem_df = pop_mem_df,
              pop_median_df = pop_median_df))
}
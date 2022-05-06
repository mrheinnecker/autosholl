
library(tiff)
library(tidyverse)
library(gridExtra)
library(cowplot)


source("/Users/Marco/git_repos/autosholl/fncts.R")

#file <- "/Users/Marco/Dropbox/Studium/Master/Praktikum_Mueller/example_data/Tier1_2_3_apical_D.tif"
files <- list.files("f:/data_sholl_analysis/single_soma", pattern = "*", full.names = T)

#file <- "f:/data_sholl_analysis/Tiffs/Apical/Deep/Tier4_5_2_apical_D.tif"

file <- files[1]
soma_xy_detection_cube_radius <- 15
soma_z_detection_radius <- 100
soma_z_detection_degree_steps <- 0.05
  
  
results <- lapply(files, function(file){
  start_time <- Sys.time()
  print(file)

  full_image <- readTIFF(file, all=T)
  nr_orig <- nrow(full_image[[1]])
  nc_orig <- ncol(full_image[[1]])
  somata <- get_somata(soma_xy_detection_cube_radius*2+1, 
                       soma_z_detection_radius, 
                       soma_z_detection_degree_steps, 
                       full_image)
  
  print("soma detected")
  SOMA <- somata[1,]
  
  apply(somata,1, function(SOMA){
    ## check soma for dendritic start sites
    main_dendrites_raw <- find_dendritic_start_sites(SOMA) 
    main_dendrites <- main_dendrites_raw[[1]] %>%
      mutate(z=SOMA[["z"]])  %>%
      arrange(h_angle)%>%
      rownames_to_column("dendrite_id")
    
    control_plot <- main_dendrites_raw[[2]]
    
    soma_radius <- main_dendrites_raw[[3]]    
    ## defining important parameters of current soma

    n_main_dendrites <- nrow(main_dendrites)
    
    soma_reg <- tibble(x=round(SOMA[["x"]]-soma_radius),
                       xend=round(SOMA[["x"]]+soma_radius),
                       y=round(SOMA[["y"]]-soma_radius),
                       yend=round(SOMA[["y"]]+soma_radius))
    
    #DENDRITE <- main_dendrites[2,]
    elongated_dendrites <- apply(main_dendrites, 1, function(DENDRITE){

      cat(paste("\ndendrite number:",DENDRITE[["dendrite_id"]]))
      
      horizontal_subdivisions <- 60
      vertical_subdivisions <- 10
      horizontal_detection_angle <- 60
      vertical_detection_angle <- 12
      intensity_cutoff <- 0.5
      steps <- 10
      
      elongate_dendrite(DENDRITE %>% as.numeric() %>% set_names(nm=names(DENDRITE)),
                        SOMA,
                        horizontal_subdivisions,
                        vertical_subdivisions,
                        horizontal_detection_angle,
                        vertical_detection_angle,
                        intensity_cutoff,
                        steps)
    })
    ## check in fiji:
    export_dendrites(elongated_dendrites, "f:/data_sholl_analysis/test/example_data/dendrites.csv")
    
    borders_vox_raw <- lapply(1:n_main_dendrites, function(n){
      
      if(n!=n_main_dendrites){
        h_angle <- 0.5*sum(main_dendrites$h_angle[n], main_dendrites$h_angle[n+1])
      } else {
        raw <- main_dendrites$h_angle[n]+0.5*(abs(360-main_dendrites$h_angle[n])+main_dendrites$h_angle[1])
        if(raw>360){
          h_angle <- raw-360
        } else {
          h_angle <- raw
        }
      }
      elongate_3d_sphere(full_image, h_angle, 0, SOMA[["x"]], SOMA[["y"]], SOMA[["z"]], round(sqrt(nc_orig^2+nr_orig^2))) %>%
        select(x,y) %>%
        return()
      
    })      
    
    selection_vector <- rep(c(1:n_main_dendrites), 3) %>% 
      set_names(c((-n_main_dendrites+1):(2*n_main_dendrites)))
    
    dendrite_segments <- lapply(1:n_main_dendrites, create_raster, avs=80) # dendrite
    
    sl <- Reduce(function(x,y)append(x,y),dendrite_segments)
    
    nosoma_image <- remove_soma(SOMA, soma_radius, full_image)
    
    #binary_image <- bi2(sl, nosoma_image, "f:/data_sholl_analysis/test/spec_segs/60_nosoma_autocut.tif")
    
    #writeTIFF(binary_image, "f:/data_sholl_analysis/test/intermediate/binarized.tif")
    
    binary_image <- readTIFF("f:/data_sholl_analysis/test/intermediate/binarized.tif", all=T)
    #filtered_binary <- lapply(binary_image, medianblur, n=2)
    
    main_vectors_raw <- lapply(elongated_dendrites, define_full_vector_info)
    main_vectors <- lapply(main_vectors_raw, nth, 1)
    main_vectors_full=lapply(main_vectors_raw, nth, 2)
    
    
    n <- 2
    
    subd_starts_all_maind <- lapply(1:n_main_dendrites, function(n){
      
      det_rad <- 30
      z_range <- 10
      
      rel_vecs_raw <- main_vectors[[n]]
      full_vecs_raw <- main_vectors_full[[n]]
      
      rel_vecs <- rel_vecs_raw[c(2:length(rel_vecs_raw))]
      full_vecs <- full_vecs_raw[c(2:length(full_vecs_raw))]
      
      full_dendrite <- bind_rows(full_vecs) %>% mutate_all(.funs = as.numeric)
      
      all_sorrounding_vox <- define_surr_layer(rel_vecs, full_vecs, det_rad, z_range)
      
      ##QC
      
      # p <- ggplot(all_sorrounding_vox %>% group_by(x,y) %>% summarize(z=max(z)),
      #             aes(x=x, y=y, color=z))+geom_point(size=0.5)+
      #   geom_segment(data=bind_rows(main_vectors[[n]]),
      #                aes(x=xs, y=ys, xend=xe, yend=ye))
      # 
      # 
      # pdf(file="f:/data_sholl_analysis/test/full_sorrounding.pdf", width=50, height=10)
      # p
      # dev.off()
      
      df_centers <- find_subd_cluster(all_sorrounding_vox)
      
      rescored_subdendrites <- apply(df_centers, 1, fit_path_to_subd_cluster) %>% 
        compact() %>% 
        bind_rows()
        
        return(rescored_subdendrites)
        
    })
      

    
    
    all_sbd <- rescored_subdendrites %>%
      arrange(dist_to_soma)
    
    ## export dendrites for Fiji check
        fiji_ctrl <- all_sbd %>%
          apply(.,1, function(IS){
            tibble(x=c(IS[["x"]], IS[["x_el"]], IS[["x"]]),
                   y=c(IS[["y"]], IS[["y_el"]], IS[["y"]]))
          }) %>%
          bind_rows() %>%
          mutate(id=factor(c(1:nrow(.)), levels=c(1:nrow(.))))
        
        t <- bind_rows(fiji_ctrl, arrange(fiji_ctrl, desc(id))) %>%
          select(-id)
        
        write_csv(t, 
                  file="f:/data_sholl_analysis/test/example_data/subd_full_surr.csv")

    
  }) # soma

}) # file  


    # 
    # vector_lengths <- lapply(1:length(ELD), function(PT){
    #   if(PT==length(ELD)){return(NULL)}
    #   sqrt(sum(abs(ELD[[PT+1]][c("x", "y", "z")]-ELD[[PT]][c("x", "y", "z")])^2)) %>%
    #     return()
    #   
    # }) 
    # 
    # 
    # 


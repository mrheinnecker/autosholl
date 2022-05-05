
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
    
    subdendrite_starts <- lapply(1:n_main_dendrites, function(n){
      
      det_rad <- 30
      z_range <- 10
      
       ## define x-y start points 
      rel_vecs_raw <- main_vectors[[n]]
      full_vecs_raw <- main_vectors_full[[n]]
      
      rel_vecs <- rel_vecs_raw[c(2:length(rel_vecs_raw))]
      full_vecs <- full_vecs_raw[c(2:length(full_vecs_raw))]
      
      v <- 1
      #all_subdendrites <- lapply(1:3, function(v){
      all_subdendrites <- lapply(1:length(rel_vecs), function(v){ 
        cat(paste("\nvec:",v))
        
        VEC <- rel_vecs[[v]]
        FVEC <- full_vecs[[v]]
        
        ortho_t <- VEC[["ha"]]+90
        ortho_b <- VEC[["ha"]]-90
      
        
        ## define starts
          
        if(v==1){
          start_t <- elongate_line(ortho_t, 0, VEC[["xs"]], VEC[["ys"]], VEC[["zs"]], det_rad)
          start_b <- elongate_line(ortho_b, 0, VEC[["xs"]], VEC[["ys"]], VEC[["zs"]], det_rad)        
        } else {
          prevVEC <- rel_vecs[[v-1]]
          elgt_start <- abs(det_rad*tan(deg2rad(0.5*abs((VEC[["ha"]]-90)-(prevVEC[["ha"]]-90)))))
          
          if(prevVEC[["ha"]]>VEC[["ha"]]){
            
            adj_t <- elongate_line(VEC[["ha"]]-180, 0, VEC[["xs"]], VEC[["ys"]], VEC[["zs"]], elgt)
            adj_b <- elongate_line(VEC[["ha"]], 0, VEC[["xs"]], VEC[["ys"]], VEC[["zs"]], elgt)
            
            start_t <- elongate_line(ortho_t, 0, adj_t[["x"]], adj_t[["y"]], adj_t[["z"]], det_rad)
            start_b <- elongate_line(ortho_b, 0, adj_b[["x"]], adj_b[["y"]], adj_b[["z"]], det_rad)
            
          } else {
            
            adj_b <- elongate_line(VEC[["ha"]]-180, 0, VEC[["xs"]], VEC[["ys"]], VEC[["zs"]], elgt)
            adj_t <- elongate_line(VEC[["ha"]], 0, VEC[["xs"]], VEC[["ys"]], VEC[["zs"]], elgt)
            
            start_t <- elongate_line(ortho_t, 0, adj_t[["x"]], adj_t[["y"]], adj_t[["z"]], det_rad)
            start_b <- elongate_line(ortho_b, 0, adj_b[["x"]], adj_b[["y"]], adj_b[["z"]], det_rad)
            
          }
           
        }
        
        ## define elongation lengths
        
        if(v!=length(rel_vecs)){
          
          nextVEC <- rel_vecs[[v+1]]
          
          elgt_end <- abs(det_rad*tan(deg2rad(0.5*abs((VEC[["ha"]]-90)-(nextVEC[["ha"]]-90)))))
          
          if(nextVEC[["ha"]]>VEC[["ha"]]){
            len_b <- VEC[["l"]]+elgt_end
            len_t <- VEC[["l"]]-elgt_end
          } else {
            len_b <- VEC[["l"]]-elgt_end
            len_t <- VEC[["l"]]+elgt_end
          }
          
        } else {
          
          len_b <- VEC[["l"]]
          len_t <- VEC[["l"]]
          
        }
        
     #   adj_start <- elongate_line(VEC[["ha"]]-180, 0, VEC[["xs"]], VEC[["ys"]], VEC[["zs"]], 2*det_rad)
        
    #    start_t <- elongate_line(ortho_t, 0, adj_start[["x"]], adj_start[["y"]], adj_start[["z"]], det_rad)
    #    start_b <- elongate_line(ortho_b, 0, adj_start[["x"]], adj_start[["y"]], adj_start[["z"]], det_rad)
        
        rel_z_layer <- seq(VEC[["zs"]]-z_range, VEC[["zs"]]+z_range, 1) %>%
          .[between(., 1, length(full_image))]
        
        
        ## TOP line
        
        df_top <- screen_subdendrite_starts(start_t, rel_z_layer, det_rad, len_t, VEC)
        df_bottom <- screen_subdendrite_starts(start_b, rel_z_layer, det_rad, len_b, VEC)
        
        cat("\n  screen finished")
        
        #### this for top and bottom
        
        
        #raw_subdendrites <- apply(df_centers, 1, function(line){line[c("x", "y", "z")]}, simplify=F)
      
        #SUBD <- df_centers[1,]
        
        df_centers <- bind_rows(list(df_top, df_bottom) %>% compact())
        
        if(nrow(df_centers)==0){return(NULL)}
        
        #SUBD <- df_centers[1,]
        
        cat(paste("\n  fitting", nrow(df_centers), "subdendrites"))
        rescored_subdendrites <- apply(df_centers, 1, function(SUBD){
        
          #print(SUBD[[1]])
          
          df_dist_raw <- FVEC %>%
            rowwise() %>%
            mutate(dist=dist_pts(x,y,z,SUBD[["x"]],SUBD[["y"]],SUBD[["z"]]),
                   ha=get_ha(x,y,SUBD[["x"]],SUBD[["y"]]),
                   va=get_va(dist, z, SUBD[["z"]])) %>%
            filter(dist<=4*det_rad) #%>%
          if(nrow(df_dist_raw)==0){return(NULL)}
          
          
          df_dist <- df_dist_raw %>%
            mutate(score=fit_subdendrite(SUBD, dist, ha, va)) %>%
            filter(score!=0)
            
          if(nrow(df_dist)==0){return(NULL)}
          
          final_intersection <- df_dist[which(df_dist$score==max(df_dist$score)),] %>%
            .[which(.$dist==min(.$dist)),] %>%
            mutate(dist_to_soma=dist_pts(x,y,z,SOMA[["x"]],SOMA[["y"]],SOMA[["z"]]),
                   x_el=SUBD[["x"]],
                   y_el=SUBD[["y"]],
                   z_el=SUBD[["z"]],) 
          
          return(final_intersection)
          
        }) %>% compact() %>% bind_rows()
        if(nrow(rescored_subdendrites)==0){return(NULL)}
        return(rescored_subdendrites)
        
      })
      
    })
    
    
    
    # all_sbd <- all_subdendrites %>% compact() %>%
    #   bind_rows()
    
    all_sbd <- all_subdendrites %>% compact() %>%
      lapply(function(x){
        arrange(x, dist_to_soma)
      }) %>% bind_rows()
    
    ## export dendrites for Fiji check
        fiji_ctrl <- all_sbd %>%
          apply(.,1, function(IS){
            
            tibble(x=c(IS[["x"]], IS[["x_el"]], IS[["x"]]),
                   y=c(IS[["y"]], IS[["y_el"]], IS[["y"]]))
            
          }) %>%
          #append(list(c(x=round(SOMA[["x"]]), y=round(SOMA[["y"]]))),.) %>%
          bind_rows() %>%
          mutate(id=factor(c(1:nrow(.)), levels=c(1:nrow(.))))
          #rownames_to_column("id") #%>%
          #bind_rows(., arrange(.,desc(id)))
        
        t <- bind_rows(fiji_ctrl, arrange(fiji_ctrl, desc(id))) %>%
          select(-id)
        
        write_csv(t, 
                  file="f:/data_sholl_analysis/test/example_data/subdendrites_4x_adj_overlaps.csv")

    
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


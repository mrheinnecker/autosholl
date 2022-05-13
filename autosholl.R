
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
  
  
#results <- lapply(files, function(file){
#  start_time <- Sys.time()
#  print(file)

  full_image <- readTIFF(file, all=T)
  nr_orig <- nrow(full_image[[1]])
  nc_orig <- ncol(full_image[[1]])
  somata <- get_somata(soma_xy_detection_cube_radius*2+1, 
                       soma_z_detection_radius, 
                       soma_z_detection_degree_steps, 
                       full_image)
  
#  print("soma detected")
  SOMA <- somata[1,] %>% mutate_all(.funs = round) 
  
##  apply(somata,1, function(SOMA){
    ## for each soma a list of all branches and knots is created at this point:
    # full_branches <- list(pos=tibble(x=round(SOMA[["x"]]),
    #                                  y=round(SOMA[["y"]]),
    #                                  z=round(SOMA[["z"]])))
    
    ## check soma for dendritic start sites
    main_dendrites_raw <- find_dendritic_start_sites(SOMA) 
    main_dendrites <- main_dendrites_raw[[1]] %>%
      mutate(z=SOMA[["z"]])  %>%
      arrange(h_angle)%>%
      rownames_to_column("dendrite_id")
    
    #control_plot <- main_dendrites_raw[[2]]
    
    soma_radius <- main_dendrites_raw[[3]]    
    ## defining important parameters of current soma

    n_main_dendrites <- nrow(main_dendrites)

    
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
    #export_dendrites(elongated_dendrites, "f:/data_sholl_analysis/test/example_data/dendrites.csv")
    
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
    main_vectors_df <- lapply(main_vectors, bind_rows)
    
    n <- 2
    
    MASTER <- lapply(1:n_main_dendrites, find_subdendritic_starts, 
                     main_vectors=main_vectors,
                     main_vectors_full=main_vectors_full,
                     SOMA=SOMA)
    

    ## tracing of subdendrite starts
    
    traced_MASTER <- lapply(1:n_main_dendrites, function(nMD){ ## main dendrites
      
      #nMD <- 3
      #nSD <- 8
      
      MV <- main_vectors_df[[nMD]]
      MAIND <- MASTER[[nMD]]
      
      if(is.null(MAIND)){
        return(NULL)
      }
      
      traced_MAIND <- lapply(1:length(MAIND), function(nND){ ## nodes
        cat(paste("\n", nND, "of", length(MAIND)))
        
        NODE <- MAIND[[nND]]
        
        SUBD_list <- NODE$subdend_full
        
        elgt_subd <- lapply(1:length(SUBD_list), function(nSD){ ## sub dendrites
          cat(paste("\n  subd:", nSD, "of", length(SUBD_list)))
          SUBD <- SUBD_list[[nSD]]
          SUBD_info <- SUBD$info
          
          det_rad <- 10
          z_range <- 5
        
          xs <- SUBD_info$x
          ys <- SUBD_info$y
          zs <- SUBD_info$z
          ha <- adj_deg(SUBD_info$ha+180)
          ha_cutoff <- MV %>% filter(level==SUBD_info$cut_vec) %>% pull(ha)
          coord_list <- list()
          c <- 1
          abort <- F
          while(abort==F){ 
            cat(paste("\n    elgt step:",c))
            if(c>2){
              ha_cutoff <- NULL
            }
          #for(i in 1:4){  
            centers <- screen_circular_new(det_rad, z_range,1, xs, ys, zs, ha, ha_cutoff) %>%
              find_subd_cluster_man(., 5, 10)
            
            if(length(centers)==0){
              ## no elongation detected
              abort=T
              cat("\n    stopped")
            } else if(length(centers)==1){
              ## subdendrite elongates
              elgt <- centers[[1]]
              ha <- get_ha(elgt[["x"]], elgt[["y"]], xs, ys)
              xs <- elgt[["x"]]
              ys <- elgt[["y"]]
              zs <- elgt[["z"]]
              coord_list[[c]] <- elgt
              c <- c+1
            } else {
              ## node detected
              cat("\n    node")
              abort <- T
            }
          }
          
          new_SUBD <- list(SUBD$info, bind_rows(SUBD$full_coords, bind_rows(coord_list)))
          
          return(new_SUBD)
        })
        
        new_NODE <- append(NODE[c(1:5)], list(subdend_full=elgt_subd))
        return(new_NODE)
    })
      # te <- bind_rows(coord_list) %>% select(x,y) %>% bind_rows(SUBD_info %>% select(x, y),.) %>% mutate(id=c(1:nrow(.)))
      # 
      #   write_csv(bind_rows(te%>% select(-id), arrange(te, desc(id)) %>% select(-id)), 
      #         "f:/data_sholl_analysis/test/dendrites/test.csv")  
      # 
        
    #  new_MASTER <- MASTER  
    
     # new_MASTER[[2]] <- traced_MAIND
    return(traced_MAIND)        
  })
    
    
    export_structure(traced_MASTER, "f:/data_sholl_analysis/test/dendrites/first.csv")
    
    
    
    
    
    n <- 2
    
    elongated_subdendrites_maind <- lapply(1:n_main_dendrites, function(n){
      
      cat(paste("\nmain dendrite:", n))
      
      det_rad <- 10
      z_range <- 5
      
      subd_list <- subd_starts_all_maind[[n]]
      traced_subd_list <- list()
      #SUBD <- subd_starts[3,]
      nsd <- 0
      while(nsd<length(subd_list)){
        nsd <- nsd+1
        cat(paste("\nsubd:", nsd, " of", length(subd_list)))
        SUBD <- subd_list[[nsd]]
        node <- first(SUBD)
        elongation_pt <- last(SUBD)
        nSUB <- 1
        for (nSUB in 1:length(elongation_pt)) {
          cat(paste("\n  trace:", nSUB, " of", length(elongation_pt)))
          traced <- trace_dendrite(elongation_pt[[nSUB]],
                                   headnode=node, 
                                   det_rad, 
                                   z_range)
          ## add elongated vector to current node of main list
          subd_list[[nsd]][["end"]][[nSUB]] <- traced$vec
          if(traced$sub){
            ## further sub node detected --> add node to main list
            subd_list <- append(subd_list, list(traced[c("pos", "end")]))
          } 
        }
      }
      #trace_dendrite(SUBD)
      traced <- apply(subd_starts, 1, trace_dendrite) ## subdendrite
    }) ## main dendrite
##  }) # soma
}) # file     
    
  


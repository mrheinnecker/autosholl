
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
    
    rem_mainobj_image <- remove_main_dendrites(nosoma_image, main_vectors, main_vectors_full)
    
    
    
    writeTIFF(rem_mainobj_image, "f:/data_sholl_analysis/test/intermediate/no_main_vecs_no_soma.tif")
    binary_image <- bi2(sl, rem_mainobj_image, "f:/data_sholl_analysis/test/spec_segs/60_nosoma_autocut.tif")
    
    writeTIFF(binary_image, "f:/data_sholl_analysis/test/intermediate/new_binary.tif")
    
    binary_image <- readTIFF("f:/data_sholl_analysis/test/intermediate/binarized.tif", all=T)
    #filtered_binary <- lapply(binary_image, medianblur, n=2)
    
    main_vectors_raw <- lapply(elongated_dendrites, define_full_vector_info)
    main_vectors <- lapply(main_vectors_raw, nth, 1)
    main_vectors_full=lapply(main_vectors_raw, nth, 2)
    main_vectors_df <- lapply(main_vectors, bind_rows)
    
    #n <- 2
    
    MASTER <- lapply(1:n_main_dendrites, find_subdendritic_starts, 
                     main_vectors=main_vectors,
                     main_vectors_full=main_vectors_full,
                     SOMA=SOMA)
    
    ## check for duplicated dendrites
    #rem_dup_MASTER <- lapply(1:n_main_dendrites, adjust_duplicated_starts)
    
    export_structure(MASTER, "f:/data_sholl_analysis/test/dendrites/subd_starts_adj2.csv")
    
    ## tracing of subdendrite starts
    
    traced_MASTER <- lapply(1:n_main_dendrites, function(nMD){ ## main dendrites
      
      nMD <- 2
      nND <- 17
      nSD <- 2
      
      MV <- main_vectors_df[[nMD]]
      MAIND <- MASTER[[nMD]]
      
      if(is.null(MAIND)){
        return(NULL)
      }
      
      ## for testing ... remove nodes at end of dendrite...
      # MAIND <- MAIND[c(1:15)]
      # MAIND[[15]]$subnodes <- list()
      # MAIND[[15]]$subnode_full <- list()
      nNODE_orig <- length(MAIND)
      
      all_vectors <- lapply(MAIND, get_all_vectors) %>%
        Reduce(function(x,y)append(x,y),.)

      all_vectors_fv <- lapply(all_vectors, get_full_vector_voxels)
      
      nND <- 0
      n_processed <- 0
      while(n_processed < length(MAIND)){ ## nodes
        nND <- nND + 1
 
        cat(paste("\n", nND, "of", length(MAIND)))
        
        NODE <- MAIND[[nND]]
        
        SUBD_list <- NODE$subdend_full
        
        new_subnodes <- NODE$subnodes
        new_SUBD_list <- list()
        new_SUBN_list <- NODE$subnode_full
        
        for(nSD in 1:length(SUBD_list)){
          
          cat(paste("\n  subd:", nSD, "of", length(SUBD_list)))
          SUBD <- SUBD_list[[nSD]]
          SUBD_info <- SUBD$info
          
          det_rad <- 12
          #z_range <- 5
        
          xs <- SUBD_info$x
          ys <- SUBD_info$y
          zs <- SUBD_info$z
          ha <- adj_deg(SUBD_info$ha)
          ## only define ha_cutoff if node intersects with main vectors

          if(SUBD_info$d_to_cv>45){
            ha_cutoff <- NULL
          } else {
            ha_cutoff <- MV %>% 
              filter(level==SUBD_info$clos_vec) %>% 
              pull(ha)
          }
          
          
          coord_list <- list()
          c <- 1
          abort <- F
          node_detected <- F
          remove_dend <- F
          
          while(abort==F){ 
            #st <- Sys.time()
            cat(paste("\n    elgt step:",c))
            if(c>4){
              ha_cutoff <- NULL
              screening_angle <- 160
              z_range <- 3
            } else {
              screening_angle <- 180
              z_range <- 5
            }
          #for(i in 1:4){  
           # t1 <- Sys.time()
            
            Z <- round(zs+tan(deg2rad(va))*det_rad)
            
            rel_z_layer <- seq(Z-z_range, Z+z_range, 1) %>%
              .[between(., 1, length(full_image))]
            
            centers1 <- screen_circular_new(det_rad, rel_z_layer,
                                            1, xs, ys, zs, ha, ha_cutoff, screening_angle) #%>%
          #  t2 <- Sys.time()
            centers <- find_subd_cluster_man(centers1, 5, 10)
            
            #centers <- find_subd_cluster_man(centers1, 5, 8)
            
            
          #  t3 <- Sys.time()
          #  cat(paste("  ", round(t1-st, 2), round(t2-t1, 2), round(t3-t2, 2)))
            if(length(centers)==0|c>30){
              ## no elongation detected
              abort=T
              cat("\n    stopped")
              if(c==1){
                remove_dend <- T
              }
            } else if(length(centers)==1){
              ## subdendrite elongates
              
              elgt <- centers[[1]]
              dist <- dist_pts(elgt[["x"]], elgt[["y"]], elgt[["z"]], xs, ys,zs)
              ha <- get_ha(elgt[["x"]], elgt[["y"]], xs,ys)
              va <- get_va(dist, 
                           elgt[["z"]], zs)

              ## check if elongation overlaps with already present dendrite
              vox <- elongate_3d_sphere(binary_image, ha, va, xs, ys, zs, round(dist)) %>%
                mutate(id=paste(x,y, sep="_")) %>%
                pull(id)
              
              overlap <- all_vectors_fv %>%
                bind_rows() %>%
                mutate(id=paste(x,y, sep="_")) %>%
                filter(id %in% vox)
              
              if(nrow(overlap)>0&c>1){
                abort <- T
                cat("\n    crossed")
              } else {
                xs <- elgt[["x"]]
                ys <- elgt[["y"]]
                zs <- elgt[["z"]]    
                coord_list[[c]] <- elgt
                c <- c+1
              } 
              
            } else {
              ## node detected
              cat("\n    node")
              node_detected <- T
              abort <- T
            }
            
          }
          
          if(node_detected==T){
            
            new_node_id <- length(MAIND)+1
            
            if(c==1){
              last_pos <- SUBD$full_coords[nrow(SUBD$full_coords),]
            } else {
              last_pos=as_tibble(last(coord_list))
            }
            ## remove dendrite from subdendrite list in current node
            ## by just not adding it to new list

            ## add subnode to subnode list of current node
            
            new_subnodes <- append(new_subnodes, new_node_id)
            
            ## add full subnode to full subnode list of current node
            
            new_SUBN <- list(node_id=new_node_id,
                             info=SUBD$info, 
                             full_coords=bind_rows(SUBD$full_coords, bind_rows(coord_list)))
            
            new_SUBN_list <- append(new_SUBN_list, list(new_SUBN))

            ## add new node to full node list
            
            proc_centers <- lapply(centers, format_new_subdendrite_starts, nMD=nMD, 
                                   main_vectors_full=main_vectors_full, 
                                   last_pos=last_pos)
            
            
            new_NODE <- list(node_id=new_node_id,
                             master_id=NODE$node_id,
                             subnodes=list(),
                             pos=last_pos,
                             subnode_full=list(),
                             subdend_full=proc_centers)

            MAIND <- append(MAIND, list(new_NODE))

            additional_vectors <- get_all_vectors(new_NODE) %>% lapply(get_full_vector_voxels)
          } else if(remove_dend==T){ 
             
            ## remove current dendrite from node
            #    if(sum(length(), length()))
            
            
          } else {

            ## update subdendrite list in current node
            traced_SUBD <- list(info=SUBD$info, full_coords=bind_rows(SUBD$full_coords, bind_rows(coord_list)))
            new_SUBD_list <- append(new_SUBD_list, list(traced_SUBD))

            additional_vectors <- list()
          }
          
          
          
        } ## end of for... over subdendrite list
        
        MAIND[[nND]] <- list(node_id=NODE$node_id,
                             master_id=NODE$master_id,
                             subnodes=new_subnodes,
                             pos=NODE$pos,
                             subnode_full=new_SUBN_list,
                             subdend_full=new_SUBD_list)
        
        all_vectors_fv <- append(all_vectors_fv, additional_vectors)
        
        n_processed <- n_processed + 1

    } ## end of while ... return full new MAIND
     
    return(traced_MAIND)        
  })
    
    
    
    traced_MASTER <- MASTER
    
    traced_MASTER <- list(MAIND)
    
    export_structure(traced_MASTER, "f:/data_sholl_analysis/test/dendrites/dr12_ang150_11.csv")
    
    
    
    
    
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
    
  


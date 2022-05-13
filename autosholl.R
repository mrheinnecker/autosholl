
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
    
    MASTER <- lapply(1:n_main_dendrites, function(n){
    
      cat(paste("\nmain dendrite:", n))
        
      det_rad <- 30
      z_range <- 10
      
      rel_vecs_raw <- main_vectors[[n]]
      full_vecs_raw <- main_vectors_full[[n]]
      
      rel_vecs <- rel_vecs_raw[c(2:length(rel_vecs_raw))]
      full_vecs <- full_vecs_raw[c(2:length(full_vecs_raw))]
      
      full_dendrite <- bind_rows(full_vecs) %>% 
        mutate_all(.funs = as.numeric)
      
      all_sorrounding_vox <- define_surr_layer(rel_vecs, full_vecs, det_rad, z_range) %>%
        filter(dist_to_soma>2.5*soma_radius)
      
      #print(1)
      
      #show_surrounding_layer(all_sorrounding_vox, "f:/data_sholl_analysis/test/full_sorrounding.pdf")

      df_centers <- find_subd_cluster(all_sorrounding_vox)
      
      #print(2)
      ## if no subdendrite starts are found... the main dendrite is just a subdendrite
      if(is.null(df_centers)){
        #full_branches
        return(NULL)
      }
      rescored_subdendrites <- apply(bind_rows(df_centers), 1, fit_path_to_subd_cluster, 
                                     full_dendrite=full_dendrite,
                                     det_rad=det_rad) %>% 
        compact() #%>% 
        #bind_rows()
        
      ## get hierarchy
      
      df_poslev <- rescored_subdendrites %>%
        bind_rows() %>%
        rowwise() %>%
        mutate(dist_to_soma=dist_pts(xs,ys,zs, round(SOMA[["x"]]),round(SOMA[["y"]]),round(SOMA[["z"]]))) %>%
        left_join(bind_rows(main_vectors[[n]]) %>% select(level, dir), by="level")
      
      
      df_sorted <- lapply(unique(df_poslev$level), function(LEV){
        rel_rows <- df_poslev %>%
          filter(level==LEV)
        
        if(unique(rel_rows$dir)==1){
          return(arrange(rel_rows, dist_to_soma))
        } else {
          return(arrange(rel_rows, desc(dist_to_soma)))
        }
        
      }) %>% bind_rows() %>% ungroup() %>%
        mutate(node=c(1:nrow(.)))
      
      
      #idNODE <- 12
      assigned_nodes <- lapply(df_sorted$node, function(idNODE){
        
        node <- df_sorted %>%
          filter(node==idNODE)
        
        #master_node <- idNODE-1
        
        sub_dend <- node %>% select(x=xe, y=ye, z=ze, ha, va)  %>%
          mutate(class="dend",
                 sub_id=NA)
        
        sub_dend_full <- bind_rows(select(node, x=xs,y=ys,z=zs),
                    select(node, x=xe,y=ye,z=ze))
        
        
        full_subdend_info <- list(sub_dend, sub_dend_full)
        
        if(idNODE==nrow(df_sorted)){
          ## no know subnodes so far ... just elongation directions
          #sub_node <- NULL
          subnode_id <- NULL
          #full_vec_pos <- NULL
          full_subnode_info <- NULL
        } else {
          next_node <- df_sorted %>% filter(node==idNODE+1)
          subnode_id <- next_node$node
          sub_node <- next_node %>%
            select(x=xs, y=ys, z=zs, sub_id=node) %>%
            rowwise() %>%
            mutate(ha=get_ha(x,y,node$xs, node$ys),
                   dist=dist_pts(x,y,z,node$xs, node$ys, node$zs),
                   va=get_va(dist, z, node$zs),
                   class="node") %>%
            select(-dist)
          
          if(next_node$level!=node$level){
            ## vector edges need to be inserted
            full_vec_pos <- main_vectors_df[[n]] %>% #rowwise() %>%
              filter(between(level, node$level, as.numeric(next_node$level-1))) %>%
              select(x=xe, y=ye, z=ze) %>%
              bind_rows(select(node, x=xs, y=ys, z=zs),
                        .,
                        select(next_node, x=xs, y=ys, z=zs))
          } else {
            
            full_vec_pos <- bind_rows(select(node, x=xs, y=ys, z=zs),
                                      select(next_node, x=xs, y=ys, z=zs))
            
          }
          
          full_subnode_info <- list(node_id=idNODE+1, 
                                    info=sub_node, 
                                    full_coords=full_vec_pos)
          
        }
        

        final_list <- list(node_id=node$node,
                           master_id=idNODE-1,
                           subnodes=list(subnode_id) %>% compact(),
                           pos=node %>% select(x=xs, y=ys, z=zs),
                           subnode_full=list(full_subnode_info) %>% compact(),
                           subdend_full=list(full_subdend_info)
                           )
        return(final_list)
          
        
      })
      
        return(assigned_nodes)
        
    })
      
    
    
    export_structure(MASTER, "f:/data_sholl_analysis/test/example_data/full_dendrites.csv")
    
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
    
    

fiji_traced <- lapply(traced, last) %>%
  lapply(function(IS){
    tibble(x=c(IS[["x"]], IS[["xe"]], IS[["x"]]),
           y=c(IS[["y"]], IS[["ye"]], IS[["y"]])) %>%
      return()
  }) %>%
  bind_rows()%>%
  mutate(id=factor(c(1:nrow(.)), levels=c(1:nrow(.))))



write_csv(bind_rows(fiji_traced, arrange(fiji_traced, desc(id))) %>%
            select(-id), 
          file="f:/data_sholl_analysis/test/example_data/md4_traced_dr10.csv")    
    
    
    all_sbd <- subd_starts %>%
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
                  file="f:/data_sholl_analysis/test/example_data/subd_comp_full_surr.csv")

    
 


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


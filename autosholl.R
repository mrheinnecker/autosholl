
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
    
    #DENDRITE <- main_dendrites[5,]
    elongated_dendrites <- apply(main_dendrites, 1, function(DENDRITE){

      cat(paste("\ndendrite number:",DENDRITE[["dendrite_id"]]))
      
      horizontal_subdivisions <- 60
      vertical_subdivisions <- 10
      horizontal_detection_angle <- 60
      vertical_detection_angle <- 12
      intensity_cutoff <- 0.5
      steps <- 10
      
      elongate_dendrite(DENDRITE %>% as.numeric() %>% set_names(nm=names(DENDRITE)),
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
      elongate_3d_sphere(h_angle, 0, SOMA[["x"]], SOMA[["y"]], SOMA[["z"]], round(sqrt(nc_orig^2+nr_orig^2))) %>%
        select(x,y) %>%
        return()
      
    })      
    
    selection_vector <- rep(c(1:n_main_dendrites), 3) %>% 
      set_names(c((-n_main_dendrites+1):(2*n_main_dendrites)))
    ### now select areas for cutoff use
    
    #nELD <- 1
    dendrite_segments <- lapply(1:n_main_dendrites, function(nELD){
      cat(paste("dendrite:", nELD))
################################################################################
######################## defining basic elements ###############################    
################################################################################       
      ELD_raw <- elongated_dendrites[[nELD]] 
      avs <- 120 
      overall_vector_raw <- ELD_raw[[length(ELD_raw)]][c("x", "y", "z")]-ELD_raw[[1]][c("x", "y", "z")]
################################################################################
## defining secondary elements accroding to vertical or horizonatal dendrite ###    
################################################################################      
      if(abs(overall_vector_raw[1])>abs(overall_vector_raw[2])){
        direction <- "h"
        ELD <- ELD_raw
        overall_vector <- overall_vector_raw
        borders_vox <- borders_vox_raw
        nr <- nr_orig
        nc <- nc_orig
      } else {
        direction <- "v"
        ELD <- lapply(ELD_raw, function(VEC){
          return(c(x=VEC[["y"]], y=VEC[["x"]], z=VEC[["z"]],  h_angle=VEC[["h_angle"]],  v_angle=VEC[["v_angle"]]))
        })
        overall_vector <- c(x=overall_vector_raw[["y"]], y=overall_vector_raw[["x"]], 
                            z=overall_vector_raw[["z"]])
        borders_vox <- lapply(borders_vox_raw, function(BORD){
          tibble(x=BORD$y, y=BORD$x) %>%
            return()
        })
        nc <- nr_orig
        nr <- nc_orig
      }
################################################################################
## defining tertiary elements accroding to v/h and direction (left/right) ######    
################################################################################ 
      overall_length <- overall_vector["x"]
      rv <- overall_vector["x"]/abs(overall_vector["x"])      
      
      if((direction=="h"&rv==1)|(direction=="v"&rv==-1)){
        top_border <- borders_vox[[selection_vector[[as.character(nELD)]]]]
        bottom_border <- borders_vox[[selection_vector[[as.character(nELD-1)]]]]
      } else {
        top_border <- borders_vox[[selection_vector[[as.character(nELD-1)]]]]
        bottom_border <- borders_vox[[selection_vector[[as.character(nELD)]]]]
      }

      if(rv==1){
        pixels_to_image_border <- nc-ELD[[1]]["x"]
        ELD[[length(ELD)+1]] <- c(x=nc, ELD[[length(ELD)]][c("y", "z")], h_angle=0, v_angle=0)
      } else {
        pixels_to_image_border <- 1-ELD[[1]]["x"]
        ELD[[length(ELD)+1]] <- c(x=1, ELD[[length(ELD)]][c("y", "z")], h_angle=180, v_angle=0)
      }
################################################################################
################### calculate required data from input #########################    
################################################################################  
      vector_pos <- create_df_of_vectors(ELD)
      vectors <- create_rv_of_vectors(ELD)
      xnorm_vector <- lapply(vectors, function(V){V[c(1:3)]/abs(V["x"])})
        
      n_segments <- abs(round(pixels_to_image_border/avs))
      use_length <- round(pixels_to_image_border/n_segments)
      
      res_list <- assign_vectors_to_segments(ELD, vector_pos, use_length, n_segments)
      
      full_vecs <- bind_rows(vectors) %>%
        pull(x) %>% cumsum()  
      
      y_coord_list <- list()
      segment_list <- list()
      segment_count <- 0
      for(n in c(1:n_segments)){
        cat(paste("\n  x_segment:",n))
################################################################################
################### set x and y variables for segment ##########################    
################################################################################ 
        xs <- ELD[[1]]["x"]+(n-1)*use_length+1
        xe <- ELD[[1]]["x"]+n*use_length 
        ys <- ifelse(n==1,
                     as.numeric(ELD[[1]]["y"]),
                     as.numeric(y_coord_list[[n-1]]))  
          ## select relevant vectors for that segment 
        rl <- combine_vectors(ELD, xs, xe, rv, full_vecs, xnorm_vector, res_list, n)  

        med_line <-   tibble(x=c(xs:xe)) %>%
          mutate(fac=c(rl, recursive=T)) %>%
            mutate(y=ceiling(ys+cumsum(fac))) %>%
            filter(between(x, 1, nc),
                   between(y, 1, nr))
          ## from here the segemt height and depth is measured until image border
          med <- round(0.5*(max(med_line$y)+min(med_line$y)))
          
          pix_to_bottom <- -med
          pix_to_top <- nr-med
          n_segments_top <- abs(round(pix_to_top/avs))
          n_segments_bottom <- abs(round(pix_to_bottom/avs))
          use_length_top <- round(pix_to_top/n_segments_top)
          use_length_bottom <- round(pix_to_bottom/n_segments_bottom)
          
          ## new loop for all segments
          
          start_list_top <- list()
          line_list_top <- list()
          
          for(n_top in 1:n_segments_top){
            
            if(n_top==1){
              start <- med
              line <- med_line
            } else {
              start <- start_list_top[[n_top-1]]+1
              line <- line_list_top[[n_top-1]]
            }
            
            top <- start+use_length_top
            if(top>nr){top <- nr}
            
            all_vox_top <- create_segment(xs, xe, top, line, nr, nc, top_border, "top",
                                          direction, nr_orig)
            
            start_list_top[[n_top]] <- top
            line_list_top[[n_top]] <- all_vox_top[[2]]
            
            segment_count <- segment_count+1
            segment_list[[segment_count]] <- all_vox_top[[1]]
            cat(paste("\n    top:", segment_count))
          }
          
          start_list_bottom <- list()
          line_list_bottom <- list()
          
          for(n_bottom in 1:n_segments_bottom){
            
            if(n_bottom==1){
              start <- med
              line <- med_line
            } else {
              start <- start_list_bottom[[n_bottom-1]]-1
              line <- line_list_bottom[[n_bottom-1]]
            }
            
            bottom <- start+use_length_bottom
            if(bottom<1){bottom <- 1} 
            
            all_vox_bottom <- create_segment(xs, xe, bottom, line, nr, nc, bottom_border, "bottom",
                                             direction, nr_orig)
            
            start_list_bottom[[n_bottom]] <- bottom
            line_list_bottom[[n_bottom]] <- all_vox_bottom[[2]]
            
            segment_count <- segment_count+1
            cat(paste("\n    bottom:", segment_count))
            segment_list[[segment_count]] <- all_vox_bottom[[1]]
          }       
              
          y_coord_list[[n]] <- med_line[[nrow(med_line), "y"]]
          

        }
    return(segment_list)  
    }) # dendrite
    
    res <- dendrite_segments
    
    test_segmentation(c(res[[1]], res[[2]], res[[3]], res[[4]]), 
                      full_image, 
                      "f:/data_sholl_analysis/test/spec_segs/testseg.tif", 
                      soma_reg)    
    
    
    normalize_regions(c(res[[1]], res[[2]], res[[3]], res[[4]]), 
                      full_image, 
                      "f:/data_sholl_analysis/test/spec_segs/final_120.tif", 
                      soma_reg)
    
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


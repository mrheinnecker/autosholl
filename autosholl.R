
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
  nr <- nrow(full_image[[1]])
  nc <- ncol(full_image[[1]])
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
      
      #return(list(knot_list, full_vector_list))
    })
    ## check in fiji:
    export_dendrites(elongated_dendrites, "f:/data_sholl_analysis/test/example_data/dendrites.csv")
    
    borders_vox <- lapply(1:nrow(main_dendrites), function(n){
      
      if(n!=nrow(main_dendrites)){
        h_angle <- 0.5*sum(main_dendrites$h_angle[n], main_dendrites$h_angle[n+1])
      } else {
        raw <- main_dendrites$h_angle[n]+0.5*(abs(360-main_dendrites$h_angle[n])+main_dendrites$h_angle[1])
        if(raw>360){
          h_angle <- raw-360
        } else {
          h_angle <- raw
        }
      }
      
      elongate_3d_sphere(h_angle, 0, SOMA[["x"]], SOMA[["y"]], SOMA[["z"]], nc) %>%
        select(x,y) %>%
        return()
      
    })      
    
    # borders_vox <- sapply(borders, function(h_angle){
    #         
    #         elongate_3d_sphere(h_angle, 0, SOMA[["x"]], SOMA[["y"]], SOMA[["z"]], nc) %>%
    #           select(x,y) %>%
    #           return()
    #         
    # }, simplify=F)
    ### now select areas for cutoff use
    
    nELD <- 2
    dendrite_segments <- lapply(1:length(elongated_dendrites), function(nELD){
    #dendrite_segments <- lapply(c(2,4), function(nELD){  
      print(paste("dendrite:", nELD))
      ELD <- elongated_dendrites[[nELD]] 
      avs <- 200 
      overall_vector <- ELD[[length(ELD)]][c("x", "y", "z")]-ELD[[1]][c("x", "y", "z")]
      

      if(abs(overall_vector[1])>abs(overall_vector[2])){
################################################################################
################################################################################
################## horizonatl dendrite orientation #############################     
################################################################################ 
################################################################################        
        overall_length <- overall_vector["x"]
        rv <- overall_vector["x"]/abs(overall_vector["x"])
        if(rv==1){
################################################################################  
################## vector direction: right #####################################          
################################################################################    
          pixels_to_image_border <- nc-ELD[[1]]["x"]
          ELD[[length(ELD)+1]] <- c(x=nc, ELD[[length(ELD)]][c("y", "z")], h_angle=0, v_angle=0)
          if(nELD==1){
            top_border <- borders_vox[[1]]
            bottom_border <- borders_vox[[length(borders_vox)]]
          } else {
            top_border <- borders_vox[[nELD]]
            bottom_border <- borders_vox[[nELD-1]]  
          }
        } else {
################################################################################  
################## vector direction: left  #####################################          
################################################################################          
          pixels_to_image_border <- 1-ELD[[1]]["x"]
          ELD[[length(ELD)+1]] <- c(x=1, ELD[[length(ELD)]][c("y", "z")], h_angle=180, v_angle=0)
          if(nELD==1){
            top_border <- borders_vox[[length(borders_vox)]]
            bottom_border <- borders_vox[[1]]
          } else {
            top_border <- borders_vox[[nELD-1]]
            bottom_border <- borders_vox[[nELD]]  
          }
        }
      } else {
################################################################################
################################################################################
################## vertical dendrite orientation ###############################        
################################################################################ 
################################################################################  
        overall_length <- overall_vector["y"]
        rv <- overall_vector["y"]/abs(overall_vector["y"])
        if(rv==1){
################################################################################  
################## vector direction: top   #####################################          
################################################################################
          pixels_to_image_border <- nr-ELD[[1]]["y"]
          ELD[[length(ELD)+1]] <- c(ELD[[length(ELD)]]["x"], y=nr, ELD[[length(ELD)]]["z"], h_angle=90, v_angle=0)
        } else {
################################################################################  
################## vector direction: bottom ####################################          
################################################################################
          pixels_to_image_border <- 1-ELD[[1]]["y"]
          ELD[[length(ELD)+1]] <- c(ELD[[length(ELD)]]["x"], y=1, ELD[[length(ELD)]]["z"], h_angle=270, v_angle=0)
        } 
      }                          
      
      vector_pos <- create_df_of_vectors(ELD)
      vectors <- create_rv_of_vectors(ELD)
      xnorm_vector <- lapply(vectors, function(V){V[c(1:3)]/abs(V["x"])})
      ynorm_vector <-  lapply(vectors, function(V){V[c(1:3)]/abs(V["y"])}) 
        
      n_segments <- abs(round(pixels_to_image_border/avs))
      use_length <- round(pixels_to_image_border/n_segments)
      
      res_list <- assign_vectors_to_segments(ELD, vector_pos, use_length, n_segments)
      
      
      full_vecs <- bind_rows(vectors) %>%
        pull(x) %>% cumsum()   
      
      final_list <- list()
      segment_list <- list()
      for(n in c(1:n_segments)){
        print(n)
          ## set x-limits of segment
          xs <- ELD[[1]]["x"]+(n-1)*use_length+1
          xe <- ELD[[1]]["x"]+n*use_length 
          ys <- ifelse(n==1,
                       as.numeric(ELD[[1]]["y"]),
                       as.numeric(final_list[[n-1]][nrow(final_list[[n-1]]),"y"]))  
          
       #   print(paste("start:", xs, "end:", xe))
          
      #}
          ## select relevant vectors for that segment
          rel_vecs <- full_vecs[res_list[[n]]]
          rl <- list()
          last_end <- xs-rv
          for(i in 1:length(rel_vecs)){
            #print(i)
            if(i==1){
              start <- xs
            } else {
              start <- last_end+rv
            } 
            
            if(i==length(rel_vecs)){
              end <- xe 
            } else {
              end <- ELD[[1]]["x"]+rel_vecs[i]
            }
            last_end <- end
            #print(paste(start, "-", end))
            rl[[i]] <- rep(xnorm_vector[[res_list[[n]][i]]]["y"], abs(end-start+rv))
          }
          

          
          ## create medium line (runs at the main dendrite)
          #tibble(x=xs+c(1:abs(use_length)*rv)) %>%
          med_line <-   tibble(x=c(xs:xe)) %>%
            
            
            mutate(fac=Reduce(function(x,y)c(x,y),rl)) %>%
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
          
          top <- med+use_length_top
          if(top>nr){top <- nr}
          
          rel_of_border_top <- top_border %>%
            filter(x %in% c(xs:xe)) %>%
            filter(y<top) %>%
            group_by(x) %>%
            summarize(y=min(y))
          
           
          top_line <- tibble(x=c(xs:xe),
                             y=top) %>%
            filter(!x %in% rel_of_border_top$x) %>%
            bind_rows(rel_of_border_top)%>%
            filter(between(x, 1, nc),
                   between(y, 1, nr))
           
          
          all_vox_top <- lapply(top_line$x, function(X){
            #print(X)
            s <- top_line %>%
              filter(x==X) %>%
              pull(y)
            e <- med_line %>%
              filter(x==X) %>%
              pull(y)
            xn <- (X-1)*nr
            return(c(s:e)+xn)
            
          }) %>% Reduce(function(x,y)c(x,y),.)
          
          
          bottom <- med+use_length_bottom
          if(bottom<1){bottom <- 1}
          
          rel_of_border_bottom <- bottom_border %>%
            
            
            filter(x %in% c(xs:xe)) %>%
            filter(y>bottom) %>%
            group_by(x) %>%
            summarize(y=max(y))
          
          
          bottom_line <- tibble(x=c(xs:xe),
                             y=bottom) %>%
            filter(!x %in% rel_of_border_bottom$x) %>%
            bind_rows(rel_of_border_bottom)%>%
            filter(between(x, 1, nc),
                   between(y, 1, nr))
          
          
          all_vox_bottom <- lapply(bottom_line$x, function(X){
            
            s <- bottom_line %>%
              filter(x==X) %>%
              pull(y)
            e <- med_line %>%
              filter(x==X) %>%
              pull(y)
            xn <- (X-1)*nr
            return(c(s:e)+xn)
            
          }) %>% Reduce(function(x,y)c(x,y),.)
          
          
          final_list[[n]] <- med_line
          segment_list[[2*n-1]] <- all_vox_top
          segment_list[[2*n]] <- all_vox_bottom
          ########################################################################
          # comp_img <- lapply(full_image, function(LAYER){
          #   
          #   LAYER[all_vox_top] <- 1
          #   LAYER[all_vox_bottom] <- 0.5
          #   return(LAYER)
          # })
          # writeTIFF(comp_img, paste0("f:/data_sholl_analysis/test/check",n,".tif"))
          ########################################################################
        }
    return(segment_list)
        # normalize_regions(segment_list, 
        #                   full_image, 
        #                   "f:/data_sholl_analysis/test/spec_segs/test_right_dendrite2.tif", 
        #                   soma_reg)
    }) # dendrite
    
    normalize_regions(tl, 
                                         full_image, 
                                         "f:/data_sholl_analysis/test/spec_segs/both.tif", 
                                         soma_reg)
    fd <- Reduce(function(x,y)c(x,y),segment_list)
    
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


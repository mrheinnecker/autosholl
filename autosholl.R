
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
      
      #return(list(knot_list, full_vector_list))
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
    dendrite_segments <- lapply(1:length(elongated_dendrites), function(nELD){
      print(paste("dendrite:", nELD))
      ELD_raw <- elongated_dendrites[[nELD]] 
      avs <- 200 
      overall_vector_raw <- ELD_raw[[length(ELD_raw)]][c("x", "y", "z")]-ELD_raw[[1]][c("x", "y", "z")]
      
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
 
      vector_pos <- create_df_of_vectors(ELD)
      vectors <- create_rv_of_vectors(ELD)
      xnorm_vector <- lapply(vectors, function(V){V[c(1:3)]/abs(V["x"])})
        
      n_segments <- abs(round(pixels_to_image_border/avs))
      use_length <- round(pixels_to_image_border/n_segments)
      
      res_list <- assign_vectors_to_segments(ELD, vector_pos, use_length, n_segments)
      
      full_vecs <- bind_rows(vectors) %>%
        pull(x) %>% cumsum()  
      
      final_list <- list()
      segment_list <- list()
      time_list <- list()
      for(n in c(1:n_segments)){
        print(n)
#t0 <- Sys.time()
        ## set x-limits of segment
        xs <- ELD[[1]]["x"]+(n-1)*use_length+1
        xe <- ELD[[1]]["x"]+n*use_length 
        ys <- ifelse(n==1,
                     as.numeric(ELD[[1]]["y"]),
                     as.numeric(final_list[[n-1]][nrow(final_list[[n-1]]),"y"]))  
          ## select relevant vectors for that segment
#t1 <- Sys.time()  
        rl <- combine_vectors(ELD, xs, xe, rv, full_vecs, xnorm_vector, res_list, n)
#t2 <- Sys.time()   

        med_line <-   tibble(x=c(xs:xe)) %>%
            #mutate(fac=Reduce(function(x,y)c(x,y),rl)) %>%
          mutate(fac=c(rl, recursive=T)) %>%
            mutate(y=ceiling(ys+cumsum(fac))) %>%
            filter(between(x, 1, nc),
                   between(y, 1, nr))
#t3 <- Sys.time()   
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
          bottom <- med+use_length_bottom
          if(bottom<1){bottom <- 1} 
          
#t4 <- Sys.time()
          
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
          
#t5 <- Sys.time()     
          all_vox_top_raw <- lapply(top_line$x, function(X){
            s <- top_line[[which(top_line$x==X), "y"]]
            e <- med_line[[which(med_line$x==X), "y"]]
            xn <- (X-1)*nr
            return(c(s:e)+xn)
            
          }) %>% c(recursive=T)#Reduce(function(x,y)c(x,y),.)
          
          all_vox_bottom_raw <- lapply(bottom_line$x, function(X){
            s <- bottom_line[[which(bottom_line$x==X), "y"]]
            e <- med_line[[which(med_line$x==X), "y"]]
            xn <- (X-1)*nr
            return(c(s:e)+xn)
            
          }) %>% c(recursive=T)#%>% Reduce(function(x,y)c(x,y),.)
#t6 <- Sys.time()   
          if(direction=="h"){
            all_vox_top <- all_vox_top_raw
            all_vox_bottom <- all_vox_bottom_raw
          } else {
            all_vox_top <- sapply(all_vox_top_raw, retransform_index, nr=nr, nr_orig=nr_orig) 
            all_vox_bottom <- sapply(all_vox_bottom_raw, retransform_index, nr=nr, nr_orig=nr_orig) 
          }
#t7 <- Sys.time()      
          final_list[[n]] <- med_line
          segment_list[[2*n-1]] <- all_vox_top
          segment_list[[2*n]] <- all_vox_bottom
          # time_list[[n]] <- tibble(t01=t1-t0,
          #                          t12=t2-t1,
          #                          t23=t3-t2,
          #                          t34=t4-t3,
          #                          t45=t5-t4,
          #                          t56=t6-t5,
          #                          t67=t7-t6,
          #                          n_seg=paste(2*n-1, "&", 2*n),
          #                          dendrite=nELD)
        }
    #return(list(segment_list, time_list))
    return(segment_list)  
    }) # dendrite
    
    # times <- lapply(dendrite_segments, last) %>% lapply(bind_rows) %>% bind_rows() %>% 
    #   mutate_at(.vars=c(names(.)[1:7]),.funs=round, digits=3)
    
    #res <- lapply(dendrite_segments, first)
    res <- dendrite_segments
    
    tl <- c(res[[1]], res[[2]], res[[3]], res[[4]])
    vl <- c(res[[1]], res[[3]])
    hl <- c(res[[2]], res[[4]])
    normalize_regions(tl, 
                      full_image, 
                      "f:/data_sholl_analysis/test/spec_segs/all.tif", 
                       soma_reg)
    normalize_regions(vl, 
                      full_image, 
                      "f:/data_sholl_analysis/test/spec_segs/vertical.tif", 
                      soma_reg)
    normalize_regions(hl, 
                      full_image, 
                      "f:/data_sholl_analysis/test/spec_segs/horizontal.tif", 
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


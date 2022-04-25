library(imager)
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
  
  print("soma detected")
  SOMA <- somata[1,]
  
#  apply(somata,1, function(SOMA){
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
      elongate_3d_sphere(h_angle, 0, SOMA[["x"]], SOMA[["y"]], SOMA[["z"]], round(sqrt(nc_orig^2+nr_orig^2)), full_image) %>%
        select(x,y) %>%
        return()
      
    })      
    
    selection_vector <- rep(c(1:n_main_dendrites), 3) %>% 
      set_names(c((-n_main_dendrites+1):(2*n_main_dendrites)))
    
    nosoma_image <- remove_soma(full_image, soma_reg, SOMA)
    
    nELD <- 2
    #DENDRITE <- main_dendrites[2,]
    elongated_dendrites <- lapply(1:n_main_dendrites, function(nELD){

      DENDRITE <- main_dendrites[nELD,]
      cat(paste("\ndendrite number:",DENDRITE[["dendrite_id"]]))
      
      
      d0 <- define_dendrite_from_raw_image(DENDRITE)
    
      ## es muss für jedes segment die average main dndrite dicke gefunden werden + die z ebene
      
      d1_raw <- segment_along_dendrite(d0, F, 150) 
      
      d1 <- first(d1_raw)
      
      bin_img <- binarise_image(lapply(first(d1), first), nosoma_image, 0.9)
      
      d2 <- lapply(first(d1), fill_segment_to_rectangle)
      
      #all_vecs <- lapply(nth(d1, 2), cart2pol_vec)
      all_vecs <- nth(d1, 2)[1:20]
      test_segments <- d2[1:6]
      soi <- d2[[2]]
      comb_sub <- lapply(test_segments, find_sub_dendrites, all_vecs=NULL)
      
    })
    ## check in fiji:
    export_dendrites(elongated_dendrites, "f:/data_sholl_analysis/test/example_data/dendrites.csv")
    

    ### now select areas for cutoff use
    

    plot(as.cimg(full_image[[30]]))
    with(bind_rows(all_vecs)[19,],
         nfline(theta, rho, col="red"))
    
    
    res <- dendrite_segments
    
    test_segmentation(c(res[[1]], res[[2]], res[[3]], res[[4]]), 
                      full_image, 
                      "f:/data_sholl_analysis/test/spec_segs/testseg.tif", 
                      soma_reg)    
    
    
    normalize_regions(c(res[[1]], res[[2]], res[[3]], res[[4]]), 
                      full_image, 
                      "f:/data_sholl_analysis/test/spec_segs/final_120.tif", 
                      soma_reg)
    
 # }) # soma

#}) # file  


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


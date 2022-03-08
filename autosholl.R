
library(tiff)
library(tidyverse)
library(gridExtra)


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
  somata <- get_somata(soma_xy_detection_cube_radius*2+1, 
                       soma_z_detection_radius, 
                       soma_z_detection_degree_steps, 
                       full_image)
  
  print("soma detected")
  SOMA <- somata[1,]
  apply(somata,1, function(SOMA){
    
    main_dendrites <- find_dendritic_start_sites(SOMA)
    
    apply(main_dendrites, 1, function(DENDRITE){
      
      
      
    })
    
  })
                          

  # return(somata %>%
  #          mutate(img=file %>% str_split("/") %>% unlist() %>% last())
  #          )
print(Sys.time()-start_time)
}) #%>% 
  #bind_rows()


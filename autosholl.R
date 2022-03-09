
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
  somata <- get_somata(soma_xy_detection_cube_radius*2+1, 
                       soma_z_detection_radius, 
                       soma_z_detection_degree_steps, 
                       full_image)
  
  print("soma detected")
  SOMA <- somata[1,]
  
  apply(somata,1, function(SOMA){
    
    main_dendrites_raw <- find_dendritic_start_sites(SOMA) 
    main_dendrites <- main_dendrites_raw[[1]] %>%
      mutate(z=SOMA[["z"]])
    
    control_plot <- main_dendrites_raw[[2]]
    
    DENDRITE <- main_dendrites[5,]
    apply(main_dendrites, 1, function(DENDRITE){
      
      knot_list <- list()
      next_pos <- c(DENDRITE[["x"]],DENDRITE[["y"]],DENDRITE[["z"]], DENDRITE[["maxima"]])
      c <- 1
      while(c<11){
        print(c)
        next_pos <- elongate_dendrite(next_pos) 
        knot_list[[c]] <- next_pos %>% set_names(c("x","y","z", "angle"))
        c <- c+1
      }
      
      
      
      
      
    })
    
  })
                          

  # return(somata %>%
  #          mutate(img=file %>% str_split("/") %>% unlist() %>% last())
  #          )
print(Sys.time()-start_time)
}) #%>% 
  #bind_rows()



SOMA <- somata[1,]
DENDRITE <- main_dendrites[5,]

pos <- knot_list[[1]]


source("/Users/Marco/git_repos/autosholl/fncts.R")
knot_list <- list()
next_pos <- c(x=DENDRITE[["x"]],y=DENDRITE[["y"]],z=DENDRITE[["z"]], h_angle=DENDRITE[["maxima"]], v_angle=0)
knot_list[[1]] <- next_pos
c <- 1
while(c<21){
  print(c)
  st <- Sys.time()
  next_pos <- elongate_dendrite(next_pos, 10, 60)
  c <- c+1 
  knot_list[[c]] <- next_pos %>% set_names(c("x","y","z", "h_angle", "v_angle"))
  print(Sys.time()-st)
}




control_data <- bind_rows(knot_list) %>% 
                         rownames_to_column("n")

p <- control_plot+geom_tile(inherit.aes=F, data=control_data, 
                            aes(x=x, y=y),fill="red")+
  geom_text(inherit.aes=F, data=bind_rows(knot_list) %>% 
              rownames_to_column("n"), aes(x=x, y=y, label=n))




pdf(file = paste0("c:/Users/Marco/Dropbox/Studium/Master/Praktikum_Mueller/",
                  file %>% str_split("/") %>% unlist() %>% last() %>% str_replace(".tif", "_dendrite_elongation.pdf"),
                  ".pdf"),
    width = 40, height=8)
grid.arrange(
  p
)
dev.off()








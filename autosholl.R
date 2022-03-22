
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
      mutate(z=SOMA[["z"]]) %>%
      rownames_to_column("dendrite_id")
    
    control_plot <- main_dendrites_raw[[2]]
    
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
    
    control_data_fiji <- lapply(elongated_dendrites, function(LO){
      
      LO[[1]] <- c(x=round(SOMA[["x"]]), y=round(SOMA[["y"]]))
      
      raw <- LO %>%
        bind_rows %>% select(x,y) %>%
        rownames_to_column("counter") %>%
        mutate(counter=as.numeric(counter)) %>%
        arrange(counter)
      
      raw_desc <- raw %>%
        arrange(desc(counter))
      
      return(bind_rows(raw, raw_desc) %>% select(-counter)
             )
      
    }) %>% bind_rows()
    
  })
  # return(somata %>%
  #          mutate(img=file %>% str_split("/") %>% unlist() %>% last())
  #          )
print(Sys.time()-start_time)
}) #%>% 
  #bind_rows()


control_new <- elongated_dendrites%>% bind_rows()


control_new <- full_vector_list %>% bind_rows()


p <- control_plot+
  geom_tile(inherit.aes=F, data=control_new, 
            aes(x=x, y=y, fill=z), width=10, height=10)


to_exp <- control_new %>% select(x,y) %>% mutate(y=nrow(full_image[[1]])-y)

write_csv(control_data_fiji, 
          file="c:/Users/Marco/Dropbox/Studium/Master/Praktikum_Mueller/dendrites.csv")


# SOMA <- somata[1,]
# DENDRITE <- main_dendrites[5,]
# 
# pos <- knot_list[[2]]
# n_vc <- 10
# n_hc <- 60
# 
# #source("/Users/Marco/git_repos/autosholl/fncts.R")
# knot_list <- list()
# full_vector_list <- list()
# next_pos <- c(x=DENDRITE[["x"]],y=DENDRITE[["y"]],z=DENDRITE[["z"]], h_angle=DENDRITE[["maxima"]], v_angle=0)
# knot_list[[1]] <- next_pos
# knot_list[[2]] <- next_pos
# c <- 2
# while(c<3|sum(knot_list[[c-1]][1:3]==knot_list[[c]][1:3])!=3){
#   print(c-1)
#   st <- Sys.time()
#   raw_dendrite <- elongate_dendrite(next_pos, n_vc, n_hc)
#   next_pos <- raw_dendrite[[1]]
#   full_vector_list[[c]] <- raw_dendrite[[2]]
#   c <- c+1 
#   knot_list[[c]] <- next_pos %>% set_names(c("x","y","z", "h_angle", "v_angle"))
#   print(Sys.time()-st)
# }


control_new <- full_vector_list %>% bind_rows()


p <- control_plot+
  geom_tile(inherit.aes=F, data=control_new, 
            aes(x=x, y=y, fill=z))






# 
# 
# control_data <- bind_rows(knot_list) %>% 
#                          rownames_to_column("n")
# 
# p <- control_plot+geom_tile(inherit.aes=F, data=control_data, 
#                             aes(x=x, y=y),fill="red")+
#   geom_text(inherit.aes=F, data=bind_rows(knot_list) %>% 
#               rownames_to_column("n"), aes(x=x, y=y, label=n))
# 



pdf(file = paste0("c:/Users/Marco/Dropbox/Studium/Master/Praktikum_Mueller/",
                  file %>% str_split("/") %>% unlist() %>% last() %>% str_replace(".tif", "_dendrite_elongation.pdf"),
                  ".pdf"),
    width = 40, height=8)
grid.arrange(
  p
)
dev.off()









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
  main_dendrites <- apply(somata,1, function(SOMA){
    
    
    deg_step <- 0.005
    #GRAD <- 18
    #SOMA <- somata[1,]
    r <- 300
    
    mn = lapply(seq(deg_step,1,deg_step)*360, function(GRAD){
      #print(GRAD)
      
      xr <- round(cos(deg2rad(GRAD))*-1, 10)
      
      yr <- round(sin(deg2rad(GRAD)),10)
      
      z_layer <- SOMA[["z"]]
      
      f <- tibble(n=1:r) %>%
        mutate(x=ceiling(SOMA[["y"]]+(n*xr)),
               y=ceiling(SOMA[["x"]]+(n*yr))) %>%
        rowwise() %>%
        mutate(i=select_intensity(x,y,full_image[[z_layer]]),
               deg=GRAD) %>%
        filter(!is.na(i))
      
      take_until <- f %>%
        filter(i<0.75) %>%
        pull(n) %>%
        min()
      
      l <- f %>% filter(n<take_until) %>% nrow()
      
      # return(c(deg=GRAD, 
      #          n=l
      #          ))

      return(f%>% filter(n<take_until))
      
    }) %>% 
      bind_rows()
    
    print("screened for dendrites")
    
    all_local_extreme <- c(mn$deg-360,mn$deg, mn$deg+360) %>%
      .[which(.>-100&.<460)] %>%
      density(bw=5) %>%
      get_minmax()
    
    ## debugging
    # all_local_extreme_debug <- c(mn$deg-360,mn$deg, mn$deg+360) %>%
    #   .[which(.>-100&.<460)] %>%
    #   density(bw=5)
    # df <- tibble(x=all_local_extreme_debug$x,
    #              y=all_local_extreme_debug$y)
    # p <- ggplot(df,
    #        aes(x=x,y=y))+
    #   geom_line()
    ## finished debugging
    
    minimum <- all_local_extreme  %>% 
      arrange(y) %>% 
      pull(x) %>% 
      .[1:2]
    
    relevant_maxima <- all_local_extreme %>%
      filter(between(x, 0, 360),
             type=="maximum") 
    
    print("maxima detected")
    
    rescored_maxima <- lapply(relevant_maxima$x, function(GRAD){
      
      xr <- round(cos(deg2rad(GRAD))*-1, 10)
      
      yr <- round(sin(deg2rad(GRAD)),10)
      
      z_layer <- SOMA[["z"]]
      
      f <- tibble(n=1:r) %>%
        mutate(x=ceiling(SOMA[["y"]]+(n*xr)),
               y=ceiling(SOMA[["x"]]+(n*yr))) %>%
        rowwise() %>%
        mutate(i=select_intensity(x,y,full_image[[z_layer]]),
               deg=GRAD) %>%
        filter(!is.na(i))
      
      take_until <- f %>%
        filter(i<0.85) %>%
        pull(n) %>%
        min()
      
      l <- f %>% filter(n<take_until) %>% pull(n) %>% max()
      
      return(f%>% filter(n==l) %>% select(x,y) %>% mutate(maxima=GRAD))
      
    }) %>% bind_rows()

    print("dendrites adjusted")
    print(names(mn))
    control_plot <- ggplot(mn, aes(x=y, y=-x, fill=factor(deg, levels=unique(deg))))+
      geom_tile(show.legend = F)+
      geom_tile(data=rescored_maxima, aes(x=y,y=-x), fill="black", height=2, width=2)  
    
    
    pdf(file = paste0("c:/Users/Marco/Dropbox/Studium/Master/Praktikum_Mueller/",
                      file %>% str_split("/") %>% unlist() %>% last() %>% str_replace(".tif", ""),
                      ".pdf"),
        height = 0.1*abs(min(mn$x)-max(mn$x)), width=0.1*abs(min(mn$y)-max(mn$y)))

    grid.arrange(control_plot)
    dev.off()
    
    #return(rescored_maxima)
    
  })
  
  # return(somata %>%
  #          mutate(img=file %>% str_split("/") %>% unlist() %>% last())
  #          )
print(Sys.time()-start_time)
}) #%>% 
  #bind_rows()


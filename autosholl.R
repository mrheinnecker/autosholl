
library(tiff)
library(tidyverse)

source("/Users/Marco/git_repos/autosholl/fncts.R")

#file <- "/Users/Marco/Dropbox/Studium/Master/Praktikum_Mueller/example_data/Tier1_2_3_apical_D.tif"
files <- list.files("f:/data_sholl_analysis/single_soma", pattern = "*", full.names = T)

#file <- "f:/data_sholl_analysis/Tiffs/Apical/Deep/Tier4_5_2_apical_D.tif"

file <- files[1]

results <- lapply(files, function(file){
  print(file)
  data <- readTIFF(file, all=T)
  cube_size <- 31

  xy_soma <- get_somata(cube_size, data)
  
  r <- 100
  
  deg_step <- 0.05
  res <- lapply(1:length(data), function(Z){
    c(z=Z,
    mn= lapply(seq(deg_step,1,deg_step)*360, function(GRAD){
      
      xr <-deg2rad(GRAD) %>% sin
      yr <- sqrt(1-xr^2)
      
      f <- tibble(n=1:r) %>%
        mutate(x=round(xy_soma$y+(n*xr)),
               y=round(xy_soma$x+(n*yr))) %>%
        rowwise() %>%
        mutate(i=select_intensity(x,y,data[[Z]]))
      
      md <- abs(diff(f$i)) %>% max()
    
      return(md)
        
    }) %>% unlist() %>% mean()) %>%
      return()
  }) %>% bind_rows()
    
  
  z_final <- res %>%
    filter(mn>quantile(res$mn, 0.95)) %>%
    pull(z) %>%
    median() %>%
    round()
  
  
  
  
  return(tibble(img=file %>% str_split("/") %>% unlist() %>% last(),
                z=z_final) %>%
           bind_cols(xy_soma %>% select(-z))
                
           )

}) %>% 
  bind_rows()












wss <- lapply(c(2:nrow(filtered)-1), function(i){
  
  return(c(wss=sum(kmeans(filtered , centers=i)$withinss),
           n_clust=i))
  
}) %>% bind_rows() %>%
  bind_rows(tibble(n_clust=1, wss=(nrow(filtered)-1)*sum(apply(filtered,2,var))))





# ggplot(filtered, aes(x,y, fill=sum))+
#   facet_wrap(~z)+
#   geom_tile()






#files <- file
#files <- rep(file, 10)
results <- lapply(files, function(file){
  
  start_time <- Sys.time()  
  print(file)
  ##############################  
  
  
  data <- readTIFF(file, all=T)
  cutoff1 <- 0.99

  cube_size <- 9
  
  (cube_size-1)/2
  
  z_vec <- c(((cube_size-1)/2):(length(data)-(cube_size-1)/2))
  
  ## loop over z axis
  lapply()
  
  
  
  })







  print("loading layers ...")
  full_data <- load_all_layers(data)
  
  
  print("adjusting x-y raster ...")
  max_intensities <- full_data %>%
    filter(intensity>=quantile(full_data$intensity, 0.99)) %>%
    mutate_at(.vars = c("x", "y", "z"),
              .funs = as.numeric)
  
  counts <- max_intensities %>% group_by(z) %>% tally() %>% filter(n>200)%>% nrow()
  input <- max_intensities
  ratio <- 0
  
  
  while(ratio<0.5*counts/length(data)){
    
      raw_input <- adjust_xy_raster(input)
      input <- raw_input[[1]]
      ratio <- raw_input[[2]]  
   # print(ratio)
  #  print(input %>% summarize(xmin=min(x),
                              # xmax=max(x),
                              # ymin=min(y),
                              # ymax=max(y)))
  }
  
  ## check perfromance
  
  ###################
  print(paste(min(input$x), min(input$y),max(input$x)-min(input$x), max(input$y)-min(input$y), sep=", "))
  ###################
  
  ## next adjust z-layer
  print("loading soma region ...")
  soma_region <- load_soma_region(data, input)
  
  int_value_of_soma <- quantile(soma_region$intensity, 0.8)
  
  
  ##############
  
  ## test
  
  ############
  
  
  # hist_all_layers <- lapply(c(1:length(data)), function(LAYER){
  #   
  #   xy <- soma_region %>% 
  #     filter(z==LAYER)
  #   
  #   d <- density(xy$intensity)
  #   
  #   e <- get_minmax(d) %>% arrange(desc(x))
  #   
  #   deri <- tibble(x=d$x,
  #                  y=d$y) %>%
  #     mutate(z=LAYER) 
  # 
  #   
  #   # if(is.bimodal(deri$y)){
  #   #   return(deri)
  #   # } else {
  #   #   return(NULL)
  #   # }
  # 
  # }) %>% compact() %>%
  #   bind_rows() %>%
  #   mutate(zz=factor(z, levels = c(1:length(data))))
  # 
  # p <-  ggplot(hist_all_layers %>% filter(z %in% c(31:33)), 
  #              aes(x, y, color=zz))+
  #   geom_line()
  # 
  # pdf(file = "c:/Users/Marco/Dropbox/Studium/Master/Praktikum_Mueller/invest_soma_detection.pdf",
  #     width = 15, height=15)
  # 
  # p
  # dev.off()
  
  
  
  print("estimating z-layer ...")
  funs <- lapply(c(1:length(data)), function(LAYER){
    
    xy <- soma_region %>% 
      filter(z==LAYER)
    
    d <- density(xy$intensity)
    
    deri <- tibble(x=d$x,
                   y=d$y,
                   y1=c(0, diff(d$y))) 
    
    e <- get_minmax(d) %>% arrange(desc(x))
    
    ## split density function at highest minimum:
  
    split_value <- e %>% filter(type=="minimum") %>% pull(x) %>% first()
    
    
    high_values <- deri %>% filter(x>split_value)
    high_integral <- get_are_under_curve(high_values$x, high_values$y)

    
    low_values <- deri %>% filter(x<split_value)
    low_integral <- get_are_under_curve(low_values$x, low_values$y)
    
    
    # area <- tibble(r=e$x[1:2], l=e$x[2:3]) %>%
    #   mutate(dist=abs(r-l))
    # score_dist <- e %>% filter(type=="maximum") %>% pull(x) %>% .[1:2] %>% diff() %>% abs()
    # 
    # score_height <- e %>% mutate(norm_y=y/max(e$y)) %>% pull(norm_y) %>% diff() %>% abs() %>% sum()
    # 
    # score_slope  <- apply(area,1, function(ROW){
    # 
    #    score_slope <- deri %>% filter(between(x, ROW[["l"]], ROW[["r"]])) %>%
    #       pull(y1) %>% abs() %>% max() %>%return()
    #   }) %>% sum()
    # 
    # 
    # score_height_ratio=e %>% filter(type=="maximum") %>% pull(y) %>% .[1:2]
    
    return(c(z=LAYER,
             x_max=e %>% filter(type=="maximum") %>% pull(x) %>% first(),
             hint=high_integral,
             lint=low_integral
             # n_extremes=nrow(e),
             # p_unimodal=dip.test(deri$y)$p.value,
             # score_height_ratio=score_height_ratio[2]/score_height_ratio[1],
             # score_height=score_height, 
             # score_slope=score_slope, 
             # score_dist=score_dist
             )
           )
    #c(mc=bimodality_coefficient(deri$y), z=LAYER) %>% return()
    
  }) %>% bind_rows() #%>% rowwise() %>%
    # mutate(score=sum(score_height, score_dist, score_slope),
    #        diff=abs(hint-lint))
    # 
  
  
  
  
  
  # sel <- funs %>% filter(score>=quantile(.$score, 0.9, na.rm=T)) %>% pull(z) %>% median(na.rm=T) %>% round()
  # 
  # print(paste("HHHHHHHHHHHHHHHHHHHHH", sel))
  
  sel <- funs %>% 
    filter(x_max>int_value_of_soma) %>% 
    arrange(desc(lint)) %>% 
    pull(z) %>%
    .[1:10]
    
  
  
  tibble(file %>% str_split("/") %>% unlist() %>% last(), 
         median(sel), 
         paste(sel, collapse = ", ")) %>% 
    set_names(nm=c("img", "layer", "best_10")) %>%
    bind_cols(tibble(time=Sys.time()-start_time)) %>%
    bind_cols(input %>% summarize(xmin=min(x),
                            xmax=max(x),
                            ymin=min(y),
                            ymax=max(y))) %>%
    return()

}) %>% bind_rows()
# dat <- tibble(x=d$x, 
#               y=d$y, 
#               y1=c(0, diff(d$y)))
#               
#               
#               y2=c(0, diff(c(0,diff(d$y)))))
# 
# ggplot(dat)+
#   geom_line(aes(x=x, y=y), color="red")+
#   geom_line(aes(x=x, y=y1), color="blue")+
#   geom_line(aes(x=x, y=y2))


# sel <- funs %>%
#   filter(!is.na(mx_left)) %>%
#   mutate(side_mn_ratio=mx_left/mn_right,
#          mx_mn_ratio=mx_right/mn_right) %>%
#   filter(mx_mn_ratio>side_mn_ratio) %>%
#   filter(side_mn_ratio>=quantile(.$side_mn_ratio, 0.9)) %>%
#   mutate(final_ratio=mx_mn_ratio/side_mn_ratio) %>%
#   filter(final_ratio==min(.$final_ratio)) %>%
#   pull(z)



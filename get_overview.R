
library(tiff)
library(tidyverse)

source("/Users/Marco/git_repos/autosholl/fncts.R")

file <- "/Users/Marco/Dropbox/Studium/Master/Praktikum_Mueller/example_data/Tier1_2_3_apical_D.tif"
files <- list.files("f:/data_sholl_analysis/single_soma", pattern = "*", full.names = T)

file <- files[4]

data <- readTIFF(file, all=T)
cutoff1 <- 0.99

full_data <- lapply(c(1:length(data)), function(LAYER){
  print(LAYER)
  
  l1 <- data[[LAYER]]
  l1[which(l1<=cutoff1)] <- NA
  
  l2 <- l1 %>% as_tibble() %>%
    rownames_to_column("y") %>%
    gather(raw_x, intensity, -y) %>%
    filter(!is.na(intensity)) %>%
  #l1 %>%
   # filter(intensity>=quantile(l1$intensity, 0.9)) %>%
    mutate(x=str_match(raw_x, "\\d+"),
           z=LAYER) %>%
    select(-raw_x) %>%
    return()
  
}) %>% bind_rows()



max_intensities <- full_data %>%
 # filter(intensity==max(full_data$intensity)) %>%
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
  
}

## next adjust z-layer
soma_region <- lapply(c(1:length(data)), function(LAYER){
  print(LAYER)
  
  l1 <- data[[LAYER]][c(min(input$y):max(input$y)), c(min(input$x):max(input$x))] %>% as_tibble() %>%
    rownames_to_column("y") %>%
    gather(raw_x, intensity, -y) %>%
    mutate(x=str_match(raw_x, "\\d+"),
           z=LAYER) %>%
    select(-raw_x) %>%
    return()
  
}) %>% bind_rows()


funs <- lapply(c(1:length(data)), function(LAYER){
  
  xy <- soma_region %>% 
    filter(z==LAYER)
  
  d <- density(xy$intensity)

  e <- get_minmax(d) %>% arrange(desc(x))
  
  # deri <- tibble(x=d$x, 
  #                y=d$y, 
  #                y1=c(0, diff(d$y)))
  # 
  # area <- tibble(r=e$x[1:2], l=e$x[2:3]) %>%
  #   apply(.,1, function(ROW){
  #     
  #    rl <- deri %>% filter(between(x, ROW[["l"]], ROW[["r"]])) %>%
  #       pull(y1) %>% abs() %>% max() %>%
  #       #tibble(ROW[["l"]], ROW[["r"]],max(rl), min(rl)) %>%
  #       return()
  #     
  #   }) %>% sum() 
  # 
  # return(c(z=LAYER, score=area))
  
  highest_max <- e %>%
    filter(type=="maximum") %>%
    filter(x==max(.$x)) %>%
    pull(y)

  highest_min <- e %>%
    filter(type=="minimum") %>%
    filter(x==max(.$x)) %>%
    pull(y)

  side_max <- e %>%
    filter(type=="maximum") %>%
    filter(x==max(.$x[-which(.$x==max(.$x))])) %>%
    pull(y)

  c(z=LAYER,
    n_events=nrow(e),
    mx_right = highest_max,
    mn_right = highest_min,
    mx_left = side_max
  )%>% return()
  
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


sel <- funs %>%
  filter(!is.na(mx_left)) %>%
  mutate(side_mn_ratio=mx_left/mn_right,
         mx_mn_ratio=mx_right/mn_right) %>%
  filter(mx_mn_ratio>side_mn_ratio) %>%
  filter(side_mn_ratio>=quantile(.$side_mn_ratio, 0.9)) %>%
  mutate(final_ratio=mx_mn_ratio/side_mn_ratio) %>%
  filter(final_ratio==min(.$final_ratio)) %>%
  pull(z)



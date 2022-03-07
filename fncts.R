
get_somata <- function(cube_size, data){
  
  z_vec <- seq(((cube_size+1)/2), (length(data)-(cube_size+1)/2), cube_size)
  y_vec <- seq(((cube_size+1)/2), (nrow(data[[1]])-(cube_size+1)/2), cube_size)
  x_vec <- seq(((cube_size+1)/2), (ncol(data[[1]])-(cube_size+1)/2), cube_size)
  ## loop over z axis
  n_its <- length(data)*length(y_vec)*length(x_vec)
  n_vox_per_cube <- cube_size^3
  
  first_test <- lapply(z_vec, function(Z){
    lapply(y_vec, function(Y){
      lapply(x_vec, function(X){
        zsum <- lapply(c((Z-(cube_size-1)/2):(Z+(cube_size-1)/2)), function(Z2){
          data[[Z2]][c((Y-(cube_size-1)/2):(Y+(cube_size-1)/2)), 
                     c((X-(cube_size-1)/2):(X+(cube_size-1)/2))] %>% 
            sum() %>%
            return()
        }) %>% 
          as.numeric() %>% 
          sum()
        return(c(x=X, y=Y, z=Z, sum=zsum))
      }) %>% 
        bind_rows() %>% 
        return()
    }) %>% 
      bind_rows() %>% 
      return()
  }) %>% 
    bind_rows()
  
  filtered <- first_test %>%
    filter(sum>0.9*n_vox_per_cube, z!=16)
  
  km <- kmeans(filtered%>% select(x,y,z), centers = 1)
  
  res <- tibble(wss=sqrt(km$withinss/nrow(filtered))) %>%
           bind_cols(km$centers)
  return(res)
}


rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}

select_intensity <- function(x,y,mat){

  return(mat[x,y])
}



# 
# load_all_layers <- function(data){
#   
#   lapply(c(1:length(data)), function(LAYER){
#     #print(LAYER)
#     
#     l1 <- data[[LAYER]] 
#     
#     colnames(l1) <- c(1:ncol(l1))
#     rownames(l1) <- c(1:nrow(l1))
#     l1[which(l1<=cutoff1|is.na(l1))] <- 0
#     
#     l2 <- l1[, which(colSums(l1) != 0)]
#     
#     if(0 %in% dim(l2)){
#       return(NULL)
#     } 
#     if(class(l2)=="numeric"){
#       l3 <- tibble(x=as.character(which(colSums(l1) != 0)),intensity=l2, y=names(l2))
#     } else {
#       l3 <- l2[which(rowSums(l2) != 0),]%>% as_tibble(rownames=NA) %>%
#         rownames_to_column("y") %>%
#         gather(raw_x, intensity, -y) %>%
#         mutate(x=str_match(raw_x, "\\d+"))%>%
#         select(-raw_x)
#     }
#     
#     
#     if(0 %in% dim(l3)){
#       return(NULL)
#     } 
#     
#     l3  %>%
#       #l1 %>%
#       # filter(intensity>=quantile(l1$intensity, 0.9)) %>%
#       mutate(z=LAYER)  %>%
#       return()
#     
#   }) %>% 
#     compact() %>%
#     bind_rows() %>%
#     return()
# }

# load_soma_region <- function(data, input){
#   lapply(c(1:length(data)), function(LAYER){ 
#   data[[LAYER]][c(min(input$y):max(input$y)), c(min(input$x):max(input$x))] %>% as_tibble() %>%
#     rownames_to_column("y") %>%
#     gather(raw_x, intensity, -y) %>%
#     mutate(x=str_match(raw_x, "\\d+"),
#            z=LAYER) %>%
#     select(-raw_x) %>%
#     return()
#   
#   }) %>% 
#   bind_rows() %>%
#   return()
#   
# }


#intensity_table <- max_intensities

# 
# adjust_xy_raster <- function(intensity_table){
#   
#   
#   n_groups <- 9
#   factor <- n_groups^(1/2)
#   x_range <- ceiling((max(intensity_table$x)-min(intensity_table$x))/factor)
#   y_range <- ceiling((max(intensity_table$y)-min(intensity_table$y))/factor)
#   
#   
#   grouped_data <- intensity_table %>%
#     mutate(x_group=1+ceiling((x-min(intensity_table$x))/(x_range)),
#            y_group=1+ceiling((y-min(intensity_table$y))/(y_range)))
# 
#     # mutate(x_group=ceiling((x)/(x_range)),
#     #        y_group=ceiling((y)/(y_range)))
# 
#   
# #  print(max(grouped_data$x_group))
# #  print(max(grouped_data$y_group))
#   test_groups <- grouped_data %>%
#     group_by(x_group, y_group) %>%
#     summarize(x_start=min(x),
#               x_end=max(x),
#               y_start=min(y),
#               y_end=max(y),
#               n=length(intensity))
#   
#   
#   x_groups <- test_groups %>%
#     group_by(x_group) %>%
#     summarize(mn=min(x_start),
#               mx=max(x_end)) %>%
#     mutate(x_length=mx-mn) %>%
#     select(-mn, -mx)
#   
#   y_groups <- test_groups %>%
#     group_by(y_group) %>%
#     summarize(mn=min(y_start),
#               mx=max(y_end)) %>%
#     mutate(y_length=mx-mn) %>%
#     select(-mn, -mx)
#   
#  # print(x_groups)
# #  print(y_groups)
#   
#   final <- test_groups %>%
#     left_join(x_groups, by="x_group") %>%
#     left_join(y_groups, by="y_group") %>%
#     mutate(max_vox=length(data)*(x_length+1)*(y_length+1),
#            det_ratio=n/max_vox) %>%
#     filter(det_ratio==max(.$det_ratio)) %>%
#     mutate(x_adj_start=x_start-0.5*x_length,
#            y_adj_start=y_start-0.5*y_length,
#            x_adj_end=x_end+0.5*x_length,
#            y_adj_end=y_end+0.5*y_length) 
#     
#     
#   rat <- round(final$det_ratio,2)
#   # grouped_data %>%
#   #   filter()
#   tab <- intensity_table %>%
#     filter(between(x, final$x_adj_start, final$x_adj_end),
#            between(y, final$y_adj_start, final$y_adj_end)) 
#   
#  # tab %>% summarize(xmin=min(x),
#  #                  xmax=max(x),
#  #                 ymin=min(y),
#  #                ymax=max(y)) %>% print()
#   
#   
#   return(list(tab, rat))
#   
# }
# 
# get_minmax <- function(d){
#   max <- which(diff(sign(diff(d$y))) < 0) + 1
#   min <- which(diff(sign(diff(d$y))) > 0) + 1
#   data.frame(x = d$x[max], y = d$y[max])
#   bind_rows(
#     tibble(x = d$x[max], 
#            y = d$y[max],
#            type="maximum"),
#     tibble(x = d$x[min], 
#            y = d$y[min],
#            type="minimum"),
#   ) %>%
#     return()
# }
# 
# 
# get_are_under_curve <- function(x,y){
#   
# 
#   id <- order(x)
#   
#   sum(diff(x[id])*rollmean(y[id],2)) %>% return()
#   
# }
# 









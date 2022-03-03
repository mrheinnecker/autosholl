

adjust_xy_raster <- function(intensity_table){
  
  
  n_groups <- 9
  factor <- n_groups^(1/2)
  x_range <- ceiling((max(intensity_table$x)-min(intensity_table$x))/factor)
  y_range <- ceiling((max(intensity_table$y)-min(intensity_table$y))/factor)
  
  
  grouped_data <- intensity_table %>%
    mutate(x_group=1+ceiling((x-min(intensity_table$x))/(x_range)),
           y_group=1+ceiling((y-min(intensity_table$y))/(y_range)))

    # mutate(x_group=ceiling((x)/(x_range)),
    #        y_group=ceiling((y)/(y_range)))

  
#  print(max(grouped_data$x_group))
#  print(max(grouped_data$y_group))
  test_groups <- grouped_data %>%
    group_by(x_group, y_group) %>%
    summarize(x_start=min(x),
              x_end=max(x),
              y_start=min(y),
              y_end=max(y),
              n=length(intensity))
  
  
  x_groups <- test_groups %>%
    group_by(x_group) %>%
    summarize(mn=min(x_start),
              mx=max(x_end)) %>%
    mutate(x_length=mx-mn) %>%
    select(-mn, -mx)
  
  y_groups <- test_groups %>%
    group_by(y_group) %>%
    summarize(mn=min(y_start),
              mx=max(y_end)) %>%
    mutate(y_length=mx-mn) %>%
    select(-mn, -mx)
  
 # print(x_groups)
#  print(y_groups)
  
  final <- test_groups %>%
    left_join(x_groups, by="x_group") %>%
    left_join(y_groups, by="y_group") %>%
    mutate(max_vox=length(data)*(x_length+1)*(y_length+1),
           det_ratio=n/max_vox) %>%
    filter(det_ratio==max(.$det_ratio)) %>%
    mutate(x_adj_start=x_start-0.5*x_length,
           y_adj_start=y_start-0.5*y_length,
           x_adj_end=x_end+0.5*x_length,
           y_adj_end=y_end+0.5*y_length) 
    
    
  rat <- round(final$det_ratio,2)
  # grouped_data %>%
  #   filter()
  tab <- intensity_table %>%
    filter(between(x, final$x_adj_start, final$x_adj_end),
           between(y, final$y_adj_start, final$y_adj_end)) 
  
 tab %>% summarize(xmin=min(x),
                  xmax=max(x),
                 ymin=min(y),
                ymax=max(y)) %>% print()
  
  
  return(list(tab, rat))
  
}

get_minmax <- function(d){
  max <- which(diff(sign(diff(d$y))) < 0) + 1
  min <- which(diff(sign(diff(d$y))) > 0) + 1
  data.frame(x = d$x[max], y = d$y[max])
  bind_rows(
    tibble(x = d$x[max], 
           y = d$y[max],
           type="maximum"),
    tibble(x = d$x[min], 
           y = d$y[min],
           type="minimum"),
  ) %>%
    return()
}




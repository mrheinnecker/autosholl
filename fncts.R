
get_somata <- function(cube_size, r, deg_step ,data){
  
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
  
  xy_soma <- tibble(wss=sqrt(km$withinss/nrow(filtered))) %>%
           bind_cols(km$centers)
  
  z_raw <- lapply(1:length(data), function(Z){
    c(z=Z,
      mn= lapply(seq(deg_step,1,deg_step)*360, function(GRAD){
        
        xr <- round(cos(deg2rad(GRAD))*-1, 10)
        
        yr <- round(sin(deg2rad(GRAD)),10)
        
        f <- tibble(n=1:r) %>%
          mutate(x=round(xy_soma$y+(n*xr)),
                 y=round(xy_soma$x+(n*yr))) %>%
          rowwise() %>%
          mutate(i=select_intensity(x,y,data[[Z]]))%>%
          filter(!is.na(i))
        
        md <- abs(diff(f$i)) %>% max()
        
        return(md)
        
      }) %>% unlist() %>% mean()) %>%
      return()
  }) %>% bind_rows()
  
  
  xy_soma$z <- z_raw %>%
    filter(mn>quantile(z_raw$mn, 0.95)) %>%
    pull(z) %>%
    median() %>%
    round()
  
  
  
  
  return(xy_soma)
}


rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}

select_intensity <- function(x,y,mat){
  if(between(x, 1, nrow(mat))&between(y, 1, ncol(mat))){
    return(mat[x,y])
  } else {
    return(NA)
  }
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

find_dendritic_start_sites <- function(SOMA){
  
  
  deg_step <- 0.005
  #GRAD <- 18
  #SOMA <- somata[1,]
  r <- 300
  
  mn = lapply(seq(deg_step,1,deg_step)*360, function(GRAD){
    #print(GRAD)
    
    yr <- round(cos(deg2rad(GRAD))*-1, 10)
    
    xr <- round(sin(deg2rad(GRAD)),10)
    
    z_layer <- SOMA[["z"]]
    
    f <- tibble(n=1:r) %>%
      mutate(y=ceiling(SOMA[["y"]]+(n*yr)),
             x=ceiling(SOMA[["x"]]+(n*xr))) %>%
      rowwise() %>%
      mutate(i=select_intensity(y,x,full_image[[z_layer]]),
             deg=GRAD) %>%
      filter(!is.na(i))
    
    take_until <- f %>%
      filter(i<0.75) %>%
      pull(n) %>%
      min()
    
    #l <- f %>% filter(n<take_until) %>% nrow()
    
    return(f%>% filter(n<take_until))
    
  }) %>% 
    bind_rows()
  
  print("screened for dendrites")
  
  dens_func <- c(mn$deg-360,mn$deg, mn$deg+360) %>%
    .[which(.>-100&.<460)] %>%
    density(bw=5) 
    
  all_local_extreme <- get_minmax(dens_func)
  
  relevant_maxima <- all_local_extreme %>%
    filter(between(x, 0, 360),
           type=="maximum") 
  
  control_plot_density <- tibble(x=dens_func$x,
                                 y=dens_func$y) %>%
    ggplot(aes(x=x, y=y))+
    geom_line()+
    geom_segment(data=relevant_maxima, linetype=5,
                 aes(x=x, xend=x, y=0, yend=y))+
    geom_text(data=relevant_maxima, aes(label=round(x), x=x, y=y),
              vjust=0, color="red")+
    coord_cartesian(xlim=c(0,360))
  print("maxima detected")
  
  rescored_maxima <- lapply(relevant_maxima$x, function(GRAD){
    
    yr <- round(cos(deg2rad(GRAD))*-1, 10)
    
    xr <- round(sin(deg2rad(GRAD)),10)
    
    z_layer <- SOMA[["z"]]
    
    f <- tibble(n=1:r) %>%
      mutate(y=ceiling(SOMA[["y"]]+(n*yr)),
             x=ceiling(SOMA[["x"]]+(n*xr))) %>%
      rowwise() %>%
      mutate(i=select_intensity(y,x,full_image[[z_layer]]),
             deg=GRAD) %>%
      filter(!is.na(i))
    
    take_until <- f %>%
      filter(i<0.85) %>%
      pull(n) %>%
      min()
    
    l <- f %>% filter(n<take_until) %>% pull(n) %>% max()
    
    return(f%>% filter(n==l) %>% select(x,y) %>% mutate(maxima=GRAD))
    
  }) %>% 
    bind_rows() #%>%
    #mutate(tier_id=)
  
  print("dendrites adjusted")
  
  control_plot <- ggplot(mn, aes(x=x, y=y, fill=deg %>% factor(levels=unique(mn$deg))))+
    geom_tile(show.legend = F)+
    geom_tile(inherit.aes=F,
              data=rescored_maxima, aes(x=x, y=y), fill="black", height=2, width=4)+
    geom_text(inherit.aes=F,
              data=rescored_maxima, aes(x=x, y=y, label=round(maxima)), color="white")
  
  
  pdf(file = paste0("c:/Users/Marco/Dropbox/Studium/Master/Praktikum_Mueller/",
                    file %>% str_split("/") %>% unlist() %>% last() %>% str_replace(".tif", ""),
                    ".pdf"),
      width = (0.1*abs(min(mn$x)-max(mn$x))), height=0.2*abs(min(mn$y)-max(mn$y)))
  grid.arrange(
  plot_grid(control_plot,
            control_plot_density,
            align="v",
            ncol=1,
            rel_heights = c(1,1))
  )
  dev.off()
  
  return(rescored_maxima)
  
}








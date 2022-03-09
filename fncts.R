
# elongate_3d_sphare_old <- function(GRAD, z_layer, r){
#   
#   yr <- round(cos(deg2rad(GRAD))*-1, 10)
#   
#   xr <- round(sin(deg2rad(GRAD)),10)
#   
#   tibble(n=1:r) %>%
#     mutate(y=ceiling(SOMA[["y"]]+(n*yr)),
#            x=ceiling(SOMA[["x"]]+(n*xr))) %>%
#     rowwise() %>%
#     mutate(i=select_intensity(y,x,full_image[[z_layer]]),
#            deg=GRAD) %>%
#     filter(!is.na(i)) %>%
#     return()
#   
# }
elongate_3d_sphere <- function(GRAD, zGRAD, x_start, y_start, z_start, r){

  zr <- round(cos(deg2rad(zGRAD)), 10)

  zo <- round(sin(deg2rad(zGRAD)),10)
  
  yr <- round(sin(deg2rad(GRAD)), 10)

  xr <- round(cos(deg2rad(GRAD)),10)

  tibble(n=1:r) %>%
    mutate(y=ceiling(y_start+(n*yr*zr)),
           x=ceiling(x_start+(n*xr*zr)),
           z=ceiling(z_start+(n*zo))) %>%
    rowwise() %>%
    mutate(i=select_intensity(x,y,z,full_image),
           deg=GRAD) %>%
    filter(!is.na(i)) %>%
    return()

}

#GRAD <- most_likely_elongation 
#zGRAD <- vertical_angle

elongate_3d_sphere_unlim <- function(GRAD, zGRAD, x_start, y_start, z_start, cutoff){
  
  steps <- 10
  
  steps_to_skip <- 12
  
  zr <- round(cos(deg2rad(zGRAD)), 10)
  
  zo <- round(sin(deg2rad(zGRAD)),10)
  
  yr <- round(sin(deg2rad(GRAD)), 10)
  
  xr <- round(cos(deg2rad(GRAD)),10)
  
  # tibble(n=1:r) %>%
  #   mutate(y=ceiling(y_start+(n*yr*zr)),
  #          x=ceiling(x_start+(n*xr*zr)),
  #          z=ceiling(z_start+(n*zo))) %>%
  #   rowwise() %>%
  #   mutate(i=select_intensity(y,x,z,full_image),
  #          deg=GRAD) %>%
  #   filter(!is.na(i)) %>%
  #   return()
  n <- steps
  i <- select_intensity(x_start,
                        y_start,
                        z_start,
                        full_image)
  
  xp <- x_start
  yp <- y_start
  zp <- z_start
  int_list <- lapply(1:12, list)
  co <- length(int_list)+1
  res_list <- list()
  
  while(median(unlist(int_list[c((co-12):co)]), na.rm=T)>cutoff&
        median(c(i,unlist(int_list)), na.rm=T)>cutoff&
        !is.na(i)){
    # print(round(c(x,
    #         y,
    #         z,
    #         i,
    #         median(unlist(int_list[c((co-12):co)]), na.rm=T),
    #         median(c(i,unlist(int_list)), na.rm=T)),3))
    
    
    x <- ceiling(x_start+(n*xr*zr))
    y <- ceiling(y_start+(n*yr*zr))
    z <- ceiling(z_start+(n*zo))
    i=select_intensity(x+seq(-2,2,1),
                       y+seq(-2,2,1),
                       z+seq(-1,1,1),
                       full_image)
    xp=x
    yp=y
    zp=z
    int_list[[co]] <- i
    res_list[[co]] <- c(x,y,z)
    n <- n+steps
    co <- co+1
  }
  
  if(co-steps_to_skip-1<=steps_to_skip){
    crds <- c(res_list[[steps_to_skip+1]], GRAD)
  } else {
    crds <- c(res_list[[co-steps_to_skip-1]], GRAD)
  }
  
  return(list(dens=tibble(n=c(1:n), deg=GRAD),
              coords=crds))
  
}




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
        
        f <- elongate_3d_sphere(GRAD, 0, xy_soma$x, xy_soma$y ,Z, 100)
        
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

select_intensity <- function(x,y,z,full_image){
  #print(between(z, 1, length(full_image)))
  
  if(between(min(x), 1, ncol(full_image[[1]]))&between(max(x), 1, ncol(full_image[[1]]))&
     between(min(y), 1, nrow(full_image[[1]]))&between(max(y), 1, nrow(full_image[[1]]))&
     between(min(z), 1, length(full_image))&between(max(z), 1, length(full_image))){
    #print(1)
    
    # if(between(x, 1, ncol(full_image[[1]]))&
    #    between(y, 1, nrow(full_image[[1]]))&
    #    between(z, 1, length(full_image))){  
    
   # return(full_image[[z]][y,x])
    return(lapply(z, function(Z){return(full_image[[Z]][y,x])}) %>% unlist() %>% median())
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
  r <- 300
  
  mn = lapply(seq(deg_step,1,deg_step)*360, function(GRAD){
    
    f <- elongate_3d_sphere(GRAD,0,SOMA[["x"]],SOMA[["y"]], SOMA[["z"]], 300)
    
    take_until <- f %>%
      filter(i<0.75) %>%
      pull(n) %>%
      min()
    
    return(f%>% filter(n<take_until))
    
  }) %>% 
    bind_rows()
  
  print("screened for dendrites")
  
  dens_func <- c(mn$deg-360,mn$deg, mn$deg+360) %>%
    .[which(.>-200&.<540)] %>%
    density(bw=5) 
    
  all_local_extreme <- get_minmax(dens_func)
  
  lowest_minimum_x <- all_local_extreme %>% filter(y==min(all_local_extreme$y)) %>% pull(x)
  if(lowest_minimum_x<180){
    relevant_maxima <- all_local_extreme %>%
    filter(between(x, lowest_minimum_x, lowest_minimum_x+360),
           type=="maximum") %>%
      mutate(x=ifelse(x<0, x+360, x))
  } else {
    relevant_maxima <- all_local_extreme %>%
      filter(between(x, lowest_minimum_x-360, lowest_minimum_x),
             type=="maximum")  %>%
      mutate(x=ifelse(x<0, x+360, x))  
    
  }
   
  
  control_plot_density <- tibble(x=dens_func$x,
                                 y=dens_func$y) %>%
    ggplot(aes(x=x, y=y))+
    geom_line()+
    geom_segment(data=relevant_maxima, linetype=5,
                 aes(x=x, xend=x, y=0, yend=y))+
    geom_text(data=relevant_maxima, aes(label=round(x), x=x, y=y),
              vjust=0, color="red")+
    geom_segment(data=tibble(x=c(0, 360)),
                 inherit.aes=F,
                 aes(x=x, xend=x), y=0, yend=1, color="red")
    #coord_cartesian(xlim=c(0,360))
  print("maxima detected")
  
  rescored_maxima <- lapply(relevant_maxima$x, function(GRAD){
    
    f <- elongate_3d_sphere(GRAD,0, SOMA[["x"]],SOMA[["y"]], SOMA[["z"]], 300)
    
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
  
  control_plot <- ggplot(mn, aes(x=x, y=y, 
                                 fill=SOMA[["z"]]
                                 #fill=deg %>% factor(levels=unique(mn$deg))
                                 ))+
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
  
  return(list(rescored_maxima, control_plot))
  
}



elongate_dendrite <- function(pos, n_vc, n_hc){
  cutoff <- 0.5
  horizontal_detection_angle <- 60 
  vertical_detection_angle <- 12
  
  horizontal_angle_input <- pos[4]
  vertical_angle_input <- pos[5]
  
  print(pos)
  # print(paste("horizontal detection range:", 
  #             round(horizontal_angle_input-0.5*horizontal_detection_angle, 2), 
  #             "to", 
  #             round(horizontal_angle_input+0.5*horizontal_detection_angle, 2)))
  # dist_from_soma_center <- 
  #   sqrt(abs(DENDRITE[["x"]]-SOMA[["x"]])^2+
  #          abs(DENDRITE[["y"]]-SOMA[["y"]])^2)

  x_start <- pos[1]
  y_start <- pos[2]
  z_start <- pos[3]  
  horizontal_screening_range <- seq(horizontal_angle_input-0.5*horizontal_detection_angle,
                 horizontal_angle_input+0.5*horizontal_detection_angle,
                 horizontal_detection_angle/n_hc)
  
  vertical_screening_range <- seq(vertical_angle_input-0.5*vertical_detection_angle,
                                  vertical_angle_input+0.5*vertical_detection_angle,
                                  vertical_detection_angle/n_vc)
  
  ft <- lapply(vertical_screening_range, function(zGRAD){
    #print(zGRAD)
    full_screen <- lapply(horizontal_screening_range, function(GRAD){
      #print(GRAD)
      f <- elongate_3d_sphere_unlim(GRAD, zGRAD, x_start, y_start, z_start, cutoff)  
      
      return(f[["dens"]])
      
    }) %>% bind_rows() %>% #unlist() %>% tibble(n=.) %>%
      mutate(zdeg=zGRAD) %>% 
      return()
    
  }) %>% 
    bind_rows()
 # print(nrow(ft))
  # control_3d_density <- lapply(unique(ft$zdeg), function(zGRAD){
  # 
  #   dens_func <- ft  %>% filter(zdeg==zGRAD) %>% pull(n) %>%
  #     density(bw=5)
  # 
  #   return(tibble(x=dens_func$x, y=dens_func$y, vd=zGRAD))
  # 
  # }) %>% bind_rows() %>%   ggplot(aes(x=as.numeric(x),
  #                                     y=as.numeric(y),
  #                                     color=as.character(vd)))+
  #   geom_line()

  dens_func <- c(ft$deg) %>%
    density(bw=5) 
  
  ### check densitxy func
  # ggplot(tibble(x=dens_func$x, y=dens_func$y),
  #      aes(x=as.numeric(x), y=as.numeric(y)))+
  # geom_line()
  ###
  
  all_local_extreme <- get_minmax(dens_func)
  
  most_likely_elongation <- all_local_extreme %>%
    filter(y==max(all_local_extreme$y)) %>%
    pull(x) %>% mean()
  
  
  adjust_z <- ft %>%
    mutate(diff=abs(deg-most_likely_elongation)) %>%
    filter(diff==min(.$diff))%>%
    pull(zdeg) %>%
    density(bw=3)
  
  vertical_angle <- get_minmax(adjust_z) %>%
    filter(y==max(.$y)) %>% pull(x) %>% mean()
  
 # ggplot(tibble(x=adjust_z$x, y=adjust_z$y), aes(x,y))+
  #  geom_line()
  f <- elongate_3d_sphere_unlim(most_likely_elongation, 
                                vertical_angle, 
                                x_start, 
                                y_start, 
                                z_start, 
                                cutoff) 
  #f <- elongate_3d_sphere(most_likely_elongation$x, 
   #                       vertical_angle, 
    #                      x_start, y_start, z_start,150) 
  fin <- f[["coords"]]
  # take_until <- f %>%
  #   filter(i<0.55) %>%
  #   pull(n) %>%
  #   min()
  # 
  # fin <- f%>% filter(n<take_until) %>%
  #   filter(n==max(.$n)-4) 
  # print(median(summed_raw$zdeg))
  # print(nrow(summed_raw))
  #print(summed_raw)
  # summed <- summed_raw[ceiling(nrow(summed_raw)/2),]
  # #print(nrow(summed))
  # #print(summed$zdeg)
  # 
  # fin <- ft %>% filter(zdeg==summed$zdeg, deg==summed$deg, n==summed$n) 
  # print(fin)
  # if(pos[2]==fin$y){
  #   new_angle <- pos[4]
  # } else {
  #   new_angle <- pos[4]-tan((pos[1]-fin$x)/(pos[2]-fin$y))
  # }
  
  #new_angle <- fin$deg
  return(c(fin, vertical_angle))
  
}





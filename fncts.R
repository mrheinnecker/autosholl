
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

elongate_3d_sphere_unlim <- function(GRAD, 
                                     zGRAD, 
                                     x_start, 
                                     y_start, 
                                     z_start, 
                                     intensity_cutoff,
                                     steps
                                     ){
  

  ## trigonometry
  zr <- round(cos(deg2rad(zGRAD)), 10)
  
  zo <- round(sin(deg2rad(zGRAD)),10)
  
  yr <- round(sin(deg2rad(GRAD)), 10)
  
  xr <- round(cos(deg2rad(GRAD)),10)
  
  ## sorrounding voxels 
  xs <- seq(-2,2,1)
  ys <- xs
  zs <- c(-1,0,1)
  
  ## intensity_cutoffs and pre assignments
  
  #steps <- 10
  
  steps_to_skip <- 12
  
  n <- steps
  i <- 1
  

  int_list <- list(i,i,i,i,i)
  #co <- length(int_list)+1
  co <- 1
  l <- co
  res_list <- list(c(x_start, y_start, z_start))
  

    
  while(!is.na(i)&
        (int_list[[co]]>intensity_cutoff|co>2)&
        (co<4|mean(unlist(int_list[(length(int_list)-4):length(int_list)]))>intensity_cutoff)
         ) { 
    #print(co)
    co <- co+1
    # print(round(c(x,
    #         y,
    #         z,
    #         i,
    #         median(unlist(int_list[c((co-12):co)]), na.rm=T),
    #         median(c(i,unlist(int_list)), na.rm=T)),3))
    
    
    x <- ceiling(x_start+(n*xr*zr))
    y <- ceiling(y_start+(n*yr*zr))
    z <- ceiling(z_start+(n*zo))
    i=select_intensity(x+xs,
                       y+ys,
                       z+zs,
                       full_image)
    int_list[[co]] <- i
    res_list[[co]] <- c(x,y,z,i=i)
    n <- n+steps
    #l <- ifelse(isTRUE(i>intensity_cutoff|!is.na(i)),co,l)
    l <- ifelse(isTRUE(i>intensity_cutoff), co, l)
  }
  
  if(co<=2){
    #crds <- c(res_list[[1]][1:3], GRAD)
    return(NULL)
    
  } else {
    return(list(dens=tibble(n=c(1:n), deg=GRAD),
              coords=c(res_list[[l]][1:3], GRAD),
              full=res_list[c(2:l)]
           #full=1
           ))
    
    #crds <-  c(res_list[[l]][1:3], GRAD)
  }
  
  
  
}

screen_for_dendrite_elongation <- function(pos, 
                              vertical_subdivisions, 
                              horizontal_subdivisions,
                              horizontal_detection_angle,
                              vertical_detection_angle,
                              intensity_cutoff,
                              steps){
  
  horizontal_angle_input <- as.numeric(pos[4])
  vertical_angle_input <- as.numeric(pos[5])
  
  #print(pos)

  x_start <- pos[1]
  y_start <- pos[2]
  z_start <- pos[3]  
  horizontal_screening_range <- seq(horizontal_angle_input-0.5*horizontal_detection_angle,
                 horizontal_angle_input+0.5*horizontal_detection_angle,
                 horizontal_detection_angle/horizontal_subdivisions)
  
  vertical_screening_range <- seq(vertical_angle_input-0.5*vertical_detection_angle,
                                  vertical_angle_input+0.5*vertical_detection_angle,
                                  vertical_detection_angle/vertical_subdivisions)
  
  ft <- lapply(vertical_screening_range, function(zGRAD){
    #print(zGRAD)
    full_screen <- lapply(horizontal_screening_range, function(GRAD){
      #print(GRAD)
      f <- elongate_3d_sphere_unlim(GRAD, zGRAD, x_start, y_start, z_start, intensity_cutoff, steps)  
      
      return(f[["dens"]])
      
    }) %>% compact() %>% bind_rows() %>% #unlist() %>% tibble(n=.) %>%
      mutate(zdeg=zGRAD) %>% 
      return()
    
  }) %>% 
    bind_rows()

  if(nrow(ft)==0){
    return(list(pos, list()))
  }
  
  dens_func <- c(ft$deg) %>%
    density(bw=5) 
  
  ### check densitxy func
 # ft %>% group_by(deg, zdeg) %>% tally() %>% View
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
                                intensity_cutoff,
                                steps) 
  fin <- f[["coords"]]
  return(list(c(fin, vertical_angle), f[["full"]]))
  
}



elongate_dendrite <- function(DENDRITE,
                  horizontal_subdivisions,
                  vertical_subdivisions,
                  horizontal_detection_angle,
                  vertical_detection_angle,
                  intensity_cutoff,
                  steps){

  knot_list <- list()
  full_vector_list <- list()
  next_pos <- c(x=DENDRITE[["x"]],
                y=DENDRITE[["y"]],
                z=DENDRITE[["z"]], 
                h_angle=DENDRITE[["h_angle"]], 
                v_angle=0)
  knot_list[[1]] <- next_pos
  knot_list[[2]] <- next_pos
  c <- 2
  while(c<3|sum(knot_list[[c-1]][1:3]==knot_list[[c]][1:3])!=3){
    cat(paste("\n  elongation step:", c-1))
    st <- Sys.time()
    raw_dendrite <- screen_for_dendrite_elongation(next_pos, 
                                                   vertical_subdivisions, 
                                                   horizontal_subdivisions, 
                                                   horizontal_detection_angle,
                                                   vertical_detection_angle,
                                                   intensity_cutoff,
                                                   steps)
    next_pos <- raw_dendrite[[1]]
    #full_vector_list[[c]] <- raw_dendrite[[2]]
    c <- c+1 
    knot_list[[c]] <- next_pos %>% set_names(c("x","y","z", "h_angle", "v_angle"))
    #print(Sys.time()-st)
  }
  return(knot_list[3:c-1])


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
  
  soma_radius <- mn %>%
    group_by(deg) %>%
    summarize(mx_raw=max(n)) %>%
    mutate(mx=ifelse(mx_raw>quantile(.$mx_raw, 0.95),
                     round(quantile(.$mx_raw, 0.95)),
                     round(mx_raw))) %>%
    pull(mx) %>%
    mean() %>%
    round()
  
  dens_func <- c(mn$deg-360,mn$deg, mn$deg+360) %>%
    .[which(.>-200&.<540)] %>%
    density(bw=6) 
    
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
    
    f <- elongate_3d_sphere(GRAD,0, SOMA[["x"]],SOMA[["y"]], SOMA[["z"]], soma_radius)
    
    # take_until <- f %>%
    #   filter(i<0.85) %>%
    #   pull(n) %>%
    #   min()
    # 
    # l <- f %>% filter(n<take_until) %>% pull(n) %>% max()
    # 
    # return(f%>% filter(n==l) %>% select(x,y) %>% mutate(maxima=GRAD))
    return(f%>% filter(n==soma_radius) %>% select(x,y) %>% mutate(h_angle=GRAD))
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
              data=rescored_maxima, aes(x=x, y=y, label=round(h_angle)), color="white")
  
  
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
  
  return(list(rescored_maxima, control_plot, soma_radius))
  
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


rad2deg <- function(rad){(rad * 180) / (pi)}
deg2rad <- function(deg){(deg * pi) / (180)}


#sl <- segment_list
#file_name <- "f:/data_sholl_analysis/test/spec_segs/ggg.tif"
normalize_regions <- function(sl, full_image, file_name, soma_reg){
  QNT <- 0.9
  ## remove soma 
  nosoma_image <- lapply(1:length(full_image), function(n){
    diff <- abs(SOMA[["z"]]-n)
    LAYER <- full_image[[n]]
    
    LAYER[c((soma_reg$y-diff*2):(soma_reg$yend+diff*2)), c((soma_reg$x-diff*2):(soma_reg$xend+diff*2))] <- NA
    return(LAYER)
    
  })
  
  
    
  full_cutoffs <- lapply(sl, function(VOX){
    
    cutoff <- lapply(nosoma_image, function(LAYER){
      
      return(LAYER[VOX])
      
    }) %>%
      Reduce(function(x,y)c(x,y),.) %>%
      quantile(QNT, na.rm=T) %>%
      return()
    
  })
  
  new_image <- nosoma_image
  
  
  for(i in 1:length(full_cutoffs)){
    cat(paste("\n  ", i))
    cutoff <- full_cutoffs[[i]]
    VOX <- sl[[i]]
    for(l in 1:length(new_image)){
      #print(l)
      new_image[[l]][VOX][new_image[[l]][VOX]>=cutoff] <- 1
      new_image[[l]][VOX][new_image[[l]][VOX]<cutoff] <- 0
      
      new_image[[l]] <- matrix(new_image[[l]], nrow = 1040)
    }
  
  }
  
  writeTIFF(new_image, file_name)
  
}



export_dendrites <- function(elongated_dendrites, file_name){
  
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
  write_csv(control_data_fiji, 
            file=file_name)
  
  
}


create_df_of_vectors <- function(ELD){
  lapply(c(1:(length(ELD)-1)), function(n){
    
    c(
      xs=ELD[[n]][["x"]],
      xe=ELD[[n+1]][["x"]],
      ys=ELD[[n]][["y"]],
      ye=ELD[[n+1]][["y"]],
      zs=ELD[[n]][["z"]],
      ze=ELD[[n+1]][["z"]]
    ) %>%
      return()
    
  }) %>% 
    bind_rows() %>% 
    rownames_to_column("id") %>%
    return()
}

create_rv_of_vectors <- function(ELD){
  lapply(1:length(ELD), function(PT){
    if(PT==length(ELD)){return(NULL)}
    return(c(ELD[[PT+1]][c("x", "y", "z")]-ELD[[PT]][c("x", "y", "z")],
             h_angle=ELD[[PT+1]]["h_angle"],
             l=sqrt(sum(abs(ELD[[PT+1]][c("x", "y")]-ELD[[PT]][c("x", "y")])^2))))
    
  }) %>%
    compact() %>%
    return()
}

assign_vectors_to_segments <- function(ELD, vector_pos, use_length, n_segments){
  
  res_list <- list()
  
  for(n in 1:n_segments){
    start <-(n-1)*use_length+ELD[[1]]["x"]
    end <-   n*use_length+ELD[[1]]["x"]
    vecs <- vector_pos %>% rowwise() %>%
      #mutate(t=ifelse(xs %in% c(start:end)|xe %in% c(start:end), T, F)) %>%
      mutate(t=ifelse(length(intersect(start:end, xs:xe))>0, T, F)) %>%
      filter(t==T) %>%
      pull(id) 
    if(length(vecs)==0){
      res_list[[n]] <- as.numeric(res_list[[n-1]])
    } else {
      res_list[[n]] <- as.numeric(vecs)
    }
  }
  return(res_list)
}









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

select_intensity <- function(x,y,z,img){
  #print(between(z, 1, length(img)))
  
  if(between(min(x), 1, ncol(img[[1]]))&between(max(x), 1, ncol(img[[1]]))&
     between(min(y), 1, nrow(img[[1]]))&between(max(y), 1, nrow(img[[1]]))&
     between(min(z), 1, length(img))&between(max(z), 1, length(img))){
    #print(1)
    
    # if(between(x, 1, ncol(img[[1]]))&
    #    between(y, 1, nrow(img[[1]]))&
    #    between(z, 1, length(img))){  
    
   # return(img[[z]][y,x])
    return(lapply(z, function(Z){return(img[[Z]][y,x])}) %>% unlist() %>% median())
  } else {
    return(NA)
  }
}


elongate_3d_sphere <- function(GRAD, zGRAD, x_start, y_start, z_start, r, image){

  zr <- round(cos(deg2rad(zGRAD)), 10)

  zo <- round(sin(deg2rad(zGRAD)),10)
  
  yr <- round(sin(deg2rad(GRAD)), 10)

  xr <- round(cos(deg2rad(GRAD)),10)

  tibble(n=1:r) %>%
    mutate(y=ceiling(y_start+(n*yr*zr)),
           x=ceiling(x_start+(n*xr*zr)),
           z=ceiling(z_start+(n*zo))) %>%
    rowwise() %>%
    mutate(i=select_intensity(x,y,z,image),
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
    
    f <- elongate_3d_sphere(GRAD,0,SOMA[["x"]],SOMA[["y"]], SOMA[["z"]], 300, full_image)
    
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
    
    f <- elongate_3d_sphere(GRAD,0, SOMA[["x"]],SOMA[["y"]], SOMA[["z"]], soma_radius, full_image)
    
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
        
        f <- elongate_3d_sphere(GRAD, 0, xy_soma$x, xy_soma$y ,Z, 100, full_image)
        
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


test_segmentation <- function(sl, full_image, file_name, soma_reg){
  QNT <- 0.9
  
  nosoma_image <- list("1"=full_image[[SOMA[["z"]]]])
  
  full_cutoffs <- lapply(sl, function(VOX){
    #print(1)
    cutoff <- lapply(nosoma_image, function(LAYER){
      
      return(LAYER[VOX])
      
    }) %>%
      Reduce(function(x,y)c(x,y),.) %>%
      quantile(QNT, na.rm=T) %>%
      return()
    
  })
  
  new_image <- nosoma_image
  
  
  for(i in 1:length(full_cutoffs)){
    #cat(paste("\n  ", i))
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
    print(1)
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


# get_single_index <- function(x,y,nr){
#   
#   nr*(x-1)+y
#   
# }

# get_xy_index <- function(i, nr){
#   y <- i%%nr
#   if(y==0){
#     y_ret <- nr
#     x <- (i-y)/nr 
#   } else {
#     y_ret <- y
#     x <- (i-y)/nr+1
#   }
#  
#   # original nr*(x-1)+y
#   return(c(x=x, y=y_ret))
#   
# }


get_single_index <- function(x,y,nr){
  
  nr*(x-1)+y
  
}

get_xy_index <- function(i, nr){
  
  y <- i%%nr
  if(y==0){
    y_ret <- nr
    x <- (i-y)/nr 
  } else {
    y_ret <- y
    x <- (i-y)/nr+1
  }
  
  return(c(x=x, y=y_ret))
  
}

retransform_index <- function(i, nr, nr_orig){
  y <- i%%nr   
  if(y==0){
    return(nr_orig*(nr-1)+(i-y)/nr)
  } else {
    return(nr_orig*(y-1)+(i-y)/nr+1)
  }
  
}



combine_vectors <- function(ELD, xs, xe, rv, full_vecs, xnorm_vector, res_list, n){
  
  rel_vecs <- full_vecs[res_list[[n]]]
  rl <- list()
  last_end <- xs-rv
  for(i in 1:length(rel_vecs)){
    #print(i)
    if(i==1){
      start <- xs
    } else {
      start <- last_end+rv
    } 
    
    if(i==length(rel_vecs)){
      end <- xe 
    } else {
      end <- ELD[[1]]["x"]+rel_vecs[i]
    }
    last_end <- end
    rl[[i]] <- rep(xnorm_vector[[res_list[[n]][i]]]["y"], abs(end-start+rv))
  }
  return(rl)
}

#sdir <- "bottom"
#cof <- bottom
#bd <- bottom_border

create_segment <- function(xs, xe, cof, line, nr, nc, bd, sdir, direction, n_orig){
  if(sdir=="bottom"){
    rob <- bd %>%
      filter(x %in% c(xs:xe)) %>%
      filter(y>cof) %>%
      group_by(x) %>%
      summarize(y=min(y))
  } else {
    rob <- bd %>%
      filter(x %in% c(xs:xe)) %>%
      filter(y<cof) %>%
      group_by(x) %>%
      summarize(y=min(y))
  }
  
  top_line <- tibble(x=c(xs:xe),
                     y=cof) %>%
    filter(!x %in% rob$x) %>%
    bind_rows(rob)%>%
    filter(between(x, 1, nc),
           between(y, 1, nr))
  
  #print(top_line)
  
    all_vox_raw <- lapply(top_line$x, function(X){
      s <- top_line[[which(top_line$x==X), "y"]]
      e <- line[[which(line$x==X), "y"]]
      xn <- (X-1)*nr
      return(c(s:e)+xn)
      
    }) %>% 
      c(recursive=T)
  
  if(direction=="h"){
    return(list(all_vox_raw, top_line))
  } else {
    sapply(all_vox_raw ,retransform_index, nr=nr, nr_orig=nr_orig) %>%
      list(., top_line) %>%
      return()
  }
  
}

define_dendrite_from_raw_image <- function(DENDRITE){
  
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
                  steps) %>%
  return()
  
  
}

segment_along_dendrite <- function(ELD_raw, seg_full_img, avs){
  
  overall_vector_raw <- ELD_raw[[length(ELD_raw)]][c("x", "y", "z")]-ELD_raw[[1]][c("x", "y", "z")]
  ################################################################################
  ## defining secondary elements accroding to vertical or horizonatal dendrite ###    
  ################################################################################      
  if(abs(overall_vector_raw[1])>abs(overall_vector_raw[2])){
    direction <- "h"
    ELD <- ELD_raw
    overall_vector <- overall_vector_raw
    borders_vox <- borders_vox_raw
    nr <- nr_orig
    nc <- nc_orig
  } else {
    direction <- "v"
    ELD <- lapply(ELD_raw, function(VEC){
      return(c(x=VEC[["y"]], y=VEC[["x"]], z=VEC[["z"]],  h_angle=VEC[["h_angle"]],  v_angle=VEC[["v_angle"]]))
    })
    overall_vector <- c(x=overall_vector_raw[["y"]], y=overall_vector_raw[["x"]], 
                        z=overall_vector_raw[["z"]])
    borders_vox <- lapply(borders_vox_raw, function(BORD){
      tibble(x=BORD$y, y=BORD$x) %>%
        return()
    })
    nc <- nr_orig
    nr <- nc_orig
  }
  ################################################################################
  ## defining tertiary elements accroding to v/h and direction (left/right) ######    
  ################################################################################ 
  overall_length <- overall_vector["x"]
  rv <- overall_vector["x"]/abs(overall_vector["x"])      
  
  if((direction=="h"&rv==1)|(direction=="v"&rv==-1)){
    top_border <- borders_vox[[selection_vector[[as.character(nELD)]]]]
    bottom_border <- borders_vox[[selection_vector[[as.character(nELD-1)]]]]
  } else {
    top_border <- borders_vox[[selection_vector[[as.character(nELD-1)]]]]
    bottom_border <- borders_vox[[selection_vector[[as.character(nELD)]]]]
  }
  
  if(rv==1){
    pixels_to_image_border <- nc-ELD[[1]]["x"]
    ELD[[length(ELD)+1]] <- c(x=nc, ELD[[length(ELD)]][c("y", "z")], h_angle=0, v_angle=0)
  } else {
    pixels_to_image_border <- 1-ELD[[1]]["x"]
    ELD[[length(ELD)+1]] <- c(x=1, ELD[[length(ELD)]][c("y", "z")], h_angle=180, v_angle=0)
  }
  ################################################################################
  ################### calculate required data from input #########################    
  ################################################################################  
  vector_pos <- create_df_of_vectors(ELD)
  vectors <- create_rv_of_vectors(ELD)
  support_vectors <- apply(vector_pos, 1, function(VEC){
    return(c(VEC["xs"], VEC["ys"], VEC["zs"]))
  }, simplify=F)
  xnorm_vector <- lapply(vectors, function(V){V[c(1:3)]/abs(V["x"])})
  
  n_segments <- abs(round(pixels_to_image_border/avs))
  use_length <- round(pixels_to_image_border/n_segments)
  
  res_list <- assign_vectors_to_segments(ELD, vector_pos, use_length, n_segments)
  
  full_vecs <- bind_rows(vectors) %>%
    pull(x) %>% cumsum()  
  
  y_coord_list <- list()
  segment_list <- list()
  segment_count <- 0
  for(n in c(1:n_segments)){
    cat(paste("\n  x_segment:",n))
    ################################################################################
    ################### set x and y variables for segment ##########################    
    ################################################################################ 
    xs <- ELD[[1]]["x"]+(n-1)*use_length+1
    xe <- ELD[[1]]["x"]+n*use_length 
    ys <- ifelse(n==1,
                 as.numeric(ELD[[1]]["y"]),
                 as.numeric(y_coord_list[[n-1]]))  
    ## select relevant vectors for that segment 
    rl <- combine_vectors(ELD, xs, xe, rv, full_vecs, xnorm_vector, res_list, n)  
    
    med_line <-   tibble(x=c(xs:xe)) %>%
      mutate(fac=c(rl, recursive=T)) %>%
      mutate(y=ceiling(ys+cumsum(fac))) %>%
      filter(between(x, 1, nc),
             between(y, 1, nr))
    ## from here the segemt height and depth is measured until image border
    med <- round(0.5*(max(med_line$y)+min(med_line$y)))
    ## mean z layer of segment
    mz <- round(0.5*sum(vector_pos[which(vector_pos$id==first(res_list[[n]])),]$zs, 
                  vector_pos[which(vector_pos$id==last(res_list[[n]])),]$ze))
    
    
    pix_to_bottom <- -med
    pix_to_top <- nr-med
    n_segments_top <- abs(round(pix_to_top/avs))
    n_segments_bottom <- abs(round(pix_to_bottom/avs))
    use_length_top <- round(pix_to_top/n_segments_top)
    use_length_bottom <- round(pix_to_bottom/n_segments_bottom)
    
    
    ## we also need to export all relevant vectors for that segment:
    full_relevant_vectors <- lapply(res_list[[n]], function(V){
      list(sv=as.numeric(support_vectors[[V]]),
           rv=vectors[[V]][1:4]) %>% return()
    })
    
    
    ## if we only want one segment from the dendrite to top and bottom
    if(seg_full_img==F){
      n_segments_top <- 1
      n_segments_bottom <- 1
    }
    
    
    
    ## new loop for all segments
    
    start_list_top <- list()
    line_list_top <- list()
    
    for(n_top in 1:n_segments_top){
      
      if(n_top==1){
        start <- med
        line <- med_line
      } else {
        start <- start_list_top[[n_top-1]]+1
        line <- line_list_top[[n_top-1]]
      }
      
      top <- start+use_length_top
      if(top>nr){top <- nr}
      
      all_vox_top <- create_segment(xs, xe, top, line, nr, nc, top_border, "top",
                                    direction, nr_orig)
      
      start_list_top[[n_top]] <- top
      line_list_top[[n_top]] <- all_vox_top[[2]]
      
      segment_count <- segment_count+1
      segment_list[[segment_count]] <- list(all_vox_top[[1]], mz, full_relevant_vectors, 
                                            c(seg_id=segment_count,
                                              dir="top",
                                              rv=rv))
      cat(paste("\n    top:", segment_count))
    }
    
    start_list_bottom <- list()
    line_list_bottom <- list()
    
    for(n_bottom in 1:n_segments_bottom){
      
      if(n_bottom==1){
        start <- med
        line <- med_line
      } else {
        start <- start_list_bottom[[n_bottom-1]]-1
        line <- line_list_bottom[[n_bottom-1]]
      }
      
      bottom <- start+use_length_bottom
      if(bottom<1){bottom <- 1} 
      
      all_vox_bottom <- create_segment(xs, xe, bottom, line, nr, nc, bottom_border, "bottom",
                                       direction, nr_orig)
      
      start_list_bottom[[n_bottom]] <- bottom
      line_list_bottom[[n_bottom]] <- all_vox_bottom[[2]]
      
      segment_count <- segment_count+1
      cat(paste("\n    bottom:", segment_count))
      segment_list[[segment_count]] <- list(all_vox_bottom[[1]], mz, full_relevant_vectors, 
                                            c(seg_id=segment_count,
                                              dir="bottom",
                                              rv=rv))
    }       
    
    y_coord_list[[n]] <- med_line[[nrow(med_line), "y"]]
    
    
  }
  
  ## for further steps also the z layer of the vectors is required
  
  # lapply(res_list, function(n){
  #   
  #   mz <- 0.5*sum(vector_pos[which(vector_pos$id==first(res_list[[n]]))], vector_pos[which(vector_pos$id==last(res_list[[n]]))])
  #   
  # })
  
  all_vectors <- lapply(1:length(vectors), function(V){
    
    list(sv=as.numeric(support_vectors[[V]]),
         rv=vectors[[V]][1:4])
    
  })
  
  return(list(segment_list, all_vectors, vector_pos))  

  
}


remove_soma <- function(full_image, soma_reg, SOMA){
  
  lapply(1:length(full_image), function(n){
    diff <- abs(SOMA[["z"]]-n)
    LAYER <- full_image[[n]]
    
    LAYER[c((soma_reg$y-diff*2):(soma_reg$yend+diff*2)), c((soma_reg$x-diff*2):(soma_reg$xend+diff*2))] <- 0
    return(LAYER)
    
  }) %>%
    return()
  
}



binarise_image <- function(sl, nosoma_image, QNT){
  #QNT <- 0.9
  ## remove soma 
  # nosoma_image <- lapply(1:length(full_image), function(n){
  #   diff <- abs(SOMA[["z"]]-n)
  #   LAYER <- full_image[[n]]
  #   
  #   LAYER[c((soma_reg$y-diff*2):(soma_reg$yend+diff*2)), c((soma_reg$x-diff*2):(soma_reg$xend+diff*2))] <- NA
  #   return(LAYER)
  #   
  # })
  
  
  full_cutoffs <- lapply(sl, function(VOX){
    print(1)
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
  
  writeTIFF(new_image, "f:/data_sholl_analysis/test/spec_segs/bintest.tif")
  return(new_image)
}

VEC <- nth(d1, 2)[[1]]

cart2pol_vec <- function(VEC){
  sv <- VEC[[1]]
  rv <- VEC[[2]]
  #d_p_sv <- c(xmin, ymin)-sv[1:2]
  # dist p zu SV
  d_p_sv <- sv[1:2]-c(1,1)
  
  cp <- d_p_sv[[1]]*rv[[2]]-d_p_sv[[2]]*rv[[1]]
  
  # betrag des RV
  d_rv <- sqrt(rv[[1]]^2+rv[[2]]^2)
  
  #rho <- abs(cp/d_rv)
  rho <- cp/d_rv
  
  if(between(rv[[4]],0, 90)){
    A <- 90
  } else if(between(rv[[4]],90, 180)){
    A <- 180
  } else if(between(rv[[4]],180, 270)){
    A <- 180
  } else {
    A <- 90
  }
  
  
  
  theta <- deg2rad(A-rv[[4]])
  
  x_coords <- c(sv[[1]], sv[[1]]+rv[[1]])
  
  return(c(rho=rho, theta=theta, xmin=min(x_coords), xmax=max(x_coords)))
}



#inp_raw <- d1[[1]]
fill_segment_to_rectangle <- function(inp){
  #inp <- first(inp_raw)
  #all_vecs <- last(inp_raw)
  segment <- first(inp)
  
  df <- lapply(segment, get_xy_index, nr=1040)
  
  xmax <- sapply(df, first) %>% max()
  xmin <- sapply(df, first) %>% min()
  ymax <- sapply(df, last) %>% max()
  ymin <- sapply(df, last) %>% min()
  
  
  ## transform vectors of segment to polar coordinates
  all_vecs <- lapply(nth(inp,3), function(VEC){
    sv <- VEC[[1]]
    rv <- VEC[[2]]

    d_p_sv <- c(xmin, ymin)-sv[1:2]

    cp <- abs(d_p_sv[[1]]*rv[[2]]-d_p_sv[[2]]*rv[[1]])

    d_rv <- sqrt(rv[[1]]^2+rv[[2]]^2)

    rho <- abs(cp/d_rv)

    theta <- deg2rad(180-rv[[4]])

    x_coords <- c(sv[[1]], sv[[1]]+rv[[1]])

     return(c(rho=rho, theta=theta, xmin=min(x_coords), xmax=max(x_coords)))
  })
  
  
  all_vox <- lapply(c(xmin:xmax), function(X){
    sapply(ymin:ymax, function(Y){
      get_single_index(X, Y, 1040)
    })
  }) %>% c(recursive=T)
  
  list(voi=segment,
       fv=all_vox[-which(all_vox %in% segment)],
       nr=ymax-ymin+1,
       nc=xmax-xmin+1,
       origin=c(x=xmin, y=ymin),
       xmin=xmin,
       xmax=xmax,
       ymin=ymin,
       ymax=ymax,
       z=nth(inp, 2),
       main_vecs=all_vecs
       ) %>% c(., inp[[4]]) %>%
    return()
  
}

pol2car <- function(r, theta){
  
  b=r/sin(theta)
  a=(-cos(theta))/sin(theta)
  return(c(a=a,b=b))
}

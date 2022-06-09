load_dependencies <- function(){
   
  library(tiff)
  library(tidyverse)
  library(gridExtra)
  library(cowplot) 
  library(Rcpp)
  library(imager)
  
}



set_options <- function(raw_image){

  
  
  list(
    ## overall info
    nr_orig = nrow(raw_image[[1]]),
    nc_orig = ncol(raw_image[[1]]),
    
    ## soma detection
    soma_xy_detection_cube_radius = 15,
    soma_z_detection_radius = 100,
    soma_z_detection_degree_steps = 0.05,
    
    ## sccreening subdendrite starts
    subd_detection_distance = 30,
    subd_detection_depth = 3,
    subd_detection_vertical_range=5,
    subd_cluster_eps=sqrt(3),
    subd_cluster_mpt=13,
    
    ## tracing subdendrites
    trace_cluster_eps=sqrt(3),
    trace_cluster_mpt=17,
    trace_detection_depth=5,
    trace_detection_distance = 12,
    trace_rescore_dist=10,
    trace_rescore_angle=40
    
  ) %>%
    return()
  
}


select_intensity <- function(x,y,z,IMG){
  if(between(min(x), 1, ncol(IMG[[1]]))&between(max(x), 1, ncol(IMG[[1]]))&
     between(min(y), 1, nrow(IMG[[1]]))&between(max(y), 1, nrow(IMG[[1]]))&
     between(min(z), 1, length(IMG))&between(max(z), 1, length(IMG))){
    return(lapply(z, function(Z){return(IMG[[Z]][y,x])}) %>% unlist() %>% median())
  } else {
    return(NA)
  }
}
# 
# elongate_line_old <- function(GRAD, zGRAD, x_start, y_start, z_start, n){
#   
#   zr <- round(cos(cppdeg2rad(zGRAD)), 10)
#   zo <- round(sin(cppdeg2rad(zGRAD)),10)
#   yr <- round(sin(cppdeg2rad(GRAD)), 10)
#   xr <- round(cos(cppdeg2rad(GRAD)),10)
#   
#   c(
#     y=ceiling(y_start+(n*yr*zr)),
#     x=ceiling(x_start+(n*xr*zr)),
#     z=ceiling(z_start+(n*zo))
#   ) %>% return()
# }



# img <- binary_image
# GRAD <- VEC[["ha"]]
# zGRAD <- VEC[["va"]]
# 
# x_start <- start["x"]
# y_start <- start["y"]
# z_start <- 10
# r <- round(VEC[["l"]])







elongate_3d_sphere <- function(img, GRAD, zGRAD, x_start, y_start, z_start, r){
  
  zr <- round(cos(cppdeg2rad(zGRAD)), 10)
  
  zo <- round(sin(cppdeg2rad(zGRAD)),10)
  
  yr <- round(sin(cppdeg2rad(GRAD)), 10)
  
  xr <- round(cos(cppdeg2rad(GRAD)),10)
  
  tibble(n=1:r) %>%
    mutate(y=ceiling(y_start+(n*yr*zr)),
           x=ceiling(x_start+(n*xr*zr)),
           z=ceiling(z_start+(n*zo))) %>%
    rowwise() %>%
    mutate(i=cppSelectIntensity(x,y,z,img),
           deg=GRAD) %>%
    .[which(!is.na(.$i)),] %>%
    return()
  
}


elongate_line_full <- function(img, GRAD, zGRAD, x_start, y_start, z_start, r){
  lapply(1:r, function(N){
    cppElongateLine(GRAD, zGRAD, x_start, y_start, z_start, N)
  }) %>%
    bind_rows() %>%
    rowwise() %>%
    mutate(i=cppSelectIntensity(x,y,z,img),
           deg=GRAD) %>%
    .[which(!is.na(.$i)),] %>%
    return()
}

elongate_line <- function(GRAD, zGRAD, x_start, y_start, z_start, n){
  
  zr <- round(cos(cppdeg2rad(zGRAD)), 10)
  zo <- round(sin(cppdeg2rad(zGRAD)),10)
  yr <- round(sin(cppdeg2rad(GRAD)), 10)
  xr <- round(cos(cppdeg2rad(GRAD)),10)
  
  c(
    y=round(y_start+(n*yr*zr)),
    x=round(x_start+(n*xr*zr)),
    z=round(z_start+(n*zo))
    
  ) %>% return()
}

#GRAD <- most_likely_elongation 
#zGRAD <- vertical_angle

elongate_3d_sphere_unlim <- function(GRAD, 
                                     zGRAD, 
                                     x_start, 
                                     y_start, 
                                     z_start, 
                                     intensity_cutoff,
                                     steps,
                                     IMG
){
  
  
  ## trigonometry
  zr <- round(cos(cppdeg2rad(zGRAD)), 10)
  
  zo <- round(sin(cppdeg2rad(zGRAD)),10)
  
  yr <- round(sin(cppdeg2rad(GRAD)), 10)
  
  xr <- round(cos(cppdeg2rad(GRAD)),10)
  
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
                       IMG)
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
                                           steps,
                                           raw_image){
  
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
      f <- elongate_3d_sphere_unlim(GRAD, zGRAD, x_start, y_start, z_start, intensity_cutoff, steps, raw_image)  
      
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
                                steps,
                                raw_image) 
  fin <- f[["coords"]]
  return(list(c(fin, vertical_angle), f[["full"]]))
  
}



elongate_dendrite <- function(DENDRITE_raw,
                              SOMA, IMG){
  
  
  horizontal_subdivisions <- 60
  vertical_subdivisions <- 10
  horizontal_detection_angle <- 60
  vertical_detection_angle <- 12
  intensity_cutoff <- 0.5
  steps <- 10
  
  
  cat(paste("\ndendrite number:",DENDRITE_raw[["dendrite_id"]]))
  DENDRITE <- DENDRITE_raw %>% as.numeric() %>% 
    set_names(nm=names(DENDRITE_raw))
  
  
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
                                                   steps,
                                                   IMG)
    next_pos <- raw_dendrite[[1]]
    #full_vector_list[[c]] <- raw_dendrite[[2]]
    c <- c+1 
    knot_list[[c]] <- next_pos %>% set_names(c("x","y","z", "h_angle", "v_angle"))
    #print(Sys.time()-st)
  }
  return(append(list(c(x=round(SOMA$x), y=round(SOMA$y), z=round(SOMA$z), h_angle=DENDRITE[["h_angle"]], v_angle=0)),knot_list[3:c-1]))
  #return(knot_list[3:c-1])
  
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

find_dendritic_start_sites <- function(SOMA, IMG){
  
  deg_step <- 0.005
  r <- 300
  
  mn = lapply(seq(deg_step,1,deg_step)*360, function(GRAD){
    
    f <- elongate_3d_sphere(IMG, GRAD,0,SOMA[["x"]],SOMA[["y"]], SOMA[["z"]], 300)
    
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
    
    f <- elongate_3d_sphere(IMG, GRAD,0, SOMA[["x"]],SOMA[["y"]], SOMA[["z"]], soma_radius)
    
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


get_somata <- function(cube_size, r, deg_step ,raw_image){
  
  z_vec <- seq(((cube_size+1)/2), (length(raw_image)-(cube_size+1)/2), cube_size)
  y_vec <- seq(((cube_size+1)/2), (nrow(raw_image[[1]])-(cube_size+1)/2), cube_size)
  x_vec <- seq(((cube_size+1)/2), (ncol(raw_image[[1]])-(cube_size+1)/2), cube_size)
  ## loop over z axis
  n_its <- length(raw_image)*length(y_vec)*length(x_vec)
  n_vox_per_cube <- cube_size^3
  
  first_test <- lapply(z_vec, function(Z){
    lapply(y_vec, function(Y){
      lapply(x_vec, function(X){
        zsum <- lapply(c((Z-(cube_size-1)/2):(Z+(cube_size-1)/2)), function(Z2){
          raw_image[[Z2]][c((Y-(cube_size-1)/2):(Y+(cube_size-1)/2)), 
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
  
  z_raw <- lapply(1:length(raw_image), function(Z){
    c(z=Z,
      mn= lapply(seq(deg_step,1,deg_step)*360, function(GRAD){
        
        f <- elongate_3d_sphere(raw_image, GRAD, 0, xy_soma$x, xy_soma$y ,Z, 100)
        
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


#rad2deg <- function(rad){(rad * 180) / (pi)}
#deg2rad <- function(deg){(deg * pi) / (180)}


#sl <- segment_list
#file_name <- "f:/data_sholl_analysis/test/spec_segs/ggg.tif"


test_segmentation <- function(sl, raw_image, file_name, soma_reg){
  QNT <- 0.9
  
  nosoma_image <- list("1"=raw_image[[SOMA[["z"]]]])
  
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


#nosoma_image <- raw_image
#VOX <- sl[[1]]
bi2 <- function(sl, nosoma_image){
  
  full_cutoffs <- lapply(sl, function(VOX){
    #print(1)
    all_vox_raw <- lapply(1:length(nosoma_image), function(L){
      
      LAYER <- nosoma_image[[L]]
      
      return(LAYER[VOX])
      
    }) %>%
      Reduce(function(x,y)c(x,y),.) 
    
    all_vox <- all_vox_raw[all_vox_raw!=0]
    
    dens_func <- density(all_vox, bw=0.001) 
    
    #ggplot(tibble(x=dens_func$x,y=dens_func$y), aes(x=x, y=y))+ geom_line()
    
    # all_local_extreme <- get_minmax(dens_func) %>%
    #   filter(x<0.95)
    # 
    # loc_ex <- all_local_extreme[which(all_local_extreme$y==max(all_local_extreme$y)),]$x
    #mn <- min(all_vox)
    #mx <- max(all_vox)
    
    df <- tibble(x=dens_func$x, y=dens_func$y, d1=c(0, diff(dens_func$y)))%>%
      filter(between(x, 0.05, 0.95))
    
    dec_row <- which(df$d1==min(df$d1))
    
    if(dec_row==nrow(df)){dec_row <- dec_row <- dec_row-1}
    
    strongest_decrease <- df[c(dec_row, dec_row+1),c(1,2)]
    
    diffs <- apply(strongest_decrease, 2, diff)
    
    cutoff <- strongest_decrease[[1,"x"]]+abs(strongest_decrease[[1,"y"]]/diffs[["y"]])*abs(diffs[["x"]])
    
    return(cutoff)
    
  })
  
  
  new_image <- nosoma_image
  
  
  for(i in 1:length(full_cutoffs)){
    #cat(paste("\n  ", i))
    cutoff <- full_cutoffs[[i]]
    VOX <- sl[[i]]
    for(l in 1:length(new_image)){
      #print(l)
      #cutoff <- cutoffs[[l]]
      new_image[[l]][VOX][new_image[[l]][VOX]>=cutoff] <- 1
      new_image[[l]][VOX][new_image[[l]][VOX]<cutoff] <- 0
      new_image[[l]] <- matrix(new_image[[l]], nrow = 1040)
    }
    
  }
  
  #writeTIFF(new_image, file_name)
  return(new_image)
}

bi3 <- function(sl, nosoma_image, QNT){
  
  full_cutoffs <- lapply(sl, function(VOX){
    #print(1)
    all_vox_raw <- lapply(1:length(nosoma_image), function(L){
      
      LAYER <- nosoma_image[[L]]
      
      return(LAYER[VOX])
      
    }) %>%
      Reduce(function(x,y)c(x,y),.) 
    
    all_vox <- all_vox_raw[all_vox_raw!=0]
    
    cutoff <- quantile(all_vox, QNT/100)
    
    # dens_func <- density(all_vox, bw=0.001) 
    # 
    # #ggplot(tibble(x=dens_func$x,y=dens_func$y), aes(x=x, y=y))+ geom_line()
    # 
    # # all_local_extreme <- get_minmax(dens_func) %>%
    # #   filter(x<0.95)
    # # 
    # # loc_ex <- all_local_extreme[which(all_local_extreme$y==max(all_local_extreme$y)),]$x
    # #mn <- min(all_vox)
    # #mx <- max(all_vox)
    # 
    # df <- tibble(x=dens_func$x, y=dens_func$y, d1=c(0, diff(dens_func$y)))%>%
    #   filter(between(x, 0.05, 0.95))
    # 
    # dec_row <- which(df$d1==min(df$d1))
    # 
    # if(dec_row==nrow(df)){dec_row <- dec_row <- dec_row-1}
    # 
    # strongest_decrease <- df[c(dec_row, dec_row+1),c(1,2)]
    # 
    # diffs <- apply(strongest_decrease, 2, diff)
    # 
    # cutoff <- strongest_decrease[[1,"x"]]+abs(strongest_decrease[[1,"y"]]/diffs[["y"]])*abs(diffs[["x"]])
    
    return(cutoff)
    
  })
  
  
  new_image <- nosoma_image
  
  
  for(i in 1:length(full_cutoffs)){
    #cat(paste("\n  ", i))
    cutoff <- full_cutoffs[[i]]
    VOX <- sl[[i]]
    for(l in 1:length(new_image)){
      #print(l)
      #cutoff <- cutoffs[[l]]
      new_image[[l]][VOX][new_image[[l]][VOX]>=cutoff] <- 1
      new_image[[l]][VOX][new_image[[l]][VOX]<cutoff] <- 0
      new_image[[l]] <- matrix(new_image[[l]], nrow = 1040)
    }
    
  }
  
  #writeTIFF(new_image, file_name)
  return(new_image)
}



normalize_regions <- function(sl, raw_image, file_name, soma_reg){
  QNT <- 0.9
  ## remove soma 
  # nosoma_image <- lapply(1:length(raw_image), function(n){
  #   diff <- abs(SOMA[["z"]]-n)
  #   LAYER <- raw_image[[n]]
  #   
  #   LAYER[c((soma_reg$y-diff*2):(soma_reg$yend+diff*2)), c((soma_reg$x-diff*2):(soma_reg$xend+diff*2))] <- NA
  #   return(LAYER)
  #   
  # })
  nosoma_image <- raw_image
  
  
  full_cutoffs <- lapply(sl, function(VOX){
    print(1)
    cutoff <- lapply(nosoma_image, function(LAYER){
      
      return(LAYER[VOX])
      
    }) %>%
      Reduce(function(x,y)c(x,y),.) %>%
      quantile(QNT, na.rm=T) %>%
      return()
    
  })
  
  new_image <- raw_image
  
  
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
get_3d_single_index <- function(x,y,z, nr){
  paste(z,get_single_index(x,y,nr), sep="_")
}

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

create_segment <- function(xs, xe, cof, line, nr, nc, bd, sdir, direction, nr_orig){
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



create_raster <- function(nELD, avs, borders_vox_raw, selection_vector, elongated_dendrites,nr_orig, nc_orig, process_full){
  cat(paste("dendrite:", nELD))
  ################################################################################
  ######################## defining basic elements ###############################    
  ################################################################################       
  ELD_raw <- elongated_dendrites[[nELD]]  
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
    
    #print(n)
    
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
    
    pix_to_bottom <- -med
    pix_to_top <- nr-med
    n_segments_top <- abs(round(pix_to_top/avs))
    n_segments_bottom <- abs(round(pix_to_bottom/avs))
    use_length_top <- round(pix_to_top/n_segments_top)
    use_length_bottom <- round(pix_to_bottom/n_segments_bottom)
    
    if(process_full==F&n_segments_top>2){
      n_segments_top <- 2
    }
    if(process_full==F&n_segments_bottom>2){
      n_segments_bottom <- 2
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
      segment_list[[segment_count]] <- all_vox_top[[1]]
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
      segment_list[[segment_count]] <- all_vox_bottom[[1]]
    }       
    #print(4)
    y_coord_list[[n]] <- med_line[[nrow(med_line), "y"]]
    
    
  }
  return(segment_list)  
}




remove_soma <- function(SOMA, soma_radius, raw_image, nr_orig, nc_orig){
  nosoma_image <- lapply(1:length(raw_image), function(n){
    r <- soma_radius+abs(round(SOMA$z)-n)
    circ <- lapply(c(-r:r), function(X){
      y=sqrt(r^2-X^2)
      bottom <- round(SOMA$y-y)
      top <- round(SOMA$y+y)
      lapply(c(bottom:top), get_single_index, x=X+round(SOMA$x), nr=nr_orig) %>%
        return()
    }) %>% unlist() 
    LAYER <- raw_image[[n]]
    LAYER[circ] <- 0
    return(LAYER)
  })  
  return(nosoma_image)
}

# IMG=images$noS_image
# main_vectors=inter$main_vectors
# main_vectors_full=inter$main_vectors_full
# rem_rad=opt$subd_detection_distance-opt$subd_detection_depth-1
# nr_orig=opt$nr_orig
# nc_orig=opt$nc_orig


remove_main_dendrites <- function(IMG, main_vectors, main_vectors_full, rem_rad, nr_orig, nc_orig){
  
  #rem_rad <- 30
  
  comb_main_dend <- lapply(1:length(main_vectors), function(n){
    full_vecs <- main_vectors_full[[n]]
    full_dendrite <- bind_rows(full_vecs) %>%
      mutate_all(.funs = as.numeric) %>%
      mutate(ind=get_single_index(x,y,nr_orig)) 
  }) %>% bind_rows()
  
  surr_vox <- lapply(1:length(main_vectors), function(n){
    print(n)
    rel_vecs <- main_vectors[[n]]
    full_vecs <- main_vectors_full[[n]]
    
    all_surrounding_vox <- define_surr_layer(rel_vecs, full_vecs, rem_rad, 0, IMG, 0) #%>%
    #mutate(ind=get_single_index(x,y,nr_orig))
    
  }) %>% bind_rows()
  
  ## get all voxel values
  
  all_vox <- lapply(c(min(surr_vox$x):max(surr_vox$x)), function(X){
    
    all_x <- surr_vox[which(surr_vox$x==X),] #%>%
    #arrange(y)
    return(tibble(x=X, y=c(min(all_x$y):max(all_x$y))))
    
    
    #return(c(first(all_x$y):last(all_x$y)))
    #return(c(min(all_x$ind):max(all_x$ind)))
    
  }) %>% bind_rows()#c(recursive=T)
  
  ## remove those which are more than 35 vox away from main vectors
  
  area_size <- 20
  n_area <- round(nc_orig/area_size)
  
  adj_vox <- lapply(1:n_area, function(N){
    cat(paste("\n", N, "of", n_area))
    start <- N*area_size
    end <- (N+1)*area_size-1
    rel_vox <- all_vox[which(between(all_vox$x, start, end)),]
    
    if(nrow(rel_vox)==0){return(NULL)}
    
    rel_MD <- comb_main_dend[which(between(comb_main_dend$x, start-(rem_rad+1), end+(rem_rad+1))),]
    
    t2 <- apply(rel_vox, 1, function(V1){
      t <- apply(rel_MD, 1, function(V2){
        cppDistPts(V1["x"], V1["y"], 0, V2["x"], V2["y"], 0)
      }) %>% min()
      
      #if(is.na(t)){t<-0}
      if(t<rem_rad){return(V1)}else{return(NULL)}
    }, simplify=F) %>%
      compact() %>%
      bind_rows() %>%
      return()
    
  }) %>% compact() %>%
    bind_rows() %>%
    mutate(ind=get_single_index(x,y,nr_orig))
  
  #ggplot(bind_rows(adj_vox), aes(x=x, y=y))+geom_point()
  
  vox_to_rem <- pull(adj_vox, ind)
  
  rem_img <- lapply(IMG, function(LAYER){
    LAYER[vox_to_rem] <- 0
    #return(matrix(LAYER, nrow = nr_orig))
    return(LAYER)
  })
  
  return(rem_img)
}






dist_pts <- function(x1, y1, z1, x2, y2, z2){
  
  
  sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)
  
}

get_va <- function(dist, z1, z2){
  ## z1 = destination, z2=start
  angle <- cpprad2deg(asin((abs(z1-z2))/dist))
  if(z2>=z1){
    fin <- 0-angle
  } else {
    fin <- 0+angle
  }
  return(fin)
}

get_ha <- function(x2,y2,x1,y1){
  ## x2 = destination... x1 = begin
  angle <- cpprad2deg(atan((abs(x1-x2))/abs(y1-y2)))
  if(x1>=x2&y1>=y2){
    fin <- 270-angle
  } else if(x1>=x2&y1<=y2){
    fin <- 90+angle
  } else if(x1<=x2&y1>=y2){
    fin <- 270+angle
  } else if(x1<=x2&y1<=y2){
    fin <- 90-angle
  } 
  return(fin)
}

fit_subdendrite <- function(SUBD, dist, ha, va, IMG){
  
  vec <- elongate_3d_sphere(IMG, ha, va, SUBD[["x"]], SUBD[["y"]], SUBD[["z"]], dist) %>%
    pull(i) %>%
    mean() %>%
    return()
  
  
}  

fit_subdendrite_new <- function(SUBD, dist, ha, va, IMG){
  
  vec <- elongate_line_full(IMG, ha, va, SUBD[["x"]], SUBD[["y"]], SUBD[["z"]], dist) #%>%
  
  max_n <- vec %>% ungroup() %>%
    mutate(cumu=cumsum(i), n=c(1:nrow(.))) %>%
    group_by(cumu) %>%
    summarize(c=length(i), mn=min(n)) %>%
    filter(c>10) 
  
  if(nrow(max_n)==0){
    vec %>%pull(i) %>%
      mean() %>%
      return()
  } else {
    return(0)
  } 
  
  
}

define_full_vector_info <- function(ELD, IMG){
  
  #ELD <- elongated_dendrites[[2]]
  
  all_vecs <- lapply(1:length(ELD), function(n){
    if(n==length(ELD)){
      return(NULL)
    } else {
      
      xs=ELD[[n]][["x"]]
      ys=ELD[[n]][["y"]]
      zs=ELD[[n]][["z"]]
      xe=ELD[[n+1]][["x"]]
      ye=ELD[[n+1]][["y"]]
      ze=ELD[[n+1]][["z"]]
      
      
      ## length of vector
      l <- sqrt((xs-xe)^2 + (ys-ye)^2 + (zs-ze)^2)
      
      
      ## check if vector diverges from soma
      dist_start_soma=cppDistPts(xs, ys, zs, 
                               round(SOMA[["x"]]), 
                               round(SOMA[["y"]]), 
                               round(SOMA[["z"]]))
      dist_end_soma=cppDistPts(xe, ye, ze, 
                             round(SOMA[["x"]]), 
                             round(SOMA[["y"]]), 
                             round(SOMA[["z"]]))
      
      direction <- ifelse(dist_start_soma<dist_end_soma,
                          1, -1)
      
      c(
        
        xs=xs,
        ys=ys,
        zs=zs,
        xe=xe,
        ye=ye,
        ze=ze,
        ha=ELD[[n+1]][["h_angle"]],
        va=ELD[[n+1]][["v_angle"]],
        l=l,
        level=n,
        dir=direction
      ) %>% return()
    } 
  }) %>% compact() #%>% 
  #bind_rows() %>%
  
  full_voxels <- lapply(all_vecs, function(VEC){
    elongate_3d_sphere(IMG, VEC[["ha"]], VEC[["va"]], VEC[["xs"]], VEC[["ys"]], VEC[["zs"]], round(VEC[["l"]]))%>% 
      select(x,y,z) %>%
      mutate(level=VEC[["level"]])
  }) 
  
  return(list(all_vecs, full_voxels))
}


# xs <- all_vectors[[1]]$xs
# ys <- all_vectors[[1]]$ys
# zs <- all_vectors[[1]]$zs
# xe <- all_vectors[[1]]$xe
# ye <- all_vectors[[1]]$ye
# ze <- all_vectors[[1]]$ze

get_full_vector_voxels <- function(V, IMG){
  
  xs <- V[["xs"]]
  ys <- V[["ys"]]
  zs <- V[["zs"]]
  xe <- V[["xe"]]
  ye <- V[["ye"]]
  ze <- V[["ze"]]
  dist <- cppDistPts(xs,ys,zs,xe,ye,ze)
  
  ha <- cppGetHA(xe, ye, xs, ys)
  
  va <- cppGetVA(dist, ze, zs)
  
  
  elongate_line_full(IMG, ha, va, xs, ys, zs, round(dist))%>% 
    select(x,y,z) %>%
    return()
  
}

adjust_duplicated_starts <- function(MAIND){
  
  all_vectors <- lapply(MAIND, get_all_subdendrite_vectors) %>% 
    Reduce(function(x,y)append(x,y),.)
  
  all_vectors_fv <- lapply(all_vectors, get_full_vector_voxels, IMG=IMG) 
  
  all_vectors_fv_labeled <- lapply(1:length(all_vectors), function(n){mutate(all_vectors_fv[[n]], vec_id=n)}) %>%
    bind_rows() %>%
    mutate(id=paste(x,y,sep="_"))
  
  dup <- lapply(1:length(all_vectors), function(n){
    
    to_check <- all_vectors_fv_labeled %>%
      filter(vec_id==n) %>%
      pull(id)
    
    to_compare <- all_vectors_fv_labeled %>%
      filter(!vec_id==n)
    
    crossed <- to_compare %>%
      filter(id %in% to_check)
    
    if(nrow(crossed)==0){
      return(NULL)
    } else {
      
      tibble(v1=min(c(n, unique(crossed$vec_id))),
             v2=max(c(n, unique(crossed$vec_id))),
             x=median(crossed$x)%>%round(),
             y=median(crossed$y)%>%round(),
             z=median(crossed$z)%>%round()) %>%
        return()
      
    }
    
  }) %>% compact() %>%
    bind_rows() %>%
    mutate(id=paste(v1, v2, sep="_")) %>%
    group_by(id) %>%
    summarize(v1=unique(v1),
              v2=unique(v2),
              x=round(median(x)),
              y=round(median(y)),
              z=round(median(z))) %>% rowwise() %>%
    mutate(v1_angle=select_angle(v1, all_vectors = all_vectors),
           v2_angle=select_angle(v2, all_vectors = all_vectors))
  
  ### angle smaller 90 deg ... set cross as new subnode
  
  
  
  ### angle larger 90 deg... do not adjust
  
  
}

select_angle <- function(n, all_vectors, col){
  
  all_vectors[[n]][[col]]
  
}



#start <- start_t
#elgt_len <- round(len_b)




screen_subdendrite_starts <- function(start, rel_z_layer, elgt_len, VEC, IMG){
  df <- lapply(rel_z_layer, function(ZS){
    # line <- elongate_3d_sphere(IMG, VEC[["ha"]], VEC[["va"]], 
    #                            start[["x"]], start[["y"]], ZS, elgt_len)
    line <- elongate_line_full(IMG, VEC[["ha"]], VEC[["va"]], 
                               start[["x"]], start[["y"]], ZS, elgt_len)
  }) %>% 
    bind_rows() %>% 
    rowwise() %>%
    mutate(dist_to_soma=cppDistPts(x, y, z, 
                                 round(SOMA[["x"]]), 
                                 round(SOMA[["y"]]), 
                                 round(SOMA[["z"]]))) %>%
    select(x,y,z,i,dist_to_soma)
  return(df)
}


#df <- bind_rows(df_comb_b, df_comb_t)

# find_subd_cluster <- function(df){
#   df_filtered <- df %>% 
#     filter(i!=0) %>%
#     select(x,y,z)
#   if(nrow(df_filtered)==0){return(NULL)} 
#   df_clustered <- df_filtered %>% 
#     ungroup() %>%
#     mutate(c=fpc::dbscan(df_filtered, eps = 3.5, MinPts = 15)$cluster) %>%
#     filter(c!=0)
#   if(nrow(df_clustered)==0){return(NULL)}
#   #ggplot(df_clustered %>% filter(c!=0), aes(x=x, y=z, color=as.character(c)))+ geom_point()
#   df_centers <- lapply(unique(df_clustered$c), function(CLUST){
#     kmeans(df_clustered %>% filter(c==CLUST) %>% select(x,y,z), centers=1)$centers %>%
#       round() %>%
#       set_names(nm=c("x","y","z"))
#   })
#   return(df_centers)
# }

#df <-centers_rem

# EPS <- 2
# MPTS <- 33

find_subd_cluster_man <- function(df, EPS, MPTS){
  df_filtered <- df %>% 
    .[which(.$i!=0),] %>%
    select(x,y,z) #
  #ggplot(df_filtered, aes(x=x, y=-y,))+geom_point()+facet_wrap(~z)
  if(nrow(df_filtered)==0){return(NULL)} 
  df_clustered <- df_filtered %>% 
    ungroup() %>%
    mutate(c=fpc::dbscan(df_filtered, eps = EPS, MinPts = MPTS)$cluster) %>%
    .[which(.$c!=0),]
  if(nrow(df_clustered)==0){return(NULL)}
  #ggplot(df_clustered, aes(x=x, y=-y,color=as.character(c)))+geom_point()+facet_wrap(~z)
  df_centers <- lapply(unique(df_clustered$c), function(CLUST){
    kmeans(df_clustered %>%  
             .[which(.$c==CLUST),] %>%
             select(x,y,z), 
           centers=1)$centers %>%
      round() %>%
      set_names(nm=c("x","y","z"))
  })
  return(df_centers)
}


find_subd_cluster_man2 <- function(Z, EPS, MPTS, xs, ys){
  
  
  
  rel_z_layer <- seq(Z-z_range, Z+z_range, 1) %>%
    .[between(., 1, length(IMG))] 
  
  centers1 <- screen_circular_new(det_rad, rel_z_layer,
                                  INC, xs, ys, zs, ha, ha_cutoff, screening_angle, IMG) %>%
    mutate(index=get_3d_single_index(x,y,z, nr_orig))

  ## remove positions of other dendrites
  centers_rem <- centers1 %>%
    .[which(!.$index %in% c(all_points, recursive=T)),]
  
  
  
  df_filtered <- df %>% 
    .[which(.$i!=0),] %>%
    select(x,y,z) #
  #ggplot(df_filtered, aes(x=x, y=-y,))+geom_point()+facet_wrap(~z)
  if(nrow(df_filtered)==0){return(NULL)} 
  
  df_clustered <- df_filtered %>% 
    ungroup() %>%
    mutate(c=fpc::dbscan(df_filtered, eps = EPS, MinPts = MPTS)$cluster) %>%
    .[which(.$c!=0),]%>%
    rowwise() %>%
    mutate(ha=cppGetHA(x,y, xs, ys))
  if(nrow(df_clustered)==0){return(NULL)}
  #ggplot(df_clustered, aes(x=x, y=-y,color=as.character(c)))+geom_point()+facet_wrap(~z)
  df_centers <- lapply(unique(df_clustered$c), function(CLUST){
    kmeans(df_clustered %>%  
             .[which(.$c==CLUST),] %>%
             select(x,y,z), 
           centers=1)$centers %>%
      round() %>%
      set_names(nm=c("x","y","z"))
  })
  return(df_centers)
}




adj_deg <- function(GRAD){
  if(GRAD>360){
    return(GRAD-360)
  } else if(GRAD<0){
    return(GRAD+360)
  } else {
    return(GRAD)
  }
}

# X <- VEC[["xe"]]
# Y <- VEC[["ye"]]
# Z <- VEC[["ze"]]
# ha <- VEC[["ha"]]
# X <- SUBD[["x_el"]]
# Y <- SUBD[["y_el"]]
# Z <- SUBD[["z_el"]]
# ha <- SUBD[["ha"]]
# ha <- 290
## es geht nur ha zwischen 90 und 270

screen_circular <- function(det_rad, z_range, X, Y, Z, ha, IMG){
  
  
  bmin=ha-90
  bmax=ha+90
  
  
  rel_z_layer <- seq(Z-z_range, Z+z_range, 1) %>%
    .[between(., 1, length(IMG))]
  
  lapply(rel_z_layer, function(ZS){
    
    df_circ_x <- tibble(x=seq(-det_rad, det_rad, 1)) %>%
      mutate(yp=sqrt(det_rad^2-x^2),
             yn=-sqrt(det_rad^2-x^2)) %>%
      gather(dir, y, -x)
    
    df_circ_y <- tibble(y=seq(-det_rad, det_rad, 1)) %>%
      mutate(xp=sqrt(det_rad^2-y^2),
             xn=-sqrt(det_rad^2-y^2)) %>%
      gather(dir, x, -y)
    
    df_circ_full <- bind_rows(df_circ_x, df_circ_y) %>%
      mutate_at(.vars = c("x", "y"), .funs = round) %>%
      mutate(id=paste(x, y, sep="_")) %>%
      filter(!duplicated(id)) %>%
      rowwise() %>%
      mutate(xa=x+X,
             ya=y+Y,
             angle=cpprad2deg(atan((y/det_rad)/(x/det_rad))),
             af=case_when(
               x<0 ~ angle+180,
               angle<0~ angle+360,
               TRUE ~ angle
             ), 
             #cs=cos(0.5*cppdeg2rad(af))
      ) #%>%
    
    if(between(bmin, 0, 360)&between(bmax, 0, 360)){
      
      df_circ <- df_circ_full %>%
        filter(between(af, bmin, bmax)) 
      
    } else if(bmin<0){
      
      df_circ <- df_circ_full %>%
        filter(between(af, bmin+360, 360)|between(af, 0, bmax))
      
    } else {
      
      df_circ <- df_circ_full %>%
        filter(between(af, bmin, 360)|between(af, 0, bmax-360))
      
    }
    
    
    df_circ %>%     
      rowwise() %>%
      mutate(i=cppSelectIntensity(xa, ya, ZS, IMG),
             z=ZS) %>%
      
      select(x=xa, y=ya, z, i) %>%
      mutate(dist_to_soma=cppDistPts(x, y, z, round(SOMA[["x"]]), round(SOMA[["y"]]), round(SOMA[["z"]]))) %>%
      return()
    # ggplot(df_circ, aes(x=x, y=y, color=cs))+geom_point()+
    #   scale_color_gradientn(#breaks=c(0, 90, 180, 270,360), 
    #                         colors = c("green" ,"red", "blue", "black"))
    
  }) %>%
    compact() %>%
    bind_rows() %>%
    return()
  
}


create_circle <- function(DET_RAD){
  df_circ_x <- tibble(x=seq(-DET_RAD, DET_RAD, 1)) %>%
    mutate(yp=sqrt(DET_RAD^2-x^2),
           yn=-sqrt(DET_RAD^2-x^2)) %>%
    gather(dir, y, -x)
  
  df_circ_y <- tibble(y=seq(-DET_RAD, DET_RAD, 1)) %>%
    mutate(xp=sqrt(DET_RAD^2-y^2),
           xn=-sqrt(DET_RAD^2-y^2)) %>%
    gather(dir, x, -y)
  
  df_circ_full <- bind_rows(df_circ_x, df_circ_y) %>%
    mutate_at(.vars = c("x", "y"), .funs = round) %>%
    mutate(id=paste(x, y, sep="_")) %>%
    filter(!duplicated(id)) %>%
    return()
}

# 
# X=xs
# Y=ys
# Z=zs




screen_circular_new <- function(det_rad, rel_z_layer,bmin, bmax,
                                INC, X, Y, Z, ha, ha_cutoff, screening_angle, IMG){
  
  #t0 <- Sys.time()
  # if(is.null(screening_angle)){
  #   bmin <- 0
  #   bmax <- 360
  # } else {
  # 
  #   bmin=ha-0.5*screening_angle
  #   bmax=ha+0.5*screening_angle
  #   if(!is.null(ha_cutoff)){
  #     if(ha<ha_cutoff){
  #       if(bmax>ha_cutoff){
  #         bmax <- round(ha_cutoff)
  #       }
  #       
  #       if(bmin<ha_cutoff-180){
  #         bmin <- round(ha_cutoff-180)
  #       }
  #     } else {
  #       if(bmax>ha_cutoff+180){
  #         bmax <- round(ha_cutoff+180)
  #       }
  #       if(bmin<ha_cutoff){
  #         bmin <- round(ha_cutoff)
  #       }
  #     }
  #   }
  # }
  # rel_z_layer <- seq(Z-z_range, Z+z_range, 1) %>%
  #   .[between(., 1, length(raw_image))]
  
  #t1 <- Sys.time()
  #lapply(c(det_rad-INC, det_rad+INC), function(DET_RAD){
  #lapply(seq(det_rad-INC, det_rad+INC, 1), function(DET_RAD){
  
    DET_RAD <- det_rad+INC
      
    fin <- lapply(rel_z_layer, function(ZS){
      
      df_circ_x <- tibble(x=seq(-DET_RAD, DET_RAD, 1)) %>%
        mutate(yp=sqrt(DET_RAD^2-x^2),
               yn=-sqrt(DET_RAD^2-x^2)) %>%
        gather(dir, y, -x)
      
      df_circ_y <- tibble(y=seq(-DET_RAD, DET_RAD, 1)) %>%
        mutate(xp=sqrt(DET_RAD^2-y^2),
               xn=-sqrt(DET_RAD^2-y^2)) %>%
        gather(dir, x, -y)
      
      df_circ_full <- bind_rows(df_circ_x, df_circ_y) %>%
        mutate_at(.vars = c("x", "y"), .funs = round) %>%
        mutate(id=paste(x, y, sep="_")) %>%
        filter(!duplicated(id)) %>%
        rowwise() %>%
        mutate(xa=x+X,
               ya=y+Y,
               angle=cpprad2deg(atan((y/DET_RAD)/(x/DET_RAD))),
               af=case_when(
                 x<0 ~ angle+180,
                 angle<0~ angle+360,
                 TRUE ~ angle
               ), 
               #cs=cos(0.5*cppdeg2rad(af))
        ) #%>%
      
      all_vox <- lapply(c(min(df_circ_full$x):max(df_circ_full$x)), function(X){
        
        all_x <- df_circ_full[which(df_circ_full$x==X),] #%>%
        #arrange(y)
        return(tibble(x=X, y=c(min(all_x$y):max(all_x$y))))
        
      }) %>% bind_rows() %>% rowwise() %>%
        mutate(dist=cppDistPts(x,y,0, 0,0,0),
               ha=cppGetHA(x,y,0,0)) %>%
        .[which(.$dist>det_rad-INC),]
      
      #ggplot(all_vox, aes(x=x, y=y))+geom_tile()
      if(between(bmin, 0, 360)&between(bmax, 0, 360)){
        
        #df_circ <- df_circ_full[which(between(df_circ_full$af, bmin, bmax)),]
        
        df_circ <- all_vox[which(between(all_vox$ha, bmin, bmax)),]
        
      } else if(bmin<0){
        
        df_circ <- all_vox[which(between(all_vox$ha, bmin+360, 360)|between(all_vox$ha, 0, bmax)),]
        
      } else {
        df_circ <- all_vox[which(between(all_vox$ha, bmin, 360)|between(all_vox$ha, 0, bmax-360)),]
        
      }
      
      #ggplot(df_circ, aes(x=x, y=y))+geom_tile()
      df_circ %>%    
        
        
        rowwise() %>%
        
        
        
        mutate(z=ZS,
               xa=x+X,
               ya=y+Y,
               i=cppSelectIntensity(xa, ya, z, IMG=IMG)) %>%
        
        select(x=xa, y=ya, z, i) %>%
        mutate(dist_to_soma=cppDistPts(x, y, z, round(SOMA[["x"]]), round(SOMA[["y"]]), round(SOMA[["z"]]))) %>%
        return()
      # ggplot(df_circ, aes(x=x, y=y, color=cs))+geom_point()+
      #   scale_color_gradientn(#breaks=c(0, 90, 180, 270,360), 
      #                         colors = c("green" ,"red", "blue", "black"))
      
    }) %>%
      compact() %>%
      bind_rows() %>%
      return() 
  #}) %>% bind_rows() %>%
    #return()
  
}



define_surr_layer <- function(rel_vecs, full_vecs, DET_RAD, z_range, IMG, INC){
  test <- lapply(c((DET_RAD-INC):(DET_RAD+INC)), function(det_rad){
    all_sorrounding_vox <- lapply(1:length(rel_vecs), function(v){ 
      cat(paste("\nvec:",v))
      
      VEC <- rel_vecs[[v]]
      FVEC <- full_vecs[[v]]
      
      ortho_t <- VEC[["ha"]]+90
      ortho_b <- VEC[["ha"]]-90
      
      rel_z_layer <- seq(VEC[["zs"]]-z_range, VEC[["zs"]]+z_range, 1) %>%
        .[between(., 1, length(IMG))]
      
      df_circ_end <- NULL
      
      ## define starts
      
      if(v==1){
        start_t <- cppElongateLine(ortho_t, 0, VEC[["xs"]], VEC[["ys"]], VEC[["zs"]], det_rad)
        start_b <- cppElongateLine(ortho_b, 0, VEC[["xs"]], VEC[["ys"]], VEC[["zs"]], det_rad)  
        elgt <- 0
        add_for_start_b <- 0
        add_for_start_t <- 0
      } else {
        prevVEC <- rel_vecs[[v-1]]
        elgt <- abs(det_rad*tan(cppdeg2rad(0.5*abs((VEC[["ha"]]-90)-(prevVEC[["ha"]]-90)))))
        
        if(prevVEC[["ha"]]>VEC[["ha"]]){
          
          adj_t <- cppElongateLine(VEC[["ha"]]-180, 0, VEC[["xs"]], VEC[["ys"]], VEC[["zs"]], elgt)
          adj_b <- cppElongateLine(VEC[["ha"]], 0, VEC[["xs"]], VEC[["ys"]], VEC[["zs"]], elgt)
          
          start_t <- cppElongateLine(ortho_t, 0, adj_t[["x"]], adj_t[["y"]], adj_t[["z"]], det_rad)
          start_b <- cppElongateLine(ortho_b, 0, adj_b[["x"]], adj_b[["y"]], adj_b[["z"]], det_rad)
          
          add_for_start_t <- elgt
          add_for_start_b <- -elgt
        } else {
          
          adj_b <- cppElongateLine(VEC[["ha"]]-180, 0, VEC[["xs"]], VEC[["ys"]], VEC[["zs"]], elgt)
          adj_t <- cppElongateLine(VEC[["ha"]], 0, VEC[["xs"]], VEC[["ys"]], VEC[["zs"]], elgt)
          
          start_t <- cppElongateLine(ortho_t, 0, adj_t[["x"]], adj_t[["y"]], adj_t[["z"]], det_rad)
          start_b <- cppElongateLine(ortho_b, 0, adj_b[["x"]], adj_b[["y"]], adj_b[["z"]], det_rad)
          
          add_for_start_t <- -elgt
          add_for_start_b <- elgt
        }
        
      }
      
      ## define elongation lengths
      
      if(v!=length(rel_vecs)){
        
        nextVEC <- rel_vecs[[v+1]]
        elgt_end <- abs(det_rad*tan(cppdeg2rad(0.5*abs((VEC[["ha"]]-90)-(nextVEC[["ha"]]-90)))))
        
        if(nextVEC[["ha"]]>VEC[["ha"]]){
          len_b <- VEC[["l"]]+elgt_end+add_for_start_b
          len_t <- VEC[["l"]]-elgt_end+add_for_start_t
        } else {
          len_b <- VEC[["l"]]-elgt_end+add_for_start_b
          len_t <- VEC[["l"]]+elgt_end+add_for_start_t
        }
      } else {
        
        len_b <- VEC[["l"]]+add_for_start_b
        len_t <- VEC[["l"]]+add_for_start_t
        
        #### performing circular screen at end of main dendrite:
        df_circ_end <- screen_circular(det_rad, z_range,
                                       VEC[["xe"]], VEC[["ye"]], VEC[["ze"]], VEC[["ha"]],
                                       IMG)
        
        # df_circ_end <- screen_circular_new(det_rad, rel_z:layer, 2, 
        #                                VEC[["xe"]], VEC[["ye"]], VEC[["ze"]], VEC[["ha"]],
        #                                ha_cutoff, screening_angle)
        #                                #IMG)
        
        
      }
      
      if(len_t>0){
        df_top <- screen_subdendrite_starts(start_t, rel_z_layer,round(len_t), VEC, IMG)
      } else {
        df_top <- NULL
      }
      if(len_b>0){
        df_bottom <- screen_subdendrite_starts(start_b, rel_z_layer,round(len_b), VEC, IMG)
      } else {
        df_bottom <- NULL
      }
      
      df_combined <- list(df_top, df_bottom, df_circ_end) %>%
        compact() %>%
        bind_rows() %>%
        return()
      
    }) %>% bind_rows() %>%
      return()    
  }) %>%
    bind_rows() %>%
    return()
}




# SUBD <- bind_rows(df_centers)[11,]
# IMG=IMG_adj


fit_path_to_subd_cluster <- function(SUBD, full_dendrite, det_rad, IMG){
  print(1)
  
  df_dist_raw <- full_dendrite %>%
    rowwise() %>%
    mutate(dist=cppDistPts(x,y,z,as.numeric(SUBD[["x"]]),as.numeric(SUBD[["y"]]),as.numeric(SUBD[["z"]])),
           ha=cppGetHA(x,y,as.numeric(SUBD[["x"]]),as.numeric(SUBD[["y"]])),
           va=cppGetVA(dist, as.numeric(SUBD[["z"]]), z)) %>%
    .[which(.$dist<=4*det_rad),]
  print(2)
  if(nrow(df_dist_raw)==0){return(NULL)}
  
  df_dist <- df_dist_raw %>%
    mutate(int_score=fit_subdendrite_new(SUBD, dist, ha, va, IMG),
           score=dist*(1/int_score^2)) %>%
    .[which(.$int_score!=0),]
  
  if(nrow(df_dist)==0){return(NULL)}
  print(3)
  # final_intersection <- df_dist[which(df_dist$score==max(df_dist$score)),] %>%
  #   .[which(.$dist==min(.$dist)),] %>%
  #   mutate(dist_to_soma=dist_pts(x,y,z,SOMA[["x"]],SOMA[["y"]],SOMA[["z"]]),
  #          x_el=SUBD[["x"]],
  #          y_el=SUBD[["y"]],
  #          z_el=SUBD[["z"]],) 
  
  
  
  final_intersection <- df_dist[which(df_dist$score==min(df_dist$score)),] %>%
    #.[which(.$dist==min(.$dist)),] %>%
    .[1,]
  print(4)
  closest_main_vector <- df_dist_raw[which(df_dist_raw$dist==min(df_dist_raw$dist)),] %>%
    ## sometimes two main vectors have exact same distance....
    .[1,]
  print(5)
  subd_start <- final_intersection %>% select(xs=x,ys=y,zs=z, level, ha, va, dist) %>%
    mutate(xe=SUBD[["x"]],
           ye=SUBD[["y"]],
           ze=SUBD[["z"]],
           clos_vec=closest_main_vector$level,
           d_to_cv=closest_main_vector$dist)
  print(6)
  return(subd_start)
  
}



show_surrounding_layer <- function(all_sorrounding_vox, filename){
  p <- ggplot(all_sorrounding_vox %>% group_by(x,y) %>% summarize(z=max(z)),
              aes(x=x, y=y, color=z))+geom_point(size=0.25)+
    geom_segment(data=bind_rows(main_vectors[[n]]),
                 inherit.aes=F, 
                 aes(x=xs, y=ys, xend=xe, yend=ye))
  
  
  pdf(file=filename, width=50, height=10)
  grid.arrange(p)
  dev.off()
}


show_surrounding_layer_all_layers <- function(all_sorrounding_vox, filename){
  lapply(unique(all_sorrounding_vox$z), function(Z){
    
    p <- ggplot(all_sorrounding_vox %>% filter(z==Z),
                aes(x=x, y=y, color=as.character(i)))+
      geom_point(size=0.25)+
      geom_segment(data=bind_rows(main_vectors[[n]]),
                   inherit.aes=F, 
                   aes(x=xs, y=ys, xend=xe, yend=ye))
    
    
    filename=paste0("f:/data_sholl_analysis/test/clustering/", Z, ".pdf")
    
    pdf(file=filename, width=50, height=10)
    grid.arrange(p)
    dev.off()
    
    
  })

}



trace_dendrite <- function(SUBD, headnode, det_rad, z_range){
  
  #print(111)
  #print(as_tibble(SUBD)) 
  segment_list <- list()
  #df_centers <- tibble(x=SUBD[["x"]], y=SUBD[["y"]], z=SUBD[["z"]])
  df_centers <- SUBD %>% select(x,y,z)
  vec_angle <- adj_deg(SUBD$ha+180)
  c <- 1
  branch_count <- 1
  segment_list[[c]] <- df_centers
  #while(!is.null(df_centers)){
  while(branch_count>0){  
    
    cat(paste("\n    elgt step:", c))
    
    df_circ_end <- screen_circular(det_rad, z_range,
                                   df_centers[["x"]], df_centers[["y"]], df_centers[["z"]], 
                                   vec_angle, binary_image)
    
    centers <- find_subd_cluster(df_circ_end)
    
    if(is.null(centers)){
      ## subdendrite finished
      branch_count <- branch_count-1
      final_vec <- tibble(x=SUBD[["x"]],y=SUBD[["y"]],z=SUBD[["z"]])  %>%
        bind_cols(df_centers %>% select(xe=x, ye=y, ze=z))
      
      ret <- list(sub=F, vec=bind_rows(segment_list))
      
    } else if(length(centers)==1){
      ## subdendrite elonagtes
      
      df_centers <- bind_rows(centers)
      vec_angle <- cppGetHA(df_centers$x, df_centers$y, segment_list[[c]]$x, segment_list[[c]]$y) 
      
      #segment_list[[c]] <- append(segment_list[[c]], list(sub=NULL))
      
      c <- c+1
      
      segment_list[[c]] <- df_centers
    } else {
      ## new intersection
      branch_count <- branch_count-1
      #final_vec <- SUBD
      
      ## calc angles to new subdendrites 
      proc_cent <- lapply(centers, function(C){
        na <- cppGetHA(C[["x"]], C[["y"]], segment_list[[c]]$x, segment_list[[c]]$y)
        return(tibble(x=C[["x"]], y=C[["y"]], z=C[["z"]], ha=na))
      })
      
      ret <- list(sub=T, vec=bind_rows(segment_list), 
                  pos=SUBD, 
                  end=proc_cent)
    }
    
    
  }
  
  return(ret)
  
}


#MAIND_raw <- safe_MAIND

#tMASTER <- MASTER
#main_vectors <- inter$main_vectors
export_structure <- function(tMASTER, main_vectors, output_file){
  
  node_order <- lapply(tMASTER, function(MAIND_raw){   ## main dendrites
    
    #print("new maind")
    
    if(is.null(MAIND_raw)){
      ## no sub node from main dendrite --> export main dendrite vectors
      return(NULL)
    }
    
    MAIND <- lapply(MAIND_raw, function(fl){fl[c("node_id", "subnodes")]})
    
    assigned <- list()
    node_id <- 1
    ids_to_go_next <- list()
    ids_to_return <- list()
    abort <- F
    finished <- list()
    
    c <- 1
    
    while(abort==F){
      #print(c)
      ids_to_go_next <- ids_to_go_next[which(!ids_to_go_next==node_id)]
      next_nodes_raw <- MAIND[[which(unlist(lapply(MAIND, function(x){x$node_id==node_id})))]]$subnodes
      next_nodes <- next_nodes_raw[which(!next_nodes_raw %in% finished)]
      
      if(length(next_nodes)==0){
        ## end of node chain
        if(length(ids_to_go_next)==0){
          abort=T
        } else {
          assigned <- append(assigned, list(tibble(from=node_id, to=first(ids_to_return))))
          finished <- append(finished, node_id)
          node_id <- first(ids_to_return)
          ids_to_return[[1]] <- NULL
          ids_to_return <- compact(ids_to_return)
        }
      } else if(length(next_nodes)==1){
        ## only one subnode
        if(length(ids_to_return)>0){
          ids_to_return <- append(node_id, ids_to_return) %>% unique()
        }
        assigned <- append(assigned, list(tibble(from=node_id, to=next_nodes[[1]])))
        finished <- append(finished, node_id) %>% unique()
        node_id <- next_nodes[[1]]
      } else {
        if(length(ids_to_return)>0){
          ids_to_return <- append(node_id, ids_to_return) %>% unique()
        }
        ## more than one subnode
        if(length(ids_to_go_next)==0){
          ## no subdendrites left to go
          assigned <- append(assigned, list(tibble(from=node_id, to=next_nodes[[1]])))
          ids_to_return <- append(ids_to_return, node_id)
          node_id <- next_nodes[[1]]
          ids_to_go_next <- append(ids_to_go_next, next_nodes[c(2:length(next_nodes))])
        } else {
          ## still subdendrites left
          if(length(ids_to_return)==0){
            assigned <- append(assigned, list(tibble(from=node_id, to=first(ids_to_go_next))))
            ids_to_return <- append(ids_to_return, node_id)
            node_id <- first(ids_to_go_next)
            ids_to_go_next[[1]] <- NULL
            ids_to_go_next <- compact(ids_to_go_next)
          } else {
            assigned <- append(assigned, list(tibble(from=node_id, to=next_nodes[[1]])))
            ids_to_return <- append(node_id, ids_to_return)
            node_id <- next_nodes[[1]]
            ids_to_go_next <- append(ids_to_go_next, next_nodes[c(2:length(next_nodes))])
          }
        }
      }
      c <- c+1
    }
    return(assigned)
  })
  
  ## add vectors according to node order
  
  full_coords <- lapply(1:length(node_order), function(MD){
    
    cat(paste0("\n  export: MD", MD))
    
    order <- node_order[[MD]]
    MAIND <- tMASTER[[MD]]
    
    if(is.null(order)){
      ## no subdenedrites --> use main dendrite
      full_vec_list <- bind_rows(main_vectors[[MD]]) %>% select(x=xe, y=ye, z=ze) %>%
        mutate(id=c(1:nrow(.)))
    } else if(length(order)==0){
      full_vec_list <- bind_rows(main_vectors[[MD]]) %>% select(x=xe, y=ye, z=ze) %>%
        mutate(id=c(1:nrow(.)))
    } else {
      
      #WAY <- tibble(from=18, to=10)
      ## damit alle subdendriten mitgenommen werden... wird der letzte order nochmal invers angefügt
      bound_order <- append(order, list(last(order) %>% mutate(x=to, to=from, from=x) %>% select(-x)))
      
      full_vec_list <- lapply(bound_order, function(WAY){
        
        #print(WAY)
        
        start_node <- MAIND[[which(unlist(lapply(MAIND, function(x){x$node_id==WAY$from})))]]
        
        start_full_subdendrites <- start_node[["subdend_full"]]
        
        ## add dendrites
        subdend <- lapply(start_full_subdendrites, function(DEND){
          
          coords <- DEND[[2]] %>% ungroup() %>% mutate(id=c(1:nrow(.)))
          
          bind_rows(coords, arrange(coords, desc(id))) %>%
            select(x,y,z) %>%
            return()
          
        }) %>% bind_rows()
        
        
        ## now look for next node to go
        
        ## first check if next node is subnode of current node
        
        start_full_subnodes <- start_node[["subnode_full"]]
        
        ## if no subnnode is assigned yet... we need to go to the master node
        if(length(start_full_subnodes)!=0){
          rel <- which(unlist(lapply(start_full_subnodes, function(x){x[[1]]==WAY$to})))
          if(length(rel)==1){
            return(bind_rows(subdend, start_full_subnodes[[rel]]$full_coords))
          }           
        } 
        ## if at this point we havent found any subnode... we probably need to go to the master node
        end_node <- MAIND[[which(unlist(lapply(MAIND, function(x){x$node_id==WAY$to})))]]
        end_full_subnodes <- end_node[["subnode_full"]]
        rel <- which(unlist(lapply(end_full_subnodes, function(x){x[[1]]==WAY$from})))
        if(length(rel)==1){
          
          ## we need to invert the sequnce of the coordinates
          inv <- end_full_subnodes[[rel]]$full_coords %>% mutate(id=c(1:nrow(.))) %>%
            arrange(desc(id)) %>% select(-id)
          return(bind_rows(subdend, inv))
        } else {
          print("no match found???")
        }
      }) %>% bind_rows() %>%
        mutate(id=c(1:nrow(.)))
      
    }
    
    
    full_maind <- bind_rows(full_vec_list, arrange(full_vec_list, desc(id))) %>%
      select(-id, -z) %>%
      bind_rows(select(SOMA, x,y), ., select(SOMA, x,y))
    
    
  }) %>%
    bind_rows()
  
  
  write_csv(full_coords, file=output_file)
  
  
}  

get_linear_equation <- function(VEC){
  
  
  if(VEC[["xe"]]==VEC[["xs"]]){
    ## infinite m
    return(list(m=NULL, b=NULL, 
                xdef=VEC[["xs"]], 
                ydef=c(min(c(VEC[["ye"]],VEC[["ys"]])), max(c(VEC[["ye"]],VEC[["ys"]])))))
    
  }             
  
  if(VEC[["xe"]]<VEC[["xs"]]){
    xs <- VEC[["xe"]]
    xe <- VEC[["xs"]]
    ys <- VEC[["ye"]]
    ye <- VEC[["ys"]]
  } else {
    xs <- VEC[["xs"]]
    xe <- VEC[["xe"]]
    ys <- VEC[["ys"]]
    ye <- VEC[["ye"]]
  }
  
  
  
  m <- (ye-ys)/(xe-xs)
  
  b <- ys-m*xs
  
  return(list(m=m, b=b, xdef=c(min(c(xs, xe)),max(c(xs, xe))),
              ydef=c(min(c(ys, ye)),max(c(ys, ye)))))
  
}

check_for_intersections <- function(n, linear_equations){
  #print(n)
  #f1 <- linear_equations[[n]]
  f <- linear_equations[[n]]
  intersections <- c(1:length(linear_equations)) %>% 
    .[which(.!=n)] %>%
    lapply(function(i){
      is <- intersect_lineq(f, linear_equations[[i]]) 
      
      if(is.null(is)){return(NULL)}
      return(mutate(is, v2=i))
    }) %>%
    compact() %>%
    bind_rows() %>%
    mutate(v1=n)
  if(nrow(intersections)==0){
    return(NULL)
  }
  in_def_raw <- intersections %>%
    ## in my vector def
    filter(between(x, min(f$xdef), max(f$xdef))&between(y, min(f$ydef), max(f$ydef))) #%>%
  
  if(nrow(in_def_raw)==0){
    return(NULL)
  }
  
  in_def <- apply(in_def_raw,1,function(INT){
    line <- linear_equations[[INT[["v2"]]]]
    if(between(INT[["x"]], min(line$xdef), max(line$xdef))&between(INT[["y"]], min(line$ydef), max(line$ydef))){
      return(INT)
    } else {
      return(NULL)
    }
  }, simplify=F) %>%
    compact() %>%
    bind_rows()
  
  # 
  # if(is.null(f$m)){
  #   in_def <- intersections %>%
  #     filter(between(y, f$def[1], f$def[2]))
  # } else {
  #   in_def <- intersections %>% rowwise() %>%
  #     filter(between(x, f$def[1], f$def[2]))
  # }
  
  if(nrow(in_def)==0){
    return(NULL)
  } else {
    return(in_def %>% mutate_at(.vars = c("x","y"), .funs = round))
  }
  
}


intersect_lineq <- function(f1, f2){
  
  
  if((is.null(f1$m)&is.null(f2$m))){
    return(NULL)
  }
  
  
  
  if(is.null(f1$m)){
    yi <- f2$m*f1$x+f2$b
    return(tibble(x=f1$x, y=yi))
  }
  
  if(is.null(f2$m)){
    yi <- f1$m*f2$x+f1$b
    return(tibble(x=f2$x, y=yi))
  }
  
  if(f1$m==f2$m){
    return(NULL)
  }
  xi <- (f1$b-f2$b)/(f2$m-f1$m)
  yi <- f1$m*xi+f1$b
  
  return(tibble(x=xi, y=yi))
}

paste_ordered <- function(v1, v2){
  paste(min(v1, v2), max(v1, v2), sep="_")
}


# main_vectors=inter$main_vectors
# main_vectors_df=inter$main_vectors_df
# main_vectors_full=inter$main_vectors_full
# SOMA=SOMA
# IMG_screen=images$med_bin_noS_noMD_image
# IMG_adj=images$bin_noS_image_p
# soma_radius=inter$soma_radius
# minPTS=opt$subd_cluster_mpt
# EPS=opt$subd_cluster_eps
# det_rad=opt$subd_detection_distance
# z_range=opt$subd_detection_vertical_range
# n <- 1
# depth=opt$subd_detection_depth

find_subdendritic_starts <- function(n, main_vectors, main_vectors_df,main_vectors_full, 
                                     SOMA, IMG_screen, IMG_adj, soma_radius, minPTS, EPS, det_rad, z_range,
                                     depth){
  
  cat(paste0("\n  process: MD:", n))
  #cat(paste0("\n    screening"))
  # det_rad <- 30
  # z_range <- 10
  
  rel_vecs_raw <- main_vectors[[n]]
  full_vecs_raw <- main_vectors_full[[n]]
  
  rel_vecs <- rel_vecs_raw[c(2:length(rel_vecs_raw))]
  full_vecs <- full_vecs_raw[c(2:length(full_vecs_raw))]
  
  full_dendrite <- bind_rows(full_vecs) %>% 
    mutate_all(.funs = as.numeric)
  
  all_sorrounding_vox <- define_surr_layer(rel_vecs, full_vecs, det_rad, z_range, IMG_screen, depth) %>%
    .[which(.$dist_to_soma>2.5*soma_radius),]
  #print(1)
  
  #show_surrounding_layer_all_layers(all_sorrounding_vox, "f:/data_sholl_analysis/test/full_sorrounding.pdf")
  ## orig... withour neighboring voxels... eps=3.5 minPTS=15
  #cat(paste0("\n    clustering"))
  df_centers <- find_subd_cluster_man(all_sorrounding_vox, EPS, minPTS)
  
  print(2)
  ## if no subdendrite starts are found... the main dendrite is just a subdendrite
  if(is.null(df_centers)){
    #full_branches
    return(NULL)
  }
  retraced_subd <- apply(bind_rows(df_centers), 1, fit_path_to_subd_cluster, 
                                     full_dendrite=full_dendrite,
                                     det_rad=det_rad,
                                     IMG=IMG_adj) %>% 
    compact() 
  print(3)
  ## remove duplicated subdendirte starts
  
  if(is.null(retraced_subd)){return(NULL)}
  rescored_subdendrites_raw <- lapply(1:length(retraced_subd), function(nDEND){
#print(nDEND)
    DEND <- retraced_subd[[nDEND]]

    sel_vec <- c(1:length(retraced_subd)) %>% .[which(.!=nDEND)]

    dist <- lapply(sel_vec, function(nV){
      V2 <- retraced_subd[[nV]]
      d <- cppDistPts(V2[["xe"]],V2[["ye"]],V2[["ze"]], DEND[["xe"]],DEND[["ye"]],DEND[["ze"]])
      return(c(v=nV, dist=d))
    }) %>% bind_rows() %>%
      .[which(.$dist<15),]
    
    
    
    if(nrow(dist)==0){
      return(DEND)
    } else {#if(nrow(dist)==1){
      ## check which one is the smaller one
       retraced_subd[c(dist$v, nDEND)] %>% 
         bind_rows() %>%
         filter(dist==min(.$dist)) %>%
         #.[which(dist==min(.$dist)),] %>% 
         return()
    } 
    
  }) %>% unique()
  
  print(4)
  #cat(paste0("\n    adjusting"))
  
  ## check intersecting subdendrites
  
  if(is.null(rescored_subdendrites_raw)){return(NULL)}
  
  linear_equations <- lapply(rescored_subdendrites_raw, get_linear_equation)
  
  intersections_raw <- lapply(1:length(linear_equations), check_for_intersections, 
                              linear_equations=linear_equations) %>%
    compact() #%>%
  #print(3)
  
  intersections_raw <- list()
  
  if(length(intersections_raw)!=0){
    
    intersections <- intersections_raw %>%
      bind_rows() %>%
      rowwise() %>%
      mutate(id=paste_ordered(v1, v2)) %>%
      #mutate(id=paste(min(c(v1, v2)),max(c(v1, v2)), sep="_")) %>%
      group_by(id) %>%
      summarize(x=round(median(x)),
                y=round(median(y)),
      ) %>% 
      mutate(v1=str_split(id, "_") %>% map_chr(.,1) %>% as.numeric(),
             v2=str_split(id, "_") %>% map_chr(.,2)%>% as.numeric()) %>%
      rowwise() %>%
      mutate(v1_angle=select_angle(v1, all_vectors = rescored_subdendrites_raw, col="ha"),
             v2_angle=select_angle(v2, all_vectors = rescored_subdendrites_raw, col="ha"),
             z1=select_angle(v1, all_vectors = rescored_subdendrites_raw, col="ze"),
             z2=select_angle(v2, all_vectors = rescored_subdendrites_raw, col="ze"),
             dist1=select_angle(v1, all_vectors = rescored_subdendrites_raw, col="dist"),
             dist2=select_angle(v2, all_vectors = rescored_subdendrites_raw, col="dist"),
             angle_diff=ifelse(abs(v1_angle-v2_angle)>180,
                               360-abs(v1_angle-v2_angle),
                               abs(v1_angle-v2_angle))) %>%
      filter(angle_diff<90,
             abs(z1-z2)<8)
    
    
    ## now create new list of duplicated subdendrites
    
    print(4)
    #rel_vec_linear_eq <- lapply(rel_vecs, get_linear_equation)
    if(nrow(intersections)!=0){
      adjusted_subdendrite_starts <- apply(intersections, 1, function(PAIR){
        
        #print(PAIR)
        v1 <- rescored_subdendrites_raw[[as.numeric(PAIR[["v1"]])]]
        v2 <- rescored_subdendrites_raw[[as.numeric(PAIR[["v2"]])]]
        
        ## if close to each other... do what i did
        
        dist_pp <- cppDistPts(v1$xe, v1$ye, v1$ze, v2$xe, v2$ye, v2$ze)
        
        if(dist_pp<25){
          ha <- adj_deg(min(as.numeric(PAIR[["v1_angle"]]), as.numeric(PAIR[["v2_angle"]]))+as.numeric(PAIR[["angle_diff"]]))
          xe <- as.numeric(PAIR[["x"]]) 
          ye <- as.numeric(PAIR[["y"]])
          ze <- round(mean(as.numeric(PAIR[["z1"]]),as.numeric(PAIR[["z2"]])))
          ## check which one has the lower distance to new node
          d1 <- cppDistPts(xe,ye,ze, v1$xs, v1$ys, v1$zs)
          d2 <- cppDistPts(xe,ye,ze, v2$xs, v2$ys, v2$zs)
          if(d1<d2){
            start <- v1
            d <- d1
          } else {
            start <- v2
            d <- d2
          }
          xs <- start$xs
          ys <- start$ys
          zs <- start$zs
          va <- cppGetVA(d, ze, zs)
          
          return(list(
            tibble(xs=xs,
                   ys=ys,
                   zs=zs,
                   level=start$level,
                   clos_vec=start$clos_vec,
                   d_to_cv=start$d_to_cv,
                   ha=ha,
                   va=va,
                   dist=d,
                   xe=xe,
                   ye=ye,
                   ze=ze),
            c(PAIR[["v1"]], PAIR[["v2"]])
          )
          )
          
        } else {
          ## distance larger 25 .. check connection
          ## load intesnsities of in between voxels
          ha_pp <- cppGetHA(v2$xe, v2$ye, v1$xe, v1$ye)
          va_pp <- cppGetVA(dist_pp, v2$ze, v1$ze)
          
          int <- elongate_line_full(IMG_adj, ha_pp, va_pp, v1$xe, v1$ye, v1$ze, dist_pp)
          
          score_mn <- mean(int$i)
          score_gap <-  int %>% ungroup() %>%
            mutate(cumu=cumsum(i)) %>%
            group_by(cumu) %>%
            summarize(c=length(i), mn=min(n)) %>%
            filter(c==max(.$c)) %>%
            pull(c)
          
          if(score_gap<10&score_mn>0.5){
            ## connected ..  probably same dendrite
            ## keep, but connect ends
            ## remove longer one
            if(v1$dist>v2$dist){
              #to_rem <- v1
              to_keep <- v1
            } else {
              #to_rem <- v2
              to_keep <- v2
            }
            
            return(list(to_keep, c(PAIR[["v1"]], PAIR[["v2"]])))
            
          } else {
            ## not connected
            ## no node adjustment
            return(NULL)
          }
          
        }
        
      }) %>% compact()
      
      
      to_keep <- c(1:length(rescored_subdendrites_raw)) %>%
        .[which(!. %in% as.numeric(lapply(adjusted_subdendrite_starts, last) %>% unlist()))]
      
      rescored_subdendrites <- append(rescored_subdendrites_raw[to_keep], lapply(adjusted_subdendrite_starts, first))  
    } else {
      rescored_subdendrites <- rescored_subdendrites_raw
    }
  } else {
    rescored_subdendrites <- rescored_subdendrites_raw
  }  
  #print(5)
  ## get hierarchy
  #cat(paste0("\n    structurize"))
  df_poslev <- rescored_subdendrites %>%
    bind_rows() %>%
    rowwise() %>%
    mutate(dist_to_soma=cppDistPts(xs,ys,zs, round(SOMA[["x"]]),round(SOMA[["y"]]),round(SOMA[["z"]]))) %>%
    left_join(bind_rows(main_vectors[[n]]) %>% select(level, dir), by="level")
  
  
  df_sorted <- lapply(sort(unique(df_poslev$level)), function(LEV){
    rel_rows <- df_poslev %>%
      filter(level==LEV)
    
    if(unique(rel_rows$dir)==1){
      return(arrange(rel_rows, dist_to_soma))
    } else {
      return(arrange(rel_rows, desc(dist_to_soma)))
    }
    
  }) %>% bind_rows() %>% ungroup() %>%
    mutate(node=c(1:nrow(.)))
  
  
  #idNODE <- 12
  #print(6)
  assigned_nodes <- lapply(df_sorted$node, function(idNODE){
    
    node <- df_sorted %>%
      filter(node==idNODE)
    
    #master_node <- idNODE-1
    
    sub_dend <- node %>% select(x=xe, y=ye, z=ze, ha, va, cut_vec=level, clos_vec, d_to_cv)  %>%
      mutate(class="dend",
             sub_id=NA,
             ha=adj_deg(ha+180),
             orient="side")
    
    sub_dend_full <- bind_rows(select(node, x=xs,y=ys,z=zs),
                               select(node, x=xe,y=ye,z=ze))
    
    
    full_subdend_info <- list(list(info=sub_dend, full_coords=sub_dend_full))
    
    if(idNODE==nrow(df_sorted)){
      ## no know subnodes so far ... just elongation directions
      #sub_node <- NULL
      subnode_id <- NULL
      #full_vec_pos <- NULL
      full_subnode_info <- NULL
      rel_vecs_df <- bind_rows(rel_vecs_raw)

      last_vec_info <-rel_vecs_df %>% 
        filter(level==max(.$level)) %>%
        select(x=xe, y=ye, z=ze, ha, va) %>%
        mutate(d_to_cv=cppDistPts(x,y,z, node$xs, node$ys, node$zs)) %>%
        bind_cols(node %>% select(cut_vec=level, clos_vec)  %>% rowwise() %>%
        mutate(class="dend",
               sub_id=NA,
               orient="end"))
      
      
      vecs_to_add <- node$level:length(rel_vecs_raw)
      
      last_vec_coords <- lapply(vecs_to_add, function(V){
        VEC <- rel_vecs_df %>% filter(level==V) %>% select(x=xe, y=ye, z=ze) %>% return()
      }) %>%
        bind_rows() %>%
        bind_rows(select(node, x=xs,y=ys,z=zs),.)
      
      full_subdend_info <- append(full_subdend_info, 
                                  list(list(info=last_vec_info, 
                                            full_coords=last_vec_coords)))
    } else {
      next_node <- df_sorted %>% filter(node==idNODE+1)
      subnode_id <- next_node$node
      sub_node <- next_node %>%
        select(x=xs, y=ys, z=zs, sub_id=node) %>%
        rowwise() %>%
        mutate(ha=cppGetHA(x,y,node$xs, node$ys),
               dist=cppDistPts(x,y,z,node$xs, node$ys, node$zs),
               va=cppGetVA(dist, node$zs, z),
               class="node") %>%
        select(-dist)
      
      if(next_node$level!=node$level){
        ## vector edges need to be inserted
        full_vec_pos <- main_vectors_df[[n]] %>% #rowwise() %>%
          filter(between(level, node$level, as.numeric(next_node$level-1))) %>%
          select(x=xe, y=ye, z=ze) %>%
          bind_rows(select(node, x=xs, y=ys, z=zs),
                    .,
                    select(next_node, x=xs, y=ys, z=zs))
      } else {
        
        full_vec_pos <- bind_rows(select(node, x=xs, y=ys, z=zs),
                                  select(next_node, x=xs, y=ys, z=zs))
        
      }
      
      full_subnode_info <- list(node_id=idNODE+1, 
                                info=sub_node, 
                                full_coords=full_vec_pos)
      
    }
    
    
    final_list <- list(node_id=node$node,
                       master_id=idNODE-1,
                       subnodes=list(subnode_id) %>% compact(),
                       pos=node %>% select(x=xs, y=ys, z=zs),
                       subnode_full=list(full_subnode_info) %>% compact(),
                       subdend_full=full_subdend_info
    )
    return(final_list)
    
    
  })
  
  return(assigned_nodes)
  
}


get_all_subdendrite_vectors <- function(D){
  
  sd <- lapply(D$subdend_full, function(S){
    vec <- lapply(2:nrow(S$full_coords), function(V){
      tibble(xe=S$full_coords[[V,"x"]],
             ye=S$full_coords[[V,"y"]],
             ze=S$full_coords[[V,"z"]],
             xs=S$full_coords[[V-1,"x"]],
             ys=S$full_coords[[V-1,"y"]],
             zs=S$full_coords[[V-1,"z"]],
             ha=S$info$ha)
    }) %>%
      #mutate(ha=S$info$ha) %>% 
      return()
  }) %>% Reduce(function(x,y)append(x,y),.) %>%
    return()
  
}


get_all_vectors <- function(D){
  
  ## for subnodes
  sn <- lapply(D$subnode_full, function(S){
    vec <- lapply(2:nrow(S$full_coords), function(V){
      
      tibble(xe=S$full_coords[[V,"x"]],
             ye=S$full_coords[[V,"y"]],
             ze=S$full_coords[[V,"z"]],
             xs=S$full_coords[[V-1,"x"]],
             ys=S$full_coords[[V-1,"y"]],
             zs=S$full_coords[[V-1,"z"]],)
      
    }) %>%
      return()
  }) %>% Reduce(function(x,y)append(x,y),.)
  
  sd <- lapply(D$subdend_full, function(S){
    vec <- lapply(2:nrow(S$full_coords), function(V){
      tibble(xe=S$full_coords[[V,"x"]],
             ye=S$full_coords[[V,"y"]],
             ze=S$full_coords[[V,"z"]],
             xs=S$full_coords[[V-1,"x"]],
             ys=S$full_coords[[V-1,"y"]],
             zs=S$full_coords[[V-1,"z"]],)
    }) %>%
      return()
  }) %>% Reduce(function(x,y)append(x,y),.)
  sf <- append(sn, sd) %>%
    return()
}


format_new_subdendrite_starts <- function(CENT_raw, nMD, main_vectors_full, last_pos, new_node_id){
  
  ## insert new ha cutoff here
  CENT <- CENT_raw %>% as_tibble()
  ## calc distance to closest main vector
  ## i just do it by taking all pts and calc all dists
  full_vecs_raw <- main_vectors_full[[nMD]]
  full_vecs <- full_vecs_raw[c(2:length(full_vecs_raw))]
  
  closest_main_vec <- bind_rows(full_vecs) %>% 
    mutate_all(.funs = as.numeric) %>%
    mutate(dist=cppDistPts(x,y,z,CENT$x, CENT$y, CENT$z)) %>%
    .[which(.$dist==min(.$dist)),] %>%
    .[1,]
  
  
  info <- CENT %>% 
    mutate(ha=cppGetHA(x,y, last_pos$x, last_pos$y),
           va=cppGetVA(cppDistPts(x,y,z, last_pos$x, last_pos$y, last_pos$z), z, last_pos$z),
           class="node",
           sub_id=new_node_id,
           clos_vec=closest_main_vec$level,
           d_to_cv=closest_main_vec$dist,
           orient="elongated"
    )
  
  full_coords <- bind_rows(last_pos,
                           CENT)
  
  
  return(list(info=info, full_coords=full_coords))
}

get_vox_of_borders <- function(n, main_dendrites, nc_orig, nr_orig, IMG){
  
  if(n!=nrow(main_dendrites)){
    h_angle <- 0.5*sum(main_dendrites$h_angle[n], main_dendrites$h_angle[n+1])
  } else {
    raw <- main_dendrites$h_angle[n]+0.5*(abs(360-main_dendrites$h_angle[n])+main_dendrites$h_angle[1])
    if(raw>360){
      h_angle <- raw-360
    } else {
      h_angle <- raw
    }
  }
  elongate_3d_sphere(IMG, h_angle, 0, SOMA[["x"]], SOMA[["y"]], SOMA[["z"]], round(sqrt(nc_orig^2+nr_orig^2))) %>%
    select(x,y) %>%
    return()
  
}


binarize_image <- function(opt, inter, n_main_dendrites, IMG, 
                           elongated_dendrites, raster_size, process_full, method){
  
  borders_vox_raw <- lapply(1:n_main_dendrites, get_vox_of_borders, 
                            main_dendrites=inter$main_dendrites,
                            nc_orig=opt$nc_orig,
                            nr_orig=opt$nr_orig,
                            IMG=images$raw_image)
  
  selection_vector <- rep(c(1:n_main_dendrites), 3) %>% 
    set_names(c((-n_main_dendrites+1):(2*n_main_dendrites)))
  
  dendrite_segments <- lapply(1:n_main_dendrites, create_raster, avs=raster_size,
                              borders_vox_raw=borders_vox_raw,
                              selection_vector=selection_vector, elongated_dendrites,
                              nr_orig=opt$nr_orig,
                              nc_orig=opt$nc_orig,
                              process_full=process_full)
  
  sl <- Reduce(function(x,y)append(x,y),dendrite_segments)
  
  if(method=="man"){
    binary_image <- bi2(sl, IMG)
  } else {
    qnt <- str_match(method, "\\d+") %>% as.numeric()
    binary_image <- bi3(sl, IMG, qnt)
  }
  
  
  return(binary_image)
}

rescore_traced_subnodes <- function(centers, last_pos, RESC_DIST, RESC_ANG){
  
  
  
  df_centers <- bind_rows(centers) %>%
    rowwise() %>%
    mutate(dist=cppDistPts(x,y,z, last_pos$x, last_pos$y, last_pos$z),
           ha=cppGetHA(x,y,last_pos$x, last_pos$y)) %>%
    ungroup() %>%
    mutate(id=c(1:nrow(.)))
  
  dups <- lapply(df_centers$id, function(nC){
    df_dist <- lapply(df_centers$id %>% .[which(!.==nC)],
                      function(C){
                        d <- cppDistPts(centers[[nC]][["x"]],
                                        centers[[nC]][["y"]],
                                        centers[[nC]][["z"]],
                                        centers[[C]][["x"]],
                                        centers[[C]][["y"]],
                                        centers[[C]][["z"]])
                        
                        diff_ha=abs(df_centers[[which(df_centers$id==nC),"ha"]]-
                                      df_centers[[which(df_centers$id==C),"ha"]])
                        
                        c(v1=min(nC, C), 
                          v2=max(nC, C), 
                          dist=d,
                          diff_ha=diff_ha) %>%
                          return()
                        
                      }) %>%
      bind_rows() %>%
      return()
  }) %>% unique() %>%
    bind_rows() %>%
    #filter(dist<RESC_DIST&(diff_ha<RESC_ANG|diff_ha>(360-RESC_ANG)))
    .[which(.$dist<RESC_DIST&(.$diff_ha<RESC_ANG|.$diff_ha>(360-RESC_ANG))),]
  if(nrow(dups)>0){
    
    
    remove <- apply(dups, 1, function(PAIR){
      
      df_centers %>%
        .[which(.$id %in% c(PAIR[["v1"]], PAIR[["v2"]])),] %>%
        .[which(!.$dist==max(.$dist)),] %>%
        pull(id) %>%
        return()
    }) %>% unlist()
    
    keep <- df_centers %>%
      filter(!id %in% remove) %>%
      pull(id)
    #.[[which(!df_centers$id %in% remove),"id"]]
    
    centers <- centers[keep]
  }
  
  return(centers)
  
}

get_circle_around_point <- function(P, RAD, nr_orig){
  
  df_circ_full <- create_circle(RAD)
  
  all_vox <- lapply(c(min(df_circ_full$x):max(df_circ_full$x)), function(X){
    
    all_x <- df_circ_full[which(df_circ_full$x==X),] #%>%
    #arrange(y)
    return(tibble(x=X, y=c(min(all_x$y):max(all_x$y))))
    
  }) %>% bind_rows() %>%
    mutate(x=x+P[["x"]],
           y=y+P[["y"]],
           short=get_single_index(x,y,nr_orig)) %>%
    pull(short) %>%
    return()
}

get_sphere_around_point <- function(P, RAD, nr_orig, nz_orig){
  
  lapply(c((P[["z"]]-RAD):(P[["z"]]+RAD)) %>% .[which(between(., 1, nz_orig))], function(Z){
    get_circle_around_point(P, RAD, nr_orig) %>% paste(Z, ., sep="_")
  }) %>% c(recursive=T)
  
}

# 
nMD=2
IMG=images$bin_noS_noMD_image
MASTER=MASTER
main_vectors_df=inter$main_vectors_df
main_vectors_full=inter$main_vectors_full
EPS_orig=EPS_orig
MPTS_orig=MPTS_orig
INC=INC_orig
det_rad=DR_orig
RESC_DIST=opt$trace_rescore_dist
RESC_ANG=opt$trace_rescore_angle
nr_orig=opt$nr_orig


trace_subdendrites <- function(nMD, MASTER, main_vectors_df, main_vectors_full, 
                               IMG, EPS_orig, MPTS_orig, INC, det_rad, RESC_DIST, RESC_ANG, nr_orig){ ## main dendrites
  
   # nMD <- 2
   # nND <- 1
   # nSD <- 1
  
  MV <- main_vectors_df[[nMD]]
  MAIND <- MASTER[[nMD]]
  
  if(is.null(MAIND)){
    return(NULL)
  }
  
  ## for testing ... remove nodes at end of dendrite...
  MAIND <- MAIND[c(1:16)]
  MAIND[[16]]$subnodes <- list()
  MAIND[[16]]$subnode_full <- list()
  
  
  
  
  nNODE_orig <- length(MAIND)
  
  all_vectors <- lapply(MAIND, get_all_vectors) %>%
    Reduce(function(x,y)append(x,y),.)
  
  all_vectors_fv <- lapply(all_vectors, get_full_vector_voxels, IMG=IMG)
  
  
  all_points <- lapply(all_vectors, function(V){
    return(list(V[1:3] %>% set_names(c("x", "y", "z")), V[4:6]%>% set_names(c("x", "y", "z"))))
  }) %>% 
    Reduce(function(x,y)append(x,y),.) %>%
    unique() %>%
    lapply(get_sphere_around_point, RAD=5, nr_orig=nr_orig, nz_orig=length(IMG))
  
  
  nND <- 0
  n_processed <- 0
  #print(1)
  #nND <- 3
  while(n_processed < length(MAIND)){ ## nodes
    
  ### nur main nodes:  
  #while(n_processed < 16){
    
      
    nND <- nND + 1
    
    cat(paste("\n", nND, "of", length(MAIND)))
    
    NODE <- MAIND[[nND]]
    
    SUBD_list <- NODE$subdend_full
    
    new_subnodes <- NODE$subnodes
    new_SUBD_list <- list()
    new_SUBN_list <- NODE$subnode_full
    #print(2)
    #nSD <- 1
    for(nSD in 1:length(SUBD_list)){
      
      cat(paste("\n  subd:", nSD, "of", length(SUBD_list)))
      SUBD <- SUBD_list[[nSD]]
      SUBD_info <- SUBD$info
      
      #det_rad <- 12
      #z_range <- 5
      
      xs <- SUBD_info$x
      ys <- SUBD_info$y
      zs <- SUBD_info$z
      ha <- adj_deg(SUBD_info$ha)
      ## only define ha_cutoff if node intersects with main vectors
      

        # ha_cutoff <- MV %>%
        #   filter(level==SUBD_info$clos_vec) %>%
        #   pull(ha)
      # }
      
      
      coord_list <- list()
      c <- 1
      abort <- F
      node_detected <- F
      remove_dend <- F
      z_cutoff <- 12
      
      
      while(abort==F){ 
      #for(i in c(1:3)){
          #print(3)
          #st <- Sys.time()
        cat(paste("\n    elgt step:",c))
        if(SUBD_info$orient=="end"){
          ha_cutoff <- NULL
          screening_angle <- 180
          z_range <- 5
          det_rad <- 20
          alt <- T
        } else if(SUBD_info$orient=="side"&c==1){
          screening_angle <- 220
          z_range <- 8
          ha_cutoff <- NULL
          alt <- T
        } else if(SUBD_info$orient=="elongated"&c==1){
          screening_angle <- 130
          z_range <- 5
          ha_cutoff <- NULL
        } else if(c<3){
          ha_cutoff <- NULL
          screening_angle <- 120
          z_range <- 3
          alt <- F
        } else {
          ha_cutoff <- NULL
          screening_angle <- 120
          z_range <- 3
          alt <- F
        }
        #for(i in 1:4){  
        # t1 <- Sys.time()
        
        Z <- round(zs+tan(cppdeg2rad(SUBD_info$va))*det_rad)
        # if(Z>z_cutoff){
        #   cat(paste("  Z:", Z))
        #   Z <- z_cutoff
        #   
        # } else if(Z<-z_cutoff){
        #   cat(paste("  Z:", Z))
        #   Z <- -z_cutoff
        #   
        # }
###########################################
        
        cluster_ok=F
        DR <- det_rad
        ZR <- z_range
        MPTS <- MPTS_orig
        EPS <- EPS_orig
        while(cluster_ok==F){
        
          rel_z_layer <- seq(Z-round(ZR), Z+round(ZR), 1) %>%
              .[between(., 1, length(IMG))]
  
          if(is.null(screening_angle)){
            bmin <- 0
            bmax <- 360
          } else {
            
            bmin=ha-0.5*screening_angle
            bmax=ha+0.5*screening_angle
            if(!is.null(ha_cutoff)){
              if(ha<ha_cutoff){
                if(bmax>ha_cutoff){
                  bmax <- round(ha_cutoff)
                }
                
                if(bmin<ha_cutoff-180){
                  bmin <- round(ha_cutoff-180)
                }
              } else {
                if(bmax>ha_cutoff+180){
                  bmax <- round(ha_cutoff+180)
                }
                if(bmin<ha_cutoff){
                  bmin <- round(ha_cutoff)
                }
              }
            }
          }
          
          
          
          
          centers1 <- screen_circular_new(DR, rel_z_layer, bmin, bmax,
                                            INC, xs, ys, zs, ha, ha_cutoff, screening_angle, IMG) #%>%
          if(nrow(centers1)==0){
            centers2 <- NULL
            cluster_ok <- T
          } else {
            
            
              ## remove positions of other dendrites
            centers_rem <- centers1 %>%
              mutate(index=get_3d_single_index(x,y,z, nr_orig)) %>%
                  .[which(!.$index %in% c(all_points, recursive=T)),] 
            if(nrow(centers_rem)==0){
              centers2 <- NULL
              cluster_ok <- T
            } else {
              
            
              df_filtered <- centers_rem %>% 
                .[which(.$i!=0),] %>%
                select(x,y,z) #
            #ggplot(df_filtered, aes(x=x, y=-y,))+geom_point()+facet_wrap(~z)
              if(nrow(df_filtered)==0){
                centers2 <- NULL
                cluster_ok <- T
              } else {
                
                df_clustered_raw <- df_filtered %>% 
                  ungroup() %>%
                  mutate(c=fpc::dbscan(df_filtered, eps = EPS, MinPts = MPTS)$cluster) %>%
                  #mutate(c=fpc::dbscan(df_filtered, eps = sqrt(3), MinPts = 27)$cluster) %>%
                  .[which(.$c!=0),] #%>%
                
                if(nrow(df_clustered_raw)==0){
                  centers2 <- NULL
                  cluster_ok <- T
                } else {
                
                  df_clustered <- df_clustered_raw %>%
                    rowwise() %>%
                    mutate(ha_man=adj_deg(cppGetHA(x,y, xs, ys)-bmin))  
                  #ggplot(df_clustered, aes(x=x, y=-y,color=as.character(c)))+geom_point()+facet_wrap(~z)
                  
                  if(nrow(df_clustered)==0){
                    centers2 <- NULL
                    cluster_ok <- T
                  } else {
                    
                    ## check if clusters are ok:
                    spreading <- lapply(unique(df_clustered$c), function(CLUST){
                      clu <- df_clustered %>%  
                        .[which(.$c==CLUST),] #%>%
                      
                      mn_ha <- min(clu$ha_man)
                      mx_ha <- max(clu$ha_man)
                      if(abs(mn_ha-mx_ha)>140){
                          return(1)
                      } else {
                          return(NULL)
                      } 
                      
                    }) %>% compact()
      
                  if(length(spreading)==0){
                    cluster_ok <- T
                  } else {
                    DR <- DR+3
                    ZR <- ZR+0.2
                    #MPTS <- MPTS+4
                  }
                  
                }
                
                  centers2 <- lapply(unique(df_clustered$c), function(CLUST){
                      kmeans(df_clustered %>%  
                             .[which(.$c==CLUST),] %>%
                             select(x,y,z), 
                           centers=1)$centers %>%
                      round() %>%
                      set_names(nm=c("x","y","z"))
                  })  
                } 
              }
              if(cluster_ok==F){
                cat("  - cluster not valid")
              }
            }
          }
        }        

        #ggplot(df_clustered, aes(x=x, y=-y,color=as.character(c)))+geom_point()+facet_wrap(~z)
      
        
        
        
        
        
        
#############################################
        
        ## we need this while loop if no good clusters are found to adjust the screening parameters
        # cluster_ok <- F
        # cc <- 1
        # while(cluster_ok==F&cc<3){
        # cc <- cc+1
          # 

          
            #centers2 <- find_subd_cluster_man(centers_rem, EPS, MPTS, alt)
            
            if(length(centers2)==0|c>20){
              ## no elongation detected
              abort=T
              cat("\n    stopped")
              if(c==1){
                remove_dend <- T
              }
              only_one=F
            } else if(length(centers2)>1){
              
              if(c==1){
                last_pos <- SUBD$full_coords[nrow(SUBD$full_coords),]
              } else {
                last_pos=as_tibble(last(coord_list))
              }
              
              
            centers <- rescore_traced_subnodes(centers2, last_pos, RESC_DIST, RESC_ANG)
              
              
              if(length(centers)==1){
                cat("\n    removed duplicated")
              } else {
                ## node detected
                
                ## check if its realy a node or just duplicated
                ## calc distance and angle to master node
                cat("\n    node")
                node_detected <- T
                abort <- T
              }
              
              only_one <- ifelse(length(centers)==1, T, F)
            } else {
              only_one <- T
              centers <- centers2
              ## elongate... only one center
            }          
          
            
            vox_to_add <- lapply(centers, get_sphere_around_point, RAD=5, nr_orig=nr_orig, nz_orig=length(IMG))
          
            
            all_points <- append(all_points, vox_to_add)
        #     if(alt==T&nrow(centers)>1){
        #       cluster_ok <- F
        #       det_rad <- det_rad
        #       z_range <- z_range+1
        #     } else {
        #       cluster_ok <- T
        #     }
        #     
        # } ## while for cluster ok
        

           
          
          
           
            
            
          if(only_one==T){
            ## subdendrite elongates
            
            elgt <- centers[[1]]
            dist <- cppDistPts(elgt[["x"]], elgt[["y"]], elgt[["z"]], xs, ys,zs)
            ha <- cppGetHA(elgt[["x"]], elgt[["y"]], xs,ys)
            va <- cppGetVA(dist, 
                         elgt[["z"]], zs)
            
            ## check if elongation overlaps with already present dendrite
            vox <- elongate_line_full(IMG, ha, va, xs, ys, zs, round(dist)) %>%
              mutate(id=paste(x,y, sep="_")) %>%
              pull(id)
            
            overlap <- all_vectors_fv %>%
              bind_rows() %>%
              mutate(id=paste(x,y, sep="_")) %>%
              filter(id %in% vox)
            
            if(nrow(overlap)>0&c>1){
              abort <- T
              cat("\n    crossed")
            } else {
              xs <- elgt[["x"]]
              ys <- elgt[["y"]]
              zs <- elgt[["z"]]    
              coord_list[[c]] <- elgt
              c <- c+1
            } 
          } 
        #}
      }
      
      if(node_detected==T){
        
        new_node_id <- length(MAIND)+1
        
        if(c==1){
          last_pos <- SUBD$full_coords[nrow(SUBD$full_coords),]
        } else {
          last_pos=as_tibble(last(coord_list))
        }
        ## remove dendrite from subdendrite list in current node
        ## by just not adding it to new list
        
        ## add subnode to subnode list of current node
        
        new_subnodes <- append(new_subnodes, new_node_id)
        
        ## add full subnode to full subnode list of current node
        
        new_SUBN <- list(node_id=new_node_id,
                         info=SUBD$info, 
                         full_coords=bind_rows(SUBD$full_coords, bind_rows(coord_list)))
        
        new_SUBN_list <- append(new_SUBN_list, list(new_SUBN))
        
        ## add new node to full node list
        
        proc_centers <- lapply(centers, format_new_subdendrite_starts, nMD=nMD, 
                               main_vectors_full=main_vectors_full, 
                               last_pos=last_pos,
                               new_node_id=new_node_id)
        
        
        new_NODE <- list(node_id=new_node_id,
                         master_id=NODE$node_id,
                         subnodes=list(),
                         pos=last_pos,
                         subnode_full=list(),
                         subdend_full=proc_centers)
        
        MAIND <- append(MAIND, list(new_NODE))
        
        additional_vectors <- get_all_vectors(new_NODE) %>% lapply(get_full_vector_voxels, IMG=IMG)
      } else if(remove_dend==T){ 
        
        ## remove current dendrite from node
        #    if(sum(length(), length()))
        
        
      } else {
        
        ## update subdendrite list in current node
        traced_SUBD <- list(info=SUBD$info, full_coords=bind_rows(SUBD$full_coords, bind_rows(coord_list)))
        new_SUBD_list <- append(new_SUBD_list, list(traced_SUBD))
        
        additional_vectors <- list()
      }
      
      
      
    } ## end of for... over subdendrite list
    
    MAIND[[nND]] <- list(node_id=NODE$node_id,
                         master_id=NODE$master_id,
                         subnodes=new_subnodes,
                         pos=NODE$pos,
                         subnode_full=new_SUBN_list,
                         subdend_full=new_SUBD_list)
    
    ## add new poits to full point list
    # snl <- lapply(new_SUBN_list, function(N){
    #   return(N[["full_coords"]])
    # }) %>%
    #   bind_rows()
    # 
    # sdl <- lapply(new_SUBD_list, function(N){
    #   return(N[["full_coords"]])
    # }) %>%
    #   bind_rows()
    # 
    # full <- bind_rows(snl, sdl)
    # 
    # if(nrow(full)!=0){
    #       vox_to_add <- apply(full, 1, 
    #       get_sphere_around_point, RAD=5, nr_orig=nr_orig, nz_orig=length(IMG),
    #       simplify=F)
    # 
    # all_points <- append(all_points, vox_to_add)
    # } 
    

    
    all_vectors_fv <- append(all_vectors_fv, additional_vectors)
    
    n_processed <- n_processed + 1
    
  } ## end of while ... return full new MAIND
  
  
  # test <- MASTER
  # test[[2]] <- MAIND
  # 
  # export_structure(test, inter$main_vectors, "f:/data_sholl_analysis/test/dendrites/traced_3.csv")
  
  return(MAIND)        
}



apply_3d_median_filter <- function(IMG, r){
  
  ni <- lapply(IMG, function(L){
    
    erg <- medianblur(as.cimg(L), 2*r) %>%
      as.matrix() 
    
    erg[erg<1] <- 0
    
    return(erg)
    
  })
  return(ni)
}





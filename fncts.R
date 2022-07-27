load_dependencies <- function(){
  
  library(tiff)
  library(tidyverse)
  library(gridExtra)
  library(cowplot) 
  library(Rcpp)
  library(imager)
  library(EBImage)
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
    cutoff_soma_region = 0.9,
    
    ## sccreening subdendrite starts
    subd_max_detection_distance = 35,
    subd_detection_depth = 8,
    subd_detection_vertical_range=6,
    subd_cluster_eps=2,
    subd_cluster_mpt=27,
    
    ## tracing subdendrites
    trace_cluster_eps=sqrt(5),
    trace_cluster_mpt=41,
    trace_detection_depth=7,
    trace_detection_distance = 15,
    trace_rescore_dist=10,
    trace_rescore_angle=40,
    trace_dend_cross_radius_vertical=9,
    trace_dend_cross_radius_horizontal=7
    
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
  test <- sapply(1:r, function(N){
    cppElongateLine(GRAD, zGRAD, x_start, y_start, z_start, N)
  }) %>% t() %>%
    cbind(
      i=apply(., 1, function(row){
        cppSelectIntensity(row[1], row[2], row[3], IMG=img)
      })) %>%
    cbind(deg=GRAD) %>%
           .[which(!is.na(.[,"i"])),] %>% as_tibble() %>% return()
}


elongate_line_full_old <- function(img, GRAD, zGRAD, x_start, y_start, z_start, r){
  #t0 <- Sys.time()
  lapply(1:r, function(N){
    cppElongateLine(GRAD, zGRAD, x_start, y_start, z_start, N)
  }) %>% bind_rows() %>%
    rowwise() %>%
    mutate(i=cppSelectIntensity(x,y,z,img),
           deg=GRAD) %>%
    .[which(!is.na(.$i)),] %>%
    return()
  #print(Sys.time()-t0)
}
# GRAD <- 93
# zGRAD <- 1
# img <- raw_image
# x_start <- 200
# y_start <- 100
# z_start <- 20
# r <- 100
elongate_line_full_mat <- function(img, GRAD, zGRAD, x_start, y_start, z_start, r){
  test <- sapply(1:r, function(N){
    cppElongateLine(GRAD, zGRAD, x_start, y_start, z_start, N)
  }) %>% t() %>%
    cbind(
  i=apply(., 1, function(row){
    cppSelectIntensity(row[1], row[2], row[3], IMG=img)
  })) %>%
    cbind(deg=GRAD) %>%
    .[which(!is.na(.[,"i"])),] %>%
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



#GRAD=most_likely_elongation 
#zGRAD=vertical_angle



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


# pos=next_pos
# raw_image=IMG
# 

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
  
  ## der part is neu:
  if(is.null(f)){return(NULL)}
  
  fin <- f[["coords"]]
  return(list(c(fin, vertical_angle), f[["full"]]))
  
}

calc_intensity_cutoff <- function(x,y,z,IMG, px_to_use, nr, nc){
  EBImage::otsu(IMG[[z]][c((y-px_to_use):(y+px_to_use)) %>% .[which(between(.,1,nr))],
                         c((x-px_to_use):(x+px_to_use)) %>% .[which(between(.,1,nc))]], 
                range=c(0,1), levels = 256)
}

# 
# DENDRITE_raw <- inter$main_dendrites[1,]
# SOMA=SOMA
# IMG=raw_image
# nr=opt$nr_orig
# nc=opt$nc_orig
elongate_dendrite <- function(DENDRITE_raw,
                              SOMA, IMG, nr, nc){
  
  
  horizontal_subdivisions <- 60
  vertical_subdivisions <- 10
  horizontal_detection_angle <- 60
  vertical_detection_angle <- 12
  steps <- 10
  
  
  #intensity_cutoff <- 0.5
  #i <- IMG[[SOMA$z]]
  
  px_to_use <- 100
  
  top_left_pix <- SOMA$x-px_to_use
  
  pix_surr_soma <- matrix()
  
  #px <- IMG[[SOMA$z]][c((SOMA$y-px_to_use):(SOMA$y+px_to_use)),c((SOMA$x-px_to_use):(SOMA$x+px_to_use))]
  intensity_cutoff <- calc_intensity_cutoff(SOMA$x, SOMA$y, SOMA$z, IMG, 200, nr, nc)
  # px %>% as_tibble() %>% rownames_to_column("y") %>% 
  #   gather(x,i,-y) %>% mutate(x=str_replace(x, "V", "") %>% as.numeric(),
  #                             y=as.numeric(y)) %>%
  # ggplot(aes(x,y, fill=i))+geom_tile()
  
  
  #intensity_cutoff <- EBImage::otsu(px, range=c(0,1), levels = 256)
  
  
  #cat(paste("\n  dendrite number:",DENDRITE_raw[["dendrite_id"]],"\n  --"))
  cat("\n  --")
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
  last_int <- SOMA
  c <- 2
  abort <- F
  #for(c in 2:7){
  while(abort==F&(c<3|sum(knot_list[[c-1]][1:3]==knot_list[[c]][1:3])!=3)){
    #cat(paste("\n  elongation step:", c-1))
    cat("-")
    #print(c)
    st <- Sys.time()
    
    if(cppDistPts(next_pos[["x"]], next_pos[["y"]], next_pos[["z"]],
                  last_int[["x"]],last_int[["y"]],last_int[["z"]])>(2*px_to_use)){
      
      intensity_cutoff <- calc_intensity_cutoff(next_pos[["x"]], next_pos[["y"]], next_pos[["z"]], 
                                                IMG, px_to_use, nr, nc)
      last_int <- next_pos
      #print(intensity_cutoff)
    }
    
    raw_dendrite <- screen_for_dendrite_elongation(next_pos, 
                                                   vertical_subdivisions, 
                                                   horizontal_subdivisions, 
                                                   horizontal_detection_angle,
                                                   vertical_detection_angle,
                                                   intensity_cutoff,
                                                   steps,
                                                   IMG)
    if(is.null(raw_dendrite)){
      abort <- T
    } else {
      next_pos <- raw_dendrite[[1]]
    #full_vector_list[[c]] <- raw_dendrite[[2]]
      c <- c+1 
      knot_list[[c]] <- next_pos %>% set_names(c("x","y","z", "h_angle", "v_angle"))
    }
    
    #print(Sys.time()-st)
    #cat(paste( "\n",round(next_pos[["x"]]), round(next_pos[["y"]]), round(next_pos[["z"]])))
  }
  cat(paste("|", c))
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

# IMG=raw_image
# cutoff_soma_region=opt$cutoff_soma_region


find_dendritic_start_sites <- function(SOMA, IMG, cutoff_soma_region){
  
  deg_step <- 0.005
  r <- 100
  #t0 <- Sys.time()
  mn = lapply(seq(deg_step,1,deg_step)*360, function(GRAD){
    #print(GRAD)
    # f <- elongate_3d_sphere(IMG, GRAD,0,SOMA[["x"]],SOMA[["y"]], SOMA[["z"]], r)
    # take_until <- f %>%
    #   #filter(i<cutoff_soma_region) %>%
    #   .[which(.$i<cutoff_soma_region), "n"] %>%
    #   pull() %>%
    #   min()  
      
    f <- elongate_line_full(IMG, GRAD,0,SOMA[["x"]],SOMA[["y"]], SOMA[["z"]], r)
    
    #take_until <- max(which(f$i>cutoff_soma_region))
    
    rel <- f[which(f$i>cutoff_soma_region),]
    
    if(nrow(rel)==0){
      return(NULL)
    } else {
      return(
       rel %>% ungroup() %>% mutate(n=c(1:nrow(.)))
      
      )
    }
  }) %>% compact() %>%
    bind_rows()
  #print(Sys.time()-t0)
  #cat("\nscreened for dendrites")
  
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
  
  
  #ggplot(tibble(x=dens_func$x, y=dens_func$y), aes(x,y))+geom_line()
  
  all_local_extreme <- get_minmax(dens_func)

  ## new one by using quantiles
  lowest_minimum_x <- all_local_extreme %>% filter(y==min(all_local_extreme$y)) %>% pull(x)
  if(lowest_minimum_x<180){
      relevant_maxima_raw <- all_local_extreme %>%
        filter(between(x, lowest_minimum_x, lowest_minimum_x+360)) %>%
        arrange(x)
      
      relevant_maxima <- lapply(which(relevant_maxima_raw$type=="maximum"), function(MX){
        if(MX==1){
          pre <- relevant_maxima_raw[[nrow(relevant_maxima_raw), "y"]]
          post <- relevant_maxima_raw[[MX+1,"y"]]
        } else if(MX==nrow(relevant_maxima_raw)){
          pre <- relevant_maxima_raw[[MX-1,"y"]]
          post <- relevant_maxima_raw[[1, "y"]]
        } else {
          pre <- relevant_maxima_raw[[MX-1,"y"]]
          post <- relevant_maxima_raw[[MX+1,"y"]]
        }
        
        rel <- relevant_maxima_raw[[MX,"y"]]
        
        if(max(diff(c(pre, rel, post)))>0.00005){
          return(relevant_maxima_raw[MX,] %>% mutate(mx=max(diff(c(pre, rel, post)))))
        }
        
      }) %>% bind_rows()%>%
         mutate(x=ifelse(x<0, x+360, x))
      
    } else {
      relevant_maxima_raw <- all_local_extreme %>%
        filter(between(x, lowest_minimum_x-360, lowest_minimum_x),) %>%
        arrange(x)
      
      relevant_maxima <- lapply(which(relevant_maxima_raw$type=="maximum"), function(MX){
        if(MX==1){
          pre <- relevant_maxima_raw[[nrow(relevant_maxima_raw), "y"]]
          post <- relevant_maxima_raw[[MX+1,"y"]]
        } else if(MX==nrow(relevant_maxima_raw)){
          pre <- relevant_maxima_raw[[MX-1,"y"]]
          post <- relevant_maxima_raw[[1, "y"]]
        } else {
          pre <- relevant_maxima_raw[[MX-1,"y"]]
          post <- relevant_maxima_raw[[MX+1,"y"]]
        }
        
        rel <- relevant_maxima_raw[[MX,"y"]]
        
        if(max(abs(diff(c(pre, rel, post))))>0.0001){
          return(relevant_maxima_raw[MX,] %>% mutate(mx=max(abs(diff(c(pre, rel, post))))))
        }
        
      }) %>% bind_rows()%>%
         mutate(x=ifelse(x<0, x+360, x)) 

    }
  
  rescored_maxima <- lapply(relevant_maxima$x, function(GRAD){
    
    f <- elongate_line(GRAD,0, SOMA[["x"]],SOMA[["y"]], SOMA[["z"]], soma_radius)
    return(f %>% c(h_angle=GRAD))
  }) %>% 
    bind_rows() %>%
    select(-z)
  
  return(list(rescored_maxima, soma_radius))
  
}
## req

# cube_size=opt$soma_xy_detection_cube_radius*2+1
# r=opt$soma_z_detection_radius
# deg_step=opt$soma_z_detection_degree_steps
# 


get_somata <- function(cube_size, r, deg_step ,raw_image){
  #t0 <- Sys.time()
  z_vec <- seq(((cube_size+1)/2), (length(raw_image)-(cube_size+1)/2), cube_size)
  y_vec <- seq(((cube_size+1)/2), (nrow(raw_image[[1]])-(cube_size+1)/2), cube_size)
  x_vec <- seq(((cube_size+1)/2), (ncol(raw_image[[1]])-(cube_size+1)/2), cube_size)
  ## loop over z axis
  n_its <- length(raw_image)*length(y_vec)*length(x_vec)
  n_vox_per_cube <- cube_size^3
  #t0 <- Sys.time()
  first_test <- lapply(z_vec, function(Z){
    lapply(y_vec, function(Y){
      lapply(x_vec, function(X){
        # sel_vec <- mapply(cppGetSingleIndex, 
        #                   c((X-(cube_size-1)/2):(X+(cube_size-1)/2)),
        #                   c((Y-(cube_size-1)/2):(Y+(cube_size-1)/2)), 
        #                   nr=nr_orig)
        zsum <- lapply(c((Z-(cube_size-1)/2):(Z+(cube_size-1)/2)), function(Z2){
          #raw_image[[Z2]][sel_vec] %>%
            
           raw_image[[Z2]][c((Y-(cube_size-1)/2):(Y+(cube_size-1)/2)),
                           c((X-(cube_size-1)/2):(X+(cube_size-1)/2))]%>% 
            
            sum() %>%
            return()
        }) %>% 
          as.numeric() %>% 
          sum()
        return(c(x=X, y=Y, z=Z, sum=zsum))
      }) %>% 
        #t() %>%
        bind_rows() %>% 
        return()
    }) %>% 
      #t() %>%
      bind_rows() %>% 
      return()
  }) %>% 
    #t()
    bind_rows()
  
  # t1 <- Sys.time()
  # print(t1-t0)
  filtered <- first_test %>%
    #.[which(.$sum>0.9*n_vox_per_cube),]
    filter(sum>0.9*n_vox_per_cube, z!=16)
  
  km <- kmeans(filtered%>% select(x,y,z), centers = 1)
  
  xy_soma <- tibble(wss=sqrt(km$withinss/nrow(filtered))) %>%
    bind_cols(km$centers)
  
  #t2 <- Sys.time()
  
  z_raw <- lapply(1:length(raw_image), function(Z){
    
      mn= sapply(seq(deg_step,1,deg_step)*360, function(GRAD){
        f <- elongate_line_full_mat(raw_image, GRAD, 0, xy_soma$x, xy_soma$y ,Z, 100) %>%
          cbind(dff=c(abs(diff(.[,"i"])),0))
        md <- max(f[,"dff"])
        rel <- which(f[,"dff"]==md) %>% min()
        md <- abs(diff(f[,"i"])) %>% max()
        return(c(dff=md, mx=rel))
      }) %>% t() %>%
        #bind_rows() %>%
        apply(.,2,mean)
        #unlist() %>% 
      #  median()
      return(c(z=Z,mn))
  }) %>% bind_rows() %>%
    mutate(score=mx/dff^0.5)
  
  if(nrow(z_raw %>%
          .[which(.$dff>0.1),])==0){
    
    xy_soma$z <- z_raw %>%
      .[which(.$dff==max(.$dff)),] %>%
      pull(z) %>%
      median()
  } else {
      xy_soma$z <- z_raw %>%
    .[which(.$dff>0.1),] %>%
    .[which(.$score==min(.$score)),] %>%
    pull(z) %>%
    median()
  }
  

  
  #t3 <- Sys.time()
  #print(t3-t2)
  # xy_soma$z <- z_raw %>%
  #   #filter(mn>quantile(z_raw$mn, 0.95)) %>%
  #   .[which(.$mn==max(.$mn)),] %>%
  #   pull(z) %>%
  #   median() %>%
  #   round()
  #t4 <- Sys.time()
  return(xy_soma)
}
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


bin_otsu2 <- function(sl, IMG){
  
  
  full_cutoffs <- lapply(1:length(sl), function(nVOX){
    #print(nVOX)
    VOX <- first(sl[[nVOX]])
  
    Z <- last(sl[[nVOX]])
    
    to_get_cutoff <- IMG[[Z]][VOX] %>% .[which(.!=0)]
    
    cutoff <- EBImage::otsu(matrix(to_get_cutoff), 
                range=c(0,1), levels = 256)
    
    #return(abs(1-cutoff))
     
    return(cutoff) 
  })
  
  new_image <- IMG
  #st1 <- Sys.time()
  for(i in 1:length(full_cutoffs)){
    #t0 <- Sys.time()
    cutoff <- full_cutoffs[[i]]
    VOX <- first(sl[[i]])
    for(l in 1:length(new_image)){
      new_image[[l]][VOX][new_image[[l]][VOX]>=cutoff] <- 1
      new_image[[l]][VOX][new_image[[l]][VOX]<cutoff] <- 0
      #new_image[[l]] <- matrix(new_image[[l]], nrow = 1040)
      #print(paste(i, l, cutoff))
      #if(class(new_image[[l]])=="numeric"){print(paste(i, l, cutoff))}
      
    }
    #cat(paste("\n  ", i, Sys.time()-t0))
  }
  # st2 <- Sys.time()
  # print(st2-st1)
  return(new_image)
  
}




#nosoma_image <- raw_image
#VOX <- sl[[1]]
#nosoma_image <- IMG

bi2 <- function(sl, nosoma_image){
  #t0 <- Sys.time()
  
  #VOX <- sl[[510]]
  
  full_cutoffs <- lapply(1:length(sl), function(nVOX){
    #print(nVOX)
    
    VOX <- sl[[nVOX]]
    
    all_vox_raw <- lapply(1:length(nosoma_image), function(L){
      
      LAYER <- nosoma_image[[L]]
      
      return(LAYER[VOX])
      
    }) %>%
      Reduce(function(x,y)c(x,y),.) 
    
    all_vox <- all_vox_raw[all_vox_raw!=0]
    
    if(length(all_vox)==0){return(1)}
    
    dens_func <- density(all_vox, bw=0.001) 
    
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
  #st1 <- Sys.time()
  for(i in 1:length(full_cutoffs)){
    #t0 <- Sys.time()
    cutoff <- full_cutoffs[[i]]
    VOX <- sl[[i]]
    for(l in 1:length(new_image)){
      new_image[[l]][VOX][new_image[[l]][VOX]>=cutoff] <- 1
      new_image[[l]][VOX][new_image[[l]][VOX]<cutoff] <- 0
      #new_image[[l]] <- matrix(new_image[[l]], nrow = 1040)
      #print(paste(i, l, cutoff))
      #if(class(new_image[[l]])=="numeric"){print(paste(i, l, cutoff))}
      
    }
    #cat(paste("\n  ", i, Sys.time()-t0))
  }
  # st2 <- Sys.time()
  # print(st2-st1)
  return(new_image)
}
#writeTIFF(new_image, "f:/data_sholl_analysis/test/runs/testbin.tif")
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

export_dendrites <- function(elongated_dendrites, file_name, SOMA){
  
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
# 
# 
# cof=top
# bd=top_border
# sdir="top"
# 



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
  
  ## hier ist das problem dass wenn ein segment aus dem bild herausragt...
  ## also x oder y jeweils unter 1 gehen... die jeweisl andere coordinate,
  ## die zur verbindung genutzt wird abschenidet... es muss also auf 1 oder max border
  ## geändert werden
  
  ################################################################################
  ## old
  # top_line <- tibble(x=c(xs:xe),
  #                    y=cof) %>%
  #   filter(!x %in% rob$x) %>%
  #   bind_rows(rob)%>%
  #   filter(between(x, 1, nc),
  #          between(y, 1, nr))
  ################################################################################
  top_line <- tibble(x=c(xs:xe),
                     y=cof) %>%
                       filter(!x %in% rob$x) %>%
                       bind_rows(rob) %>%
    mutate(x=ifelse(x<1, 
                    1, 
                    ifelse(x>nc, 
                           nc, 
                           x)),
           y=ifelse(y<1,
                    1,
                    ifelse(y>nr,
                           nr,
                           y))
           )
  #print(top_line)
  
  all_vox_raw <- lapply(top_line$x, function(X){
    #print(X)
    s <- top_line[which(top_line$x==X), "y"] %>% pull %>% unique()
    e <- line[which(line$x==X), "y"] %>% pull %>% unique()
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
# 
# nELD <- 3
# avs=raster_size
# borders_vox_raw=borders_vox_raw
# selection_vector=selection_vector
# nr_orig=nr_orig
# nc_orig=nc_orig
# process_full=process_full


create_raster <- function(nELD, avs, borders_vox_raw, selection_vector, elongated_dendrites,nr_orig, nc_orig, process_full){
   # cat(paste("dendrite:", nELD))
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
  
  
  #ggplot(bind_rows(bottom_border), aes(x,-y))+geom_tile()
  
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
  
  ### here we need to check if backward segments need to be added... according to
  ### bordes which excedd the current quadrant/half
  
  
  
  use_length <- round(pixels_to_image_border/n_segments)
  
  res_list <- assign_vectors_to_segments(ELD, vector_pos, use_length, n_segments)
  
  full_vecs <- bind_rows(vectors) %>%
    pull(x) %>% cumsum()  
  
  y_coord_list <- list()
  segment_list <- list()
  segment_count <- 0
  #for(n in c(1:28)){
  for(n in c(1:n_segments)){
    
    #print(n)
    
    # cat(paste("\n  x_segment:",n))
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
    
    ## old
    # med_line <-   tibble(x=c(xs:xe)) %>%
    #   mutate(fac=c(rl, recursive=T)) %>%
    #   mutate(y=ceiling(ys+cumsum(fac))) %>%
    #   filter(between(x, 1, nc),
    #          between(y, 1, nr))
    
    med_line <-   tibble(x=c(xs:xe)) %>%
      mutate(fac=c(rl, recursive=T)) %>%
      mutate(y=ceiling(ys+cumsum(fac))) %>%
      mutate(x=ifelse(x<1, 
                      1, 
                      ifelse(x>nc, 
                             nc, 
                             x)),
             y=ifelse(y<1,
                      1,
                      ifelse(y>nr,
                             nr,
                             y))
      )
    
    
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
    
    if(n_segments_bottom<1){
      n_segments_bottom <- 1
    }
    if(n_segments_top<1){
      n_segments_top <- 1
    }
    ## new loop for all segments
    
    start_list_top <- list()
    line_list_top <- list()
    
    for(n_top in 1:n_segments_top){
      #print(n_top)
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
      #cat(paste("\n    top:", segment_count))
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
      #cat(paste("\n    bottom:", segment_count))
      segment_list[[segment_count]] <- all_vox_bottom[[1]]
    }       
    #print(4)
    y_coord_list[[n]] <- med_line[[nrow(med_line), "y"]]
    
    
  }
  return(segment_list)  
}


#  
# soma_radius=inter$soma_radius
# nr_orig=opt$nr_orig
# nc_orig=opt$nc_orig



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
    .[which(.$c>10),]
  #filter(c>10) 
  
  if(nrow(max_n)==0){
    vec %>%pull(i) %>%
      mean() %>%
      return()
  } else {
    return(0)
  } 
  
  
}

define_full_vector_info <- function(ELD, IMG, SOMA){
  
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
   #print(2)
  elongate_line_full_old(IMG, ha, va, xs, ys, zs, round(dist)) %>% 
  
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
# n <- 2
# main_vectors=inter$main_vectors
# main_vectors_df=inter$main_vectors_df
# main_vectors_full=inter$main_vectors_full
# SOMA=SOMA
# IMG_screen=med_bin_noS_noMD_image
# IMG_adj=bin_noS_image_p
# soma_radius=inter$soma_radius
# MPTS=opt$subd_cluster_mpt
# EPS=opt$subd_cluster_eps
# det_rad=opt$subd_max_detection_distance
# z_range=opt$subd_detection_vertical_range
# nr_orig=opt$nr_orig
# screen_subd_starts_circle=inter$subd_max_detection_circle
find_subdendritic_starts_new <- function(n, main_vectors, main_vectors_df,main_vectors_full, 
                                         SOMA, IMG_screen, IMG_adj, soma_radius, MPTS, EPS, det_rad, z_range,
                                         depth, nr_orig, screen_subd_starts_circle){
  
  #find_subdendritic_starts <- function(IMG, main_vectors, main_vectors_full, rem_rad, nr_orig, nc_orig)  
  
  #t0 <- Sys.time()
  cat(paste0("\n  process: MD:", n))
  
  rel_vecs_raw <- main_vectors[[n]]
  full_vecs_raw <- main_vectors_full[[n]]
  
  rel_vecs <- rel_vecs_raw[c(2:length(rel_vecs_raw))]
  full_vecs <- full_vecs_raw[c(2:length(full_vecs_raw))]
  
  full_dendrite <- bind_rows(full_vecs) %>% 
    mutate_all(.funs = as.numeric)
  
  comb_main_dend <- full_vecs_raw %>% bind_rows() %>%
    mutate_all(.funs = as.numeric) %>%
    mutate(ind=get_single_index(x,y,nr_orig))
  #t1 <- Sys.time()
  all_surrounding_vox <- get_surr_layer(comb_main_dend, SOMA, nr_orig, 
                                        IMG_screen, soma_radius, z_range, det_rad)
  #print(nrow(all_surrounding_vox))
  if(nrow(all_surrounding_vox)==0){
    return(NULL)
  }
  #print(t2-t1)
  #ggplot(all_surrounding_vox, aes(x=x, y=y))+facet_wrap(~z)+geom_tile()
  #show_surrounding_layer_all_layers(all_sorrounding_vox, "f:/data_sholl_analysis/test/full_sorrounding.pdf")
  ## orig... withour neighboring voxels... eps=3.5 minPTS=15
  #cat(paste0("\n    clustering"))
  df_centers <- find_subd_cluster_man(all_surrounding_vox, EPS, MPTS)
  #t2 <- Sys.time()
  #print(2)
  ## if no subdendrite starts are found... the main dendrite is just a subdendrite
  if(is.null(df_centers)){
    #print(2)
    #full_branches
    return(NULL)
  }
  retraced_subd <- apply(bind_rows(df_centers), 1, fit_path_to_subd_cluster, 
                         full_dendrite=full_dendrite,
                         det_rad=det_rad,
                         IMG=IMG_adj) %>% 
    compact() 
  #print(3)
  ## remove duplicated subdendirte starts
  #t3 <- Sys.time()
  if(is.null(retraced_subd)){#print(3)
    return(NULL)}
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
  
  #print(4)
  #cat(paste0("\n    adjusting"))
  
  ## check intersecting subdendrites
  #t4 <- Sys.time()
  if(is.null(rescored_subdendrites_raw)){#print(4)
    return(NULL)}
  
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
    
    #print(4)
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
            #print(5)
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
  
  #rescored_subdendrites <- retraced_subd
  #t5 <- Sys.time()
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
  
  
  #idNODE <- 3
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
      
      dist_to_main <- cppDistPts(last_vec_info$x, last_vec_info$y, last_vec_info$z,
                                 node$xs, node$ys, node$zs)
      
      
      if(dist_to_main<4){
        ## remove dendrite
      } else {
        vecs_to_add <- node$level:length(rel_vecs_raw)
        
        last_vec_coords <- lapply(vecs_to_add, function(V){
          VEC <- rel_vecs_df %>% filter(level==V) %>% select(x=xe, y=ye, z=ze) %>% return()
        }) %>%
          bind_rows() %>%
          bind_rows(select(node, x=xs,y=ys,z=zs),.)
        
        ## remove those which are almost the same
        full_subdend_info <- append(full_subdend_info, 
                                    list(list(info=last_vec_info, 
                                              full_coords=last_vec_coords)))        
      }
      

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
  
  nodes_to_merge <- df_sorted %>%
    mutate(id=paste(xs,ys,zs, sep="_")) %>%
    group_by(id) %>%
    summarize(to_merge=paste(node, collapse = ", "))
  
  
  final_nodes <- lapply(nodes_to_merge$to_merge, function(NODES){
    #print(1)
    if(str_count(NODES, ",")==0){
      #print(2)
      return(assigned_nodes[[as.numeric(NODES)]])
    } else {
      
      nodes_to_merge <- as.numeric(unlist(str_split(NODES, ", ")))
      
      new_node_id <- min(nodes_to_merge)
      master_id <- assigned_nodes[[new_node_id]]$master_id
      
      pos <- assigned_nodes[[new_node_id]]$pos
      
      subnodes_raw <- lapply(assigned_nodes[nodes_to_merge], function(LO){
        LO[["subnodes"]]
      }) %>% unlist() %>% .[which(!. %in% nodes_to_merge)]
      
      if(length(subnodes_raw)==0){
        subnodes <- list()
      }
      
      subnode_full <- assigned_nodes[[max(nodes_to_merge)]]$subnode_full
      
      subdend_full <- lapply(assigned_nodes[nodes_to_merge], function(LO){
        LO[["subdend_full"]]
      }) %>% Reduce(function(x,y)append(x,y),.)
      
      list(node_id=new_node_id,
           master_id=master_id,
           subnodes=subnodes,
           pos=pos,
           subnode_full=subnode_full,
           subdend_full=subdend_full) %>%
        return()
      
    }
    
  })
  
  return(final_nodes)
  #t6 <- Sys.time()
}
#df <-centers_rem
#df <- all_surrounding_vox
# EPS <- 2
# MPTS <- 33
find_subd_cluster_man <- function(df, EPS, MPTS){
  df_filtered <- df %>% 
    #   #.[which(.$i!=0),] %>%
    select(x,y,z) #
  #ggplot(df_filtered, aes(x=x, y=-y,))+geom_point()+facet_wrap(~z)
  #if(nrow(df_filtered)==0){return(NULL)} 
  df_clustered <- df_filtered %>% 
    #select(x,y,z) %>%
    ungroup() %>%
    mutate(c=fpc::dbscan(df_filtered, eps = EPS, MinPts = MPTS)$cluster) %>%
    .[which(.$c!=0),]
  if(nrow(df_clustered)==0){return(NULL)}
  ggplot(df_clustered, aes(x=x, y=-y,color=as.character(c)))+geom_point()+facet_wrap(~z)
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

create_circle_mat <- function(DET_RAD){
  if(is.integer(DET_RAD/2)){
    ss <- 2*pi/(DET_RAD*12)
  } else {
    ss <- 2*pi/(DET_RAD*2*(5+(DET_RAD-2)))
  }
  sapply(seq(ss, 2*pi, ss), function(a){
    return(c(x=cos(a)*DET_RAD,
             y=sin(a)*DET_RAD) %>% round())
  }) %>% t() %>%
    return()
  
}

create_circle <- function(DET_RAD){
  if(is.integer(DET_RAD/2)){
    ss <- 2*pi/(DET_RAD*12)
  } else {
    ss <- 2*pi/(DET_RAD*2*(5+(DET_RAD-2)))
  }
  
  lapply(seq(ss, 2*pi, ss), function(a){
    return(c(x=cos(a)*DET_RAD,
             y=sin(a)*DET_RAD) %>% round())
  }) %>% unique() %>% bind_rows() %>%
    return()
}
# X=xs
# Y=ys
# Z=zs
screen_circular_new <- function(det_rad, rel_z_layer,bmin, bmax,
                                INC, X, Y, Z, ha, ha_cutoff, screening_angle, IMG,
                                circle, SOMA){
  
  DET_RAD <- det_rad+INC
  somx <- round(SOMA[["x"]])
  somy <- round(SOMA[["y"]])
  somz <- round(SOMA[["z"]])
  fin <- lapply(rel_z_layer, function(ZS){
    #sst0 <- Sys.time()
    # df_circ_x <- tibble(x=seq(-DET_RAD, DET_RAD, 1)) %>%
    #   mutate(yp=sqrt(DET_RAD^2-x^2),
    #          yn=-sqrt(DET_RAD^2-x^2)) %>%
    #   gather(dir, y, -x)
    # 
    # df_circ_y <- tibble(y=seq(-DET_RAD, DET_RAD, 1)) %>%
    #   mutate(xp=sqrt(DET_RAD^2-y^2),
    #          xn=-sqrt(DET_RAD^2-y^2)) %>%
    #   gather(dir, x, -y)
    # 
    # df_circ_full <- bind_rows(df_circ_x, df_circ_y) %>%
    #   mutate_at(.vars = c("x", "y"), .funs = round) %>%
    #   mutate(id=paste(x, y, sep="_")) %>%
    #   filter(!duplicated(id)) %>%
    #   rowwise() %>%
    #   mutate(xa=x+X,
    #          ya=y+Y,
    #          angle=cpprad2deg(atan((y/DET_RAD)/(x/DET_RAD))),
    #          af=case_when(
    #            x<0 ~ angle+180,
    #            angle<0~ angle+360,
    #            TRUE ~ angle
    #          ),
    #   )
    #sst0 <- Sys.time()
    df_circ_full <- circle %>%
      rowwise() %>%
      mutate(xa=x+X,
             ya=y+Y,
             angle=cpprad2deg(atan((y/DET_RAD)/(x/DET_RAD))),
             # af=case_when(
             #   x<0 ~ angle+180,
             #   angle<0~ angle+360,
             #   TRUE ~ angle
             # ),
             af=ifelse(x<0,
                       angle+180,
                       ifelse(angle<0, 
                              angle+360, 
                              angle))
      )
    #sst1 <- Sys.time()
    #print(sst1-sst0)
    all_vox <- lapply(c(min(df_circ_full$x):max(df_circ_full$x)), function(X){
      
      all_x <- df_circ_full[which(df_circ_full$x==X),] #%>%
      #arrange(y)
      #return(tibble(x=X, y=c(min(all_x$y):max(all_x$y))))
      y=c(min(all_x$y):max(all_x$y))
      x=rep(X, length(y))
      return(cbind(x,y) %>% set_names(nm=c("x", "y")))
    }) %>% 
      bind_rows() %>%
      #Reduce(function(x,y)rbind(x,y),.) #%>%
      # apply(.,1,function(ROW){
      #   dist=cppDistPts(ROW["x"],ROW["y"],0, 0,0,0)
      #   ha=cppGetHA(ROW["x"],ROW["y"],0,0)
      #   return(c(ROW, dist=dist, ha=ha))
      #  }) #%>% 
      # t() %>% 
      # as_tibble() %>%
      # #bind_rows() %>% 
      rowwise() %>%
      mutate(dist=cppDistPts(x,y,0, 0,0,0),
             ha=cppGetHA(x,y,0,0)) %>%
      .[which(.$dist>det_rad-INC),]
    #sst2 <- Sys.time()
    #print(sst2-sst1)
    #ggplot(all_vox, aes(x=x, y=y))+geom_tile()
    if(between(bmin, 0, 360)&between(bmax, 0, 360)){
      df_circ <- all_vox[which(between(all_vox$ha, bmin, bmax)),]
    } else if(bmin<0){
      df_circ <- all_vox[which(between(all_vox$ha, bmin+360, 360)|between(all_vox$ha, 0, bmax)),]
    } else {
      df_circ <- all_vox[which(between(all_vox$ha, bmin, 360)|between(all_vox$ha, 0, bmax-360)),]
    }
    
    #ggplot(df_circ, aes(x=x, y=y))+geom_tile()
    ret <- df_circ %>%    
      rowwise() %>%
      mutate(z=ZS,
             x=x+X,
             y=y+Y,
             i=cppSelectIntensity(x, y, z, IMG=IMG)) %>%
      
      #select(x=xa, y=ya, z, i) %>%
      mutate(dist_to_soma=cppDistPts(x, y, z, somx, somy, somz)) #%>%
    #sst3 <- Sys.time()
    #print(sst3-sst2)
    return(ret[c("x","y","z","i","dist_to_soma")])
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
# SUBD <- bind_rows(df_centers)[8,]
# full_dendrite=full_dendrite
# det_rad=det_rad
# IMG=IMG_adj
fit_path_to_subd_cluster <- function(SUBD, full_dendrite, det_rad, IMG){
  #print(1)
  # t0 <- Sys.time()
  df_dist_raw <- full_dendrite %>%
    rowwise() %>%
    mutate(dist=cppDistPts(x,y,z,as.numeric(SUBD[["x"]]),as.numeric(SUBD[["y"]]),as.numeric(SUBD[["z"]])),
           ha=cppGetHA(x,y,as.numeric(SUBD[["x"]]),as.numeric(SUBD[["y"]])),
           va=cppGetVA(dist, as.numeric(SUBD[["z"]]), z)) %>%
    .[which(.$dist<=4*det_rad),]
  #print(2)
  # t1 <- Sys.time()
  ############################################################### make this faster (2.1 secs)
  if(nrow(df_dist_raw)==0){return(NULL)}
  
  df_dist <- df_dist_raw %>%
    mutate(int_score=fit_subdendrite_new(SUBD, dist, ha, va, IMG),
           score=dist*(1/int_score^2)) %>%
    .[which(.$int_score>0.7),]
  ###############################################################
  #  t2 <- Sys.time()
  if(nrow(df_dist)==0){return(NULL)}
  
  
  final_intersection <- df_dist[which(df_dist$int_score==max(df_dist$int_score)),] %>%
    .[which(.$score==min(.$score)),] %>%
    
    #.[which(.$dist==min(.$dist)),] %>%
    .[1,]
  #print(4)
  #  t3 <- Sys.time()
  closest_main_vector <- df_dist_raw[which(df_dist_raw$dist==min(df_dist_raw$dist)),] %>%
    ## sometimes two main vectors have exact same distance....
    .[1,]
  #print(5)
  #  t4 <- Sys.time()
  subd_start <- final_intersection %>% select(xs=x,ys=y,zs=z, level, ha, va, dist) %>%
    mutate(xe=SUBD[["x"]],
           ye=SUBD[["y"]],
           ze=SUBD[["z"]],
           clos_vec=closest_main_vector$level,
           d_to_cv=closest_main_vector$dist)
  #print(6)
  #  t5 <- Sys.time()
  return(subd_start)
  
}
#MAIND_raw <- safe_MAIND
#tMASTER <- MASTER
#main_vectors <- inter$main_vectors
export_structure <- function(tMASTER, main_vectors, SOMA, output_file){
  
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

format_new_subdendrite_starts <- function(CENT_raw, full_vecs_raw, last_pos, new_node_id){
  
  ## insert new ha cutoff here
  CENT <- CENT_raw %>% as_tibble()
  ## calc distance to closest main vector
  ## i just do it by taking all pts and calc all dists
  #full_vecs_raw <- main_vectors_full[[nMD]]
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

get_vox_of_borders <- function(n, main_dendrites, nc_orig, nr_orig, IMG, SOMA){
  
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

get_circle_around_point <- function(P, df_circ_full, nr_orig){
  
  #df_circ_full <- create_circle(RAD)
  
  #print(min(df_circ_full$x))
  #print(max(df_circ_full$x))
  
  all_vox <- lapply(c(min(df_circ_full$x):max(df_circ_full$x)), function(X){
    
    all_x <- df_circ_full[which(df_circ_full$x==X),] #%>%
    #arrange(y)
    #return(tibble(x=X, y=c(min(all_x$y):max(all_x$y))))
    y=c(min(all_x$y):max(all_x$y))
    x=rep(X, length(y))
    return(cbind(x,y) %>% set_names(nm=c("x", "y")))
  }) %>% 
    bind_rows() %>%
    mutate(x=x+P[["x"]],
           y=y+P[["y"]],
           short=get_single_index(x,y,nr_orig)) %>%
    pull(short) %>%
    return()
}

get_circle_around_point_xy <- function(P, df_circ_full, nr_orig){
  
  #df_circ_full <- create_circle(RAD)
  all_vox <- lapply(c(min(df_circ_full$x):max(df_circ_full$x)), function(X){
    
    all_x <- df_circ_full[which(df_circ_full$x==X),] #%>%
    #arrange(y)
    #return(tibble(x=X, y=c(min(all_x$y):max(all_x$y))))
    y=c(min(all_x$y):max(all_x$y))+P[["y"]]
    x=rep(X, length(y))+P[["x"]]
    return(cbind(x,y) %>% set_names(nm=c("x", "y")))
  }) %>% 
    bind_rows() %>%
    # mutate(x=x+P[["x"]],
    #        y=y+P[["y"]],
    #        #short=get_single_index(x,y,nr_orig)
    #        ) %>%
    #pull(short) %>%
    return()
}

get_sphere_around_point <- function(P, RADv, circle, nr_orig, nz_orig){
  
  lapply(c((P[["z"]]-RADv):(P[["z"]]+RADv)) %>% .[which(between(., 1, nz_orig))], function(Z){
    get_circle_around_point(P, circle, nr_orig) %>% paste(Z, ., sep="_")
  }) %>% c(recursive=T)
  
}

get_all_vector_speheres <- function(all_vectors_fv_comb, RADv, RADh, circle, nr_orig, IMG){
  all_points <- list()
  
  nra <- nrow(all_vectors_fv_comb)
  
  
  
  steps <- seq(1, nra, round(RADh/2)) %>% c(.,nra) %>% unique()
  
#  cat(paste("       assigning", nra, "pts (", length(steps), "iterations)"))
  for(nPT in steps){
    
    #print(nPT)
    
    PT <- all_vectors_fv_comb[nPT,]
    apt <- get_sphere_around_point(PT, RADv, circle, nr_orig, length(IMG))
    
    all_points <- c(all_points, apt) %>% unique()
    
  }
  return(all_points)
} 
# n=2
# IMG=images$bin_noS_noMD_image
# MAIND=MASTER[[n]]
# MV=inter$main_vectors_df[[n]]
# main_vectors_full=inter$main_vectors_full[[n]]
# 
# all_vectors = inter$all_vectors[[n]]
# all_vectors_fv =inter$all_vectors_fv[[n]]
# all_points =inter$all_points[[n]]
# EPS_orig=opt$trace_cluster_eps
# MPTS_orig=opt$trace_cluster_mpt
# INC=opt$trace_detection_depth
# 
# RADv=opt$trace_dend_cross_radius_vertical
# RADh=opt$trace_dend_cross_radius_horizontal
# 
# det_rad=opt$trace_detection_distance
# RESC_DIST=opt$trace_rescore_dist
# RESC_ANG=opt$trace_rescore_angle
# nr_orig=opt$nr_orig
trace_subdendrites <- function(MAIND, MV, main_vectors_full,
                               all_vectors, all_points,
                               IMG, EPS_orig, MPTS_orig, INC, 
                               RADv, RADh, cross_circle,
                               det_rad, RESC_DIST, RESC_ANG, nr_orig,
                               circle, SOMA){  
  
  red_vec <- lapply(all_vectors, function(V) {
    return(list(c(x=V[["xs"]],y=V[["ys"]],z=V[["zs"]], id=V[["id"]], ha=V[["ha"]]),
                c(x=V[["xe"]],y=V[["ye"]],z=V[["ze"]], id=V[["id"]], ha=V[["ha"]])
    ))
  }) %>% Reduce(function(x,y)append(x,y),.) %>% unique() %>% bind_rows()
  
  if(is.null(MAIND)){
    return(NULL)
  }
  nND <- 0
  n_processed <- 0
  #nND <- 25
  while(n_processed < length(MAIND)){ ## nodes
    
    nND <- nND + 1
    
    #cat(paste("\n  node:", nND, "of", length(MAIND)))
    
    NODE <- MAIND[[nND]]
    
    SUBD_list <- NODE$subdend_full
    
    new_subnodes <- NODE$subnodes
    new_SUBD_list <- list()
    new_SUBN_list <- NODE$subnode_full
    
    cat(paste0("\n|--", nND, "/", length(MAIND), ": ", length(SUBD_list)))
    #nSD <- 1
    for(nSD in 1:length(SUBD_list)){
      
      #cat(paste("\n    subd:", nSD, "of", length(SUBD_list), "|"))
      cat("\n|  ")
      
      SUBD <- SUBD_list[[nSD]]
      SUBD_info <- SUBD$info
      
      xs <- SUBD_info$x
      ys <- SUBD_info$y
      zs <- SUBD_info$z
      ha <- adj_deg(SUBD_info$ha)
      
      coord_list <- list()
      ha_list <- list()
      c <- 1
      abort <- F
      node_detected <- F
      remove_dend <- F
      z_cutoff <- 12
      #for(i in 1){
      while(abort==F){
        #t0 <- Sys.time()
        #cat(paste("\n      elgt step:",c))
        cat("-")
        if(SUBD_info$orient=="end"){
          ha_cutoff <- NULL
          screening_angle <- 180
          z_range <- 6
          det_rad <- 20
          alt <- T
        } else if(SUBD_info$orient=="side"&c==1){
          screening_angle <- 220
          z_range <- 8
          ha_cutoff <- NULL
          alt <- T
        } else if(SUBD_info$orient=="elongated"&c==1){
          screening_angle <- 130
          z_range <- 6
          ha_cutoff <- NULL
        } else if(c<3){
          ha_cutoff <- NULL
          screening_angle <- 120
          z_range <- 6
          alt <- F
        } else {
          ha_cutoff <- NULL
          screening_angle <- 120
          z_range <- 6
          alt <- F
        }
        
        Z <- round(zs+tan(cppdeg2rad(SUBD_info$va))*det_rad)
        
        ###########################################
        
        cluster_ok=F
        DR <- det_rad
        ZR <- z_range
        MPTS <- MPTS_orig
        EPS <- EPS_orig
        #t1 <- Sys.time()
        while(cluster_ok==F){
          st1 <- Sys.time()
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
          #st2 <- Sys.time()
          
          
          
          centers1 <- screen_circular_new(DR, rel_z_layer, bmin, bmax,
                                          INC, xs, ys, zs, ha, ha_cutoff, screening_angle, IMG, circle, SOMA) #%>%
          
          # st3 <- Sys.time()
          # print(st3-st2)
          #ggplot(centers1, aes(x=x, y=-y))+geom_point()+facet_wrap(~z)
          
          if(nrow(centers1)==0){
            centers2 <- NULL
            cluster_ok <- T
          } else {
            
            
            ## remove positions of other dendrites
            # centers_rem <- centers1 %>%
            #   mutate(index=get_3d_single_index(x,y,z, nr_orig)) %>%
            #       .[which(!.$index %in% c(all_points, recursive=T)),]
            centers_rem <- centers1 %>%
              mutate(index=get_3d_single_index(x,y,z, nr_orig)) %>%
              .[which(!.$index %in% all_points),]
            #st4 <- Sys.time()
            #ggplot(centers_rem, aes(x=x, y=-y))+geom_point()+facet_wrap(~z)
            
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
                  #st5 <- Sys.time()
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
                #cat("  - cluster not valid")
              }
            }
          }
          #st6 <- Sys.time()
        }  
        #t2 <- Sys.time()    
        if(length(centers2)==0){
          ## no elongation detected
          abort=T
          #cat("\n      stopped")
          
          if(c==1){
            remove_dend <- T
            cat(paste("X"))
          } else {
            cat(paste("|"))
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
            #cat("\n      removed duplicated")
          } else {
            cat("<")
            node_detected <- T
            abort <- T
          }
          
          only_one <- ifelse(length(centers)==1, T, F)
        } else {
          only_one <- T
          centers <- centers2
          ## elongate... only one center
        }          
        #t3 <- Sys.time()
        
        if(only_one==T){
          ## subdendrite elongates
          
          elgt <- centers[[1]]
          dist <- cppDistPts(elgt[["x"]], elgt[["y"]], elgt[["z"]], xs, ys,zs)
          ha <- cppGetHA(elgt[["x"]], elgt[["y"]], xs,ys)
          va <- cppGetVA(dist, 
                         elgt[["z"]], zs)
          nha <- ha
          # overlap <- red_vec %>% rowwise() %>%
          #   #bind_rows() %>%
          #   mutate(dist=cppDistPts(x,y,z,elgt[["x"]],elgt[["y"]],elgt[["z"]])) %>%
          #   .[which(.$dist<10),]
          overlap <- red_vec %>% rowwise() %>%
            #bind_rows() %>%
            mutate(dist=cppDistPts(x,y,z,elgt[["x"]],elgt[["y"]],elgt[["z"]])) %>%
            #.[which(.$dist<100),] #%>%
            mutate(xyd=cppDistPts(x,y,0,elgt[["x"]],elgt[["y"]],0),
                   zd=abs(z-elgt[["z"]]),
                   score=sqrt(xyd^2+(0.125*zd^2))) %>%
            .[which(.$score<10),]
          
          #overlap <- lapply(all_vectors, function()
          
          if(nrow(overlap)>0){
            ovfin <- overlap %>% #rename(vha=ha) %>% rowwise() %>%
              mutate(dha=abs(ha-nha),
                     dha2=ifelse(dha>180, 360-dha, dha),
                     diff_ha=ifelse(dha2>90, abs(180-dha2), dha2)) %>%
              .[which(.$diff_ha<45),]
          } else {
            ovfin <- overlap
          }
          
          if(nrow(ovfin)>0&c>1){
            abort <- T
            #cat("\n      crossed")
          } else {
            xs <- elgt[["x"]]
            ys <- elgt[["y"]]
            zs <- elgt[["z"]]    
            coord_list[[c]] <- elgt
            ha_list[[c]] <- ha
            c <- c+1
          } 
        } 
        #t4 <- Sys.time()
      } ## end of elongations (main while loop)
      
      #cat(paste("|", c))
      
      
      ################################################################################     
      if(node_detected==T){
        
        new_node_id <- length(MAIND)+1
        
        if(c==1){
          last_pos <- SUBD$full_coords[nrow(SUBD$full_coords),]
        } else {
          last_pos=as_tibble(last(coord_list)) #%>% select(-ha)
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
        
        proc_centers <- lapply(centers, format_new_subdendrite_starts, 
                               full_vecs_raw=main_vectors_full, 
                               last_pos=last_pos,
                               new_node_id=new_node_id)
        
        
        new_NODE <- list(node_id=new_node_id,
                         master_id=NODE$node_id,
                         subnodes=list(),
                         pos=last_pos,
                         subnode_full=list(),
                         subdend_full=proc_centers)
        
        MAIND <- append(MAIND, list(new_NODE))
        
        #additional_vectors <- get_all_vectors(new_NODE) %>% lapply(get_full_vector_voxels, IMG=IMG)
        
        #new_vectors <- get_all_vectors(new_NODE)
        
        
        if(length(coord_list)>1){
          new_vectors <- lapply(2:length(coord_list), function(V){
            
            tibble(xe=coord_list[[V]][["x"]],
                   ye=coord_list[[V]][["y"]],
                   ze=coord_list[[V]][["z"]],
                   xs=coord_list[[V-1]][["x"]],
                   ys=coord_list[[V-1]][["y"]],
                   zs=coord_list[[V-1]][["z"]])
            
          }) %>% append(get_all_vectors(new_NODE))
        } else {
          new_vectors <- get_all_vectors(new_NODE)
        }
        
        
        new_red_vec <- new_vectors %>%
          lapply(function(V) {
            ha <- cppGetHA(V[["xe"]], V[["ye"]],V[["xs"]], V[["ys"]])
            return(list(c(x=V[["xs"]],y=V[["ys"]],z=V[["zs"]], id=NA, ha=ha),
                        c(x=V[["xe"]],y=V[["ye"]],z=V[["ze"]], id=NA, ha=ha)
            ))
          }) %>% Reduce(function(x,y)append(x,y),.) %>% unique() %>% bind_rows()
        
        #new_all_points <- assign_vector_voxels(new_vectors)
        
        ##make them right
        #circ2 <- create_circle(RADh)
        
        new_all_points <- assign_vector_voxels(new_vectors,
                                             IMG, cross_circle, RADv, RADh, nr_orig)
        
      } else if(remove_dend==T){ 
        
        ## remove current dendrite from node
        #    if(sum(length(), length()))
        
        
      } else {
        
        ## update subdendrite list in current node
        traced_SUBD <- list(info=SUBD$info, full_coords=bind_rows(SUBD$full_coords, bind_rows(coord_list)))
        new_SUBD_list <- append(new_SUBD_list, list(traced_SUBD))
        
        new_red_vec <- bind_rows(coord_list) %>% mutate(ha=unlist(ha_list))
        
        new_all_points <-  lapply(2:nrow(traced_SUBD$full_coords), function(V){
          tibble(xe=traced_SUBD$full_coords[[V,"x"]],
                 ye=traced_SUBD$full_coords[[V,"y"]],
                 ze=traced_SUBD$full_coords[[V,"z"]],
                 xs=traced_SUBD$full_coords[[V-1,"x"]],
                 ys=traced_SUBD$full_coords[[V-1,"y"]],
                 zs=traced_SUBD$full_coords[[V-1,"z"]],)
        }) %>% assign_vector_voxels(.,
                                    IMG, cross_circle, RADv, RADh, nr_orig)
      }
      
      red_vec <- bind_rows(red_vec, new_red_vec) 
      all_points <- append(all_points, new_all_points) %>% unique()
    } ## end of for... over subdendrite list
    
    MAIND[[nND]] <- list(node_id=NODE$node_id,
                         master_id=NODE$master_id,
                         subnodes=new_subnodes,
                         pos=NODE$pos,
                         subnode_full=new_SUBN_list,
                         subdend_full=new_SUBD_list)
    
    n_processed <- n_processed + 1
    
  } ## end of while ... return full new MAIND
  
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

bin_otsu <- function(sl, IMG){
  
  #VOX <- sl[[55]]
  
  full_cutoffs <- lapply(sl, function(VOX){
    
    all_vox_raw <- lapply(1:length(IMG), function(L){
      
      LAYER <- IMG[[L]]
      
      return(LAYER[VOX])
      
    }) %>%
      Reduce(function(x,y)c(x,y),.) %>%
      .[which(!.==0)]
    
    #t <- as.matrix(all_vox_raw)
    cutoff <- otsu_own(all_vox_raw, c(0,1), levels=512)
    # o <- otsu(as.matrix(all_vox_raw), c(0,1), levels=256)
    # oo <- otsu_own(as.matrix(all_vox_raw), c(0,1), levels=256)
    return(cutoff)
    
    #print(paste(o, "::", oo, o==oo, abs(o-oo)))
    
  })
  
  new_image <- IMG
  
  
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


# 
# n_main_dendrites = inter$n_main_dendrites
# main_dendrites=inter$main_dendrites
# IMG=noS_noMD_image
# elongated_dendrites=inter$elongated_dendrites
# raster_size=80
# process_full=T
# method="man"
# SOMA=SOMA
# nr_orig=opt$nr_orig
# nc_orig=opt$nc_orig
# main_vectors_full <- inter$main_vectors_full



binarize_image <- function(n_main_dendrites, main_dendrites, IMG,
                           elongated_dendrites, raster_size, process_full, 
                           method, SOMA, nr_orig, nc_orig, main_vectors_full){
  cat("\nmake binary image")
  #t0 <- Sys.time()
  borders_vox_raw <- lapply(1:n_main_dendrites, get_vox_of_borders, 
                            main_dendrites=main_dendrites,
                            nc_orig=nc_orig,
                            nr_orig=nr_orig,
                            IMG=IMG,
                            SOMA=SOMA)
  
  ## plot borders
  #ggplot(bind_rows(borders_vox_raw), aes(x, -y))+geom_tile()
  
  
  
  selection_vector <- rep(c(1:n_main_dendrites), 3) %>% 
    set_names(c((-n_main_dendrites+1):(2*n_main_dendrites)))
  cat("\n  create raster image")
  
  if(process_full==T){
    
    sl <- get_circular_raster(SOMA, 
                              main_vectors_full, 
                              nr=nr_orig, 
                              nc=nc_orig, 
                              n_main_dendrites, 
                              raster_size) %>%
      lapply(.,first)
    
  } else {
    
      dendrite_segments <- lapply(1:n_main_dendrites, create_raster, 
                              avs=raster_size,
                              borders_vox_raw=borders_vox_raw,
                              selection_vector=selection_vector, 
                              elongated_dendrites,
                              nr_orig=nr_orig,
                              nc_orig=nc_orig,
                              process_full=process_full)
      
      sl <- Reduce(function(x,y)append(x,y),dendrite_segments)
  }
  

  
  
  #t1 <- Sys.time()
  cat("\n  binarize")
  if(method=="man"){
    #binary_image <- bi2(sl, IMG)
    binary_image <- bi2(sl, IMG)
  } else if(method=="otsu"){
    binary_image <- bin_otsu(sl, IMG)
  } else if(method=="otsu2"){
    
    binary_image <- bin_otsu2(sl, IMG)
    
  } else {
    qnt <- str_match(method, "\\d+") %>% as.numeric()
    binary_image <- bi3(sl, IMG, qnt)
  }
  writeTIFF(binary_image, "f:/data_sholl_analysis/run/Tier11_3_3_apical_D/test_circular2.tif")
  
  return(binary_image)
}
#nr_orig=1040
recompile_segments <- function(sl, nr_orig){
  
  
  
  retsl <- lapply(sl, function(VOX){
    print(1)
    a <- lapply(VOX, get_xy_index, nr=nr_orig) %>% bind_rows()
    return(c(mnx=min(a$x), mxx=max(a$x), mny=min(a$y), mxy=max(a$y)))
    
  })
  
  df <- bind_rows(retsl) %>% rownames_to_column("name") %>% rowwise() %>%
    mutate(px=0.5*(mnx+mxx), py=0.5*(mny+mxy))
  ggplot(df, aes(x=px, y=py, label=name))+
    geom_text()+
    geom_segment(inherit.aes = F,
                 data=df,
                 aes(x=mnx, xend=mxx, y=mny, yend=mny))
  
}

extract_all_vectors <- function(MAIND){
  #print(1)
  raw <- lapply(MAIND, get_all_vectors) %>%
    Reduce(function(x,y)append(x,y),.)%>%
    bind_rows() #%>%
  if(nrow(raw)==0){return(NULL)}
  r2 <- raw %>%  mutate(id=c(1:nrow(.)))
  lapply(1:nrow(r2),function(nV){
    V <- r2[nV,]
    ha <- cppGetHA(V[["xe"]], V[["ye"]],V[["xs"]], V[["ys"]])
    return(c(V, ha=ha))
  }) %>% return()
}
# IMG=bin_noS_noMD_image
# circle=inter$trace_dend_cross_circle
# RADv=opt$trace_dend_cross_radius_vertical
# RADh=opt$trace_dend_cross_radius_horizontal
# nr_orig=opt$nr_orig
# AV <- inter$all_vectors[[2]]
assign_vector_voxels <- function(AV, IMG, circle, RADv, RADh, nr_orig){
  #print("a")
  AVF <-  lapply(AV, get_full_vector_voxels, IMG=IMG)
  #print("b")
  if(length(AVF)==0){return(NULL)}
  get_all_vector_speheres(bind_rows(AVF), 
                          RADv, 
                          RADh, 
                          #inter$trace_dend_cross_circle,
                          circle,
                          nr_orig, 
                          IMG)
}

get_surr_layer_doesntwork <- function(comb_main_dend, SOMA, nr_orig, IMG_screen, soma_radius, z_range, screen_subd_starts_circle){
  # t0 <- Sys.time()
  for(CT in 1:nrow(comb_main_dend)){
    VP <- comb_main_dend[CT,]
    sc <- get_circle_around_point(VP, screen_subd_starts_circle, nr_orig)
    if(CT==1){
      circ <- sc
    } else {
      circ <- c(circ, sc) %>% unique()
    }
  } 
  
  z_raw <- unique(comb_main_dend$z)
  res_vec <- lapply((z_raw-z_range):(z_raw+z_range), function(Z){
    
    all_sorrounding_vox <- circ %>% 
      lapply(function(IND){
        raw <- cppGetXYIndex(IND, nr_orig)
        
        i <- cppSelectIntensity(raw[["x"]], raw[["y"]],Z, IMG_screen)
        if(i==0){
          return(NULL)
        }
        d <- cppDistPts(raw[["x"]], raw[["y"]],Z, SOMA[["x"]], SOMA[["y"]], SOMA[["z"]])
        return(c(raw, z=Z, i=i, dist_to_soma=d))
      }) %>% compact() %>%
      bind_rows() %>%
      .[which(.$dist_to_soma>3*soma_radius),] %>%
      return()    
  }) %>%
    bind_rows()
  
}
#r <- det_rad
get_surr_layer_new <- function(comb_main_dend, SOMA, nr_orig, IMG_screen, soma_radius, z_range, r){
  # t0 <- Sys.time()
  for(CT in 1:nrow(comb_main_dend)){
    VP <- comb_main_dend[CT,]
    sc <- lapply(c(-r:r), function(X){
      y=sqrt(r^2-X^2)
      bottom <- round(VP[["y"]]-y)
      top <- round(VP[["y"]]+y)
      lapply(c(bottom:top), cppGetSingleIndex, x=X+round(VP[["x"]]), nr=nr_orig) %>%
        return()
    }) %>% 
      c(recursive=T)
    if(CT==1){
      circ <- sc
    } else {
      circ <- c(circ, sc) %>% unique()
    }
  } 
  
  z_raw <- unique(comb_main_dend$z)
  res_vec <- lapply((z_raw-z_range):(z_raw+z_range), function(Z){
    
    all_sorrounding_vox <- circ %>% 
      lapply(function(IND){
        raw <- cppGetXYIndex(IND, nr_orig)
        
        i <- cppSelectIntensity(raw[["x"]], raw[["y"]],Z, IMG_screen)
        if(i==0){
          return(NULL)
        }
        d <- cppDistPts(raw[["x"]], raw[["y"]],Z, SOMA[["x"]], SOMA[["y"]], SOMA[["z"]])
        return(c(raw, z=Z, i=i, dist_to_soma=d))
      }) %>% compact() %>%
      bind_rows() %>%
      .[which(.$dist_to_soma>3*soma_radius),] %>%
      return()    
  }) %>%
    bind_rows()
  
}


get_surr_layer <- function(comb_main_dend, SOMA, nr_orig, IMG_screen, soma_radius, z_range, det_rad){
  
  # t0 <- Sys.time()
  
  sequence <- seq(1, nrow(comb_main_dend), round(det_rad/3)) %>% c(.,nrow(comb_main_dend)) %>% unique()
  
  for(CT in sequence){
    #print(CT)
    VP <- comb_main_dend[CT,]
    sc <-lapply(c((VP[["z"]]-z_range):(VP[["z"]]+z_range)) %>% .[.>0], function(Z){
      r <- det_rad
      circ <- lapply(c(-r:r), function(X){
        y=sqrt(r^2-X^2)
        bottom <- round(VP[["y"]]-y)
        top <- round(VP[["y"]]+y)
        lapply(c(bottom:top), get_single_index, x=X+round(VP[["x"]]), nr=nr_orig) %>%
          return()
      }) %>% 
        c(recursive=T) %>%
        paste(Z, ., sep="_") %>%
        return()
    }) %>% 
      c(recursive=T)
    if(CT==1){
      res_vec <- sc
    } else {
      res_vec <- c(res_vec, sc) %>% unique()
    }
  }
  
  
  all_sorrounding_vox <- res_vec %>% 
    lapply(function(IND){
      #print(IND)
      s <- str_split(IND, "_") 
      z <- as.numeric(map_chr(s,1))
      raw <- cppGetXYIndex(as.numeric(map_chr(s,2)), nr_orig)
      i <- cppSelectIntensity(raw[["x"]], raw[["y"]],z, IMG_screen)
      if(i==0){
        #print(IND)
        #cat("    ----- 0")
        return(NULL)
      }
      d <- cppDistPts(raw[["x"]], raw[["y"]],z, SOMA[["x"]], SOMA[["y"]], SOMA[["z"]])
      return(c(raw, z=z, i=i, dist_to_soma=d))
    }) %>% compact() %>%
    bind_rows() %>%
    .[which(.$dist_to_soma>2*soma_radius),] %>%
    return()
  
  #print(Sys.time()-t0)
  
}

# 
# IMG=noS_image
# main_vectors =inter$main_vectors
# main_vectors_full=inter$main_vectors_full
# rem_rad=opt$subd_max_detection_distance-opt$subd_detection_depth-1
# nr_orig=opt$nr_orig
# nc_orig=opt$nc_orig
# 

remove_main_dendrites <- function(IMG, main_vectors, main_vectors_full, rem_rad, nr_orig, nc_orig){
  
  ## time improvement at its best... wont get faster
  
  comb_main_dend <- lapply(1:length(main_vectors), function(n){
    full_vecs <- main_vectors_full[[n]]
    full_dendrite <- bind_rows(full_vecs) %>%
      mutate_all(.funs = as.numeric) %>%
      mutate(ind=get_single_index(x,y,nr_orig)) 
    
    selector <- seq(1, nrow(full_dendrite), rem_rad/2) %>% c(.,nrow(full_dendrite)) %>% unique()
    
    return(full_dendrite[selector,])
    
  }) %>% bind_rows() #%>% as.matrix()
  
  cat("\nremoving main dendrites; about", nrow(comb_main_dend)*0.31, "s")
  
  rem_IMG <- lapply(1:length(IMG), function(n){
    
    #print(n)
    
    all_vox <- apply(comb_main_dend, 1, function(VP){
      r <- round(rem_rad+0.0175*(abs(round(VP[["z"]])-n))^2)
      circ <- lapply(c(-r:r), function(X){
        y=sqrt(r^2-X^2)
        bottom <- round(VP[["y"]]-y)
        top <- round(VP[["y"]]+y)
        lapply(c(bottom:top), get_single_index, x=X+round(VP[["x"]]), nr=nr_orig) %>%
          return()
      }) %>% unlist() %>%
        .[which(.<(nc_orig*nr_orig)&.>0)]#%>%
        #return()
      
      
      # if(max(circ)>nc_orig*nr_orig){print(VP)
      #   print(max(circ))}
      return(circ)
      
    }) %>%
      c(recursive=T) %>%
      unique() 
    
    
    
    LAYER <- IMG[[n]]
    LAYER[all_vox] <- 0
    return(LAYER)
  })
  return(rem_IMG)
  
}

main <- function(file, run_dir){
  rt <- list()
  rt[[1]] <- Sys.time()
  
  tmp_dir <- file.path(run_dir, "tmp")
  dir.create(tmp_dir)
  
  #images <- list(
  raw_image=readTIFF(file, all=T)
  #  )
  # run_dir <- "f:/data_sholl_analysis/test/parabolar"
  # load(file.path(run_dir, "images/images_new.RData"))
  opt <- set_options(raw_image)
  inter <- list()
  ## find all somata in the image
  somata <- get_somata(opt$soma_xy_detection_cube_radius*2+1, 
                       opt$soma_z_detection_radius, 
                       opt$soma_z_detection_degree_steps, 
                       raw_image)
  rt[[2]] <- Sys.time()
  ## loop over somata ... now just do it for one
  SOMA <- somata[1,] %>% mutate_all(.funs = round) 
  cat(paste("\n ", nrow(somata), "somata detected"))
  ## find start sites of main dendrites from soma
  inter$main_dendrites_raw <- find_dendritic_start_sites(SOMA, raw_image, 
                                                         #0.98
                                                         opt$cutoff_soma_region
                                                         )
  
  ### ab hier könnte raw image weg
  
  rt[[3]] <- Sys.time()
  inter$main_dendrites <- inter$main_dendrites_raw[[1]] %>% 
    mutate(z=SOMA[["z"]])  %>%
    arrange(h_angle) %>%
    rownames_to_column("dendrite_id")
  
  
  
  inter$soma_radius <- inter$main_dendrites_raw[[2]]    
  inter$n_main_dendrites <- nrow(inter$main_dendrites)
  
  cat(paste("\n ", inter$n_main_dendrites, "main dendrites detected\n  elongating main dendrites:"))
  
  ## elongate main dendrites
  inter$elongated_dendrites <- apply(inter$main_dendrites, 1, elongate_dendrite,
                                     SOMA=SOMA, IMG=raw_image, nr=opt$nr_orig, nc=opt$nc_orig)
  rt[[4]] <- Sys.time()
  export_dendrites(inter$elongated_dendrites, file.path(run_dir, "elongated_dendrites.csv"), SOMA=SOMA)
  
  write_tsv(tibble(), file=file.path(run_dir, paste0("soma_z_", SOMA$z, ".tsv")))
  
# }  
#     #print(1)
# 
# secmain <- function(){
  inter$main_vectors_raw <- lapply(inter$elongated_dendrites, define_full_vector_info,
                                   IMG=raw_image, SOMA=SOMA)
   #print(2)
  inter$main_vectors <- lapply(inter$main_vectors_raw, nth, 1)
  inter$main_vectors_full=lapply(inter$main_vectors_raw, nth, 2)
  inter$main_vectors_df <- lapply(inter$main_vectors, bind_rows)
  #print(3)
  ## 34 sec
  noS_image <- remove_soma(SOMA, inter$soma_radius, raw_image, opt$nr_orig, opt$nc_orig)
  rt[[5]] <- Sys.time()
  writeTIFF(raw_image, file.path(run_dir, file %>% str_split("/") %>% unlist() %>% last()))
  rm(raw_image)
  #print(4)

  noS_noMD_image <- remove_main_dendrites(noS_image,
                                          inter$main_vectors,
                                          inter$main_vectors_full,
                                          opt$subd_max_detection_distance-opt$subd_detection_depth-1,
                                          opt$nr_orig,
                                          opt$nc_orig)


  writeTIFF(noS_image, file.path(tmp_dir, "noS_image.tif"))
  rm(noS_image)
  #writeTIFF(noS_noMD_image, file.path(tmp_dir, "noS_noMD_image.tif"))

  ## 2 mins until here

  rt[[6]] <- Sys.time()

  #print(1)



  ## remove noSOMA image
  #images <- images[which(!names(images %in% c("noS_image")))]

  #print(2)
  bin_noS_noMD_image <- binarize_image(inter$n_main_dendrites,
                                       inter$main_dendrites,
                                       noS_noMD_image,
                                       inter$elongated_dendrites,
                                       80,
                                       T,
                                       "man",
                                       SOMA=SOMA,
                                       opt$nr_orig,
                                       opt$nc_orig,
                                       inter$main_vectors_full)
  rm(noS_noMD_image)
  writeTIFF(bin_noS_noMD_image, file.path(run_dir, "bin_noS_noMD_image.tif"))
  rt[[7]] <- Sys.time()
  
  #test <- lapply(bin_noS_noMD_image, as.matrix)
  
#### hier rest wieder rein machen
  



  #writeTIFF(bin_noS_noMD_image, file.path(tmp_dir, "bin_noS_noMD_image.tif"))
  
  #med_bin_noS_image_p <- apply_3d_median_filter(bin_noS_image_p, 2)

  inter$subd_max_detection_circle <- create_circle(opt$subd_max_detection_distance)

  # export all images
  ## 4 mins until here

  med_bin_noS_noMD_image <- apply_3d_median_filter(bin_noS_noMD_image, 2)
  writeTIFF(med_bin_noS_noMD_image, file.path(run_dir, "med_bin_noS_noMD_image.tif"))
  rm(bin_noS_noMD_image)


  noS_image <- readTIFF(file.path(tmp_dir, "noS_image.tif"), all=T)
  bin_noS_image_p <- binarize_image(inter$n_main_dendrites,
                                    inter$main_dendrites,
                                    noS_image,
                                    inter$elongated_dendrites,
                                    80,
                                    F,
                                    "man",
                                    SOMA=SOMA,
                                    opt$nr_orig,
                                    opt$nc_orig)

  rm(noS_image)
  #
  rt[[8]] <- Sys.time() 
  MASTER <- lapply(1:inter$n_main_dendrites, find_subdendritic_starts_new,
                   main_vectors=inter$main_vectors,
                   main_vectors_df=inter$main_vectors_df,
                   main_vectors_full=inter$main_vectors_full,
                   SOMA=SOMA,
                   IMG_screen=med_bin_noS_noMD_image,
                   IMG_adj=bin_noS_image_p,
                   soma_radius=inter$soma_radius,
                   MPTS=opt$subd_cluster_mpt,
                   EPS=opt$subd_cluster_eps,
                   det_rad=opt$subd_max_detection_distance,
                   z_range=opt$subd_detection_vertical_range,
                   nr_orig=opt$nr_orig,
                   screen_subd_starts_circle=inter$subd_max_detection_circle
  )
  # print("Master created")
  #
  rt[[9]] <- Sys.time()
  
  time <- tibble(step=paste(c(1:(length(rt)-1)),c(2:length(rt)), sep="-"), time=diff(unlist(rt)))
  write_tsv(time, file=file.path(run_dir, "required_time_00.tsv"))
  #### test
  
  # test <- find_subdendritic_starts_new(1,
  # main_vectors=inter$main_vectors,
  # main_vectors_df=inter$main_vectors_df,
  # main_vectors_full=inter$main_vectors_full,
  # SOMA=SOMA,
  # IMG_screen=med_bin_noS_noMD_image,
  # IMG_adj=bin_noS_image_p,
  # soma_radius=inter$soma_radius,
  # MPTS=opt$subd_cluster_mpt,
  # EPS=opt$subd_cluster_eps,
  # det_rad=opt$subd_max_detection_distance,
  # z_range=opt$subd_detection_vertical_range,
  # nr_orig=opt$nr_orig,
  # screen_subd_starts_circle=inter$subd_max_detection_circle)
  
  
  ####
  
  export_structure(MASTER, inter$main_vectors, SOMA,
                   file.path(run_dir, "subdendrite_starts.csv"))
  
  rm(med_bin_noS_noMD_image)
  rm(bin_noS_image_p)

  bin_noS_noMD_image <- readTIFF(file.path(run_dir, "bin_noS_noMD_image.tif"), all=T)
  #images <- images[which(!names(images %in% c("med_bin_noS_noMD_image", "bin_noS_image_p")))]
  rt6 <- Sys.time()

  #save(images, file=file.path(run_dir, "images/images_new.RData"))

  inter$all_vectors <- lapply(MASTER, extract_all_vectors)
  #print(1)

  inter$trace_dend_cross_circle <- create_circle(opt$trace_dend_cross_radius_horizontal)
  #print(2)
  inter$trace_dend_detection_circle <- create_circle(opt$trace_detection_distance+opt$trace_detection_depth)

  inter$all_points <- lapply(inter$all_vectors, assign_vector_voxels,
                             IMG=bin_noS_noMD_image,
                             circle=inter$trace_dend_cross_circle,
                             RADv=opt$trace_dend_cross_radius_vertical,
                             RADh=opt$trace_dend_cross_radius_horizontal,
                             nr_orig=opt$nr_orig)



  #print(3)



  #print("prepared for tracing")




  traced_MASTER <- lapply(1:inter$n_main_dendrites, function(n){
    cat(paste("\nmain-dendrite:", n))
    trace_subdendrites(
      IMG=bin_noS_noMD_image,
      MAIND=MASTER[[n]],
      MV=inter$main_vectors_df[[n]],
      main_vectors_full=inter$main_vectors_full[[n]],

      all_vectors = inter$all_vectors[[n]],
      all_points =inter$all_points[[n]],

      EPS_orig=opt$trace_cluster_eps,
      MPTS_orig=opt$trace_cluster_mpt,
      INC=opt$trace_detection_depth,

      RADv=opt$trace_dend_cross_radius_vertical,
      RADh=opt$trace_dend_cross_radius_horizontal,
      cross_circle=inter$trace_dend_cross_circle,

      det_rad=opt$trace_detection_distance,
      RESC_DIST=opt$trace_rescore_dist,
      RESC_ANG=opt$trace_rescore_angle,
      nr_orig=opt$nr_orig,
      circle=inter$trace_dend_detection_circle,

      SOMA=SOMA)


  })
  rt7 <- Sys.time()

  export_structure(traced_MASTER, inter$main_vectors, SOMA,
                   file.path(run_dir, "traced_dendrites.csv"))

  #writeTIFF(bin_noS_noMD_image, file.path(run_dir, "main_binary.tif"))
  
}





# main_vectors_full=main_vectors_full
# SOMA=SOMA
# start=start
# end=end



get_main_vector_angles_at_segment_range <- function(n, main_vectors_full, SOMA, start, end){
  
  #print(n)
  
  mvf <- bind_rows(main_vectors_full[[n]]) %>%
    mutate(dts=cppDistPts(x,y,1,SOMA$x, SOMA$y, 1)) %>%
    filter(between(dts, start, end))
  
  
  
  if(nrow(mvf)>=2){
    
    adj_ha <- cppGetHA(mvf[[nrow(mvf), "x"]],mvf[[nrow(mvf), "y"]],
                    #mvf[[1, "x"]], mvf[[1, "y"]])
                    SOMA$x, SOMA$y)
    
   return(list(adj_ha, round(mean(mvf$z))))        
  } else {
    
    adj_ha <- cppGetHA(main_vectors_full[[n]] %>% .[[length(.)]] %>% .[[nrow(.), "x"]],
                    main_vectors_full[[n]] %>% .[[length(.)]] %>% .[[nrow(.), "y"]],
                    #main_vectors_full[[n]] %>% .[[length(.)]] %>% .[[1, "x"]],
                    #main_vectors_full[[n]] %>% .[[length(.)]] %>% .[[1, "y"]]
                    SOMA$x, SOMA$y
                    )
    return(list(adj_ha, main_vectors_full[[n]] %>% .[[length(.)]] %>% .[[nrow(.), "z"]]))
  }
  
  
}  
get_circle_cuts_old <- function(n, ajd_ha_main_vecs_raw, theoretical_band_size){  
  ajd_ha_main_vecs <- lapply(ajd_ha_main_vecs_raw, first) %>% unlist()
  ## get all cuts though circle
  if(n==1){
    h_angle_p <- adj_deg(0.5*sum(ajd_ha_main_vecs[n], ajd_ha_main_vecs[n+1]))
    h_angle_n <- adj_deg(ajd_ha_main_vecs[n]+0.5*(abs(360-ajd_ha_main_vecs[n])+ajd_ha_main_vecs[nrow(main_dendrites)]))
  } else if(n==nrow(main_dendrites)){
    h_angle_p <- adj_deg(ajd_ha_main_vecs[n]+0.5*(abs(360-ajd_ha_main_vecs[n])+ajd_ha_main_vecs[1]))
    h_angle_n <- adj_deg(0.5*sum(ajd_ha_main_vecs[n], ajd_ha_main_vecs[n-1]))
  } else {
    h_angle_p <- 0.5*sum(ajd_ha_main_vecs[n], ajd_ha_main_vecs[n+1])
    h_angle_n <- 0.5*sum(ajd_ha_main_vecs[n], ajd_ha_main_vecs[n-1])
  }
  return(list(c(h_angle_n, ajd_ha_main_vecs[n], ajd_ha_main_vecs_raw[[n]][[2]]),
              c(ajd_ha_main_vecs[n], h_angle_p, ajd_ha_main_vecs_raw[[n]][[2]])))
  
}


get_circle_cuts <- function(ajd_ha_main_vecs_raw, theoretical_band_size){
  
  max_pix <- 10000
  
  ajd_ha_main_vecs <- lapply(ajd_ha_main_vecs_raw, first) %>% unlist()
  
  diff_vec <- diff(ajd_ha_main_vecs) %>% c(., 360-ajd_ha_main_vecs[length(ajd_ha_main_vecs)]+ajd_ha_main_vecs[1])
  
  main_cuts <- c(ajd_ha_main_vecs, ajd_ha_main_vecs+0.5*diff_vec) %>% sapply(., adj_deg) %>% sort()
  
  diff_vec2 <- diff(main_cuts) %>% c(., 360-main_cuts[length(main_cuts)]+main_cuts[1])
  
  ts <- theoretical_band_size*(diff_vec2/360)
  
  lapply(1:length(main_cuts), function(N){
    if(N!=length(main_cuts)){
      a1 <- main_cuts[N]
      a2 <- main_cuts[N+1]
      npx <- ts[N]
      if(npx>max_pix){
        subs <- ceiling(npx/max_pix)
        deg_step <- (a2-a1)/subs
        return(seq(a1, a2, deg_step))
      } else {
        return(c(a1,a2))
      }      
    } else {
      a1 <- main_cuts[N]
      a2 <- main_cuts[1]
      npx <- ts[N]
      if(npx>max_pix){
        subs <- ceiling(npx/max_pix)
        deg_step <- (a2+(360-a1))/subs
        return(seq(a1, a2+360, deg_step) %>% lapply(., adj_deg))
      } else {
        return(c(a1,a2))
      }
    }

  }) %>% c(recursive=T) %>% round(1) %>% unique() %>% sort() %>%
    return()
  
}




#main_vectors_full <- inter$main_vectors_full

# #
# nr=nr_orig
# nc=nc_orig
# avs=raster_size


get_circular_raster <- function(SOMA, 
                                main_vectors_full, 
                                nr, 
                                nc, 
                                n_main_dendrites, 
                                avs){
  
  max_distance_to_any_corner <- lapply(list(c(1,1), c(1,nr), c(nc, nr), c(nc, 1)), function(C){
    cppDistPts(C[[1]], C[[2]], 1, SOMA$x, SOMA$y, 1)
  }) %>% unlist() %>% max() %>% ceiling()
  
  n_segments <- ceiling(max_distance_to_any_corner/avs)
  t0 <- Sys.time()
  all_pix_of_image <- t(sapply(c(1:(nr*nc)), cppGetXYIndex, nr=nr)) %>%
    cbind(mapply(cppGetHA, .[,"x"], .[,"y"], xs=SOMA$x, ys=SOMA$y)) %>%
    cbind(mapply(cppDistPts, .[,"x"], .[,"y"],zs=1, xe=SOMA$x, ye=SOMA$y, ze=1)) %>%
    cbind(1:(nr*nc))
  #print(Sys.time()-t0)
  
  sl_list <- list()
  
  for(MS in c(1:n_segments)){
  #for(MS in 1:5){
    #print(MS)
    
    #cat(paste("\nMAIn SEG", MS, "\n"))
    start <- MS*avs-avs
    end <- MS*avs
    
    theoretical_band_size <- (pi*end^2)-(pi*start^2)
    
    
    ajd_ha_main_vecs_raw <- lapply(1:n_main_dendrites, get_main_vector_angles_at_segment_range, 
                               main_vectors_full=main_vectors_full,
                               SOMA=SOMA,
                               start=start,
                               end=end) #%>% unlist()
    
    #print(unlist(ajd_ha_main_vecs_raw))
    
    cc_raw <- get_circle_cuts( 
                 ajd_ha_main_vecs_raw=ajd_ha_main_vecs_raw,
                 theoretical_band_size=theoretical_band_size)# %>%
      #Reduce(function(x,y)append(x,y),.)
    cc <- matrix(c(cc_raw[1:length(cc_raw)],
                   c(cc_raw[2:length(cc_raw)], cc_raw[1])), ncol=2)
    
    
    sl_per_dist <- lapply(1:nrow(cc), function(nCP){
      #print(nCP)
      min_angle <- cc[nCP,1]
      max_angle <- cc[nCP,2]
      if(nCP==nrow(cc)){
        rel <- all_pix_of_image[which(
          (between(all_pix_of_image[,3], min_angle, 360)|
             between(all_pix_of_image[,3], 0, max_angle))&
            between(all_pix_of_image[,4], start, end)
        ),] 
      } else {
        rel <- all_pix_of_image[which(
          between(all_pix_of_image[,3], min_angle, max_angle)&between(all_pix_of_image[,4], start, end)
        ),] 
        
      }
      if(class(rel)=="numeric"){
        rel <- rel %>% as.matrix() %>% t()
      }
      if(nrow(rel)==0){
        #cat("zero")
        return(NULL)
        #print(1)
      } else {
        return(list(pix=rel[,5], s=MS))
        ## for controling (use ggplot below)
        #return(rel[,1:2]  %>% as_tibble() %>% mutate(n=nCP))
      }
    }) %>% compact()
    #print(length(sl_per_dist))

    sl_list <- append(sl_list, sl_per_dist)
    
  }
  
  return(sl_list)
  
}














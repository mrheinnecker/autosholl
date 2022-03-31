
library(tiff)
library(tidyverse)
library(gridExtra)
library(cowplot)


source("/Users/Marco/git_repos/autosholl/fncts.R")

#file <- "/Users/Marco/Dropbox/Studium/Master/Praktikum_Mueller/example_data/Tier1_2_3_apical_D.tif"
files <- list.files("f:/data_sholl_analysis/single_soma", pattern = "*", full.names = T)

#file <- "f:/data_sholl_analysis/Tiffs/Apical/Deep/Tier4_5_2_apical_D.tif"

file <- files[1]
soma_xy_detection_cube_radius <- 15
soma_z_detection_radius <- 100
soma_z_detection_degree_steps <- 0.05
  
  
results <- lapply(files, function(file){
  start_time <- Sys.time()
  print(file)

  full_image <- readTIFF(file, all=T)
  nr <- nrow(full_image[[1]])
  somata <- get_somata(soma_xy_detection_cube_radius*2+1, 
                       soma_z_detection_radius, 
                       soma_z_detection_degree_steps, 
                       full_image)
  
  print("soma detected")
  SOMA <- somata[1,]
  
  apply(somata,1, function(SOMA){
    
    main_dendrites_raw <- find_dendritic_start_sites(SOMA) 
    main_dendrites <- main_dendrites_raw[[1]] %>%
      mutate(z=SOMA[["z"]]) %>%
      rownames_to_column("dendrite_id")
    
    control_plot <- main_dendrites_raw[[2]]
    
    #DENDRITE <- main_dendrites[5,]
    elongated_dendrites <- apply(main_dendrites, 1, function(DENDRITE){

      cat(paste("\ndendrite number:",DENDRITE[["dendrite_id"]]))
      
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
                        steps)
      
      #return(list(knot_list, full_vector_list))
    })
    
    ### now select areas for cutoff use
    ELD <- elongated_dendrites[[4]] 
    
    vector_pos <- lapply(c(1:(length(ELD)-1)), function(n){
      
      c(
        xs=ELD[[n]][["x"]],
        xe=ELD[[n+1]][["x"]],
        ys=ELD[[n]][["y"]],
        ye=ELD[[n+1]][["y"]],
        zs=ELD[[n]][["z"]],
        ze=ELD[[n+1]][["z"]]
      ) %>%
        return()
      
    }) %>% bind_rows() %>% rownames_to_column("id")
    
    avs <- 200    
    vectors <- lapply(1:length(ELD), function(PT){
      if(PT==length(ELD)){return(NULL)}
      return(c(ELD[[PT+1]][c("x", "y", "z")]-ELD[[PT]][c("x", "y", "z")],
               h_angle=ELD[[PT+1]]["h_angle"],
               l=sqrt(sum(abs(ELD[[PT+1]][c("x", "y")]-ELD[[PT]][c("x", "y")])^2))))
      
    }) %>%
      compact()
    xnorm_vector <- lapply(vectors, function(V){
      return(V[c(1:3)]/abs(V["x"]))
    })
    ynorm_vector <-  lapply(vectors, function(V){
      return(V[c(1:3)]/abs(V["y"]))
    }) 
    overall_vector <- ELD[[length(ELD)]][c("x", "y", "z")]-ELD[[1]][c("x", "y", "z")]
    if(abs(overall_vector[1])>abs(overall_vector[2])){
        ## horizontal dendrite oriantation
      overall_length <- overall_vector["x"]
      n_segments <- abs(round(overall_length/avs))
      use_length <- round(overall_length/n_segments)
      full_vecs <- bind_rows(vectors) %>%
        pull(x) %>% cumsum()      
      #left_vecs <- full_vecs
      res_list <- list()
      
      ######### new
      
      for(n in 1:n_segments){
        #vec <- c(((n-1)*use_length):(n*use_length))+ELD[[1]]["x"]
        
        start <-(n-1)*use_length+ELD[[1]]["x"]
        end <-   n*use_length+ELD[[1]]["x"]
        vecs <- vector_pos %>%
          mutate(t=ifelse(xs %in% c(start:end)|xe %in% c(start:end), T, F)) %>%
          filter(t==T) %>%
          pull(id) 
        if(length(vecs)==0){
          res_list[[n]] <- as.numeric(res_list[[n-1]])
        } else {
          res_list[[n]] <- as.numeric(vecs)
        }
        
      }
      
      
      ######### old
      # for(n in rev(1:n_segments)){
      #   which_vecs <- which(ceiling(left_vecs/use_length)==n)
      #   if(length(which_vecs)==0){
      #     res_list[[n]] <- length(left_vecs)
      #   } else {
      #     res_list[[n]] <- c(max(which_vecs)+1, which_vecs) %>% .[which(.<=length(full_vecs))]
      #     left_vecs <- left_vecs[-c((which_vecs[1]+1):length(full_vecs))]
      #   }
      # }
      
      final_list <- list()
      for(n in c(1:n_segments)){
        print(n)
        xs <- ELD[[1]]["x"]+n*use_length-use_length+1
        xe <- ELD[[1]]["x"]+n*use_length 
        ys <- ifelse(n==1,
                     as.numeric(ELD[[1]]["y"]),
                     as.numeric(final_list[[n-1]][nrow(final_list[[n-1]]),"y"]))  
        
        rv <- use_length/abs(use_length)
        
        rel_vecs <- full_vecs[res_list[[n]]]
        rl <- list()
        last_end <- xs-rv
        for(i in 1:length(rel_vecs)){
          #print(i)
          if(i==1){
            start <- xs
          } else {
            start <- last_end+1
          } 
          
          if(i==length(rel_vecs)){
            end <- xe 
          } else {
            end <- ELD[[1]]["x"]+rel_vecs[i]
          }
          last_end <- end
          rl[[i]] <- rep(xnorm_vector[[res_list[[n]][i]]]["y"], abs(end-start+1))
        }
        
        line <- tibble(x=xs+c(1:abs(use_length)*rv)) %>%
          mutate(fac=Reduce(function(x,y)c(x,y),rl)) %>%
          mutate(y=ceiling(ys+cumsum(fac)))
        
        final_list[[n]] <- line
        
        med <- round(0.5*(max(line$y)+min(line$y)))
        
        top <- med+avs
        
        all_vox <- apply(line, 1, function(C){
          
          xn <- (C[["x"]]-1)*nr
          yn <- 0
          return(c((C[["y"]]):(top))+xn+yn)
          
        }) %>% Reduce(function(x,y)c(x,y),.)
        
        ########################################################################
        comp_img <- lapply(full_image, function(LAYER){
          
          LAYER[all_vox] <- 1
          return(LAYER)
        })
        writeTIFF(comp_img, paste0("f:/data_sholl_analysis/test/check",n,".tif"))
        ########################################################################
      }
        
        
        
        
        
        
        
        # lapply(c(min(line$x):max(line$x)), function(x){
        #   print(x)
        # })
        
        
        # tl <- tibble(x=c(min(line$x):max(line$x)),
        #              y=med+avs) %>%
        #   left_join(line, by="x")
        
        
        
        # start_first_vec <- vectors[[min(res_list[[n]])]]
        # 
        # 
        # med <- 
        # 
        # 
        # tl <- tibble(x=)
        # bl <-
        # ll <-
        # rl <- 
        # 
        # cl <-   
          
      
 
      
        hl <- ELD[[length(ELD)]]["x"]-ELD[[1]]["x"] 
        
    } else {
        ## vertical dedrite oriantation
      overall_length <- overall_vector["y"]    
      vl <- ELD[[length(ELD)]]["y"]-ELD[[1]]["y"]
    }
    
    

    

    
    
    
    vector_lengths <- lapply(1:length(ELD), function(PT){
      if(PT==length(ELD)){return(NULL)}
      sqrt(sum(abs(ELD[[PT+1]][c("x", "y", "z")]-ELD[[PT]][c("x", "y", "z")])^2)) %>%
        return()
      
    }) 
    
    
    
    
    
    
    
    
    
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
    
  })
  # return(somata %>%
  #          mutate(img=file %>% str_split("/") %>% unlist() %>% last())
  #          )
  write_csv(control_data_fiji, 
            file="f:/data_sholl_analysis/test/example_data/dendrites.csv")
  
print(Sys.time()-start_time)
}) #%>% 
  #bind_rows()


control_new <- elongated_dendrites%>% bind_rows()


control_new <- full_vector_list %>% bind_rows()


p <- control_plot+
  geom_tile(inherit.aes=F, data=control_new, 
            aes(x=x, y=y, fill=z), width=10, height=10)


to_exp <- control_new %>% select(x,y) %>% mutate(y=nrow(full_image[[1]])-y)

write_csv(control_data_fiji, 
          file="c:/Users/Marco/Dropbox/Studium/Master/Praktikum_Mueller/dendrites.csv")


# SOMA <- somata[1,]
# DENDRITE <- main_dendrites[5,]
# 
# pos <- knot_list[[2]]
# n_vc <- 10
# n_hc <- 60
# 
# #source("/Users/Marco/git_repos/autosholl/fncts.R")
# knot_list <- list()
# full_vector_list <- list()
# next_pos <- c(x=DENDRITE[["x"]],y=DENDRITE[["y"]],z=DENDRITE[["z"]], h_angle=DENDRITE[["maxima"]], v_angle=0)
# knot_list[[1]] <- next_pos
# knot_list[[2]] <- next_pos
# c <- 2
# while(c<3|sum(knot_list[[c-1]][1:3]==knot_list[[c]][1:3])!=3){
#   print(c-1)
#   st <- Sys.time()
#   raw_dendrite <- elongate_dendrite(next_pos, n_vc, n_hc)
#   next_pos <- raw_dendrite[[1]]
#   full_vector_list[[c]] <- raw_dendrite[[2]]
#   c <- c+1 
#   knot_list[[c]] <- next_pos %>% set_names(c("x","y","z", "h_angle", "v_angle"))
#   print(Sys.time()-st)
# }


control_new <- full_vector_list %>% bind_rows()


p <- control_plot+
  geom_tile(inherit.aes=F, data=control_new, 
            aes(x=x, y=y, fill=z))






# 
# 
# control_data <- bind_rows(knot_list) %>% 
#                          rownames_to_column("n")
# 
# p <- control_plot+geom_tile(inherit.aes=F, data=control_data, 
#                             aes(x=x, y=y),fill="red")+
#   geom_text(inherit.aes=F, data=bind_rows(knot_list) %>% 
#               rownames_to_column("n"), aes(x=x, y=y, label=n))
# 



pdf(file = paste0("c:/Users/Marco/Dropbox/Studium/Master/Praktikum_Mueller/",
                  file %>% str_split("/") %>% unlist() %>% last() %>% str_replace(".tif", "_dendrite_elongation.pdf"),
                  ".pdf"),
    width = 40, height=8)
grid.arrange(
  p
)
dev.off()








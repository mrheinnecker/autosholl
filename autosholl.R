### this is the main script of autosholl
## load functions and dependencies
source("/Users/Marco/git_repos/autosholl/fncts.R")
load_dependencies()
Rcpp::sourceCpp(file="c:/Users/Marco/git_repos/autosholl/cpp_fncts.cpp")
## load input files
files <- list.files("f:/data_sholl_analysis/single_soma", pattern = "*", full.names = T)
file <- files[1]
## three main objects are used... list of all images; global options; intermediate data
#images <- list(raw_image=readTIFF(file, all=T))

load("f:/data_sholl_analysis/test/images/images.RData")


opt <- set_options(images$raw_image)
inter <- list()
## find all somata in the image
somata <- get_somata(opt$soma_xy_detection_cube_radius*2+1, 
                     opt$soma_z_detection_radius, 
                     opt$soma_z_detection_degree_steps, 
                     images$raw_image)
## loop over somata ... now just do it for one
SOMA <- somata[1,] %>% mutate_all(.funs = round) 

## find start sites of main dendrites from soma
inter$main_dendrites_raw <- find_dendritic_start_sites(SOMA, images$raw_image)
inter$main_dendrites <- inter$main_dendrites_raw[[1]] %>% 
  mutate(z=SOMA[["z"]])  %>%
  arrange(h_angle) %>%
  rownames_to_column("dendrite_id")

inter$soma_radius <- inter$main_dendrites_raw[[3]]    
inter$n_main_dendrites <- nrow(inter$main_dendrites)

## elongate main dendrites
inter$elongated_dendrites <- apply(inter$main_dendrites, 1, elongate_dendrite,
                             SOMA=SOMA, IMG=images$raw_image)

# export_dendrites(elongated_dendrites, "f:/data_sholl_analysis/test/dendrites/main_dend.tsv")

inter$main_vectors_raw <- lapply(inter$elongated_dendrites, define_full_vector_info, IMG=images$raw_image)
inter$main_vectors <- lapply(inter$main_vectors_raw, nth, 1)
inter$main_vectors_full=lapply(inter$main_vectors_raw, nth, 2)
inter$main_vectors_df <- lapply(inter$main_vectors, bind_rows) 


images$noS_image <- remove_soma(SOMA, inter$soma_radius, images$raw_image, opt$nr_orig, opt$nc_orig)
images$noS_noMD_image <- remove_main_dendrites(images$noS_image, 
                                               inter$main_vectors, 
                                               inter$main_vectors_full, 
                                               opt$subd_max_detection_distance-opt$subd_detection_depth-1, 
                                               opt$nr_orig, 
                                               opt$nc_orig)


images$bin_noS_image_p <- binarize_image(opt, inter, 
                                          inter$n_main_dendrites, 
                                          images$noS_image, 
                                          inter$elongated_dendrites,80, F, "man")


images$bin_noS_noMD_image <- binarize_image(opt, inter, 
                                            inter$n_main_dendrites, 
                                            images$noS_noMD_image, 
                                            inter$elongated_dendrites, 80, T, "man")

images$med_bin_noS_noMD_image <- apply_3d_median_filter(images$bin_noS_noMD_image, 2)


## export all images
# 
# run_dir <- "f:/data_sholl_analysis/test/parabolar"
# writeTIFF(images$raw_image, file.path(run_dir, "images/raw.tif"))
# writeTIFF(images$noS_image, file.path(run_dir, "images/noS.tif"))
# writeTIFF(images$noS_noMD_image, file.path(run_dir, "images/noS_noMD.tif"))
# writeTIFF(images$bin_noS_image_p, file.path(run_dir, "images/bin_noS_partly.tif"))
# writeTIFF(images$bin_noS_noMD_image, file.path(run_dir, "images/bin_noS_noMD.tif"))
# writeTIFF(images$med_bin_noS_noMD_image, file.path(run_dir, "images/bin_med_noS_noMD.tif"))
# 


# MASTER <- lapply(1:inter$n_main_dendrites, find_subdendritic_starts_old, 
#                  main_vectors=inter$main_vectors,
#                  main_vectors_df=inter$main_vectors_df,
#                  main_vectors_full=inter$main_vectors_full,
#                  SOMA=SOMA,
#                  IMG_screen=images$med_bin_noS_noMD_image,
#                  IMG_adj=images$bin_noS_image_p,
#                  soma_radius=inter$soma_radius,
#                  minPTS=opt$subd_cluster_mpt,
#                  EPS=opt$subd_cluster_eps,
#                  det_rad=opt$subd_detection_distance,
#                  z_range=opt$subd_detection_vertical_range,
#                  depth=opt$subd_detection_depth)


MASTER <- lapply(1:inter$n_main_dendrites, find_subdendritic_starts_new, 
                 main_vectors=inter$main_vectors,
                 main_vectors_df=inter$main_vectors_df,
                 main_vectors_full=inter$main_vectors_full,
                 SOMA=SOMA,
                 IMG_screen=images$med_bin_noS_noMD_image,
                 IMG_adj=images$bin_noS_image_p,
                 soma_radius=inter$soma_radius,
                 minPTS=opt$subd_cluster_mpt,
                 #minPTS=17,
                 EPS=opt$subd_cluster_eps,
                 det_rad=opt$subd_max_detection_distance,
                 z_range=opt$subd_detection_vertical_range,
                 #depth=opt$subd_detection_depth
                 nr_orig=opt$nr_orig
                 )



#save(images, file=file.path(run_dir, "images/images_new.RData"))




export_structure(MASTER, inter$main_vectors, file.path(run_dir, "dendrites/subdendrite_starts_new_surr.csv"))

traced_MASTER <- lapply(1:inter$n_main_dendrites, trace_subdendrites, 
                        IMG=images$bin_noS_noMD_image,
                        MASTER=MASTER,
                        main_vectors_df=inter$main_vectors_df,
                        main_vectors_full=inter$main_vectors_full,
                        EPS=opt$trace_cluster_eps,
                        MPTS=opt$trace_cluster_mpt,
                        INC=opt$trace_detection_depth,
                        det_rad=opt$trace_detection_distance)




















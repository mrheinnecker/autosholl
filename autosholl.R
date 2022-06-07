### this is the main script of autosholl
## load functions and dependencies
source("/Users/Marco/git_repos/autosholl/fncts.R")
load_dependencies()
Rcpp::sourceCpp(file="c:/Users/Marco/git_repos/autosholl/cpp_fncts.cpp")
## load input files
files <- list.files("f:/data_sholl_analysis/single_soma", pattern = "*", full.names = T)
file <- files[1]
## three main objects are used... list of all images; global options; intermediate data
images <- list(raw_image=readTIFF(file, all=T))
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
                                               opt$subd_detection_distance-opt$subd_detection_depth-1, 
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


 # writeTIFF(images$bin_noS_image_p, "f:/data_sholl_analysis/test/images/DR30_depth3/bin_noS.tif")
 # writeTIFF(images$bin_noS_noMD_image, "f:/data_sholl_analysis/test/images/DR30_depth3/bin_noS_noMD.tif")
 # writeTIFF(images$med_bin_noS_noMD_image, "f:/data_sholl_analysis/test/images/DR30_depth3/med_bin_noS_noMD.tif")
# images$bin_nosoma_image <- readTIFF("f:/data_sholl_analysis/test/intermediate/binarized.tif", all=T)
# images$bin_rem <- readTIFF(paste0("f:/data_sholl_analysis/test/binary/bin_90_25.tif"), all=T)
# 

MASTER <- lapply(1:inter$n_main_dendrites, find_subdendritic_starts, 
                 main_vectors=inter$main_vectors,
                 main_vectors_df=inter$main_vectors_df,
                 main_vectors_full=inter$main_vectors_full,
                 SOMA=SOMA,
                 IMG_screen=images$med_bin_noS_noMD_image,
                 IMG_adj=images$bin_noS_image_p,
                 soma_radius=inter$soma_radius,
                 minPTS=opt$subd_cluster_mpt,
                 EPS=opt$subd_cluster_eps,
                 det_rad=opt$subd_detection_distance,
                 z_range=opt$subd_detection_vertical_range,
                 depth=opt$subd_detection_depth)


save(images, file="f:/data_sholl_analysis/test/images/images.RData")

export_structure(MASTER, inter$main_vectors, "f:/data_sholl_analysis/test/dendrites/final_sectry.csv")

traced_MASTER <- lapply(1:inter$n_main_dendrites, trace_subdendrites, 
                        IMG=images$bin_noS_noMD_image,
                        MASTER=MASTER,
                        main_vectors_df=inter$main_vectors_df,
                        main_vectors_full=inter$main_vectors_full,
                        EPS=opt$trace_cluster_eps,
                        MPTS=opt$trace_cluster_mpt,
                        INC=opt$trace_detection_depth,
                        det_rad=opt$trace_detection_distance)




















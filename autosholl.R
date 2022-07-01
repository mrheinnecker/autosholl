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
# run_dir <- "f:/data_sholl_analysis/test/parabolar"
# load(file.path(run_dir, "images/images_new.RData"))


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

images$med_bin_noS_image_p <- apply_3d_median_filter(images$bin_noS_image_p, 2)
# export all images

MASTER <- lapply(1:inter$n_main_dendrites, find_subdendritic_starts_new, 
                 main_vectors=inter$main_vectors,
                 main_vectors_df=inter$main_vectors_df,
                 main_vectors_full=inter$main_vectors_full,
                 SOMA=SOMA,
                 IMG_screen=images$med_bin_noS_noMD_image,
                 IMG_adj=images$bin_noS_image_p,
                 soma_radius=inter$soma_radius,
                 MPTS=opt$subd_cluster_mpt,
                 EPS=opt$subd_cluster_eps,
                 det_rad=opt$subd_max_detection_distance,
                 z_range=opt$subd_detection_vertical_range,
                 nr_orig=opt$nr_orig
)






#save(images, file=file.path(run_dir, "images/images_new.RData"))






inter$all_vectors <- lapply(MASTER, function(MAIND){
  lapply(MAIND, get_all_vectors) %>%
  Reduce(function(x,y)append(x,y),.)
})
inter$all_vectors_fv <- lapply(inter$all_vectors, function(AV){
  print(1)
  lapply(AV, get_full_vector_voxels, IMG=images$bin_noS_noMD_image)
})
  
  
#AVF <- inter$all_vectors_fv[3]

inter$all_points <- lapply(inter$all_vectors_fv, function(AVF){
  #print(bind_rows(AVF))
  print(2)
  if(length(AVF)==0){return(NULL)}
    get_all_vector_speheres(bind_rows(AVF), opt$trace_dend_cross_radius_vertical, 
                                            opt$trace_dend_cross_radius_horizontal, 
                                            opt$nr_orig, 
                                      images$bin_noS_noMD_image)
})
  




traced_MASTER <- lapply(1:inter$n_main_dendrites, function(n){
  
     trace_subdendrites( 
                        IMG=images$bin_noS_noMD_image,
                        MAIND=MASTER[[n]],
                        MV=inter$main_vectors_df[[n]],
                        main_vectors_full=inter$main_vectors_full[[n]],
                        
                        all_vectors = inter$all_vectors[[n]], 
                        all_vectors_fv =inter$all_vectors_fv[[n]] %>% bind_rows(), 
                        all_points =inter$all_points[[n]],
                        
                        EPS_orig=opt$trace_cluster_eps,
                        MPTS_orig=opt$trace_cluster_mpt,
                        INC=opt$trace_detection_depth,
                        
                        RADv=opt$trace_dend_cross_radius_vertical,
                        RADh=opt$trace_dend_cross_radius_horizontal,
                        
                        det_rad=opt$trace_detection_distance,
                        RESC_DIST=opt$trace_rescore_dist,
                        RESC_ANG=opt$trace_rescore_angle,
                        nr_orig=opt$nr_orig)
  
  
})
                        
                     
## exporting
run_dir <- "f:/data_sholl_analysis/test/runs/010722_01"
dir.create(run_dir)
export_structure(traced_MASTER, inter$main_vectors,
                 file.path(run_dir, "dendrites_full_new.csv"))



export_structure(MASTER, inter$main_vectors, file.path(run_dir, "subdendrite_starts_final_rep.csv"))



# 
# writeTIFF(images$raw_image, file.path(run_dir, "images/raw.tif"))
# writeTIFF(images$noS_image, file.path(run_dir, "images/noS.tif"))
# writeTIFF(images$noS_noMD_image, file.path(run_dir, "images/noS_noMD.tif"))
# writeTIFF(images$bin_noS_image_p, file.path(run_dir, "images/bin_noS_partly_otsumod5.tif"))
 writeTIFF(images$bin_noS_noMD_image, file.path(run_dir, "bin_noS_noMD.tif"))
# writeTIFF(images$med_bin_noS_noMD_image, file.path(run_dir, "images/bin_med_noS_noMD.tif"))
# writeTIFF(images$med_bin_noS_image_p, file.path(run_dir, "images/bin_med_noS_partly.tif"))
# save(images, file=file.path(run_dir, "images/images_new.RData"))















traced_MASTER_old <- lapply(1:inter$n_main_dendrites, trace_subdendrites, 
                        IMG=images$bin_noS_noMD_image,
                        MASTER=MASTER,
                        main_vectors_df=inter$main_vectors_df,
                        main_vectors_full=inter$main_vectors_full,
                        
                        all_vectors = inter$all_vectors, 
                        all_vectors_fv =inter$all_vectors_fv, 
                        all_points =inter$all_points,
                        
                        EPS_orig=opt$trace_cluster_eps,
                        MPTS_orig=opt$trace_cluster_mpt,
                        INC=opt$trace_detection_depth,
                        det_rad=opt$trace_detection_distance,
                        RESC_DIST=opt$trace_rescore_dist,
                        RESC_ANG=opt$trace_rescore_angle,
                        nr_orig=opt$nr_orig)



export_structure(traced_MASTER, inter$main_vectors,
                 paste0("f:/data_sholl_analysis/test/parabolar/soa/dendrites_new_rem.csv"))



# 
# 
# traced_MASTER <- trace_subdendrites(nMD=2,  
#                                     IMG=images$bin_noS_noMD_image,
#                                     MASTER=MASTER,
#                                     main_vectors_df=inter$main_vectors_df,
#                                     main_vectors_full=inter$main_vectors_full,
#                                     EPS=opt$trace_cluster_eps,
#                                     MPTS=opt$trace_cluster_mpt,
#                                     INC=opt$trace_detection_depth,
#                                     det_rad=opt$trace_detection_distance,
#                                     RESC_DIST=opt$trace_rescore_dist,
#                                     RESC_ANG=opt$trace_rescore_angle,
#                                     nr_orig=opt$nr_orig)










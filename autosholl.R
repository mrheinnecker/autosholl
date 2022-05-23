### this is the main script of autosholl
## load functions and dependencies
source("/Users/Marco/git_repos/autosholl/new_fncts.R")
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


images$nosoma_image <- remove_soma(SOMA, inter$soma_radius, images$raw_image, opt$nr_orig, opt$nc_orig)
images$rem_mainobj_image <- remove_main_dendrites(images$nosoma_image, inter$main_vectors, 
                                                  inter$main_vectors_full, 30, opt$nr_orig, opt$nc_orig)


images$bin_nosoma_image <- binarize_image(opt, inter, 
                                          inter$n_main_dendrites, images$nosoma_image, 
                                          inter$elongated_dendrites,80, F)
images$bin_rem_mainobj_image <- binarize_image(opt, inter$n_main_dendrites, images$rem_mainobj_image, 80, T)


# writeTIFF(images$bin_nosoma_image, "f:/data_sholl_analysis/test/intermediate/screen_binary.tif")

images$bin_nosoma_image <- readTIFF("f:/data_sholl_analysis/test/intermediate/binarized.tif", all=T)


MASTER <- lapply(1:inter$n_main_dendrites, find_subdendritic_starts, 
                 main_vectors=inter$main_vectors,
                 main_vectors_df=inter$main_vectors_df,
                 main_vectors_full=inter$main_vectors_full,
                 SOMA=SOMA,
                 IMG=images$bin_nosoma_image,
                 soma_radius=inter$soma_radius,
                 minPTS=15,
                 EPS=3.5)

#export_structure(MASTER, inter$main_vectors, "f:/data_sholl_analysis/test/dendrites/cpp_try.csv")

traced_MASTER <- lapply(1:inter$n_main_dendrites, trace_subdendrites, 
                        MASTER=MASTER,
                        main_vectors_df=main_vectors_df)





### this is the main script of autosholl
## load functions and dependencies
source("/Users/Marco/git_repos/autosholl/fncts.R")
load_dependencies()
Rcpp::sourceCpp(file="c:/Users/Marco/git_repos/autosholl/cpp_fncts.cpp")




## load input files
files <- list.files("f:/data_sholl_analysis/single_soma", pattern = "*", full.names = T)
#files <- files[3:length(files)]
file <- files[1]
#file <- "f:/data_sholl_analysis/single_soma/Tier6_5_2_apical_D.tif"
main_dir <- "f:/data_sholl_analysis/run"
#dir.create(main_dir)
lapply(files, function(file){
  
  #print(file)
  
  run_dir <- file %>% str_split("/") %>% unlist() %>% last() %>% str_replace(".tif", "") %>%
    file.path(main_dir, .)
  
  dir.create(run_dir)
  
  cat(paste0("\n", file))
  
  main(file, run_dir)
  
  
})



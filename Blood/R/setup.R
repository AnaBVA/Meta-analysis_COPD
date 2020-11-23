library(stringr)

TODAY = format(Sys.time(),'%Y-%m-%d')

DATA <- function(x){
  str_c(file.path(here::here("Blood/data/")),paste0("Blood ",x,collapse = ""))
}

OUTPUT <- function(x){
  str_c(file.path(here::here("Blood/output_data/")),paste0("Blood ",x,collapse = ""))
}

DOWNLOAD <- function(x){
  str_c(file.path(here::here("Blood/download/")),paste0("Blood ",x,collapse = ""))
}

FIG <- function(x){
  str_c(file.path(here::here("Blood/fig/")),paste0("Blood ",x,collapse = ""))
}


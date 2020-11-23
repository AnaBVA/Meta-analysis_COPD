library(stringr)

TODAY = format(Sys.time(),'%Y-%m-%d')

DATA <- function(x){
  str_c(file.path(here::here("Blood/data/")),paste0(x,collapse = ""))
}

OUTPUT <- function(x){
  str_c(file.path(here::here("Blood/output_data/")),"Blood_",paste0(x,collapse = ""))
}

DOWNLOAD <- function(x){
  str_c(file.path(here::here("Blood/download/")),"Blood_",paste0(x,collapse = ""))
}

FIG <- function(x){
  str_c(file.path(here::here("Blood/fig/")),"Blood_",paste0(x,collapse = ""))
}


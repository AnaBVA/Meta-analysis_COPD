library(stringr)

TODAY = format(Sys.time(),'%Y-%m-%d')

DATA <- function(x){
  str_c(file.path(here::here("Tissue/data/")),paste0(x,collapse = ""))
}

OUTPUT <- function(x){
  str_c(file.path(here::here("Tissue/output_data/")),paste0(x,collapse = ""))
}

DOWNLOAD <- function(x){
  str_c(file.path(here::here("Tissue/download/")),paste0(x,collapse = ""))
}

FIG <- function(x){
  str_c(file.path(here::here("Tissue/fig/")),paste0(x,collapse = ""))
}


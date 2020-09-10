
TODAY = format(Sys.time(),'%Y-%m-%d')

DATA <- function(x){
  str_c(file.path(here::here("data/")),paste0(x,collapse = ""))
}

OUTPUT <- function(x){
  str_c(file.path(here::here("output_data/")),paste0(x,collapse = ""))
}

DOWNLOAD <- function(x){
  str_c(file.path(here::here("download/")),paste0(x,collapse = ""))
}

FIG <- function(x){
  str_c(file.path(here::here("fig/")),paste0(x,collapse = ""))
}


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

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

myBreaks <- function(svn_df, paletteLength = 30){
  dataZ <- t(apply(svn_df, 1, cal_z_score))
  dataZ <- na.omit(dataZ)
  #paletteLength = paletteLength - length(seq(min(dataZ),quantile(dataZ)[2],1)) - length(seq(quantile(dataZ)[4], max(dataZ),1))
  myBreaks <- c(seq(min(dataZ),quantile(dataZ)[2], length.out=ceiling(paletteLength/4)),
                seq(quantile(dataZ)[2], 0, length.out=ceiling(paletteLength/4)),
                seq(max(dataZ)/paletteLength, quantile(dataZ)[4], length.out=floor(paletteLength/4)+1),
                seq(quantile(dataZ)[4], max(dataZ), length.out=floor(paletteLength/4))
  )
  myBreaks <- myBreaks[!duplicated(myBreaks)]
  # myBreaks <- c(seq(min(dataZ), 0, length.out=ceiling(paletteLength/2) + 1), 
  #               seq(max(dataZ)/paletteLength, max(dataZ), length.out=floor(paletteLength/2))
  # )
  return(myBreaks)
}


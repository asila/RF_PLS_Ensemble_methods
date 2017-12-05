gg_spectra <- function(raw, wd, save = TRUE){

suppressWarnings(library(readr))
suppressWarnings(library(dplyr))
suppressWarnings(library(prospectr))
suppressWarnings(library(DT))
suppressWarnings(library(reshape2))
suppressWarnings(library(ggplot2))
suppressWarnings(library(soil.spec))

# set working directory with data
wd <- wd

#Read MIR spectra
raw <- raw

#Average
colnames(raw) <- c("SSN", colnames(raw[,-1]))

raw0 <- raw %>%

  group_by(SSN) %>%
  	
  	summarise_all(mean)

raw0 <- as.data.frame(raw0)

#get wavenumbers from column of raw0

wavenumbers<-as.numeric(substr(colnames(raw0[,-1]),2,19))

colnames(raw0) <- c("SSN", wavenumbers)

spec.m <- melt(raw0, id = "SSN")

p <- ggplot(data = spec.m, aes(x = as.numeric(as.vector(variable)),y = value,group = SSN)) +
    geom_line(size = 0.14,aes(col = "brown"), alpha = 0.3) +

      ggtitle("Raw MIR spectra") +
  
  xlim(rev(range(wavenumbers)))+
  
  #xlim(c(4002,610)) +
  
    ylim(range(spec.m$value)) + 
  
    xlab(expression("Wavenumbers cm"^-1)) +
  
    ylab("Absorbance") + 
    #theme with white background
    theme_bw() +
    #eliminates background, gridlines, and chart border
    theme(
        plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
    )
p <- p + theme(plot.title = element_text(hjust =0.5))

p <- p + theme(legend.position = "none")

return(p)

if (save == TRUE){
wd <- wd
ggsave(file = paste0(wd,"/raw_spectra_signatures.png"),p,height = 4, width = 7)
}

}

# Use gg_spectra illustration

raw <- read.csv("~/Downloads/MIR Technoserve Soil.csv")[,-c(2:6)]

wd <- "~/Downloads"

gg_spectra(raw, wd, save = TRUE)
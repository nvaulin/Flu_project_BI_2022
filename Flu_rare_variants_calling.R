rm(list=setdiff(ls(),c()))
gc()

packages <- c("tidyverse", 
              "data.table")

install.packages(setdiff(packages, rownames(installed.packages())))

library("tidyverse")
library("data.table")

workdir <- getwd()
datadir <-paste0(workdir,"/variant_calling")

filenames <-  "^variants_list_*"
files <- list.files(path = datadir, pattern = filenames)

df<-data.frame( position = c(1),
                prev_nucl = c(1),
                new_nucl = c(1),
                freq = c("1"),
                sample = c("1"),
                rep = c(1)
)[-1,]

for(file in files){
  filesplit <- strsplit(file ,split='_', fixed=TRUE)
  sample <-  filesplit[[1]][3]
  if(sample == "control"){
    rep <- filesplit[[1]][4]
    rep <- strsplit(rep ,split='.', fixed=TRUE)[[1]][1]
  }else{
    rep <- 0
    sample <- strsplit(sample ,split='.', fixed=TRUE)[[1]][1]
  }
  data <- read.table(paste(datadir, file, sep = "/"), sep = ' ',  
             header = FALSE, 
             stringsAsFactors = TRUE)
  data <- cbind(data, sample, rep)
  colnames(data) <- c("position", "prev_nucl", "new_nucl", "freq",
                      "sample", "rep")
  df <- rbind(df,data)
}

rm(list=c("data", "filesplit", "file", "filenames", "files", "rep", "sample"))

df$freq <- lapply(df$freq, str_replace_all, "%", "")

df$freq <- as.numeric(df$freq)
df$rep <-  as.numeric(df$rep)
df$sample <- as.factor(df$sample)

control_freq_borders <- df  %>% subset(sample == "control") %>%  
                                select(freq, sample, rep) %>% 
       aggregate(. ~ sample*rep, 
                 function(x) c(mean = mean(x), sd = sd(x), 
                               upper_border = mean(x) + 3*sd(x)))

freq_treshold <- max(control_freq_borders$freq[,"upper_border"])

person_variants <- df %>% subset(sample == "person") %>% 
                          select(position, prev_nucl, new_nucl, freq) %>% 
                          filter(freq > freq_treshold & freq < 90)

colnames(person_variants)[colnames(person_variants) == 'freq'] <- 'freq, %'

write.csv(person_variants, "person_variants_selected.csv")


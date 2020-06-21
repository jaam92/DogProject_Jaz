library(tidyverse)
df = read.table("~/Documents/DogProject_Jaz/DidNotUse_SamsData/FinalSamsData_autoOnly.fam", stringsAsFactors= F)
y = expand(df, V1, V2) %>% filter(V1 != V2) 
z = y[!duplicated(t(apply(y[c("V1", "V2")], 1, sort))), ] %>%
  mutate(V1 = paste0(V1, ".bed"),
         V2 = paste0(V2, ".bed"))

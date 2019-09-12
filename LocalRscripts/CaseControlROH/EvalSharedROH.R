library(data.table)
library(dplyr)

#Load Files
IW_pheno = read.delim("~/Documents/DogProject_Jaz/CaseControlROH/IrishWolfhounds_Epilepsy.txt")
ROH = fread("~/Documents/DogProject_Jaz/CaseControlROH/MergedFile_CornellCanineFitak_UnrelatedsOnly.hom.overlap", fill = T)

#separate ROH by case and control
IW_control = IW_pheno %>% filter(epilepsy_irishWolfhounds == 1)
IW_case = IW_pheno %>% filter(epilepsy_irishWolfhounds == 2)

IW_control_ROH = IW_control_ROH = ROH %>% filter(FID %in% IW_control$dogID) %>% group_by(POOL) %>% count()  %>% rename(countControl = n) %>% as.data.frame() 
IW_case_ROH = ROH %>% filter(FID %in% IW_case$dogID) %>% group_by(POOL) %>% count() %>% rename(countCase = n) %>% as.data.frame() 

#union ROH
ROH_consensus = ROH %>% filter(FID == "CON") %>% mutate(controlCount = IW_control_ROH$countControl[match(POOL, IW_control_ROH$POOL)], caseCount = IW_case_ROH$countCase[match(POOL, IW_case_ROH$POOL)]) %>% filter(controlCount >= 2 | caseCount >= 2)

ROH_union = ROH %>% filter(FID == "UNION") %>% mutate(controlCount = IW_control_ROH$countControl[match(POOL, IW_control_ROH$POOL)], caseCount = IW_case_ROH$countCase[match(POOL, IW_case_ROH$POOL)]) %>% filter(controlCount >= 2 | caseCount >= 2)

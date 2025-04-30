library(tidyverse)
library(magrittr)
library(here)


gbif_pow_exactjoin <- function(gbif_data = "Cacti_data_send_to_Santiago.csv",
                               checklist = "wcvp_names.csv",
                               taxon_status_cats = c("Accepted","Orthographic","Synonym","Unplaced"),
                               output_file="Cacti_step1.Rdata") { 
  aa <- data.table::fread(here("data",gbif_data),encoding = "Latin-1") %>% as_tibble()
  aa %<>% mutate(ID=1:nrow(.)) %>% relocate(ID,.before=1) %>% separate(scientificName,into=c("G","S","extra"),extra = "merge",sep=" ") %>% unite("binomial",c(G,S),sep="_",remove=F) %>% unite("scientificName",c(G,S,extra),sep=" ",remove=T) %>% filter(species!="")
  
  list_names <- data.table::fread(here("~/Documents/7.MAPOTECA/wcvp/",checklist),sep="|",quote="") %>% as_tibble() # %>% ## data provided by Kew (we are unable to provide it here)
  
  list_names %<>% unite("scientificName",c(taxon_name,taxon_authors),sep=" ",remove = F) %>% filter(species!="") %>% mutate(Accepted_Name = scientificName[match(accepted_plant_name_id,plant_name_id,nomatch=NA)]) %>% as_tibble() %>%  filter(taxon_status %in% all_of(taxon_status_cats)) %>% unite("binomial",c(genus,species),sep="_",remove = F) %>% filter(accepted_plant_name_id!="")
  
## EXACT JOIN 
  joined <- aa %>% #rename("taxon_name"=scientificName) %>% 
    left_join(.,list_names,by = "scientificName")
  cat("Searching for duplicate IDs due to multiple matches!","\n")
  joined %>% pull(ID) %>% duplicated() -> dups
  joined$ID[which(dups)] -> id_dups
  cat("Found",length(id_dups),"duplicated IDs","\n")
  if(length(id_dups) != 0){ 
  joined %>% filter(ID %in% all_of(id_dups)) -> to_correct
  joined <- joined[which(!joined$ID %in% id_dups),]
  pow_distributions="wcvp_distribution.txt"
  pow_dist <- data.table::fread(here("wcvp_2022",pow_distributions),sep="|",quote="") %>% as_tibble() %>% #left_join(.,list_names,by="accepted_db_id") %>% 
    filter(introduced==0)
  joined %<>% mutate(Dist_correction=T)
  to_correct %<>% mutate(Dist_correction=T)
  cat("First step: resolving by distribution","\n")
  for(i in 1:nrow(to_correct)){
    cat("processing",i,"\r")
    pow_dist %>% filter(plant_name_id==to_correct$accepted_plant_name_id[i]) %>% pull(continent_code_l1) -> tar_dist
    tar_dist %in% c(1:6,9) -> res
    if(sum(res)!=0) {to_correct$Dist_correction[i] <- FALSE}
  }
  to_correct %<>% filter(Dist_correction)
  rbind(to_correct,joined) -> joined
  joined %>% pull(ID) %>% duplicated() -> dups
  joined$ID[which(dups)] -> id_dups
  cat("Remaining",length(id_dups),"duplicated IDs","\n")
  joined %>% filter(ID %in% all_of(id_dups)) -> to_correct
  joined <- joined[which(!joined$ID %in% id_dups),]
  to_correct %>% filter(homotypic_synonym=="") %>% group_split(ID) -> to_correct
  cat("Second step: resolve manually","\n")
  for(i in 1:length(to_correct)){
    cat("processing",i,"\n")
    if(nrow(to_correct[[i]])==1) next
    to_correct[[i]] %>% select(scientificName,Accepted_Name,homotypic_synonym,accepted_plant_name_id,taxon_status,geographic_area) %>% print()
    to_select <- as.numeric(readline("Which one to keep?"))
    if(is.na(to_select)) next
    to_correct[[i]] %<>% slice(-all_of(to_select))
  }
  data.table::rbindlist(to_correct) %>% as_tibble() -> to_correct
  rbind(to_correct,joined) -> joined
  joined %>% pull(ID) %>% duplicated() -> dups
  joined$ID[which(dups)] -> id_dups
  cat("Still",length(id_dups),"some unsolved records!","\n")
  } 
  joined %>% filter(is.na(plant_name_id)) %>% distinct(scientificName) %>% nrow() -> spp
  cat(spp,"unresolved names left!","\n")
  save(joined,file=here("interim/",output_file))
}
gbif_pow_first_fuzzy <- function(data="Cacti_step1.Rdata",
                                 checklist = "wcvp_names.csv",
                                 taxon_status_cats = c("Accepted","Orthographic","Synonym","Unplaced"),
                                 output_file="Cacti_step2.Rdata"){ 
  load(here("interim",data))
  joined %>% filter(is.na(plant_name_id)) %>% distinct(scientificName) -> spp
  spp %>% filter(grepl("×",scientificName)) -> hybrids
  spp %>% filter(grepl("\\?",scientificName)) -> noidea
  spp %>% filter(!grepl("×",scientificName),!grepl("\\?",scientificName)) -> species
  
  list_names <- data.table::fread(here("~/Documents/7.MAPOTECA/wcvp/",checklist),sep="|",quote="") %>% as_tibble() # %>% ## data provided by Kew (we are unable to provide it here)
  
  list_names %<>% unite("scientificName",c(taxon_name,taxon_authors),sep=" ",remove = F)
  list_names %<>% filter(species!="") %>% mutate(Accepted_Name = scientificName[match(accepted_plant_name_id,plant_name_id,nomatch=NA)]) %>% as_tibble() %>% filter(taxon_status %in% all_of(taxon_status_cats)) %>% unite("binomial",c(genus,species),sep="_",remove = F) %>% filter(accepted_plant_name_id!="")
  
  joined_nana <- joined %>% filter(is.na(plant_name_id))
  joined %<>% filter(!is.na(plant_name_id))
  
  resolved_matches <- list()
  non_matches <- list()
  multiple_matches <- list()
  max.dist=0.15
  for (i in 1:nrow(species)){
    species %>% slice(i) %>% pull(scientificName) -> quien
    cat(i,"Attempting match on:", quien,"\n")
    joined_nana %>% filter(scientificName %in% all_of(quien)) -> target
    target %>% slice(1) %>% select(binomial.x,scientificName,accepted_plant_name_id,Accepted_Name) %>% rename(binomial=binomial.x) -> sole_target
    
    list_names %>% filter(binomial== sole_target$binomial) -> names_tomatch
    fuzzyjoin::stringdist_left_join(species %>% slice(i),names_tomatch,by="scientificName",method="jw",distance_col="distance",max_dist=max.dist) -> result
    if(nrow(result)>1){ 
      if(length(unique(result$accepted_plant_name_id))==1){
        cat("        Single match!","\n")
        result[1,] -> transfer
        target$accepted_plant_name_id <- transfer$accepted_plant_name_id
        target$Accepted_Name <- transfer$Accepted_Name
        resolved_matches[[i]] <- target
      } 
      if(length(unique(result$accepted_plant_name_id))>1) {
        if(min(result$distance)<0.01){
          cat("        Single match!","\n")
          result[which.min(result$distance),] -> transfer
          target$accepted_plant_name_id <- transfer$accepted_plant_name_id
          target$Accepted_Name <- transfer$Accepted_Name
          resolved_matches[[i]] <- target
        }
        if(min(result$distance)>0.01){
          cat("        Multiple matches!","\n")
          target$accepted_plant_name_id <- NA
          target$Accepted_Name <- NA
          multiple_matches[[i]] <- target
        }
      }
      next
    }
    if(is.na(result$plant_name_id)){ cat("       No match!","\n")
      target$accepted_plant_name_id <- NA
      target$Accepted_Name <- NA
      non_matches[[i]] <- target
      next}
    if(!is.na(result$plant_name_id)){
      cat("        Single match!","\n")
      result %>% select(scientificName.x,accepted_plant_name_id,Accepted_Name) -> transfer
      target$accepted_plant_name_id <- transfer$accepted_plant_name_id
      target$Accepted_Name <- transfer$Accepted_Name
      resolved_matches[[i]] <- target
    }
  }
  
  resolved_matches %>% data.table::rbindlist(.) %>% as_tibble() -> joined_fuzzy
  joined_fuzzy %<>% bind_rows(joined,.)
  joined_fuzzy %>% distinct(Accepted_Name)
  
  save(joined_fuzzy,file=here("interim",output_file))
  save(non_matches,file=here("interim","non_matches_fuzzy.Rdata"))
  save(multiple_matches,file=here("interim","multiples_matches_fuzzy.Rdata"))
  
}

gbif_pow_second_fuzzy <- function(output_file="Cacti_final.v1.Rdata"){ 
  load(here("interim/non_matches_fuzzy.Rdata"))
  load("interim/Cacti_step2.Rdata")
  non_matches %<>% data.table::rbindlist() %>% as_tibble() 
  non_matches %>% distinct(scientificName) -> species
  resolved_matches <- list()
  non_matches_final <- list()
  multiple_matches_final <- list()
  max.dist=0.15
  for (i in 1:nrow(species)){
    species %>% slice(i) %>% pull(scientificName) -> quien
    cat(i,"Attempting match on:", quien,"\n")
    non_matches %>% filter(scientificName %in% all_of(quien)) -> target
    target %>% slice(1) %>% select(binomial.x,scientificName,accepted_plant_name_id,Accepted_Name,genus.x) %>% rename(binomial=binomial.x,genus=genus.x) -> sole_target
    
    list_names %>% filter(genus == sole_target$genus) -> names_tomatch
    fuzzyjoin::stringdist_left_join(species %>% slice(i),names_tomatch,by="scientificName",method="jw",distance_col="distance",max_dist=max.dist) -> result
    if(nrow(result)>1){ 
      if(length(unique(result$accepted_plant_name_id))==1){
        cat("        Single match!","\n")
        result[1,] -> transfer
        target$accepted_plant_name_id <- transfer$accepted_plant_name_id
        target$Accepted_Name <- transfer$Accepted_Name
        resolved_matches[[i]] <- target
      } 
      if(length(unique(result$accepted_plant_name_id))>1) {
        if(min(result$distance)<0.01){
          cat("        Single match!","\n")
          result[which.min(result$distance),] -> transfer
          target$accepted_plant_name_id <- transfer$accepted_plant_name_id
          target$Accepted_Name <- transfer$Accepted_Name
          resolved_matches[[i]] <- target
        }
        if(min(result$distance)>0.01){
          cat("        Multiple matches!","\n")
          target$accepted_plant_name_id <- NA
          target$Accepted_Name <- NA
          multiple_matches_final[[i]] <- target
        }
      }
      next
    }
    if(is.na(result$plant_name_id)){ cat("       No match!","\n")
      target$accepted_plant_name_id <- NA
      target$Accepted_Name <- NA
      non_matches_final[[i]] <- target
      next}
    if(!is.na(result$plant_name_id)){
      cat("        Single match!","\n")
      result %>% select(scientificName.x,accepted_plant_name_id,Accepted_Name) -> transfer
      target$accepted_plant_name_id <- transfer$accepted_plant_name_id
      target$Accepted_Name <- transfer$Accepted_Name
      resolved_matches[[i]] <- target
    }
  }
  non_matches_final %>% data.table::rbindlist(.) %>% as_tibble() 
  multiple_matches %>% data.table::rbindlist(.) %>% as_tibble() %>% distinct(scientificName)
  
  resolved_matches %>% data.table::rbindlist(.) %>% as_tibble() -> joined_fuzzy_f
  joined_fuzzy %<>% bind_rows(.,joined_fuzzy_f)
  joined_fuzzy %>% distinct(Accepted_Name)
  save(joined_fuzzy,file=here("output",output_file))
  save(non_matches_final,file="output/no_matches_final.v1.Rdata")
  save(multiple_matches_final,file="output/multiples_matches_final.v1.Rdata")
  
}
geographic_filter <- function(data="Cacti_final.v1.Rdata",
                              output_file = "Cacti_finalPOW.v1.Rdata",
                              perform_tests=c("centroids","institutions", "equal", "gbif","capitals", "zeros","seas")) {
  load(here("output",data))
  ### APPLY FILTERS. NOT USING THE OUTLIER TEST, BECAUSE WE BASE THIS ON KEW'S DISTRIBUTIONS.
  joined_fuzzy %<>% filter(!is.na(decimalLongitude))
  joined_fuzzy %<>% filter(!is.na(decimalLatitude))
  
  cat("Filtering using CoordinateCleaner","\r")
  to_filter <- joined_fuzzy %>% select(decimalLongitude,decimalLatitude) %>% as.data.frame() 
  f <- CoordinateCleaner::clean_coordinates(to_filter, lon = "decimalLongitude",lat = "decimalLatitude",centroids_rad = 1000, centroids_detail = "both",tests=perform_tests,species=NULL,value="flagged",seas_scale = 110)
  joined_fuzzy %<>% mutate(Outlier_Test = f)
  save(joined_fuzzy,file=here("output",output_file))
}
powdist_filter <- function(data="Cacti_finalPOW.v1.Rdata",
                           output_file = "Cacti_finalPOWdist.v1.Rdata",
                           pow_distributions="wcvp_distribution.csv",
                           wgsrpd = "level3/level3.shp") { 
  load(here("output",data))
  pow_dist <- data.table::fread(here("~/Documents/7.MAPOTECA/wcvp/",pow_distributions),sep="|",quote="") %>% as_tibble() %>% #left_join(.,list_names,by="accepted_db_id") %>% 
    filter(introduced==0)
  cat("....... reading polygons from WGSRPD","\r")
  poly <- terra::vect(here("~/Documents/7.MAPOTECA/wgsrpd-master/",wgsrpd))
  points <- terra::vect(joined_fuzzy, geom=c("decimalLongitude", "decimalLatitude"), crs=poly, keepgeom=TRUE)
  points_powo <- terra::extract(poly,points)
  joined_fuzzy %<>% bind_cols(.,points_powo)
  cat("Filtering using KEW's distributions","\n")
  joined_fuzzy %>% group_split(Accepted_Name) -> aver
  some <- list()
  for (i in 1:length(aver)){
    cat(i,"--","Getting POW data for:", unique(aver[[i]]$Accepted_Name),"               ","\r")
    temp <- aver[[i]] %>% mutate(POW_distribution = FALSE)
    if(is.na(unique(aver[[i]]$Accepted_Name))) {
      some[[i]] <- temp %>% mutate(POW_distribution = TRUE)
    }
    matching <- !is.na(match(temp$LEVEL3_COD,pow_dist %>% filter(plant_name_id==unique(temp$accepted_plant_name_id)) %>% pull(area_code_l3),nomatch = NA))
    if(length(which(matching)) == 0) { 
      some[[i]] <- temp %>% mutate(POW_distribution = TRUE)
      next
    }
    if(length(which(matching)) != 0) some[[i]] <- temp %>% mutate(POW_distribution = matching)
  }
  some %<>% data.table::rbindlist(.) %>% as_tibble()
  cat("-------------- DONE --------------","\r")
  save(some,file=here("output",output_file))
}
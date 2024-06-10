# YOU WILL NEED AN ALLELE TABLE (IMPORTED FROM AN allele_data.txt) AND MASKED REFERENCES
# THE ALLELE TABLE IS CALLED allele.data IN THIS CODE. CHANGE IF NECESSARY

# COLUMN NAMES IN ALLELE TABLE HERE ARE SampleID, Locus, ASV, Reads, Allele, PseudoCIGAR
# RENAME IF NECESSARY

# THE OUTPUT IS allele.data.fixed 
# CHECK THE LAST LINES IN THE CODE THAT COLLAPSE READS FROM ASV WITH SAME PSEUDOCIGAR 
# AND REMOVES ALLELE NAME (AS IT IS NOT CONSISTENT ACROSS RUNS)

# CHANGE PATH TO MASKED REFERENCES 
masked_reference = read.delim("G:/My Drive/EPPIcenter/Data/SensitivitySpecificity/masking/masked_references_pseudocigar.txt") %>% 
  select(-ASV_ref_homopol,-ASV_ref_tr)

fix_PseudoCIGAR <- function(alleletable){
  
  parse_PseudoCIGAR <- function(PseudoCIGAR) {
    str_extract_all(PseudoCIGAR, "\\d+\\+\\d+N|\\d+I=[A-Z]+|\\d+D=[A-Z]+|\\d+[A-Z]")[[1]]
  }
  
  
  replace_elements <- function(pseudo_list, ref_list, idx_list) {
    asv_list <- rep("", length(pseudo_list))
    asv_list[idx_list] <- ref_list
    return(asv_list)
  }
  
  
  pseudocigar.tofix = alleletable %>% 
    ungroup() %>% 
    select(Locus,PseudoCIGAR) %>% 
    distinct() %>% 
    mutate(MaskPseudoCIGAR = map_chr(PseudoCIGAR, 
                                     ~ str_extract_all(.x, "[0-9]+\\+[0-9]+N") %>% 
                                       unlist() %>% 
                                       paste(collapse = ""))) %>% 
    mutate(MaskPseudoCIGAR = ifelse(MaskPseudoCIGAR=="",".",MaskPseudoCIGAR)) %>% 
    left_join(masked_reference %>%
                select(-ASV_ref) %>% 
                rename(RefPseudoCIGAR = PseudoCIGAR),
              by = "Locus") %>% 
    distinct() %>% 
    filter(MaskPseudoCIGAR != RefPseudoCIGAR)%>% 
    rowwise() %>% 
    mutate(PseudoCIGAR_list = list(parse_PseudoCIGAR(PseudoCIGAR)),
           RefPseudoCIGAR_list = list(parse_PseudoCIGAR(RefPseudoCIGAR)),
           PseudoCIGAR_list_indexes = list(as.numeric(str_remove_all(PseudoCIGAR_list,"I=|D=|[ATCG]+|\\+\\d+N"))),
           PseudoCIGAR_list_N = list(map_dbl(as.numeric(str_remove_all(PseudoCIGAR_list,"\\d+\\+|\\d+I=[ATGC]+|\\d+[ATGC]|\\d+D=[ATGC]+|N")), ~ replace_na(.x, 0))),
           PseudoCIGAR_list_Nidx = list(which(PseudoCIGAR_list_N!=0)),
           RefPseudoCIGAR_asASV = list(replace_elements(PseudoCIGAR_list, RefPseudoCIGAR_list, PseudoCIGAR_list_Nidx)),
           RefPseudoCIGAR_asASV_N = list(map_dbl(as.numeric(str_remove_all(RefPseudoCIGAR_asASV,"\\d+\\+|N")), ~ replace_na(.x, 0))),
           Ndiff = list(PseudoCIGAR_list_N - RefPseudoCIGAR_asASV_N),
           idxdiff = list(lag(cumsum(unlist(Ndiff)),default = 0)),
           PseudoCIGAR_list_N = list(as.character(PseudoCIGAR_list_N - Ndiff)),
           PseudoCIGAR_list_N = list(map_chr(PseudoCIGAR_list_N, ~ ifelse(.x == "0", "", .x))),
           PseudoCIGAR_list_indexes = list(PseudoCIGAR_list_indexes - idxdiff),
           PseudoCIGAR_list_other = list(str_remove_all(PseudoCIGAR_list,"\\d+|N")),
           PseudoCIGAR_list_whereN = list(map_chr(str_extract(PseudoCIGAR_list,"N"),~ replace_na(.x, ""))),
           newPseudoCIGAR = list(map2_chr(PseudoCIGAR_list_indexes,PseudoCIGAR_list_other,~paste0(.x,.y))),
           newPseudoCIGAR = list(map2_chr(newPseudoCIGAR,PseudoCIGAR_list_N,~paste0(.x,.y))),
           newPseudoCIGAR = list(map2_chr(newPseudoCIGAR,PseudoCIGAR_list_whereN,~paste0(.x,.y))),
           newPseudoCIGAR = paste(newPseudoCIGAR,collapse="")
    )%>% 
    select(Locus,RefPseudoCIGAR,PseudoCIGAR,newPseudoCIGAR)
  
  alleletable.fixed = alleletable %>% 
      left_join(pseudocigar.tofix %>% 
                  select(Locus,PseudoCIGAR,newPseudoCIGAR),
                by = c("Locus","PseudoCIGAR")) %>% 
    mutate(PseudoCIGAR = ifelse(is.na(newPseudoCIGAR),PseudoCIGAR,newPseudoCIGAR)) %>% 
    select(-newPseudoCIGAR) 
  
  return(alleletable.fixed)
}



# THIS IS YOUR NEW TABLE WITH PSEUDO CIGARS FIXED
allele.data.fixed = fix_PseudoCIGAR(allele.data)

rm(list = c("masked_reference","fix_PseudoCIGAR"))

# RECOMMENDED CODE TO COLLAPSE ALL READS WITH SAME PSEUDOCIGAR AND DIFFERENT ASV:
allele.data.fixed = allele.data.fixed %>% 
  group_by(SampleID,Locus,PseudoCIGAR) %>% 
  summarise(Reads=sum(Reads))%>% 
  ungroup()


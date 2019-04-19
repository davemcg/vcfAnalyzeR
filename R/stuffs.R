vcf_import <- function(vcf){
  vcf <- vcfR::read.vcfR(vcf)
  #merge <- cbind(vcf@fix, vcf@gt) %>% as_tibble()
  vcf
}



filler <- function(info, regex = 'phyloP_100way=', filler_value = ';phyloP_100way=-1'){
  info[!grepl(regex, info)] <- paste0(info[!grepl(regex, info)], filler_value)
  info
}

appender <- function(tibble, info_field, numeric = TRUE){
  info <- filler(tibble$INFO, regex = paste0(info_field,"\\="),
                 filler_value = paste0(";", info_field, "=-1"))
  if (numeric){
    tibble[,info_field] <- str_split(info, ';') %>%
      unlist() %>%
      grep(paste0('^',info_field, '\\='),., value =T) %>%
      gsub(paste0(info_field,'='),'',.) %>%
      as.numeric()
  } else {
    tibble[,info_field] <- str_split(info, ';') %>%
      unlist() %>%
      grep(paste0('^', info_field, '\\='),., value =T) %>%
      gsub(paste0(info_field,'='),'',.)
  }

  tibble
}

vcf <- vcf_import('~/Desktop/EGAD00001002656_2019.GATK.ANNO.novel_exons.vcf.gz')
merge <- cbind(vcf@fix, vcf@gt) %>% as_tibble()

info2 <- filler(merge$INFO)
merge$phyloP_100way <- str_split(info2, ';') %>% unlist() %>% grep('^phyloP_100way\\=',., value =T) %>% gsub('phyloP_100way=','',.) %>% as.numeric()

merge$AC <- str_split(info2, ';') %>% unlist() %>% grep('^AC\\=',., value =T) %>% gsub('AC=','',.) %>% as.numeric()

info2 <- filler(merge$INFO, regex = "gno_af_all\\=", filler_value = ";gno_af_all=-1")
merge$gno_af_all <- str_split(info2, ';') %>% unlist() %>% grep('^gno_af_all\\=',., value =T) %>% gsub('gno_af_all=','',.) %>% as.numeric()


info2 <- filler(merge$INFO, regex = "gno_hom\\=", filler_value = ";gno_hom=-1")
merge$gno_hom <- str_split(info2, ';') %>% unlist() %>% grep('^gno_hom\\=',., value =T) %>% gsub('gno_hom=','',.) %>% as.numeric()

info2 <- filler(merge$INFO, regex = "gno_ac_popmax\\=", filler_value = ";gno_ac_popmax=-1")
merge[,'gno_ac_popmax'] <- str_split(info2, ';') %>% unlist() %>% grep('^gno_ac_popmax\\=',., value =T) %>% gsub('gno_ac_popmax=','',.) %>% as.numeric()


samples <- grep('^[B|W|G]\\d', colnames(merge), value =T)


HET_COUNT <- numeric(nrow(merge))
HOM_COUNT <- numeric(nrow(merge))
for (i in seq(1,nrow(merge))){
  HET_COUNT[i] <- ((merge[i,samples] %>% map(~ str_extract(.x, pattern = '^[\\d|\\.][\\/|\\|][\\d\\.]')) %>% unlist()) == "0/1") %>% sum()
  HOM_COUNT[i] <- ((merge[i,samples] %>% map(~ str_extract(.x, pattern = '^[\\d|\\.][\\/|\\|][\\d\\.]')) %>% unlist()) == "1/1") %>% sum()
}
merge$HET_COUNT <- HET_COUNT
merge$HOM_COUNT <- HOM_COUNT

merge <- appender(merge, 'CSQ', numeric = FALSE)
merge <- merge %>% rowwise() %>% mutate(Gene = str_split(CSQ, pattern = '\\|')[[1]][6])


merge %>% filter(HOM_COUNT > 0, gno_hom < 1, gno_af_all < 0.001, HOM_COUNT < 10, HET_COUNT < 10) %>% select(-contains('B2'), -contains('G0'), -contains('W0'), -CSQ, -INFO, -FILTER, -FORMAT) %>% arrange(-phyloP_100way)

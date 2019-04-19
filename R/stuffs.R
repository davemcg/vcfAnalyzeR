vcf_import <- function(vcf){
  vcf <- vcfR::read.vcfR(vcf)
  #merge <- cbind(vcf@fix, vcf@gt) %>% as_tibble()
  vcf
}



filler <- function(info, regex = 'phyloP_100way=', filler_value = ';phyloP_100way=-1'){
  info[!grepl(regex, info)] <- paste0(info[!grepl(regex, info)], filler_value)
  info
}

appender <- function(tibble, regex, filler, new_col_name){

  tibble[,new_col_name] <- str_split(info2, ';') %>% unlist() %>% grep('^phyloP_100way\\=',., value =T) %>% gsub('phyloP_100way=','',.) %>% as.numeric()
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

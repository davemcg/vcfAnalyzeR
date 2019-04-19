vcf_import <- function(vcf){
  vcf <- vcfR::read.vcfR(vcf)
  #merge <- cbind(vcf@fix, vcf@gt) %>% as_tibble()
  vcf
}



filler <- function(info, regex = 'phyloP_100way=', filler_value = ';phyloP_100way=-1'){
  info[!grepl(regex, info)] <- paste0(info[!grepl(regex, info)], filler_value)
  info
}

vcf <- vcf_import('~/Desktop/EGAD00001002656_2019.GATK.ANNO.novel_exons.vcf.gz')
merge <- cbind(vcf@fix, vcf@gt) %>% as_tibble()
merge$phyloP_100way <- str_split(info2, ';') %>% unlist() %>% grep('^phyloP_100way\\=',., value =T) %>% gsub('phyloP_100way=','',.) %>% as.numeric()

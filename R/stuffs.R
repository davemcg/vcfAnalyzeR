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


retnet <- read_tsv('~/git/variant_prioritization/data/retnet_hgncIDs_2017-03-28.txt', col_names = FALSE)


un_solved <- read_delim('~/git/EGA_EGAD00001002656_NGS_reanalyze/EGAD00001002656_2.ped', delim = ' ')
samples <- grep('^[B|W|G]\\d', colnames(merge), value =T)
unsolved_samples <- samples[samples %in% (un_solved %>% filter(Status == 'Unsolved') %>% pull(Patient))]
partially_samples <- samples[samples %in% (un_solved %>% filter(Status == 'PartiallySolved') %>% pull(Patient))]
solved_samples <- samples[!samples %in% c(unsolved_samples, partially_samples)]

HET_COUNT_SOLVED <- numeric(nrow(merge))
HOM_COUNT_SOLVED <- numeric(nrow(merge))
HET_COUNT_PSOLVED <- numeric(nrow(merge))
HOM_COUNT_PSOLVED <- numeric(nrow(merge))
HET_COUNT_UNSOLVED <- numeric(nrow(merge))
HOM_COUNT_UNSOLVED <- numeric(nrow(merge))
HET_SAMPLES <-
for (i in seq(1,nrow(merge))){
  HET_COUNT_SOLVED[i] <- ((merge[i,solved_samples] %>% map(~ str_extract(.x, pattern = '^[\\d|\\.][\\/|\\|][\\d\\.]')) %>% unlist()) == "0/1") %>% sum()
  HET_COUNT_PSOLVED[i] <- ((merge[i,partially_samples] %>% map(~ str_extract(.x, pattern = '^[\\d|\\.][\\/|\\|][\\d\\.]')) %>% unlist()) == "0/1") %>% sum()
  HET_COUNT_UNSOLVED[i] <- ((merge[i,unsolved_samples] %>% map(~ str_extract(.x, pattern = '^[\\d|\\.][\\/|\\|][\\d\\.]')) %>% unlist()) == "0/1") %>% sum()
  HOM_COUNT_SOLVED[i] <- ((merge[i,solved_samples] %>% map(~ str_extract(.x, pattern = '^[\\d|\\.][\\/|\\|][\\d\\.]')) %>% unlist()) == "1/1") %>% sum()
  HOM_COUNT_PSOLVED[i] <- ((merge[i,partially_samples] %>% map(~ str_extract(.x, pattern = '^[\\d|\\.][\\/|\\|][\\d\\.]')) %>% unlist()) == "1/1") %>% sum()
  HOM_COUNT_UNSOLVED[i] <- ((merge[i,unsolved_samples] %>% map(~ str_extract(.x, pattern = '^[\\d|\\.][\\/|\\|][\\d\\.]')) %>% unlist()) == "1/1") %>% sum()
}
merge$HET_COUNT_SOLVED <- HET_COUNT_SOLVED
merge$HET_COUNT_PSOLVED <- HET_COUNT_PSOLVED
merge$HET_COUNT_UNSOLVED <- HET_COUNT_UNSOLVED
merge$HOM_COUNT_SOLVED <- HOM_COUNT_SOLVED
merge$HOM_COUNT_PSOLVED <- HOM_COUNT_PSOLVED
merge$HOM_COUNT_UNSOLVED <- HOM_COUNT_UNSOLVED

merge <- appender(merge, 'CSQ', numeric = FALSE)
merge <- merge %>% rowwise() %>% mutate(Gene = str_split(CSQ, pattern = '\\|')[[1]][6])

# Try to ID hom alt variants
merge %>%
  dplyr::filter(HOM_COUNT_UNSOLVED > 0, gno_hom < 1, gno_af_all < 0.001, HOM_COUNT_UNSOLVED < 20, HOM_COUNT_SOLVED < 1, HET_COUNT_SOLVED < 5, HET_COUNT_UNSOLVED < 5) %>%
  mutate(RetNet = case_when(Gene %in% retnet$X1 ~ 'Yes',
                            TRUE ~ 'No'),
         region = case_when(grepl('UTR', CSQ) ~ 'UTR',
                            grepl('intron', CSQ) ~ 'Intron',
                            TRUE ~ 'Other')) %>%
  select(-contains('B2'), -contains('G0'), -contains('W0'), -contains('HOME'), -CSQ, -INFO, -FILTER, -FORMAT) %>%
  arrange(-phyloP_100way) %>%
  unique() %>%
  data.frame()


# Try to ID het variants
merge %>% dplyr::filter(gno_af_all == -1, HET_COUNT_SOLVED < 5, HOM_COUNT_SOLVED < 1) %>%
  mutate(RetNet = case_when(Gene %in% retnet$X1 ~ 'Yes',
                            TRUE ~ 'No'),
         Delta = HET_COUNT_UNSOLVED - HET_COUNT_SOLVED,
         region = case_when(grepl('UTR', CSQ) ~ 'UTR',
                            grepl('intron', CSQ) ~ 'Intron',
                            TRUE ~ 'Other')) %>%
  filter(Delta > 0) %>%
  select(-contains('B2'), -contains('G0'), -contains('W0'), -contains('HOME'), -CSQ, -INFO, -FILTER, -FORMAT) %>%
  arrange(-Delta) %>% unique() %>%
  mutate(Sample_Het = names(ID(CHROM, POS, REF, ALT, "0/1")) %>% substr(., 1, 7) %>% na.omit() %>% unique() %>% paste(., collapse = ', ')) %>%
  left_join(un_solved, by = c('Sample_Het' = 'Patient')) %>%
  data.frame()

# hom_samples
ID <- function(chrom, pos, ref, alt, genotype){
  ((merge[,c('CHROM', 'POS', 'REF', 'ALT', unsolved_samples)] %>% filter(CHROM==chrom, POS==pos, REF==ref, ALT==alt) %>% select(-POS) %>%  map(~ str_extract(.x, pattern = '^[\\d|\\.][\\/|\\|][\\d\\.]')) %>% unlist()) == genotype)[((merge[,c('CHROM','POS', 'REF', 'ALT', unsolved_samples)] %>% filter(CHROM==chrom, POS==pos, REF==ref, ALT==alt) %>% select(-POS) %>%  map(~ str_extract(.x, pattern = '^[\\d|\\.][\\/|\\|][\\d\\.]')) %>% unlist()) == genotype) == TRUE]
}


# interesting hom variants, the sample, and the phenotype
merge %>%
  dplyr::filter(HOM_COUNT_UNSOLVED > 0, gno_hom < 1, gno_af_all < 0.001, HOM_COUNT_UNSOLVED < 20, HOM_COUNT_SOLVED < 1, HET_COUNT_SOLVED < 5, HET_COUNT_UNSOLVED < 5) %>%
  mutate(RetNet = case_when(Gene %in% retnet$X1 ~ 'Yes',
                            TRUE ~ 'No'),
         region = case_when(grepl('UTR', CSQ) ~ 'UTR',
                            grepl('intron', CSQ) ~ 'Intron',
                            TRUE ~ 'Other')) %>%
  select(-contains('B2'), -contains('G0'), -contains('W0'), -contains('HOME'), -CSQ, -INFO, -FILTER, -FORMAT) %>%
  arrange(-phyloP_100way) %>%
  unique() %>% mutate(Sample_Hom = names(ID(CHROM, POS, REF, ALT, "1/1")) %>% substr(., 1, 7) %>% na.omit() %>% unique() %>% paste(., collapse = ', ')) %>%
  data.frame() %>%
  left_join(un_solved, by = c('Sample_Hom' = 'Patient'))

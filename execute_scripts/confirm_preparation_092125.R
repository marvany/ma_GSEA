library(data.table)


## actual input TW coding
main.dir = "/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/Resources/TW/laptop_to_minerva/TW_mod_processCols/coding" # for non-bulk tissues
p.bulk = "/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/Resources/TW/laptop_to_minerva/TW_mod_bulk_processCols/coding"


# prepare dt1.main
dt1.main = rbindlist(lapply(list.files(main.dir, recursive = TRUE, full.names = TRUE, pattern = '.tsv'), fread))
head(dt1.main[, .(tissue)])
setnames(dt1.main, 'p', 'pvalue')
names(dt1.main)


# prepare dt1.bulk
dt1.bulk = rbindlist(lapply(list.files(p.bulk, recursive = TRUE, full.names = TRUE, pattern = '.csv'), fread))
names(dt1.bulk)
dt1.bulk[, tissue := 'Bulk']
setnames(dt1.bulk, 'TRAIT', 'gwas')
setnames(dt1.bulk, 'effect_size', 'beta')
names(dt1.bulk)


# merge dt1.bulk and dt1.main
to_keep = names(dt1.main)[names(dt1.main) %in% names(dt1.bulk)]
dt1 = rbind (dt1.main[, ..to_keep], dt1.bulk[, ..to_keep])
names(dt1)


## original
dt2.main = fread('/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/Resources/TW/coding/psychad_TWFile.csv', sep = ',')
dt2.bulk = fread('/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/Resources/TW/coding/bulk_TW.tsv', sep = '\t')


# prepare dt2.bulk
setnames(dt2.bulk, 'TRAIT', 'gwas')
dt2.bulk[, tissue := 'Bulk']
names(dt2.bulk)


# prepare dt2.main
names(dt2.main)
setnames(dt2.main, c('TRAIT'), c('gwas'))
dt2.main[, tissue := NULL]
dt2.main[, tissue := TIM]
dt2.main[, TIM := NULL]


# merge dt2 main and bulk
to_keep = names(dt2.main)[names(dt2.main) %in% names(dt2.bulk)]
dt2.main = dt2.main[, ..to_keep]
dt2.bulk = dt2.bulk[, ..to_keep]
dt2 = rbind(dt2.main, dt2.bulk)

# merge dt1 and dt2
names(dt2)
names(dt1)
to_keep2 = names(dt1)[names(dt1) %in% names(dt2)]

dt1 = dt1[, ..to_keep2]
dt2 = dt2[, ..to_keep2]

names(dt1)
names(dt2)

unique(dt1$gwas) %in% unique(dt2$gwas)
unique(dt1$tissue) %in% unique(dt2$tissue)
all(unique(dt1$gene) %in% unique(dt2$gene))

dt = merge(dt1, dt2, by = c('gwas', 'tissue', 'gene'), suffix = c('_dt1', '_dt2'))
cor(dt[, zscore_dt1], dt[, zscore_dt2])
all(dt[, zscore_dt1] == dt[, zscore_dt2])
all(dt[, zscore_dt1] == dt[, zscore_dt2])
all(dt[, pvalue_dt1] == dt[, pvalue_dt2])

head(dt)


##########  TW2 ----------------------------------------------
### ACTUAL INPUT
main.dir = "/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/Resources/TW/laptop_to_minerva/TW_mod_processCols/coding" # since we are focusing on p-value it doesnt matter if choose original TW_mod
p.bulk = "/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/Resources/TW/laptop_to_minerva/TW2_mod_bulk_processedCols/coding/bulk"

dt3.main = rbindlist(lapply(list.files(main.dir, recursive = TRUE, full.names = TRUE, pattern = '.tsv'), function(x){
   fread(x)[, gwas := basename(dirname(dirname(x)))]
   })) 
head(dt3.main[, gwas])

dt3.bulk = rbindlist(lapply(list.files(p.bulk, recursive = TRUE, full.names = TRUE, pattern = '.tsv'), function(x){
   fread(x)[, gwas := basename(dirname(dirname(x)))]
   })) 
head(dt3.bulk[, gwas])


names(dt3.main)
setnames(dt3.main, 'p', 'pvalue')
head(dt3.main[, .(tissue)])
names(dt3.bulk)
head(dt3.bulk[, .(TIM, tissue)])
dt3.bulk[, TIM := NULL]

tokeep = names(dt3.bulk)[names(dt3.bulk) %in% names(dt3.main)]

dt3.bulk = dt3.bulk[, ..tokeep]
dt3.main = dt3.main[, ..tokeep]

dt3 = rbind(dt3.main, dt3.bulk)
names(dt3)
head(dt3[, .(tissue, gwas)])

dt3[, gwas := gsub('_TW_coding*', '', gwas)]

### ORIGINAL
dt4.bulk.main = fread('/sc/arion/projects/va-biobank/PROJECTS/ma_GSEA/Resources/TW/tw2_coding/TW_Bulk_PsychAD.tsv')

names(dt4.bulk.main)
head(dt4.bulk.main[, .(TIM, TRAIT, tissue)])
setnames(dt4.bulk.main, c('TRAIT'), c('gwas'))
dt4.bulk.main[, tissue := NULL]
dt4.bulk.main[, tissue := TIM]
dt4.bulk.main[, TIM := NULL]
names(dt4.bulk.main)
dt4 = dt4.bulk.main


all(names(dt3), names(dt4))
dt = merge(dt3, dt4, by = c('gwas', 'tissue', 'gene'), suffix = c('_dt1', '_dt2'))
cor(dt[, zscore_dt1], dt[, zscore_dt2])
all(dt[, zscore_dt1] == dt[, zscore_dt2])
all(dt[, zscore_dt1] == dt[, zscore_dt2])
all(dt[, pvalue_dt1] == dt[, pvalue_dt2])

library('edgeR')
library('plyr')

expression <-read.table(file = '/Users/pwangel/Downloads/pluripotent_rnaseq_data.tsv', sep = '\t', header = TRUE, row.names=1)
annotations <-read.delim(file = '/Users/pwangel/Downloads/pluripotent_rnaseq_metadata.tsv', sep = '\t', header = TRUE, row.names=1, fill=TRUE)
ensembl_to_symbol <- read.table(file= '/Users/pwangel/Data/ensembl_hg38.91/gene_to_symbol_ensembl91_human.tsv', sep='\t', header=FALSE, col.names = c('ensembl', 'symbol'))

state_group <- factor(annotations$LM_Group_COLOR)
dataset <- factor(annotations$dataset)

design <- model.matrix(~0+state_group+dataset)
colnames(design) <- gsub("state_group", "", colnames(design))
colnames(design) <- gsub("dataset", "", colnames(design))

contr.matrix <- makeContrasts(
   NaivevsPrimed = naive-primed, 
   levels = colnames(design))

y <- DGEList(counts=as.matrix(expression), group=as.matrix(state_group))
y <- calcNormFactors(y, method = "TMM")
y$samples$norm.factors

v <- voom(y, design)

fit <- lmFit(v, design)
fit <- contrasts.fit(fit, contrasts=contr.matrix)
efit <- eBayes(fit)

topTable(efit, n=500)
write.csv(merge(topTable(efit,n=500),ensembl_to_symbol,by.x=0,by.y="ensembl", all.x=TRUE, sort=FALSE), file = sprintf("/Users/pwangel/Downloads/de_primed_vs_naive.tsv"),sep='\t')

#design <- model.matrix(~sc_state_group)
#y <- estimateDisp(y, design)
#y <- estimateDisp(y)

#fit <- glmQLFit(y, design)
#contrast = c(1,-1)
#qlf <- glmQLFTest(fit, contrast=contrast)
#merge(topTags(qlf,n=500),ensembl_to_symbol,by.x=0,by.y="ensembl", all.x=TRUE, sort=FALSE)
#write.csv(merge(topTags(qlf,n=500),ensembl_to_symbol,by.x=0,by.y="ensembl", all.x=TRUE, sort=FALSE), file = sprintf("/Users/pwangel/Downloads/de_primed_vs_naive.tsv"),sep='\t')

#et <- exactTest(y)
#merge(topTags(et,n=500),ensembl_to_symbol,by.x=0,by.y="ensembl", all.x=TRUE, sort=FALSE)
#write.csv(merge(topTags(et,n=500),ensembl_to_symbol,by.x=0,by.y="ensembl", all.x=TRUE, sort=FALSE), file = sprintf("/Users/pwangel/Downloads/de_primed_vs_naive.tsv"),sep='\t')

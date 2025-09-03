gsea.path = "/sc/arion/projects/va-biobank/marios/Clinical.Significance/final_revision/fisher_results_V3"

multi3_wrapper(
  twas.list = ready_for_Fisher_analysis_TWAS_list[3], #for NO_focus   #####all.TWAS.forGSEA.fisher #for FOCUS,
  backgrounds = list(NULL), #for NO_focus ####list(Imputable.Genes.NO.BULK, NULL) #for FOCUS,   # please do not name the elements included in this list
  gsea.path = gsea.path, # all_combinations_TWAS_GSEA_fisher_V3 uses protein specific background for Bulk and has improved aesthetics
  validation = Validation.Gene.List,
  defined.genesets =  4, # NULL for generalised run, 
  tissues_to_examine = c("/Bulk", "/sn-Pseudohomogenate", "/All EUR Class", "/All EUR Subclass", "/All EUR Class + Subclass"),
  show.title = F
)

# vale to neo hotfix sto gsea_orplot function
# alla3e to multi.enrich.dfs gia na kaneis ta modifications kalytera


# Sanan, please start here.
#multi.enrich.dfs <- fread("/sc/arion/projects/va-biobank/marios/Clinical.Significance/final_revision/multi.enrich.dfs.csv")

View(multi.enrich.dfs)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
saveplot = T
title.of.plot = "TWAS_Fisher_input_p_0.01_PROTEIN_GreatGenes"
gsea.outdir.or = "/sc/arion/projects/va-biobank/marios/Clinical.Significance/final_revision/fisher_results_V3/TWAS_Fisher_input_p_0.01_PROTEIN_GreatGenes/ORplots/all.traits.1.graph"
show.title = F
plot.width = 8.5
plot.height = 8
returnplot = T


multi.enrich.dfs <- multi.enrich.dfs[!grepl("All EUR Class \\+ Subclass", multi.enrich.dfs$name_full),]

multi.enrich.dfs$name_full[multi.enrich.dfs$name_full == "All EUR Class"] <- "Class-Aggregate"
multi.enrich.dfs$name_full[multi.enrich.dfs$name_full == "All EUR Subclass"] <- "Subclass-Aggregate"
multi.enrich.dfs$name_full[multi.enrich.dfs$name_full == "Bulk"] <- "Bulk" 
multi.enrich.dfs$name_full[multi.enrich.dfs$name_full == "sn-Pseudohomogenate"] <- "Pseudobulk Bulk" 
multi.enrich.dfs$name_full <- factor(multi.enrich.dfs$name_full,levels=c("Subclass-Aggregate","Class-Aggregate","Pseudobulk Bulk", "Bulk"))

multi.enrich.dfs$legend <- factor(multi.enrich.dfs$name_full,levels=c("Subclass-Aggregate","Class-Aggregate","Pseudobulk Bulk", "Bulk"))
multi.enrich.dfs$Set <- unlist(lapply(multi.enrich.dfs$Set, function(x) strsplit(x,split="_")[[1]][1]))
multi.enrich.dfs$Set[multi.enrich.dfs$Set == "Alcoholism"] <- "AUD"
multi.enrich.dfs$Set <- factor(multi.enrich.dfs$Set,levels=rev(c("BD","SCZ","AD","MS","AUD","ALS","PD","Anorexia","Migraines","MDD","ADHD","Insomnia")))


# this is used in scale_color_manual in values = ... This is to more easily manipulate the colors and order in your plot.
colors_to_include <- rep(c(cbbPalette, cbPalette[1]), 2)[1:(length(unique(multi.enrich.dfs$name_full))) + 1]



#multi.enrich.dfs <- df

all.traits.p.OR.CADsig <-  ggplot( # reference line is a dashed line corresponds to OR = 1 (by geom_hline)                                                      
  
  multi.enrich.dfs,          # text labels (by geom_text) will be added in each point for the stars
  aes(                      # there will be a legend to explain the colors and shapes, at the top right of the plot.
    x=Set,       # there will be multiple panels (by facet_grid) 
    y=odds.ratio, # the plot will be using a classic theme (theme_classic) with a base font size of 10
    label=star)) +
  geom_hline(aes(yintercept = 1), size = .25, linetype ="dashed") +
  geom_errorbar(aes(ymax = OR.ci95.up, ymin = OR.ci95.down, group = name_full, color = name_full), size = 0.5,
                width = 0.25, position = position_dodge(width = 0.65)) +
  
  geom_point(aes(
    # color = Set,
    group = name_full, 
    # shape = Set),
    size = intersection), 
    position = position_dodge(width = 0.65)) +
  scale_size_continuous(range = c(0.5,2))+
  #geom_point(aes(color = Set), size = 2.5, position = position_dodge(width = 0.65)) +
  geom_text(
    aes(group = name_full, y = OR.ci95.up), # to align the labels
    vjust = 0.7,
    hjust = -0.5,
    position = position_dodge(width = 0.65)) +
  # The following is optional to paint the error bars according to the pallete
  scale_color_manual( # color in the middle of the errorbar (modifies ge)
    name = "Set", # name of the legend
    values = colors_to_include, # edw praktika enwnei thn cbbPalette me to prwto element tou cbPalette kai kanei repeat 2 fores to neo palette
    labels = rev(levels(multi.enrich.dfs$name_full)),                                              # sth synexeia xrhsimopoiei tous deiktes gia na vrei ta xrwmata pou tha valei
    breaks = rev(levels(multi.enrich.dfs$name_full))) +
  scale_shape_manual( # shape in middle errorbar
    name = "Set",  # define the name in the legend
    values=c(rep(16, 9), rep(17,9))[1:length(unique(multi.enrich.dfs$Set))],  # defines the various shapes
    labels = unique(multi.enrich.dfs$Set),  # specifies the labels in the legend
    breaks = unique(multi.enrich.dfs$Set)) + # 
  theme_classic(base_size = 10) +
  xlab("")+
  ylab(expression(paste("Enrichment odds ratio for pathways"))) +
  labs(size = "Intersection"
  )+
  # scale_x_discrete(limits = levels(multi.enrich.dfs$name_full)) +
  scale_y_continuous(trans='log10',
                     expand = c(0.2, 0)) +
  theme(legend.position = c(0.18, 0.82),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) + # c(0.77, 0.27)) +
  guides(color=guide_legend(title = paste0("TWAS"),
                            ncol=1, byrow=F, keywidth=0.1, keyheight=0.1,
                            default.unit="inch"
  ),
  shape=guide_legend(title = paste0("TWAS"),
                     ncol=1, byrow=F, keywidth=0.1, keyheight=0.1,
                     default.unit="inch")) +
  coord_flip()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

if(length(unique(multi.enrich.dfs$group)) > 1) all.traits.p.OR.CADsig <- all.traits.p.OR.CADsig + facet_wrap(~group) # this should be included if you are doing comparisons accross validation genesets
if(show.title) all.traits.p.OR.CADsig <- all.traits.p.OR.CADsig + labs(title = title.of.plot)
#remove bracket!!
#   list.files(getwd())
#list.files(paste0(getwd(), "/Scripts/GSEA_ORplot_Functions.R"))
if (saveplot == T) {
  if(!is.null(title.of.plot)) {
    file.name <- gsub(" ", "_", title.of.plot)
    ggsave(paste0(gsea.outdir.or, "/", file.name, ".pdf"), all.traits.p.OR.CADsig, width = plot.width, height = plot.height, useDingbats=FALSE)
    ggsave(paste0(gsea.outdir.or, "/", file.name, ".png"), all.traits.p.OR.CADsig, width = plot.width, height = plot.height)
    # all.plots.path <- gsub(basename(gsea.path), "", gsea.path)
    if(!is.null(all.plots.path)) {
      ggsave(paste0(all.plots.path, "/", title.of.plot, ".pdf"), all.traits.p.OR.CADsig, width = plot.width, height = plot.height, useDingbats = FALSE)
      ggsave(paste0(all.plots.path, "/", title.of.plot, ".png"), all.traits.p.OR.CADsig, width = plot.width, height = plot.height)
      }
  }else{
    ggsave(paste0(gsea.outdir.or, "/top", top, "_Sets_OR.pdf"), all.traits.p.OR.CADsig, width = plot.width, height = plot.height, useDingbats=FALSE)
    ggsave(paste0(gsea.outdir.or, "/top", top, "_Sets_OR.png"), all.traits.p.OR.CADsig, width = plot.width, height = plot.height)
  }
}
if (returnplot == T) return(all.traits.p.OR.CADsig)        #This excludes the subsequent line, fix this with a list

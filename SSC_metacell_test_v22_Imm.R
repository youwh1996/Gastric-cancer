library(foreach)
library(pheatmap)
library(metacell)
library(tgconfig)
source('~/Repos/scTools/scTools_meta.R')
setwd('~/Project_SSC_shared/')

set_param("scm_spike_regexp","^ERCC-","metacell")
set_param("scm_mc_mark_k_per_clust",100,"metacell") #default: 5
set_param("scm_mc_mark_min_gene_cov",0.3,"metacell") # default: 0.25
set_param("scm_mc_mark_min_gene_fold",2,"metacell") # default: 1.5

set_param("mcell_mc2d_K",30,"metacell") # default: 20
set_param("mcell_mc2d_T_edge",0.02,"metacell") # default: 0.05
set_param("mcell_mc2d_max_confu_deg",4,"metacell") # default: 5
set_param("mcell_mc2d_edge_asym",FALSE,"metacell") # default: TRUE
set_param("mcell_mc2d_proj_blur",0.02,"metacell") # default: 0.02

pjname <- 'SSC_test_v22_Imm'
dbpath <- 'data_test_v22_Imm'
fgpath <- 'figs_test_v22_Imm'
mat_id <- 'SSC_test_v22_Imm'
mtpath <- 'metadata_SSC_test_v22.txt'

unlink(dbpath, recursive=TRUE)
if(!dir.exists(dbpath)) dir.create(dbpath)
scdb_init(dbpath, force_reinit=T)
mcell_import_multi_mars(mat_id,dataset_table_fn = mtpath, base_dir = "scdata/")
mat <- scdb_mat(mat_id)
print(dim(mat@mat))
if(!dir.exists(fgpath)) dir.create(fgpath)
scfigs_init(fgpath)
mcell_plot_umis_per_cell(mat_id)
mat <- scdb_mat(mat_id)

nms <- unique(c(rownames(mat@mat), rownames(mat@ignore_gmat)))
pre_nr_term <- c("^ERCC-","^RPS","^RPL","^MT-","^MTMR","^MTND","^MTRN","^MTCO","^MRPL","^MRPS","^HBA","^HBB","^MTATP")
pre_nr_genes <- foreach(i=pre_nr_term, .combine = c) %do% grep(i, nms, v=T)
pre_ex_genes <- c("MALAT1", "XIST", "XIST_intron")
pre_bad_genes <- unique(c(pre_nr_genes, pre_ex_genes))
pre_bad_genes

mcell_mat_ignore_genes(new_mat_id=mat_id, mat_id=mat_id, pre_bad_genes, reverse=F)
load('test_v18_Imm.Rda')
load('SB382_Imm.Rda')
load('SB383_Imm.Rda')
load('SB385_Imm.Rda')
load('SB388_Imm.Rda')
load('SB391_Imm.Rda')
load('SB392_Imm.Rda')
load('SB400_Imm.Rda')
good_cells <- c(test_v18_Imm, SB382_Imm, SB383_Imm, SB385_Imm, SB388_Imm, SB391_Imm, SB392_Imm, SB400_Imm)

bad_cells <- setdiff(colnames(mat@mat), good_cells)
length(good_cells);length(bad_cells)
mcell_mat_ignore_cells(new_mat_id = mat_id, mat_id = mat_id, ig_cells = bad_cells)
mat <- scdb_mat(mat_id)
mcell_mat_ignore_small_cells2(new_mat_id=mat_id, mat_id=mat_id, 250, 100)

genes_anchors = c('FOS','FOSB','NFKBIA','NFKBIZ','JUN','ZFP36','ISG15','HMGB2','STMN1','TOP2A','MKI67','MX1','RSAD2')
tab_fn = paste(fgpath, "lateral_gmods.txt", sep="/")
mcell_mat_rpt_cor_anchors(mat_id=mat_id, gene_anchors = genes_anchors, cor_thresh = 0.1,
                          gene_anti = c(), tab_fn = tab_fn, sz_cor_thresh = 0.2)
gcor_mat = read.table(file.path(fgpath, 'lateral_gmods.txt'), header=T)
foc_genes = apply(gcor_mat[, genes_anchors], 1, which.max)
foc_genes

mcell_add_gene_stat(gstat_id="test", mat_id=mat_id, force=T)
mcell_gset_filter_varmean(gset_id="test_feats", gstat_id="test", T_vm=0.08, force_new=T)
mcell_gset_filter_cov(gset_id = "test_feats", gstat_id="test", T_tot=100, T_top3=2)
mcell_plot_gstats(gset_id="test_feats", gstat_id="test")
gset <- scdb_gset("test_feats")

pst_genes <- names(gset@gene_set)
pst_nr_term <- c("^AC[0-9]+\\.", "^AP[0-9]+\\.", "^AL[0-9]+\\.", "-AS[0-9]*$", "^MIR[0-9]", "^LINC[0-9]", "^SNOR", "^sno", "^SNR",
                 "^TMSB", "^HIST", "^HSP", "^IFI", "^HLA-", "^ATP", "-", ":", "\\.", '^KIAA',
                 "^IGJ", "^IGH", "^IGK", "^IGL", "^DNAJ", "^GZM", "^CCL", "^XCL", '^FTH', '^FTL', '^LGALS')
pst_nr_genes <- foreach(i=pst_nr_term, .combine = c) %do% grep(i, pst_genes, v=T)
pst_ex_genes <- c()
pst_bad_genes <- unique(c(pst_nr_genes, pst_ex_genes, names(foc_genes)))
pst_add_genes <- c()
final_genes <- unique(setdiff(pst_genes, pst_bad_genes), pst_add_genes)
final_genes
genes_to_gset("test_feats",final_genes,"test_feats")
gset <- scdb_gset("test_feats")
length(gset@gene_set)

mcell_add_cgraph_from_mat_bknn(
  mat_id=mat_id,
  gset_id = "test_feats",
  graph_id="test_graph",
  K=240,
  dsamp=T)
mcell_coclust_from_graph_resamp(
  coc_id="test_coc500",
  graph_id="test_graph",
  min_mc_size=30,
  p_resamp=0.75, n_resamp=500)
mcell_mc_from_coclust_balanced(
  coc_id="test_coc500",
  mat_id= mat_id,
  mc_id= "test_mc",
  K=45, min_mc_size=45, alpha=2)
mc<-scdb_mc("test_mc")
pct<-cal_mc_pct(mat,mc)

expressed_genes <- rownames(mc@mc_fp)
expressed_genes_avg_umi <- rowSums(mat@mat[expressed_genes,])/length(mc@mc)

meta <- read.table(mtpath, header=T, row.names = 1, sep="\t")
meta <- meta[intersect(rownames(meta), rownames(mc@n_bc[rowSums(mc@n_bc) > 0,])),]
meta <- meta[order(meta$Patient.ID),]
Total <- rowSums(mc@n_bc[rownames(meta),])
write.table(cbind(meta,mc@n_bc[rownames(meta),],Total), file.path(fgpath,paste0(pjname,"_nbc.csv")), sep=",", col.names=NA)
write.table(round(cbind(mc@mc_fp, expressed_genes_avg_umi),3), file.path(fgpath,paste0(pjname,"_mean.csv")), sep=",", col.names=NA)
means<-as.matrix(log(mc@mc_fp[names(gset@gene_set),],2))
png(file.path(fgpath,paste0(pjname,"_mean_cor.png")),height=20,width=20,units="in",res=300)
pheatmap(cor(means),border_color = NA)
dev.off()

mcell_gset_from_mc_markers(gset_id="test_markers", mc_id="test_mc")
mcell_mc_plot_marks(mc_id="test_mc", gset_id="test_markers", mat_id=mat_id)
mcell_mc2d_force_knn(mc2d_id="test_mc_2dproj",mc_id="test_mc", graph_id="test_graph")
mcell_mc2d_plot(mc2d_id="test_mc_2dproj")
mc2d <- scdb_mc2d("test_mc_2dproj")

mc_hc <- mcell_mc_hclust_confu(mc_id="test_mc", graph_id="test_graph")
mc_sup <- mcell_mc_hierarchy(mc_id="test_mc", mc_hc=mc_hc, T_gap=0.04)
mcell_mc_plot_hierarchy(mc_id="test_mc", 
                        graph_id="test_graph", 
                        mc_order=mc_hc$order, 
                        sup_mc = mc_sup, 
                        width=3000, heigh=2000, 
                        min_nmc=2, show_mc_ids = T)

# # mcell_plot_outlier_heatmap(mc_id="test_mc", mat_id = mat_id, T_lfc=3)
# mcell_mc_split_filt(
#   new_mc_id="test_mc_f",
#   mc_id="test_mc",
#   mat_id=mat_id,
#   T_lfc=4, plot_mats=F)
# mc_f<-scdb_mc("test_mc_f")
# 
# mcell_gset_from_mc_markers(gset_id="test_markers_f", mc_id="test_mc_f")
# mcell_mc_plot_marks(mc_id="test_mc_f", gset_id="test_markers_f", mat_id=mat_id)
# 
# # marks_colors = read.table(system.file("extdata", "pbmc_mc_colorize.txt", package="metacell"), sep="\t", h=T, stringsAsFactors=F)
# # mc_colorize("test_mc_f", marker_colors=marks_colors)
# mcell_mc2d_force_knn(mc2d_id="test_mc_f_2dproj",mc_id="test_mc_f", graph_id="test_graph")
# mcell_mc2d_plot(mc2d_id="test_mc_f_2dproj")

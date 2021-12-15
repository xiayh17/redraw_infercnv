library(stringr)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
# library(phylogram)
# library(dendextend) 


## 构建表型信息表 ----
if (!file.exists("phe.Rdata")) {
  load('recluster_sce.Rdata')
  colnames(sce@meta.data)
  table(sce$seurat_clusters)
  phe_recluster_sce = sce@meta.data
  
  load('epi_sce.Rdata')
  colnames(sce@meta.data)
  table(sce$seurat_clusters)
  phe_raw_sce = sce@meta.data
  
  load('phe_annoCNV.Rdata') 
  head(phe)
  phe_with_cnv = phe
  
  load('phe-by-markers.Rdata')
  head(phe)
  head(phe_with_cnv) 
  identical(rownames(phe_with_cnv),
            rownames(phe_raw_sce))
  identical(rownames(phe_recluster_sce),
            rownames(phe_raw_sce))
  
  celltype = phe[match( rownames(phe_with_cnv) ,rownames(phe)),
                 'celltype']
  table( celltype )
  
  
  phe=data.frame(
    bar=rownames( phe_with_cnv ),
    sample=phe_with_cnv$orig.ident,
    celltype=celltype,
    raw_cluster=phe_with_cnv$seurat_clusters,
    new_clutser = phe_recluster_sce$seurat_clusters,
    cnv_scores= log10(phe_with_cnv$cnv_scores)
  )
  head(phe)
  table(phe$sample) 
  table(phe$celltype )
  ## 保存初步整理好的表型信息表
  save(file = "phe.Rdata",phe)
}

load("phe.Rdata")

# load_infercnv_obj -------------------------------------------------------

infercnv_obj = readRDS('infercnv_output/preliminary.infercnv_obj')

# 提取cnv表达矩阵 
expr <- infercnv_obj@expr.data
# dim(expr)
# expr[1:4,1:4]
message(glue::glue("There are {ncol(expr)} cells and {nrow(expr)} genes in infercnv_obj"))

## only keep Observations(Cells) epi 
cells_expr <- as.vector(attr(expr, "dimnames")[[2]])
my_cell_indices <- match(phe$bar,cells_expr)

plot_data <- expr[,my_cell_indices]

message(glue::glue("There are {ncol(plot_data)} cells and {nrow(plot_data)} genes after filter with phe"))
# scale_data --------------------------------------------------------------
# mean
x.center=mean(expr)
# examine distribution of data that's off-center, since much of the center could
# correspond to a mass of data that has been wiped out during noise reduction
quantiles = quantile(expr[expr != x.center], c(0.01, 0.99))

# determine max distance from the center.
delta = max( abs( c(x.center - quantiles[1],  quantiles[2] - x.center) ) )
low_threshold = x.center - delta
high_threshold = x.center + delta

plot_data[plot_data < low_threshold] <- low_threshold
plot_data[plot_data > high_threshold] <- high_threshold

# chr_data ----------------------------------------------------------------
chr = infercnv_obj@gene_order$chr
table(infercnv_obj@gene_order$chr)

chr_d = data.frame(
  
  chr = infercnv_obj@gene_order$chr,
  start = infercnv_obj@gene_order$start,
  stop = infercnv_obj@gene_order$stop,
  gene = as.vector(attr(expr, "dimnames")[[1]])
  
)

chr_d$chr <- gsub("chr","",chr_d$chr)
chr_d$chr <- factor(chr_d$chr,levels = stringr::str_sort(unique(chr_d$chr),numeric = T))

# 按细胞类型分组聚类 ---------------------------------------------------------------
hclust_method="ward.D"
threads = ceiling(parallel::detectCores()-1) ## 聚类使用的核心数
# 来自infercnv包的代码
# method that returns an artificial dendrogram with only 1 branch
# since at least 2 elements need to be present for building one normally
# but we need to be able to merge this lone branch to the bigger dendrogram for plotting

.single_element_dendrogram <- function(unique_label) {
  dfake = list()
  dfake[[1]] = 1L
  attr(dfake, "class") = "dendrogram"
  attr(dfake, "height") = 1
  attr(dfake, "members") = 1L
  attr(dfake[[1]], "class") = "dendrogram"
  attr(dfake[[1]], "leaf") = TRUE
  attr(dfake[[1]], "label") = unique_label
  attr(dfake[[1]], "height") = 0
  return(dfake)
}


.pairwise_dendrogram <- function(labels) {
  dfake = list()
  dfake[[1]] = 1L
  dfake[[2]] = 1L
  attr(dfake, "class") = "dendrogram"
  attr(dfake, "height") = 1
  attr(dfake, "members") = 2L
  attr(dfake[[1]], "class") = "dendrogram"
  attr(dfake[[1]], "leaf") = TRUE
  attr(dfake[[1]], "label") = labels[1]
  attr(dfake[[1]], "height") = 0
  attr(dfake[[2]], "class") = "dendrogram"
  attr(dfake[[2]], "leaf") = TRUE
  attr(dfake[[2]], "label") = labels[2]
  attr(dfake[[2]], "height") = 0
  return(dfake)
}

## 按细胞类型对细胞分组 (需要细胞类型的列名) ----
cell_l <- split(phe$bar,phe$celltype)

# 对表达矩阵分组
plot_data_l <- lapply(seq_along(cell_l), function(x){
  
  plot_data[,cell_l[[x]]]
  # print(length(cells))
  
})

names(plot_data_l) <- names(cell_l)

# 配色代码 --------------------------------------------------------------------

## 染色体的颜色框
# 来自函数：infercnv::plot_cnv
get_group_color_palette <- function () {
  return(colorRampPalette(RColorBrewer::brewer.pal(12, "Set3")))
}
color <- get_group_color_palette()(length(unique(chr_d$chr)))
top_color <- HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = color), # 设置填充色
                       labels = levels(chr_d$chr), 
                       labels_gp = gpar(cex = 0.9, col = "black"))) 

## 表型的颜色(按表型信息的实际情况更改) ----
sample_color <- RColorBrewer::brewer.pal(8,name = "Set2")[1:length(unique(phe$sample))]
names(sample_color) <- unique(phe$sample)

cell_color <- RColorBrewer::brewer.pal(8,name = "Set2")[8-c(1:length(unique(phe$celltype))-1)]
names(cell_color) <- unique(phe$celltype)

## cluster 类型数量多的配色 ----
cluster_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
raw_cluster_color <- cluster_color_hue(length(unique(phe$raw_cluster)))
names(raw_cluster_color) <- unique(phe$raw_cluster)
new_cluster_color <- cluster_color_hue(length(unique(phe$new_clutser)))
names(new_cluster_color) <- unique(phe$new_clutser)

## 连续型变量表型配色 CNVscore (需要按实际情况更改) ----
col_score = colorRamp2(c(1, 2, 3, 4), c('#ffffff','#bae4b3','#74c476','#238b45'))

## 表达量的颜色
col_fun = colorRampPalette(c("darkblue", "white", "darkred"))(15)

## 重新计算聚类 较为耗时 可以考虑增大核心数
hc_l <- lapply(seq_along(plot_data_l), function(x){
  
  # 只对三个细胞及以上的进行聚类, 其他的直接转换成dendrogram
  if(any(is.matrix(plot_data_l[[x]]) & ncol(plot_data_l[[x]]) >= 3)){
    hclust(parallelDist::parallelDist(t(plot_data_l[[x]]),threads = threads),method = hclust_method)
  } else if (any(is.matrix(plot_data_l[[x]]) & ncol(plot_data_l[[x]]) <= 2)) {
    .pairwise_dendrogram(colnames(plot_data_l[[x]]))
  } else {
    key = names(plot_data_l)[[x]]
    .single_element_dendrogram(phe[grep(key,phe[,grep(key,phe)]),]$bar)
  }
  
})

## hclust 转换成dendrogram
dend_l <- lapply(seq_along(hc_l), function(x){
  
  if(class(hc_l[[x]]) == "hclust"){
    as.dendrogram(hc_l[[x]])
  } else {
    hc_l[[x]]
  }
  
})

## 合并分组的聚类，必须用do.call
if (length(dend_l) >=2) {
  
  rc_o <- do.call(merge,dend_l)
  
} else {
  
  # 当只有一种聚类时
  rc_o <- dend_l[[1]]
  
}

## 根据聚类模块数判断水平分割情况
if (length(hc_l)>=2) {
  split = 2 # 分隔聚类
} else {
  split = NULL 
}

# order by cluster labels
# 把表型按聚类排序
# 提取所有聚类的标签
labels_l <-lapply(dend_l, function(x){
  
  labels(x)
  
})

## 按聚类结果对表型信息排序 （需要bar code的列名） ----
phe_o <- phe[match(unlist(labels_l),phe$bar),]

## 根据聚类结果重排序的表型信息构建表型注释信息(需要根据实际情况更改) ----
cell_anno <- rowAnnotation(
  sample = phe_o$sample,
  celltype = phe_o$celltype,
  raw_cluster = phe_o$raw_cluster,
  new_cluster = phe_o$new_clutser,
  log10CNVscore = phe_o$cnv_scores,
  col = list(sample = sample_color,
             celltype =cell_color,
             raw_cluster = raw_cluster_color,
             new_cluster = new_cluster_color,
             log10CNVscore = col_score
  ))


# 画图代码 --------------------------------------------------------------------

# raster 的问题
#https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#heatmap-as-raster-image
png_units = "in"
png_res = 300
png_width =20
png_height = 10
png(filename = paste0("new_v2_reclust——test.png"),
    height = png_height,width = png_width,res = png_res,units = png_units)
ht_tumor = Heatmap(t(plot_data[,match(phe_o$bar,colnames(plot_data))]),
                   cluster_rows = rc_o, # do.call 产生的
                   split = split, # 聚类分隔
                   row_dend_width = unit(4, "cm"),
                   cluster_columns = F,
                   show_column_names = F,
                   show_row_names = F,
                   col = col_fun,
                   row_dend_reorder = F,
                   column_split = chr_d$chr,
                   left_annotation = cell_anno,
                   heatmap_legend_param = list(
                     title = "Modified Expression",
                     title_position = "leftcenter-rot", # 图例标题位置
                     at=c(0.4,1.6), #图例范围
                     legend_height = unit(3, "cm") #图例长度
                   ),
                   top_annotation = top_color,
                   use_raster = T, 
                   raster_by_magick = T,
                   raster_resize_mat= F,
                   raster_magick_filter = "Point",
                   row_title = "Observations (Cells)",
                   row_title_side = c("right"),
                   column_title = "Genomic Region",
                   column_title_side = c("bottom"))
draw(ht_tumor, heatmap_legend_side = "right") # 图例位置
dev.off()


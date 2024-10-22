library(rio)
library(plotrix)

## 2.读取数据
venn_data_index <- list.files(path = "./", pattern = "^venn_data")
venn_data <- import(venn_data_index)
sample_id <- colnames(venn_data)
otu_id <- unique(venn_data[,1])
otu_id <- otu_id[otu_id != '']
core_otu_id <- otu_id
otu_num <- length(otu_id)
## 3.遍历数据文件获取数据
for (i in 2:ncol(venn_data)) {
  otu_id <- unique(venn_data[,i])
  otu_id <- otu_id[otu_id != '']
  core_otu_id <- intersect(core_otu_id, otu_id)
  otu_num <- c(otu_num, length(otu_id))
}

core_num <- length(core_otu_id)
## 4.定义备选颜色
ellipse_col <- c('#6181BD4E','#F348004E','#64A10E4E')
## 5.创建绘图函数：flower_plot()
flower_plot <- function(sample, otu_num, core_otu, start, a, b, r, ellipse_col, circle_col) {
  par( bty = 'n', ann = F, xaxt = 'n', yaxt = 'n', mar = c(1,1,1,1))
  plot(c(0,10),c(0,10),type='n')
  n  <- length(sample)
  deg <- 360 / n
  res <- lapply(1:n, function(t){
    draw.ellipse(x = 5 + cos((start + deg * (t - 1)) * pi / 180),
                 y = 5 + sin((start + deg * (t - 1)) * pi / 180),
                 col = ellipse_col[t],
                 border = ellipse_col[t],
                 a = a, b = b, angle = deg * (t - 1))
    text(x = 5 + 2.5 * cos((start + deg * (t - 1)) * pi / 180),
         y = 5 + 2.5 * sin((start + deg * (t - 1)) * pi / 180),
         otu_num[t])
    if (deg * (t - 1) < 180 && deg * (t - 1) > 0 ) {
      text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           srt = deg * (t - 1) - start,
           adj = 1,
           cex = 1
      )
    } else {
      text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           srt = deg * (t - 1) + start,
           adj = 0,
           cex = 1
      )
    }
  })
  draw.circle(x = 5, y = 5, r = r, col = circle_col, border = NA)
  text(x = 5, y = 5, paste('Core:', core_otu))
}

## 6.调用上述函数作图（保存图片到本地文件夹中）
png('flower.png', width = 1500, height = 1500, res = 200, units = 'px')
flower_plot(sample = sample_id, otu_num = otu_num, core_otu = core_num,
            start = 90, a = 0.5, b = 2, r = 1, ellipse_col = ellipse_col, circle_col = 'white')
dev.off()
pdf('flower.pdf', width = 15, height = 15)
flower_plot(sample = sample_id, otu_num = otu_num, core_otu = core_num,
            start = 90, a = 0.5, b = 2, r = 1, ellipse_col = ellipse_col, circle_col = 'white')
dev.off()

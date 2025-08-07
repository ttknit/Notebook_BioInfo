library(ggplot2)
library(dplyr)
library(tools)
file_paths <- c(
    "./dp_file/ctrl_hifi",
    "./dp_file/rDNA_hifi",
    "./dp_file/ctrl_ngs",
    "./dp_file/rDNA_ngs"
)

all_data_list <- list()
for (i in 1:length(file_paths)) {
    current_file_path <- file_paths[i]
    if (!file.exists(current_file_path)) {
        stop(paste("错误：文件未找到，请检查路径：", current_file_path))
    }
    data <- read.table(current_file_path, header = FALSE, col.names = c("Value"))
    sample_name <- file_path_sans_ext(basename(current_file_path))
    df_current <- data.frame(
        Value = as.numeric(data$Value), # 确保数值列是数字类型
        Group = sample_name
    )

    all_data_list[[i]] <- df_current
}

combined_data <- bind_rows(all_data_list)
desired_order <- c("ctrl_hifi","rDNA_hifi","ctrl_ngs","rDNA_ngs") 
combined_data$Group <- factor(combined_data$Group, levels = desired_order)



boxplot_plot <- ggplot(combined_data, aes(x = Group, y = Value, fill = Group)) +
    geom_boxplot() + # 绘制箱线图
    # geom_boxplot(outlier.shape = NA) + # 可选：如果数据点很多，可以隐藏离群点，让图更简洁
    scale_fill_brewer(palette = "Set3") + # 为每个箱线图设置不同的颜色
    labs(
        title = "", # 图的标题
        x = "",                    # X轴标签
        y = "copy number"                         # Y轴标签
    ) +
    theme_minimal() + # 使用简洁的图表主题
    theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), # 标题居中、加粗、大小
        axis.title = element_text(size = 14),                             # 轴标题字体大小
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),     # X轴标签倾斜45度，防止重叠
        axis.text.y = element_text(size = 12),                             # Y轴标签字体大小
        legend.position = "none" # 如果 Group 已经是 X 轴标签，图例通常可以隐藏
    )
# ggsave("combined_boxplot.pdf", plot = boxplot_plot, width = 10, height = 7)
ggsave("combined_boxplot.png", plot = boxplot_plot, width = 10, height = 7, dpi = 300)

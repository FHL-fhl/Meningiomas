######5_Survival and recurrence analysis#####


library(survival)
library(survminer)
library(RColorBrewer)
library(tibble)
library(ggpp)

####### 生存曲线分析######
#导入生存信息#
b<-as.data.frame(Cluster_DRR[,c(1,5,8,9)])
rownames(b)<-b$ID
b<-b[,-1]
##将生存时间转换为月份
b$LFFR <- b$LFFR * 12 

fitd <- survdiff(Surv(LFFR, LF) ~ Cluster,
                 data      = b,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(LFFR, LF) ~ Cluster,
               data      =  b,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
# 配对生存分析
ps <- pairwise_survdiff(Surv(LFFR, LF) ~ Cluster,
                        data            = b,
                        p.adjust.method = "BH") # 这里不使用矫正，若需要矫正可以将none替换为BH
# 设置颜色
mycol <- c("#2EC4B6","#E71D36","#FF9F1C","#BDD5EA","#FFA5AB")
# 绘制基础图形
## 隐藏类标记
names(fit$strata) <- gsub("group=", "", names(fit$strata))
## 生存曲线图
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE, # 不绘制置信区间
                risk.table        = TRUE, # 生存风险表
                risk.table.col    = "strata",
                palette           = mycol, # KM曲线颜色
                data              = b,
                xlim              = c(0,120), # 时间轴，一般考虑5年（原文）或者10年长度
                size              = 1,
                break.time.by     = 12, # 时间轴的刻度（每年）
                legend.title      = "",
                xlab              = "Time (months)",
                #ylab              = "Overall Survival",
                ylab              = "Local Freedom From Recurrence",
                risk.table.y.text = FALSE,
                tables.height     = 0.3) # 风险表的高度
## 添加overall pvalue
p.lab <- paste0("log-rank test P",
                ifelse(p.val < 0.001, " < 0.001", # 若P值<0.001则标记为“<0.001”
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 0, y = 0.4, # 在y=0.55处打印overall p值
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
## 添加配对表格
addTab <- as.data.frame(as.matrix(ifelse(round(ps$p.value, 3) < 0.001, "<0.001",
                                         round(ps$p.value, 3))))
addTab[is.na(addTab)] <- "-"
df <- tibble(x = 0, y = 0, tb = list(addTab))
p$plot <- p$plot + 
  geom_table(data = df, 
             aes(x = x, y = y, label = tb), 
             table.rownames = TRUE)
p



###复发曲线分析####
#导入复发信息#
b<-as.data.frame(Cluster_DRR[,c(1,5,10,11)])
rownames(b)<-b$ID
b<-b[,-1]
##将生存时间转换为月份
b$OSb <- b$OS * 12 
# 生存分析
fitd <- survdiff(Surv(OSb, VS) ~ Cluster,
                 data      = b,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(OSb, VS) ~ Cluster,
               data      =  b,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)
# 配对生存分析
ps <- pairwise_survdiff(Surv(OSb, VS) ~ Cluster,
                        data            = b,
                        p.adjust.method = "BH") # 这里不使用矫正，若需要矫正可以将none替换为BH
# 设置颜色
mycol <- c("#2EC4B6","#E71D36","#FF9F1C","#BDD5EA","#FFA5AB")
# 绘制基础图形
## 隐藏类标记
names(fit$strata) <- gsub("group=", "", names(fit$strata))
## 生存曲线图
p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE, # 不绘制置信区间
                risk.table        = TRUE, # 生存风险表
                risk.table.col    = "strata",
                palette           = mycol, # KM曲线颜色
                data              = b,
                xlim              = c(0,120), # 时间轴，一般考虑5年（原文）或者10年长度
                size              = 1,
                break.time.by     = 12, # 时间轴的刻度（每年）
                legend.title      = "",
                xlab              = "Time (months)",
                ylab              = "Overall Survival",
                
                risk.table.y.text = FALSE,
                tables.height     = 0.3) # 风险表的高度
## 添加overall pvalue
p.lab <- paste0("log-rank test P",
                ifelse(p.val < 0.001, " < 0.001", # 若P值<0.001则标记为“<0.001”
                       paste0(" = ",round(p.val, 3))))
p$plot <- p$plot + annotate("text",
                            x = 0, y = 0.4, # 在y=0.55处打印overall p值
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)
## 添加配对表格
addTab <- as.data.frame(as.matrix(ifelse(round(ps$p.value, 3) < 0.001, "<0.001",
                                         round(ps$p.value, 3))))
addTab[is.na(addTab)] <- "-"
df <- tibble(x = 0, y = 0, tb = list(addTab))
p$plot <- p$plot + 
  geom_table(data = df, 
             aes(x = x, y = y, label = tb), 
             table.rownames = TRUE)
p



# 假设 p 是一个 ggsurvplot 对象
p_converted <- p$plot
p1_converted <- p1$plot
# 使用 cowplot 拼接两个图
library(cowplot)
combined_plot <- plot_grid(p_converted, p1_converted, ncol = 2)

# 显示拼接后的图
combined_plot

p

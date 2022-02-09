library("tools")
library("dplyr")
library("ggplot2")

data_dir <- file_path_as_absolute(paste0(getwd(), "/../data"))


# vmPFC Diff LM ----
df_amgL_vmpfc <- read.csv(paste0(data_dir, "/df_amgL-lvmPFC.csv"))
df_amgL_vmpfc$dx <- factor(df_amgL_vmpfc$dx)
df_amgL_vmpfc <- na.omit(df_amgL_vmpfc)

# find outliers of SnegLF
h_iqr <- IQR(df_amgL_vmpfc$SnegLF, na.rm = TRUE)
h_quant <- quantile(df_amgL_vmpfc$SnegLF, na.rm = TRUE, names = FALSE)
h_min <- h_quant[2] - (1.5 * h_iqr)
h_max <- h_quant[4] + (1.5 * h_iqr)

# replace outliers with min/max
ind_out <- which(df_amgL_vmpfc$SnegLF < h_min | df_amgL_vmpfc$SnegLF > h_max)
for(ind in ind_out){
  if(df_amgL_vmpfc[ind,]$SnegLF > h_max){
    df_amgL_vmpfc[ind,]$SnegLF <- h_max
  }else if(df_amgL_vmpfc[ind,]$SnegLF < h_min){
    df_amgL_vmpfc[ind,]$SnegLF <- h_min
  }
}

# find outliers of SneuLF
h_iqr <- IQR(df_amgL_vmpfc$SneuLF, na.rm = TRUE)
h_quant <- quantile(df_amgL_vmpfc$SneuLF, na.rm = TRUE, names = FALSE)
h_min <- h_quant[2] - (1.5 * h_iqr)
h_max <- h_quant[4] + (1.5 * h_iqr)

# replace outliers with min/max
ind_out <- which(df_amgL_vmpfc$SneuLF < h_min | df_amgL_vmpfc$SneuLF > h_max)
for(ind in ind_out){
  if(df_amgL_vmpfc[ind,]$SneuLF > h_max){
    df_amgL_vmpfc[ind,]$SneuLF <- h_max
  }else if(df_amgL_vmpfc[ind,]$SneuLF < h_min){
    df_amgL_vmpfc[ind,]$SneuLF <- h_min
  }
}

df_amgL_vmpfc$diff <- df_amgL_vmpfc$SnegLF - df_amgL_vmpfc$SneuLF

hist(df_amgL_vmpfc$SneuLF)
hist(df_amgL_vmpfc$SnegLF)

fit <- lm(diff ~ pscared*dx, data = df_amgL_vmpfc)
anova(fit)
summary(fit)
plot(fit)

# fit_coef <- summary(fit)$coef

# plot(df_amgL_vmpfc$pscared, df_amgL_vmpfc$diff)
# abline(a=fit_coef[1],b=fit_coef[2])
# abline(a=fit_coef[1]+fit_coef[3],b=fit_coef[2]+fit_coef[5])
# abline(a=fit_coef[1]+fit_coef[4],b=fit_coef[2]+fit_coef[6])

ggplot(data = df_amgL_vmpfc, aes(x = pscared, y = diff, color = dx)) +
  geom_point() + 
  geom_ribbon(
    stat='smooth', method = "lm", se=TRUE, alpha=0.1, 
    aes(color = NULL, group = factor(dx))) +
  geom_line(stat='smooth', method = "lm", alpha=1) +
  ggtitle("L. Amg -- vmPFC") +
  ylab("Diff: negLF-neuLF")


# vmPFC MLM -----
neg <- lm(cbind(SnegLF, SneuLF) ~ pscared*dx, data = df_amgL_vmpfc)
anova(neg)
summary(neg)

df_long <- data.frame(
  subj = rep(df_amgL_vmpfc$subj, 2),
  dx = rep(df_amgL_vmpfc$dx, 2),
  pscared = rep(df_amgL_vmpfc$pscared, 2),
  beh = rep(c("negLF", "neuLF"), each = dim(df_amgL_vmpfc)[1])
)
df_long$coef <- NA
df_long[which(df_long$beh == "negLF"),]$coef <- df_amgL_vmpfc$SnegLF
df_long[which(df_long$beh == "neuLF"),]$coef <- df_amgL_vmpfc$SneuLF

df_long %>%
  group_by(dx) %>%
  do(broom::tidy(lm(coef ~ pscared*beh)))


# dmPFC Diff LM-----
df_amgL_dmpfc <- read.csv(paste0(data_dir, "/df_amgL-ldmPFC.csv"))
df_amgL_dmpfc$dx <- factor(df_amgL_dmpfc$dx)
df_amgL_dmpfc <- na.omit(df_amgL_dmpfc)

# find outliers of SnegLF
h_iqr <- IQR(df_amgL_dmpfc$SnegLF, na.rm = TRUE)
h_quant <- quantile(df_amgL_dmpfc$SnegLF, na.rm = TRUE, names = FALSE)
h_min <- h_quant[2] - (1.5 * h_iqr)
h_max <- h_quant[4] + (1.5 * h_iqr)

# replace outliers with min/max
ind_out <- which(df_amgL_dmpfc$SnegLF < h_min | df_amgL_dmpfc$SnegLF > h_max)
for(ind in ind_out){
  if(df_amgL_dmpfc[ind,]$SnegLF > h_max){
    df_amgL_dmpfc[ind,]$SnegLF <- h_max
  }else if(df_amgL_dmpfc[ind,]$SnegLF < h_min){
    df_amgL_dmpfc[ind,]$SnegLF <- h_min
  }
}

# find outliers of SneuLF
h_iqr <- IQR(df_amgL_dmpfc$SneuLF, na.rm = TRUE)
h_quant <- quantile(df_amgL_dmpfc$SneuLF, na.rm = TRUE, names = FALSE)
h_min <- h_quant[2] - (1.5 * h_iqr)
h_max <- h_quant[4] + (1.5 * h_iqr)

# replace outliers with min/max
ind_out <- which(df_amgL_dmpfc$SneuLF < h_min | df_amgL_dmpfc$SneuLF > h_max)
for(ind in ind_out){
  if(df_amgL_dmpfc[ind,]$SneuLF > h_max){
    df_amgL_dmpfc[ind,]$SneuLF <- h_max
  }else if(df_amgL_dmpfc[ind,]$SneuLF < h_min){
    df_amgL_dmpfc[ind,]$SneuLF <- h_min
  }
}

df_amgL_dmpfc$diff <- df_amgL_dmpfc$SnegLF - df_amgL_dmpfc$SneuLF

hist(df_amgL_dmpfc$SneuLF)
hist(df_amgL_dmpfc$SnegLF)

fit <- lm(diff ~ pscared*dx, data = df_amgL_dmpfc)
anova(fit)
summary(fit)
plot(fit)

# fit_coef <- summary(fit)$coef

# plot(df_amgL_dmpfc$pscared, df_amgL_dmpfc$diff)
# abline(a=fit_coef[1],b=fit_coef[2])
# abline(a=fit_coef[1]+fit_coef[3],b=fit_coef[2]+fit_coef[5])
# abline(a=fit_coef[1]+fit_coef[4],b=fit_coef[2]+fit_coef[6])

ggplot(data = df_amgL_dmpfc, aes(x = pscared, y = diff, color = dx)) +
  geom_point() + 
  geom_smooth(method = "lm") + 
  ggtitle("L. Amg -- dmPFC") +
  ylab("Diff: negLF-neuLF")
  

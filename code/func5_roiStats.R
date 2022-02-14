library("tools")
library("dplyr")
library("ggplot2")

data_dir <- file_path_as_absolute(paste0(getwd(), "/../data"))


# amgL vmPFC LMs ----
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

hist(df_amgL_vmpfc$SneuLF)
hist(df_amgL_vmpfc$SnegLF)

fit_neg <- lm(SnegLF ~ pscared * dx, data = df_amgL_vmpfc)
summary(fit_neg)
anova(fit_neg)

fit_neu <- lm(SneuLF ~ pscared * dx, data = df_amgL_vmpfc)
summary(fit_neu)
anova(fit_neu)

ggplot(data = df_amgL_vmpfc, aes(x = pscared, y = SnegLF, color = dx)) +
  geom_point() +
  geom_ribbon(stat = "smooth", method = "lm", se = T, alpha=0.1,
              aes(color = NULL, group = factor(dx))) +
  geom_line(stat = "smooth", method = "lm", alpha = 1) +
  ggtitle("L. Amg -- vmPFC")

ggplot(data = df_amgL_vmpfc, aes(x = pscared, y = SneuLF, color = dx)) +
  geom_point() +
  geom_ribbon(stat = "smooth", method = "lm", se = T, alpha=0.1,
              aes(color = NULL, group = factor(dx))) +
  geom_line(stat = "smooth", method = "lm", alpha = 1) +
  ggtitle("L. Amg -- vmPFC")



# amgL dmPFC LMs ----
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

hist(df_amgL_dmpfc$SneuLF)
hist(df_amgL_dmpfc$SnegLF)

fit_neg <- lm(SnegLF ~ pscared * dx, data = df_amgL_dmpfc)
summary(fit_neg)
anova(fit_neg)

fit_neu <- lm(SneuLF ~ pscared * dx, data = df_amgL_dmpfc)
summary(fit_neu)
anova(fit_neu)

ggplot(data = df_amgL_dmpfc, aes(x = pscared, y = SnegLF, color = dx)) +
  geom_point() +
  geom_ribbon(stat = "smooth", method = "lm", se = T, alpha=0.1,
              aes(color = NULL, group = factor(dx))) +
  geom_line(stat = "smooth", method = "lm", alpha = 1) +
  ggtitle("L. Amg -- dmPFC")

ggplot(data = df_amgL_dmpfc, aes(x = pscared, y = SneuLF, color = dx)) +
  geom_point() +
  geom_ribbon(stat = "smooth", method = "lm", se = T, alpha=0.1,
              aes(color = NULL, group = factor(dx))) +
  geom_line(stat = "smooth", method = "lm", alpha = 1) +
  ggtitle("L. Amg -- dmPFC")



# blaL vmPFC LMs ----
df_blaL_vmpfc <- read.csv(paste0(data_dir, "/df_blaL-lvmPFC.csv"))
df_blaL_vmpfc$dx <- factor(df_blaL_vmpfc$dx)
df_blaL_vmpfc <- na.omit(df_blaL_vmpfc)

# find outliers of SnegLF
h_iqr <- IQR(df_blaL_vmpfc$SnegLF, na.rm = TRUE)
h_quant <- quantile(df_blaL_vmpfc$SnegLF, na.rm = TRUE, names = FALSE)
h_min <- h_quant[2] - (1.5 * h_iqr)
h_max <- h_quant[4] + (1.5 * h_iqr)

# replace outliers with min/max
ind_out <- which(df_blaL_vmpfc$SnegLF < h_min | df_blaL_vmpfc$SnegLF > h_max)
for(ind in ind_out){
  if(df_blaL_vmpfc[ind,]$SnegLF > h_max){
    df_blaL_vmpfc[ind,]$SnegLF <- h_max
  }else if(df_blaL_vmpfc[ind,]$SnegLF < h_min){
    df_blaL_vmpfc[ind,]$SnegLF <- h_min
  }
}

# find outliers of SneuLF
h_iqr <- IQR(df_blaL_vmpfc$SneuLF, na.rm = TRUE)
h_quant <- quantile(df_blaL_vmpfc$SneuLF, na.rm = TRUE, names = FALSE)
h_min <- h_quant[2] - (1.5 * h_iqr)
h_max <- h_quant[4] + (1.5 * h_iqr)

# replace outliers with min/max
ind_out <- which(df_blaL_vmpfc$SneuLF < h_min | df_blaL_vmpfc$SneuLF > h_max)
for(ind in ind_out){
  if(df_blaL_vmpfc[ind,]$SneuLF > h_max){
    df_blaL_vmpfc[ind,]$SneuLF <- h_max
  }else if(df_blaL_vmpfc[ind,]$SneuLF < h_min){
    df_blaL_vmpfc[ind,]$SneuLF <- h_min
  }
}

hist(df_blaL_vmpfc$SneuLF)
hist(df_blaL_vmpfc$SnegLF)

fit_neg <- lm(SnegLF ~ pscared * dx, data = df_blaL_vmpfc)
summary(fit_neg)
anova(fit_neg)

fit_neu <- lm(SneuLF ~ pscared * dx, data = df_blaL_vmpfc)
summary(fit_neu)
anova(fit_neu)

ggplot(data = df_blaL_vmpfc, aes(x = pscared, y = SnegLF, color = dx)) +
  geom_point() +
  geom_ribbon(stat = "smooth", method = "lm", se = T, alpha=0.1,
              aes(color = NULL, group = factor(dx))) +
  geom_line(stat = "smooth", method = "lm", alpha = 1) +
  ggtitle("L. BLA -- vmPFC")

ggplot(data = df_blaL_vmpfc, aes(x = pscared, y = SneuLF, color = dx)) +
  geom_point() +
  geom_ribbon(stat = "smooth", method = "lm", se = T, alpha=0.1,
              aes(color = NULL, group = factor(dx))) +
  geom_line(stat = "smooth", method = "lm", alpha = 1) +
  ggtitle("L. BLA -- vmPFC")



# blaL dmPFC LMs ----
df_blaL_dmpfc <- read.csv(paste0(data_dir, "/df_blaL-ldmPFC.csv"))
df_blaL_dmpfc$dx <- factor(df_blaL_dmpfc$dx)
df_blaL_dmpfc <- na.omit(df_blaL_dmpfc)

# find outliers of SnegLF
h_iqr <- IQR(df_blaL_dmpfc$SnegLF, na.rm = TRUE)
h_quant <- quantile(df_blaL_dmpfc$SnegLF, na.rm = TRUE, names = FALSE)
h_min <- h_quant[2] - (1.5 * h_iqr)
h_max <- h_quant[4] + (1.5 * h_iqr)

# replace outliers with min/max
ind_out <- which(df_blaL_dmpfc$SnegLF < h_min | df_blaL_dmpfc$SnegLF > h_max)
for(ind in ind_out){
  if(df_blaL_dmpfc[ind,]$SnegLF > h_max){
    df_blaL_dmpfc[ind,]$SnegLF <- h_max
  }else if(df_blaL_dmpfc[ind,]$SnegLF < h_min){
    df_blaL_dmpfc[ind,]$SnegLF <- h_min
  }
}

# find outliers of SneuLF
h_iqr <- IQR(df_blaL_dmpfc$SneuLF, na.rm = TRUE)
h_quant <- quantile(df_blaL_dmpfc$SneuLF, na.rm = TRUE, names = FALSE)
h_min <- h_quant[2] - (1.5 * h_iqr)
h_max <- h_quant[4] + (1.5 * h_iqr)

# replace outliers with min/max
ind_out <- which(df_blaL_dmpfc$SneuLF < h_min | df_blaL_dmpfc$SneuLF > h_max)
for(ind in ind_out){
  if(df_blaL_dmpfc[ind,]$SneuLF > h_max){
    df_blaL_dmpfc[ind,]$SneuLF <- h_max
  }else if(df_blaL_dmpfc[ind,]$SneuLF < h_min){
    df_blaL_dmpfc[ind,]$SneuLF <- h_min
  }
}

hist(df_blaL_dmpfc$SneuLF)
hist(df_blaL_dmpfc$SnegLF)

fit_neg <- lm(SnegLF ~ pscared * dx, data = df_blaL_dmpfc)
summary(fit_neg)
anova(fit_neg)

fit_neu <- lm(SneuLF ~ pscared * dx, data = df_blaL_dmpfc)
summary(fit_neu)
anova(fit_neu)

ggplot(data = df_blaL_dmpfc, aes(x = pscared, y = SnegLF, color = dx)) +
  geom_point() +
  geom_ribbon(stat = "smooth", method = "lm", se = T, alpha=0.1,
              aes(color = NULL, group = factor(dx))) +
  geom_line(stat = "smooth", method = "lm", alpha = 1) +
  ggtitle("L. BLA -- dmPFC")

ggplot(data = df_blaL_dmpfc, aes(x = pscared, y = SneuLF, color = dx)) +
  geom_point() +
  geom_ribbon(stat = "smooth", method = "lm", se = T, alpha=0.1,
              aes(color = NULL, group = factor(dx))) +
  geom_line(stat = "smooth", method = "lm", alpha = 1) +
  ggtitle("L. BLA -- dmPFC")

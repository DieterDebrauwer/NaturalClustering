
# R script accompanying Section 7.2 (Example 2: nested Archimedean copula) of 
# "Portfolio optimization under natural clustering" by Debrauwer and Gijbels

#In this script we illustrate how to construct figure 2 containing boxplots.
#We do assume that previous computations have been done already and have been saved locally on your PC.

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggh4x)


# Load first dataset (Setting 1: extremile)
load("simulationStudyHACextremile14SampleSize50.Rdata")
load("simulationStudyHACextremile14SampleSize50Opt1To5.Rdata")
trueExtrsingle=1.566232
trueExtrOpt1=1.5903
trueExtrOpt2=1.596338
trueExtrOpt3Representative=1.633811 
trueExtrOpt4=1.568372
trueExtrOpt5=1.590797


df1 <- data.frame(
  Method = rep(c("SingleStepES", "Opt1ES", "Opt2ES", "Opt3ES","Opt4ES", "Opt5ES"), each = 500),
  EffectSize = c(weightEstimateMatrix[,5], weightEstimateMatrixOpt1[,5], weightEstimateMatrixOpt2[,5], weightEstimateMatrixOpt3Representative[,5],
                 weightEstimateMatrixOpt4[,5], weightEstimateMatrixOpt5[,5]),
  Setting = "Setting 1 Extremile"
)

trueES1 <- data.frame(
  Method = c("SingleStepES", "Opt1ES", "Opt2ES", "Opt3ES","Opt4ES", "Opt5ES"),
  TrueValue = c(trueExtrsingle, trueExtrOpt1, trueExtrOpt2, trueExtrOpt3Representative,trueExtrOpt4,trueExtrOpt5),
  Setting = "Setting 1 Extremile"
)







# Load first dataset (Setting 2: extremile)
load("simulationStudyHACextremile23SampleSize50.Rdata")
load("simulationStudyHACextremile23SampleSize50Opt1To5.Rdata")
trueExtrsingle=1.584063
trueExtrOpt1=1.573776
trueExtrOpt2=1.637306
trueExtrOpt3Representative=1.619796
trueExtrOpt4=1.831826
trueExtrOpt5=1.616228

df2 <- data.frame(
  Method = rep(c("SingleStepES", "Opt1ES", "Opt2ES", "Opt3ES","Opt4ES", "Opt5ES"), each = 500),
  EffectSize = c(weightEstimateMatrix[,5], weightEstimateMatrixOpt1[,5], weightEstimateMatrixOpt2[,5], weightEstimateMatrixOpt3Representative[,5],
                 weightEstimateMatrixOpt4[,5], weightEstimateMatrixOpt5[,5]),
  Setting = "Setting 2 Extremile"
)

trueES2 <- data.frame(
  Method = c("SingleStepES", "Opt1ES", "Opt2ES", "Opt3ES","Opt4ES", "Opt5ES"),
  TrueValue = c(trueExtrsingle, trueExtrOpt1, trueExtrOpt2, trueExtrOpt3Representative,trueExtrOpt4,trueExtrOpt5),
  Setting = "Setting 2 Extremile"
)





# Load first dataset (Setting 3: extremile)
load("simulationHACsetting3ExtremileSampleSize50.RData")
trueExtrsingle=1.634145
trueExtrOpt1=1.624993
trueExtrOpt2=1.671036
trueExtrOpt3Representative=1.619796
trueExtrOpt4=1.964674
trueExtrOpt5=1.653633

df3 <- data.frame(
  Method = rep(c("SingleStepES", "Opt1ES", "Opt2ES", "Opt3ES","Opt4ES", "Opt5ES"), each = 500),
  EffectSize = c(weightEstimateMatrix[,5], weightEstimateMatrixOpt1[,5], weightEstimateMatrixOpt2[,5], weightEstimateMatrixOpt3Representative[,5],
                 weightEstimateMatrixOpt4[,5], weightEstimateMatrixOpt5[,5]),
  Setting = "Setting 3 Extremile"
)

trueES3 <- data.frame(
  Method = c("SingleStepES", "Opt1ES", "Opt2ES", "Opt3ES","Opt4ES", "Opt5ES"),
  TrueValue = c(trueExtrsingle, trueExtrOpt1, trueExtrOpt2, trueExtrOpt3Representative,trueExtrOpt4, trueExtrOpt5),
  Setting = "Setting 3 Extremile"
)










# Load first dataset (Setting 1: ES)
load("simulationStudyHACES14SampleSize50.Rdata")
load("simulationStudyHACES14SampleSize50Opt1To5.Rdata")
trueExtrsingle=1.903419
trueExtrOpt1=1.906987
trueExtrOpt2=1.913347
trueExtrOpt3Representative=1.951604
trueExtrOpt4=1.906898
trueExtrOpt5=1.908894

df4 <- data.frame(
  Method = rep(c("SingleStepES", "Opt1ES", "Opt2ES", "Opt3ES","Opt4ES", "Opt5ES"), each=500),
  EffectSize = c(weightEstimateMatrix[,5], weightEstimateMatrixOpt1[,5], weightEstimateMatrixOpt2[,5], weightEstimateMatrixOpt3Representative[,5],
                 weightEstimateMatrixOpt4[,5], weightEstimateMatrixOpt5[,5]),  Setting = "Setting 1 ES"
)

trueES4 <- data.frame(
  Method = c("SingleStepES", "Opt1ES", "Opt2ES", "Opt3ES","Opt4ES", "Opt5ES"),
  TrueValue = c(trueExtrsingle, trueExtrOpt1, trueExtrOpt2, trueExtrOpt3Representative,trueExtrOpt4,trueExtrOpt5),
  Setting = "Setting 1 ES"
)







# Load first dataset (Setting 2: ES)
load("simulationStudyHACES23SampleSize50.Rdata")
load("simulationStudyHACES23SampleSize50Opt1To5.Rdata")
trueExtrsingle=1.943854
trueExtrOpt1=1.932433
trueExtrOpt2=2.00597
trueExtrOpt3Representative=1.979807
trueExtrOpt4=2.221211
trueExtrOpt5=1.980336

df5 <- data.frame(
  Method = rep(c("SingleStepES", "Opt1ES", "Opt2ES", "Opt3ES","Opt4ES", "Opt5ES"), each=500),
  EffectSize = c(weightEstimateMatrix[,5], weightEstimateMatrixOpt1[,5], weightEstimateMatrixOpt2[,5], weightEstimateMatrixOpt3Representative[,5],
                 weightEstimateMatrixOpt4[,5], weightEstimateMatrixOpt5[,5]),  Setting = "Setting 2 ES"
)

trueES5 <- data.frame(
  Method = c("SingleStepES", "Opt1ES", "Opt2ES", "Opt3ES","Opt4ES", "Opt5ES"),
  TrueValue = c(trueExtrsingle, trueExtrOpt1, trueExtrOpt2, trueExtrOpt3Representative,trueExtrOpt4,trueExtrOpt5),
  Setting = "Setting 2 ES"
)





# Load first dataset (Setting 3: ES)
load("simulationHACsetting3ESSampleSize50.RData")
load("simulationHACsetting3ESSingleStepSampleSize50.RData")
trueExtrsingle=1.998166 
trueExtrOpt1=1.987322
trueExtrOpt2=2.034554
trueExtrOpt3Representative=1.979807
trueExtrOpt4=2.380255
trueExtrOpt5=2.022879


df6 <- data.frame(
  Method = rep(c("SingleStepES", "Opt1ES", "Opt2ES", "Opt3ES","Opt4ES", "Opt5ES"), each = 500),
  EffectSize = c(weightEstimateMatrix[,5], weightEstimateMatrixOpt1[,5], weightEstimateMatrixOpt2[,5], weightEstimateMatrixOpt3Representative[,5],
                 weightEstimateMatrixOpt4[,5], weightEstimateMatrixOpt5[,5]),
  Setting = "Setting 3 ES"
)

trueES6 <- data.frame(
  Method = c("SingleStepES", "Opt1ES", "Opt2ES", "Opt3ES","Opt4ES", "Opt5ES"),
  TrueValue = c(trueExtrsingle, trueExtrOpt1, trueExtrOpt2, trueExtrOpt3Representative,trueExtrOpt4, trueExtrOpt5),
  Setting = "Setting 3 ES"
)









# Combine datasets
df <- rbind(df1, df2, df3, df4,df5,df6)
trueES <- rbind(trueES1, trueES2, trueES3, trueES4,trueES5,trueES6)

# Ensure correct ordering: First all Setting 1, then Setting 2, etc.
df$Group <- factor(paste(df$Setting, df$Method, sep = " - "), 
                   levels = unique(paste(df$Setting, df$Method, sep = " - ")))
trueES$Group <- factor(paste(trueES$Setting, trueES$Method, sep = " - "), 
                       levels = levels(df$Group))

# Assign distinct colors for each setting
distinct_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3","#FF7F00","#A65628")  # More contrast
df$ColorGroup <- factor(df$Setting, levels = unique(df$Setting))
trueES$ColorGroup <- factor(trueES$Setting, levels = unique(trueES$Setting))

manual_xpos <- c(1:6, 9:14, 17:22, 25:30,33:38,41:46)  # Match your manual spacing
# Map trueES$Group to manual_xpos
trueES$xpos <- manual_xpos[match(trueES$Group, levels(df$Group))]


# Darken the colors for the lines
dark_colors <- adjustcolor(distinct_colors, red.f = 0.7, green.f = 0.7, blue.f = 0.7)  # Darker shades

ggplot(df, aes(x = Group, y = EffectSize, fill = ColorGroup)) +
  geom_boxplot(outlier.shape = 1, width = 0.5, alpha = 0.5, color = "black") +  # Lighter boxplots with black borders
  scale_x_manual(manual_xpos)+
  geom_segment(data = trueES, 
               aes(x = xpos - 0.27, xend = xpos + 0.27, 
                   y = TrueValue, yend = TrueValue, color = ColorGroup), 
               size = 1.5,lty="11") +
  scale_fill_manual(values = adjustcolor(distinct_colors, alpha.f = 0.5)) +  # More transparent fill colors
  scale_color_manual(values = dark_colors) +  # Darker lines matching boxplot colors
  scale_y_continuous(limits = c(0.5, 3.25), breaks = seq(0.5, 3.25, by = 0.25)) +
  theme_minimal() +
  labs(title = "n=50", x = "", y = " ") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),
        legend.position = "none") +  
  annotate("text", x = c(3.5,11.5,19.5,27.5,35.5,43.5),
           y = 0.5,
           label = unique(df$Setting), size = 4, fontface = "bold")  # Add setting labels below

























































































# Load first dataset (Setting 1: extremile)
load("simulationStudyHACextremile14SampleSize400.Rdata")
load("simulationStudyHACextremile14SampleSize400Opt1To5.Rdata")
trueExtrsingle=1.566232
trueExtrOpt1=1.5903
trueExtrOpt2=1.596338
trueExtrOpt3Representative=1.633811 
trueExtrOpt4=1.568372
trueExtrOpt5=1.590797


df1 <- data.frame(
  Method = rep(c("SingleStepES", "Opt1ES", "Opt2ES", "Opt3ES","Opt4ES", "Opt5ES"), each = 500),
  EffectSize = c(weightEstimateMatrix[,5], weightEstimateMatrixOpt1[,5], weightEstimateMatrixOpt2[,5], weightEstimateMatrixOpt3Representative[,5],
                 weightEstimateMatrixOpt4[,5], weightEstimateMatrixOpt5[,5]),
  Setting = "Setting 1 Extremile"
)

trueES1 <- data.frame(
  Method = c("SingleStepES", "Opt1ES", "Opt2ES", "Opt3ES","Opt4ES", "Opt5ES"),
  TrueValue = c(trueExtrsingle, trueExtrOpt1, trueExtrOpt2, trueExtrOpt3Representative,trueExtrOpt4,trueExtrOpt5),
  Setting = "Setting 1 Extremile"
)







# Load first dataset (Setting 2: extremile)
load("simulationStudyHACextremile23SampleSize400.Rdata")
load("simulationStudyHACextremile23SampleSize400Opt1To5.Rdata")
trueExtrsingle=1.584063
trueExtrOpt1=1.573776
trueExtrOpt2=1.637306
trueExtrOpt3Representative=1.619796
trueExtrOpt4=1.831826
trueExtrOpt5=1.616228

df2 <- data.frame(
  Method = rep(c("SingleStepES", "Opt1ES", "Opt2ES", "Opt3ES","Opt4ES", "Opt5ES"), each = 500),
  EffectSize = c(weightEstimateMatrix[,5], weightEstimateMatrixOpt1[,5], weightEstimateMatrixOpt2[,5], weightEstimateMatrixOpt3Representative[,5],
                 weightEstimateMatrixOpt4[,5], weightEstimateMatrixOpt5[,5]),
  Setting = "Setting 2 Extremile"
)

trueES2 <- data.frame(
  Method = c("SingleStepES", "Opt1ES", "Opt2ES", "Opt3ES","Opt4ES", "Opt5ES"),
  TrueValue = c(trueExtrsingle, trueExtrOpt1, trueExtrOpt2, trueExtrOpt3Representative,trueExtrOpt4,trueExtrOpt5),
  Setting = "Setting 2 Extremile"
)





# Load first dataset (Setting 3: extremile)
load("simulationHACsetting3ExtremileSampleSize400.RData")
trueExtrsingle=1.634145
trueExtrOpt1=1.624993
trueExtrOpt2=1.671036
trueExtrOpt3Representative=1.619796
trueExtrOpt4=1.964674
trueExtrOpt5=1.653633

df3 <- data.frame(
  Method = rep(c("SingleStepES", "Opt1ES", "Opt2ES", "Opt3ES","Opt4ES", "Opt5ES"), each = 500),
  EffectSize = c(weightEstimateMatrix[,5], weightEstimateMatrixOpt1[,5], weightEstimateMatrixOpt2[,5], weightEstimateMatrixOpt3Representative[,5],
                 weightEstimateMatrixOpt4[,5], weightEstimateMatrixOpt5[,5]),
  Setting = "Setting 3 Extremile"
)

trueES3 <- data.frame(
  Method = c("SingleStepES", "Opt1ES", "Opt2ES", "Opt3ES","Opt4ES", "Opt5ES"),
  TrueValue = c(trueExtrsingle, trueExtrOpt1, trueExtrOpt2, trueExtrOpt3Representative,trueExtrOpt4, trueExtrOpt5),
  Setting = "Setting 3 Extremile"
)










# Load first dataset (Setting 1: ES)
load("simulationStudyHACES14SampleSize400.Rdata")
load("simulationStudyHACES14SampleSize400Opt1To5.Rdata")
trueExtrsingle=1.903419
trueExtrOpt1=1.906987
trueExtrOpt2=1.913347
trueExtrOpt3Representative=1.951604
trueExtrOpt4=1.906898
trueExtrOpt5=1.908894

df4 <- data.frame(
  Method = rep(c("SingleStepES", "Opt1ES", "Opt2ES", "Opt3ES","Opt4ES", "Opt5ES"), each=500),
  EffectSize = c(weightEstimateMatrix[,5], weightEstimateMatrixOpt1[,5], weightEstimateMatrixOpt2[,5], weightEstimateMatrixOpt3Representative[,5],
                 weightEstimateMatrixOpt4[,5], weightEstimateMatrixOpt5[,5]),  Setting = "Setting 1 ES"
)

trueES4 <- data.frame(
  Method = c("SingleStepES", "Opt1ES", "Opt2ES", "Opt3ES","Opt4ES", "Opt5ES"),
  TrueValue = c(trueExtrsingle, trueExtrOpt1, trueExtrOpt2, trueExtrOpt3Representative,trueExtrOpt4,trueExtrOpt5),
  Setting = "Setting 1 ES"
)







# Load first dataset (Setting 2: ES)
load("simulationStudyHACES23SampleSize400.Rdata")
load("simulationStudyHACES23SampleSize400Opt1To5.Rdata")
trueExtrsingle=1.943854
trueExtrOpt1=1.932433
trueExtrOpt2=2.00597
trueExtrOpt3Representative=1.979807
trueExtrOpt4=2.221211
trueExtrOpt5=1.980336

df5 <- data.frame(
  Method = rep(c("SingleStepES", "Opt1ES", "Opt2ES", "Opt3ES","Opt4ES", "Opt5ES"), each=500),
  EffectSize = c(weightEstimateMatrix[,5], weightEstimateMatrixOpt1[,5], weightEstimateMatrixOpt2[,5], weightEstimateMatrixOpt3Representative[,5],
                 weightEstimateMatrixOpt4[,5], weightEstimateMatrixOpt5[,5]),  Setting = "Setting 2 ES"
)

trueES5 <- data.frame(
  Method = c("SingleStepES", "Opt1ES", "Opt2ES", "Opt3ES","Opt4ES", "Opt5ES"),
  TrueValue = c(trueExtrsingle, trueExtrOpt1, trueExtrOpt2, trueExtrOpt3Representative,trueExtrOpt4,trueExtrOpt5),
  Setting = "Setting 2 ES"
)





# Load first dataset (Setting 3: ES)
load("simulationHACsetting3ESSampleSize400.RData")
load("simulationHACsetting3ESSingleStepSampleSize400.RData")

trueExtrsingle=1.998166 
trueExtrOpt1=1.987322
trueExtrOpt2=2.034554
trueExtrOpt3Representative=1.979807
trueExtrOpt4=2.380255
trueExtrOpt5=2.022879

df6 <- data.frame(
  Method = rep(c("SingleStepES", "Opt1ES", "Opt2ES", "Opt3ES","Opt4ES", "Opt5ES"), each = 500),
  EffectSize = c(weightEstimateMatrix[,5], weightEstimateMatrixOpt1[,5], weightEstimateMatrixOpt2[,5], weightEstimateMatrixOpt3Representative[,5],
                 weightEstimateMatrixOpt4[,5], weightEstimateMatrixOpt5[,5]),
  Setting = "Setting 3 ES"
)

trueES6 <- data.frame(
  Method = c("SingleStepES", "Opt1ES", "Opt2ES", "Opt3ES","Opt4ES", "Opt5ES"),
  TrueValue = c(trueExtrsingle, trueExtrOpt1, trueExtrOpt2, trueExtrOpt3Representative,trueExtrOpt4, trueExtrOpt5),
  Setting = "Setting 3 ES"
)









# Combine datasets
df <- rbind(df1, df2, df3, df4,df5,df6)
trueES <- rbind(trueES1, trueES2, trueES3, trueES4,trueES5,trueES6)

# Ensure correct ordering: First all Setting 1, then Setting 2, etc.
df$Group <- factor(paste(df$Setting, df$Method, sep = " - "), 
                   levels = unique(paste(df$Setting, df$Method, sep = " - ")))
trueES$Group <- factor(paste(trueES$Setting, trueES$Method, sep = " - "), 
                       levels = levels(df$Group))

# Assign distinct colors for each setting
distinct_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3","#FF7F00","#A65628")  # More contrast
df$ColorGroup <- factor(df$Setting, levels = unique(df$Setting))
trueES$ColorGroup <- factor(trueES$Setting, levels = unique(trueES$Setting))

manual_xpos <- c(1:6, 9:14, 17:22, 25:30,33:38,41:46)  # Match your manual spacing
# Map trueES$Group to manual_xpos
trueES$xpos <- manual_xpos[match(trueES$Group, levels(df$Group))]


# Darken the colors for the lines
dark_colors <- adjustcolor(distinct_colors, red.f = 0.7, green.f = 0.7, blue.f = 0.7)  # Darker shades

ggplot(df, aes(x = Group, y = EffectSize, fill = ColorGroup)) +
  geom_boxplot(outlier.shape = 1, width = 0.5, alpha = 0.5, color = "black") +  # Lighter boxplots with black borders
  scale_x_manual(manual_xpos)+
  geom_segment(data = trueES, 
               aes(x = xpos - 0.27, xend = xpos + 0.27, 
                   y = TrueValue, yend = TrueValue, color = ColorGroup), 
               size = 1.5,lty="11") +
  scale_fill_manual(values = adjustcolor(distinct_colors, alpha.f = 0.5)) +  # More transparent fill colors
  scale_color_manual(values = dark_colors) +  # Darker lines matching boxplot colors
  scale_y_continuous(limits = c(0.5, 3.25), breaks = seq(0.5, 3.25, by = 0.25)) +
  theme_minimal() +
  labs(title = "n=400", x = "", y = " ") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(),
        legend.position = "none") +  
  annotate("text", x = c(3.5,11.5,19.5,27.5,35.5,43.5),
           y = 0.5,
           label = unique(df$Setting), size = 4, fontface = "bold")  # Add setting labels below








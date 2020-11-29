library(googlesheets4)
library(tidyverse)
library(reshape2)
library(cowplot)
theme_set(theme_cowplot())
library(ggsignif)
library(patchwork)
library(wesanderson)
library(gt)
library(pwr)

# Power Tests -------------------------------------------------------------
# Power to detect large effect (f=0.4) with 10 reps per construct
pwr.anova.test(k=5, n=13, sig.level = 0.05, f=0.4, power=NULL)
# power to detect smell effect (f=0.1) with 10 reps per construct
pwr.anova.test(k=5, n=13, sig.level = 0.05, f=0.1, power=NULL)
# reps needs to detect a small effect with a good power (p=0.8)
pwr.anova.test(k=5, n=NULL, sig.level = 0.05, f=0.1, power=0.8)
# reps needs to detect a large effect with a good power (p=0.8)
pwr.anova.test(k=5, n=NULL, sig.level = 0.05, f=0.4, power=0.8)
# power test for branching
# power to detect large effect
pwr.chisq.test(w=0.5, N=135, df=4, sig.level = 0.05, power=NULL)
# power to detect medium effect
pwr.chisq.test(w=0.3, N=135, df=4, sig.level = 0.05, power=NULL)
# power to detect small effect
pwr.chisq.test(w=0.1, N=135, df=4, sig.level = 0.05, power=NULL)
# reps need to detect small branching effect
pwr.chisq.test(w=0.1, N=NULL, df=4, sig.level = 0.05, power=0.8)


# Read and Clean data -----------------------------------------------------
# Read in Data
phenos <- read_sheet("https://docs.google.com/spreadsheets/d/1skhTjnfVESJ_eSaFTbt7_D1w33bDZeJ4WndzfEznxPc/edit?usp=sharing",
                     sheet="Sorted")
# Massage design factor levels
phenos$Construct <- dplyr::recode(phenos$Construct, Abeer="All")
phenos$Construct <- factor(phenos$Construct, 
                           levels = c("Wild-Type", "Empty Vector", "euFULI", "euFULII", "All"))
phenos$InfiltrationLocation <- factor(phenos$InfiltrationLocation,
                                      levels = c("Rosette", "Cauline"))
# Massage observation columns data
PhenoCols <- 4:8 # Column indicies containing growth phenotypes
RelDays <- c(11,14,17,19,27) # Days since infiltration for each column in PhenoCols
phenotypes <- c("Vegetative",
                "Bolted",
                "Flower Buds",
                "Flower1",
                "Flower2",
                "Flower3",
                "Fruit") # In cases of overlap the highest level on the plant wins
phenos[,PhenoCols] <-lapply(phenos[,PhenoCols], 
                            FUN=function(i) factor(i, levels=phenotypes))

# throw out bad plant that was mislabeled
phenos <- phenos[-c(phenos$EUID==64),]


# Create Summary Columns --------------------------------------------------
# Create Days to Bolting Column
phenos$D2Bolt <- RelDays[apply(phenos[,PhenoCols],1,  function(x) min(grep("Bolted|Flower|Fruit",x)))]
# Code date for plants we can't see
phenos$D2Bolt[phenos[max(PhenoCols)]=="Vegetative"] <- 40
#phenos$D2Bolt[apply(phenos[,c(max(PhenoCols))],1,function(x) grepl("Vegetative", x))] <- 40
phenos$D2Bolt[phenos[min(PhenoCols)]=="Flower Buds"] <- 3

# Create a Days to First Flower Column
phenos$D2FFlower <- RelDays[apply(phenos[,PhenoCols],1,function(x) min(grep("Flower\\d+|Fruit",x)))]
# Code date for plants we can't see
phenos$D2FFlower[apply(phenos[,max(PhenoCols)],1,function(x) !grepl("Flower\\d+|Fruit", x))] <- 40
phenos$D2FFlower[apply(phenos[,min(PhenoCols)],1,function(x) grepl("Flower\\d+|Fruit", x))] <- 3

# Create a Days to First Flower Column
phenos$D2FFruit <- RelDays[apply(phenos[,PhenoCols],1,function(x) min(grep("Fruit",x)))]
# Code date for plants we can't see
phenos$D2FFruit[apply(phenos[,max(PhenoCols)],1,function(x) !grepl("Fruit", x))] <- 40


# Melt and Ccst data for summary ------------------------------------------
# Reshape Data for plots
phenos.melt <- melt(phenos[,c(1:3,PhenoCols)], 
                    id.vars = c("EUID", "Construct", "InfiltrationLocation"))
phenos.melt$value <- dplyr::recode(phenos.melt$value, Flower1="Flower", Flower2="Flower", Flower3="Flower")
phenos.melt$value <- factor(phenos.melt$value, 
                            levels=c("Fruit", "Flower", "Flower Buds", "Bolted", "Vegetative"))

colnames(phenos.melt) <- c("EUID", "Construct", "InfiltrationLocation", "Date", "Stage")
# Reshape again to summarize
phenos.cast <- phenos.melt %>% count(Construct, InfiltrationLocation, Date, Stage)

# Linear Models -----------------------------------------------------------
lm.phenos1 <- lm(as.numeric(D2Bolt) ~ Construct,
                 data = phenos,
                 subset = InfiltrationLocation=="Rosette")
lm.anova1 <- anova(lm.phenos1)

lm.phenos2 <- lapply(unique(phenos$InfiltrationLocation),
                     function(i) lm(as.numeric(D2FFlower) ~ Construct,
                                    data = phenos[phenos$InfiltrationLocation==i,]))
names(lm.phenos2) <- unique(phenos$InfiltrationLocation)
lm.anova2 <- lapply(lm.phenos2, anova)
#lapply(lapply(lm.phenos2,aov), TukeyHSD) # ANOVA is ns

lm.phenos3 <- lapply(unique(phenos$InfiltrationLocation),
                        function(i) lm(as.numeric(D2FFruit) ~ Construct,
                                       data = phenos[phenos$InfiltrationLocation==i,]))
names(lm.phenos3) <- unique(phenos$InfiltrationLocation)
lm.anova3 <- lapply(lm.phenos3, anova)
#lapply(lapply(lm.phenos3,aov), TukeyHSD) # ANOVA is ns


# Neg Binomial Models -----------------------------------------------------
library(MASS)
library(car)
library(emmeans)
#Create the model for each infiltration location separately
glm.phenos1 <- lapply(as.character(unique(phenos$InfiltrationLocation)),
                     function(i) glm.nb(Branches_Nov22 ~ Construct,
                                        control=glm.control(maxit = 1000),
                                        data=phenos[phenos$InfiltrationLocation==i,]))
glm.phenos2 <- lapply(as.character(unique(phenos$InfiltrationLocation)),
                      function(i) glm(Branches_Nov22 ~ Construct,
                                         family="poisson",
                                         data=phenos[phenos$InfiltrationLocation==i,]))
names(glm.phenos2) <- unique(phenos$InfiltrationLocation)
# Show ANOVA tables
nb.anova2 <- lapply(glm.phenos2, Anova, type="II", test="LR")

# Show mean separation letters
lapply(glm.phenos2,
       function(i) {
         marginal = emmeans(i, ~Construct)
         #pwpm(marginal, adjust="tukey") # Redundant, but useful sometimes
         CLD(marginal,
             alpha=0.05,
             Letters=letters,
             type="response",
             adjust="sidak")
       })


# Plots -------------------------------------------------------------------
# Proportion in each stage at given date, separated by construct
Plot_1 <- lapply(unique(phenos.cast$Construct),
      function(i) ggplot(phenos.cast[phenos.cast$Construct==i,],
                         aes(fill=Stage, y=n, x=Date)) +
        geom_bar(position="fill", stat="identity") +
        scale_x_discrete(labels=RelDays) +
        scale_y_continuous(labels = scales::percent) +
        scale_fill_manual(drop=FALSE,
                          values=wes_palette("Zissou1", 5)) +
        ggtitle(i) +
        labs(y=element_blank(),
             x="Days Post Infiltration") +
        theme(plot.title = element_text(hjust = 0.5)))
names(Plot_1) <- unique(phenos.cast$Construct)

# Number of branches on Nov22 by construct, separated by infiltration location
pal <- c("#000000FF",
         "#808080FF",
         "#2A788EFF",
         "#FFCC00FF", 
         "#7AD151FF")
Plot_2 <- lapply(as.character(unique(phenos$InfiltrationLocation)),
                 function(i) ggplot(phenos[phenos$InfiltrationLocation==i,],
                                    aes(x=Construct, y=Branches_Nov22, fill=Construct)) +
                   geom_violin(alpha=0.6) +
                   scale_fill_manual(values=pal) +
                   ggtitle(paste(i, "Leaf Infiltration"),
                           subtitle = bquote("Likelihood Ratio Chi"^2*"="*.(round(nb.anova2[[i]]$`LR Chisq`,2))*", DF="*.(nb.anova2[[i]]$Df)*", p="*.(round(nb.anova2[[i]]$`Pr(>Chisq)`,3)))) +
                   labs(y="Branch Number") +
                   theme(plot.title = element_text(hjust = 0.5),
                         plot.subtitle = element_text(hjust=0.5),
                         legend.position = "none"))
names(Plot_2) <- unique(phenos$InfiltrationLocation)

# Days to bolting by construct for Rosette leaf infiltration
Plot_3 <- ggplot(phenos[phenos$InfiltrationLocation=="Rosette",], 
                 aes(x=Construct, y=D2Bolt, fill=Construct)) +
  labs(y="Days to Bolting") +
  geom_violin(alpha=0.6) +
  scale_fill_manual(values=pal) +
  ggtitle("Rosette Leaf Infiltration",
          subtitle=bquote("F"[.(paste(lm.anova1[,1], collapse = ","))]*"="*.(round(lm.anova1[1,4],3))*", p="*.(round(lm.anova1[1,5],3)))) +
  theme(plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5),
        legend.position = "none")

# Days to first flower by construct separated by infiltration location
Plot_4 <- lapply(as.character(unique(phenos$InfiltrationLocation)),
                 function(i) ggplot(phenos[phenos$InfiltrationLocation==i,],
                                    aes(x=Construct, y=D2FFlower, fill=Construct)) +
                   labs(y="Days to First Flower") +
                   geom_violin(alpha=0.6) +
                   scale_fill_manual(values=pal) +
                   ggtitle(paste(i, "Leaf Infiltration"),
                           subtitle=bquote("F"[.(paste(lm.anova2[[i]][,1], collapse = ","))]*"="*.(round(lm.anova2[[i]][1,4],3))*", p="*.(round(lm.anova2[[i]][1,5],3)))) +
                   theme(plot.title = element_text(hjust=0.5),
                         plot.subtitle = element_text(hjust=0.5),
                         legend.position = "none"))
names(Plot_4) <- unique(phenos$InfiltrationLocation)

# Days to first fruit by construct separated by infiltration location
Plot_5 <- lapply(as.character(unique(phenos$InfiltrationLocation)),
                 function(i) ggplot(phenos[phenos$InfiltrationLocation==i,],
                                    aes(x=Construct, y=D2FFruit, fill=Construct)) +
                   labs(y="Days to First Fruit") +
                   geom_violin(alpha=0.6) +
                   scale_fill_manual(values=pal) +
                   ggtitle(paste(i, "Leaf Infiltration"),
                           subtitle=bquote("F"[.(paste(lm.anova3[[i]][,1], collapse = ","))]*"="*.(round(lm.anova3[[i]][1,4],3))*", p="*.(round(lm.anova3[[i]][1,5],3)))) +
                   theme(plot.title = element_text(hjust=0.5),
                         plot.subtitle = element_text(hjust=0.5),
                         legend.position = "none"))
names(Plot_5) <- unique(phenos$InfiltrationLocation)

# Assemble Figures --------------------------------------------------------
(Plot_1[['euFULI']] + 
   Plot_1[['euFULII']] + 
   Plot_1[['All']] + 
   Plot_1[['Wild-Type']] + 
   Plot_1[['Empty Vector']] +
   guide_area()) +
  plot_layout(guides="collect",
              nrow=2) +
  plot_annotation(tag_levels = "A")
ggsave2("Stage_By_Construct.pdf",
       height=6,
       width=9)

(Plot_2[['Rosette']] + Plot_2[['Cauline']]) +
  plot_annotation(tag_levels = "A")
ggsave2("Branch_By_ConstructLocation.pdf",
        height=4,
        width=12)

Plot_3
ggsave2("Bolt_By_Construct.pdf",
        height=4,
        width=6)

(Plot_4[['Rosette']] + Plot_4[['Cauline']]) +
  plot_annotation(tag_levels = "A")
ggsave2("Flower_By_ConstructLocation.pdf",
        height=4,
        width=12)

(Plot_5[['Rosette']] + Plot_5[['Cauline']]) +
  plot_annotation(tag_levels = "A")
ggsave2("Fruit_By_ConstructLocation.pdf",
        height=4,
        width=12)

(Plot_5[['Rosette']] + Plot_5[['Cauline']] +
    Plot_4[['Rosette']] + Plot_4[['Cauline']] +
    Plot_3) +
  plot_layout(nrow=3) +
  plot_annotation(tag_levels = "A")
ggsave2("Days_By_ConstructLocation.pdf",
        height=12,
        width=12)


# Replicates Table --------------------------------------------------------
phenos %>%
  count(InfiltrationLocation, Construct) %>%
  spread(Construct, n) %>%
  mutate('\u03A3'=rowSums(.[2:6])) %>%
  column_to_rownames("InfiltrationLocation") %>%
  gt(rownames_to_stub = TRUE) %>%
  tab_header(title=md("**Summary of Replicates**")) %>%
  grand_summary_rows(columns = 2:7,
                     fns = list(
                       '\u03A3' = "sum"),
                     decimals=0) %>%
  tab_style(locations = list(cells_grand_summary(),
                             cells_column_labels(columns = 6),
                             cells_body(columns = 7)),
            style = list(cell_text(weight = "bold"))) %>%
  gtsave("Replicates_Table.png")




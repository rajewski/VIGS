library(googlesheets4)
library(tidyverse)
library(reshape2)
library(cowplot)
theme_set(theme_cowplot())
library(ggsignif)

# Read in Data
phenos <- read_sheet("https://docs.google.com/spreadsheets/d/1skhTjnfVESJ_eSaFTbt7_D1w33bDZeJ4WndzfEznxPc/edit?usp=sharing",
                     sheet="Sorted")
# Massage design factor levels
phenos$Construct <- factor(phenos$Construct, 
                           levels = c("Wild-Type", "Empty Vector", "euFULI", "euFULII", "Abeer"))
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

# Create Days to Bolting Column
phenos$D2Bolt <- RelDays[apply(phenos[,PhenoCols],1,  function(x) min(grep("Bolted",x)))]
# Code date for plants we can't see
phenos$D2Bolt[phenos[max(PhenoCols)]=="Vegetative"] <- 40
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

# Reshape Data for plots
phenos.melt <- melt(phenos[,c(1:3,PhenoCols)], 
                    id.vars = c("EUID", "Construct", "InfiltrationLocation"))
phenos.melt$value <- factor(phenos.melt$value, levels=phenotypes)
colnames(phenos.melt) <- c("EUID", "Construct", "InfiltrationLocation", "Date", "Stage")
# Reshape again to summarize
phenos.cast <- phenos.melt %>% count(Construct, InfiltrationLocation, Date, Stage)

# Make a linear model of the data
lm.phenos1 <- lm(as.numeric(D2Bolt) ~ Construct*InfiltrationLocation,
                data = phenos)
lm.phenos2 <- lm(as.numeric(D2Bolt) ~ Construct,
                 data = phenos,
                 subset = InfiltrationLocation=="Rosette")
lm.phenos3 <- lm(n ~ Date*Stage,
                 data = phenos.cast,
                 subset = Construct == "euFULII")


summary(lm.phenos1)
summary(lm.phenos2)

anova(lm.phenos1)
anova(lm.phenos2)
anova(lm.phenos3)

TukeyHSD(aov(lm.phenos1))
TukeyHSD(aov(lm.phenos2))

model.tables(aov(lm.phenos1))


# Neg Binomial Models -----------------------------------------------------
library(MASS)
library(car)
library(emmeans)
#Create the model for each infiltration location separately
glm.phenos1 <- lapply(unique(phenos$InfiltrationLocation),
                     function(i) glm.nb(Branches_Nov22 ~ Construct,
                                        control=glm.control(maxit = 1000),
                                        data=phenos[phenos$InfiltrationLocation==i,]))
names(glm.phenos1) <- unique(phenos$InfiltrationLocation)
# Show ANOVA tables
nb.anova1 <- lapply(glm.phenos1, Anova, type="II", test="LR")

# Show mean separation letters
lapply(glm.phenos1,
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
        ggtitle(i) +
        labs(y="Percent") +
        theme(plot.title = element_text(hjust = 0.5)))

# Number of branches on Nov22 by construct, separated by infiltration location
Plot_2 <- lapply(as.character(unique(phenos$InfiltrationLocation)),
                 function(i) ggplot(phenos[phenos$InfiltrationLocation==i,],
                                    aes(x=Construct, y=Branches_Nov22)) +
                   geom_violin() +
                   geom_boxplot(width=0.1) +
                   ggtitle(paste(i, "Leaf Infiltration"),
                           subtitle = bquote("Likelihood Ratio Chi"^2*"="*.(round(nb.anova1[[i]]$`LR Chisq`,2))*", DF="*.(nb.anova1[[i]]$Df)*", p="*.(round(nb.anova1[[i]]$`Pr(>Chisq)`,3)))) +
                   labs(y="Branch Number") +
                   theme(plot.title = element_text(hjust = 0.5),
                         plot.subtitle = element_text(hjust=0.5)))
names(Plot_2) <- unique(phenos$InfiltrationLocation)

# Days to bolting by construct for Rosette leaf infiltration
Plot_3 <- ggplot(phenos[phenos$InfiltrationLocation=="Rosette",], 
                 aes(x=Construct, y=D2Bolt)) +
  labs(y="Days to Bolting") +
  geom_violin() +
  ggtitle("Rosette Leaf Infiltration",
          subtitle=bquote("F"[.(paste(anova(lm.phenos2)[,1], collapse = ","))]*"="*.(round(anova(lm.phenos2)[1,4],3))*", p="*.(round(anova(lm.phenos2)[1,5],3)))) +
  theme(plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5))
      


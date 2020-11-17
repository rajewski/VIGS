library(googlesheets4)
# Read in Data
phenos <- read_sheet("https://docs.google.com/spreadsheets/d/1skhTjnfVESJ_eSaFTbt7_D1w33bDZeJ4WndzfEznxPc/edit?usp=sharing",
                     sheet="Sorted")
# Change order 
phenotypes <- c("Vegetative",
                "Bolted",
                "Flower Buds",
                "Flower1",
                "Flower2",
                "Flower3",
                "Fruit")
phenos$Construct <- factor(phenos$Construct, 
                           levels = c("Wild-Type", "Empty Vector", "euFULI", "euFULII", "Abeer"))
phenos$InfiltrationLocation <- factor(phenos$InfiltrationLocation,
                                       levels = c("Rosette", "Cauline"))
phenos[,4:7] <-lapply(phenos[,4:7], FUN=function(i) factor(i, levels=phenotypes))

max(factor(phenos$Nov14, ordered=TRUE))
phenos$Nov14Mod <- ifelse(grepl("Flower*|Fruit", phenos$Nov14),"3", phenos$Nov14)

# Make a linear model of the data
lm.phenos1 <- lm(as.numeric(Nov14Mod) ~ Construct*InfiltrationLocation,
                data = phenos)
lm.phenos2 <- lm(as.numeric(Nov14) ~ Construct,
                 data = phenos,
                 subset = InfiltrationLocation == "Rosette")

summary(lm.phenos1)
summary(lm.phenos2)

anova(lm.phenos1)
anova(lm.phenos2)

TukeyHSD(aov(lm.phenos1))
TukeyHSD(aov(lm.phenos2))

model.tables(aov(lm.phenos1))

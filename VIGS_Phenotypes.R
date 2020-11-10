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
                "Flower3")
phenos$Construct <- factor(phenos$Construct, 
                           levels = c("Wild-Type", "Empty Vector", "euFULI", "euFULII", "Abeer"))
phenos$InfiltrationLocation <- factor(phenos$InfiltrationLocation,
                                       levels = c("Rosette", "Cauline"))
phenos$Nov6 <- factor(phenos$Nov6,
                             levels = phenotypes)
phenos$Nov9 <- factor(phenos$Nov9,
                             levels = phenotypes)

max(factor(phenos$Nov9, ordered=TRUE))
phenos$Nov9Mod <- ifelse(grepl("Flower*", phenos$Nov9),"3", phenos$Nov9)

# Make a linear model of the data
lm.phenos1 <- lm(as.numeric(Nov9Mod) ~ Construct*InfiltrationLocation,
                data = phenos)
lm.phenos2 <- lm(as.numeric(Nov9) ~ Construct,
                 data = phenos,
                 subset = InfiltrationLocation == "Rosette")

summary(lm.phenos1)
summary(lm.phenos2)

TukeyHSD(aov(lm.phenos1))
TukeyHSD(aov(lm.phenos2))

model.tables(aov(lm.phenos1))

# This code tests a hypothesis for the Imagia technical challenge
# Kim Phan - 31 August 2019


# Load libraries
library(bestglm)
library(car)
library(caret)
library(corrplot)
library(cowplot)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(glmnet)
library(glmulti)
library(Hmisc)
library(InformationValue)
library(LogisticDx)
library(MASS)
library(mice)
library(pROC)
library(randomForest)
library(scales)
library(tableone)
library(tidyr)


# Load data
rawData <- read.csv("./hcc-data.txt",
                    stringsAsFactors = FALSE,
                    na.string = "?",
                    header = FALSE)

seed <- 123
set.seed(seed)
sample_n(rawData, 5)


# Rename columns
names(rawData) <- c("male", "symptoms", "alcohol",
                    "hepBsAg", "hepBeAg", "hepBcAb", "hepCvAb",
                    "cirrhosis", "endemic", "smoking", "diabetes", "obesity",
                    "hemochromatosis", "hypertension", "renalInsufficiency", "hiv",
                    "nash", "varices", "splenomegaly", "portalHypertension", "pvt",
                    "metastasis", "hallmark",
                    "age",
                    "alcohol_g", "smoking_pack",
                    "performance", "encephalopathy", "ascites",
                    "inr", "afp", "hb", "mcv", "leuk", "plt", "albumin",
                    "bilirubin_total", "alt", "ast", "ggt", "alp", "protein_total",
                    "creatinine",
                    "nodule_no",
                    "nodule_dim", "bilirubin_direct", "iron", "o2sat", "ferritin",
                    "lives")

# Set variable data types
data <- rawData %>%
    mutate_all(as.numeric) %>%
    mutate_at(vars(male:hallmark, nodule_no), as.factor) %>%
    mutate_at(vars(age), as.integer) %>%
    mutate(performance = factor(performance, levels = c(0, 1, 2, 3, 4), ordered = TRUE),
           encephalopathy = factor(encephalopathy, levels = c(1, 2, 3), ordered = TRUE),
           ascites = factor(ascites, levels = c(1, 2, 3), ordered = TRUE)) %>%
    mutate_at(vars(lives), as.factor)

str(data)



# MISSING DATA HANDLING ----

# Survey level of missing data
na.patterns <- naclus(data)  # hierarchical cluster of NA
print(na.patterns)
plot(na.patterns)

featureMissing <- apply(data, 2, function(x) sum(is.na(x))/length(x))
featureMissing[featureMissing > 0.25]  # vars with > 25% missing data

data.slimF <- data %>%
    dplyr::select(-c(varices, alcohol_g, smoking_pack, bilirubin_direct,
                     iron, o2sat, ferritin))

obsMissing <- apply(data.slimF, 1, function(x) sum(is.na(x))/length(x))
obsMissing[obsMissing > 0.25]  # obs with > 25% missing data
drop <- which(obsMissing > 0.25)

data.slimFO <-  data.slimF %>%
    slice(-drop)

na.patterns <- naclus(data.slimFO)
print(na.patterns)
plot(na.patterns)

# Multiple imputation of missing data
ghost.imp <- mice(data = data.slimFO, maxit = 0, print = FALSE)
pred <- quickpred(data.slimFO, mincor = 0.35, minpuc = 0.35)

toImpute <- names(data.slimFO)[apply(data.slimFO,
                                     2,
                                     function(x) sum(is.na(x))/length(x)) > 0]
pred_modified <- pred
pred_modified[toImpute, "lives"] <- 1

data.mids <- mice(data.slimFO, m = 5, maxit = 50,
                  seed = seed,
                  predictorMatrix = pred_modified)
data.imp <- mice::complete(data.mids)

# Check convergence over iterations
plot(data.mids)

# Check imputed values in continuous variables
xyplot(data.mids, inr + afp + albumin + nodule_dim ~ lives, pch = 18, cex = 1)
xyplot(data.mids, bilirubin_total + protein_total + creatinine ~ lives, pch = 18, cex = 1)



# EXPLORATORY DATA ANALYSIS ----

# Peek at the data
contents(data.imp)
describe(data.imp)

# Survival rate
table(data.imp$lives)
prop.table(table(data.imp$lives))

ggplot(data.imp, aes(lives)) +
    geom_bar(aes(y = (..count..)/sum(..count..), fill = lives)) +
    theme_light() +
    theme(legend.position = "nobe") +
    ggthemes::scale_fill_tableau(name = "Lives") +
    labs(x = "Lives",
         y = "Proportion")

# Review variable distriutions
varsPred <- names(data.imp)[-which(names(data.imp) == "lives")]
plots.univar <- list()
plots.bivar <- list()

for (i in seq_along(varsPred)) {  # create distribution plots
    if (is.numeric(data.imp[, varsPred[i]])) {  # coninuous variables
        plots.univar[[varsPred[i]]] <- ggplot(data.imp, aes_string(y = varsPred[i])) +
            geom_boxplot(fill = "#BAB0AC") +
            coord_flip() +
            theme_light() +
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank()) +
            labs(x = "",
                 y = varsPred[i])
        plots.bivar[[varsPred[i]]] <- ggplot(data.imp, aes_string(x = "lives", y = varsPred[i],
                                                                  group = "lives", fill = "lives")) +
            geom_boxplot() +
            coord_flip() +
            stat_summary(fun.y = mean, geom = "point",
                         shape = 5, size = 2, color = "black") +
            theme_light() +
            theme(legend.position = "none") +
            ggthemes::scale_fill_tableau() +
            labs(x = "Lives",
                 y = varsPred[i])
    } else {  # categorical variables
        plots.univar[[varsPred[i]]] <- ggplot(data.imp, aes_string(x = varsPred[i])) +
            geom_bar(aes(y = ..prop.., group = 1, fill = factor(..x..))) +
            scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
            theme_light() +
            theme(legend.position = "none") +
            ggthemes::scale_fill_tableau(palette = "Classic Gray 5") +
            labs(y = "Proportion",
                 x = varsPred[i])
        plots.bivar[[varsPred[i]]] <- ggplot(data.imp, aes_string(x = varsPred[i],
                                                                  fill = "lives")) +
            geom_bar(position = "fill") +
            scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
            theme_light() +
            ggthemes::scale_fill_tableau(name = "Lives") +
            labs(y = "Proportion",
                 x = varsPred[i])
    }
}

# Plot clinical risk factors
cowplot::plot_grid(plots.univar[["age"]], plots.bivar[["age"]],
                   plots.univar[["male"]], plots.bivar[["male"]],
                   plots.univar[["endemic"]], plots.bivar[["endemic"]],
                   ncol = 2, axis = "tb", align = "v")

cowplot::plot_grid(plots.univar[["diabetes"]], plots.bivar[["diabetes"]],
                   plots.univar[["obesity"]], plots.bivar[["obesity"]],
                   plots.univar[["hypertension"]], plots.bivar[["hypertension"]],
                   ncol = 2, axis = "tb", align = "v")

cowplot::plot_grid(plots.univar[["alcohol"]], plots.bivar[["alcohol"]],
                   plots.univar[["smoking"]], plots.bivar[["smoking"]],
                   ncol = 2, axis = "tb", align = "v")

cowplot::plot_grid(plots.univar[["hepBsAg"]], plots.bivar[["hepBsAg"]],
                   plots.univar[["hepBeAg"]], plots.bivar[["hepBeAg"]],
                   plots.univar[["hepBcAb"]], plots.bivar[["hepBcAb"]],
                   plots.univar[["hepCvAb"]], plots.bivar[["hepCvAb"]],
                   ncol = 2, axis = "tb", align = "v")

cowplot::plot_grid(plots.univar[["symptoms"]], plots.bivar[["symptoms"]],
                   plots.univar[["cirrhosis"]], plots.bivar[["cirrhosis"]],
                   plots.univar[["nash"]], plots.bivar[["nash"]],
                   plots.univar[["hiv"]], plots.bivar[["hiv"]],
                   ncol = 2, axis = "tb", align = "v")

cowplot::plot_grid(plots.univar[["hemochromatosis"]], plots.bivar[["hemochromatosis"]],
                   plots.univar[["protein_total"]], plots.bivar[["protein_total"]],
                   plots.univar[["creatinine"]], plots.bivar[["creatinine"]],
                   plots.univar[["renalInsufficiency"]], plots.bivar[["renalInsufficiency"]],
                   ncol = 2, axis = "tb", align = "v")

cowplot::plot_grid(plots.univar[["mcv"]], plots.bivar[["mcv"]],
                   plots.univar[["leuk"]], plots.bivar[["leuk"]],
                   ncol = 2, axis = "tb", align = "v")

# Plot tumor status markers
cowplot::plot_grid(plots.univar[["nodule_no"]], plots.bivar[["nodule_no"]],
                   plots.univar[["nodule_dim"]], plots.bivar[["nodule_dim"]],
                   plots.univar[["metastasis"]], plots.bivar[["metastasis"]],
                   plots.univar[["hallmark"]], plots.bivar[["hallmark"]],
                   ncol = 2, axis = "tb", align = "v")


# Plot liver function markers
cowplot::plot_grid(plots.univar[["bilirubin_total"]], plots.bivar[["bilirubin_total"]],
                   plots.univar[["albumin"]], plots.bivar[["albumin"]],
                   plots.univar[["inr"]], plots.bivar[["inr"]],
                   ncol = 2, axis = "tb", align = "v")

cowplot::plot_grid(plots.univar[["ascites"]], plots.bivar[["ascites"]],
                   plots.univar[["encephalopathy"]], plots.bivar[["encephalopathy"]],
                   plots.univar[["hb"]], plots.bivar[["hb"]],
                   plots.univar[["plt"]], plots.bivar[["plt"]],
                   ncol = 2, axis = "tb", align = "v")

cowplot::plot_grid(plots.univar[["alt"]], plots.bivar[["alt"]],
                   plots.univar[["ast"]], plots.bivar[["ast"]],
                   plots.univar[["alp"]], plots.bivar[["alp"]],
                   ncol = 2, axis = "tb", align = "v")

cowplot::plot_grid(plots.univar[["splenomegaly"]], plots.bivar[["splenomegaly"]],
                   plots.univar[["portalHypertension"]], plots.bivar[["portalHypertension"]],
                   plots.univar[["pvt"]], plots.bivar[["pvt"]],
                   ncol = 2, axis = "tb", align = "v")

# Plot performance status
cowplot::plot_grid(plots.univar[["performance"]], plots.bivar[["performance"]],
                   ncol = 2, axis = "tb", align = "v")


# Plot molecular biomarkers
cowplot::plot_grid(plots.univar[["afp"]], plots.bivar[["afp"]],
                   plots.univar[["ggt"]], plots.bivar[["ggt"]],
                   ncol = 2, axis = "tb", align = "v")



# GROUP CHARACTERISTICS ----
varsClinical <- c("age", "male", "endemic",
                  "diabetes", "obesity", "hypertension",
                  "alcohol", "smoking",
                  "hepBsAg", "hepBeAg", "hepBcAb", "hepCvAb",
                  "symptoms", "cirrhosis", "nash", "hiv",
                  "hemochromatosis",
                  "protein_total", "creatinine", "renalInsufficiency",
                  "mcv", "leuk")
varsBCLC <- c("nodule_no", "nodule_dim", "metastasis", "hallmark",
              "bilirubin_total", "albumin", "inr",
              "hb", "plt", "ascites", "encephalopathy",
              "alt", "ast","alp",
              "splenomegaly", "portalHypertension", "pvt",
              "performance")
varsMarker <- c("afp", "ggt")

groupData <- CreateTableOne(vars = c(varsClinical,
                                     varsBCLC,
                                     varsMarker),
                            data = data.imp,
                            factorVars = names(data.imp)[sapply(data.imp, is.factor)])

groupTable <- print(groupData,
                    catDigits = 1, contDigits = 1, pDigits = 3, missing = TRUE,
                    nonnormal = names(data)[sapply(data, is.numeric)])

comparisonData <- CreateTableOne(vars = c(varsClinical,
                                          varsBCLC,
                                          varsMarker),
                                 strata = "lives",
                                 data = data,
                                 factorVars = names(data)[sapply(data, is.factor)])

comparisonTable <- print(comparisonData,
                         catDigits = 1, contDigits = 1, pDigits = 3, missing = TRUE,
                         nonnormal = names(data)[sapply(data, is.numeric)])



# SPLIT DATASET ----

set.seed(seed)
trainIndex <- createDataPartition(data.imp$lives,
                                  p = .80,
                                  list = FALSE)
dataTrain <- data.imp[trainIndex, ]
dataTest <- data.imp[-trainIndex, ]

table(dataTrain$lives)
prop.table(table(dataTrain$lives))

table(dataTest$lives)
prop.table(table(dataTest$lives))



# BUILD THE CLINICAL MODEL ----

# Correlation matrix for clinical variables
corrplot(cor(dataTrain %>%
                 dplyr::select(one_of(varsClinical)) %>%
                 dplyr::select_if(is.numeric),
             use = "pairwise.complete.obs"),
         method = "number",
         type = "lower", order = "FPC", tl.col = "black", tl.srt = 45)

# Univariate logistic regressions
glmUnivarClinical <- cbind(varsClinical,
                           setNames(data.frame(matrix(nrow = length(varsClinical),
                                                      ncol = 3)),
                                    c("OR_", "CI_", "pval")))

for (i in seq_along(varsClinical)) {
    fit <- glm(as.formula(paste0("lives ~ ", varsClinical[i])),
               data = dataTrain,
               family = binomial)
    beta <- coef(fit)
    ci <- confint(fit)
    pval <- coef(summary(fit))[2, 4]
    glmUnivarClinical[i, 2] <- round(exp(beta[2]), 2)
    glmUnivarClinical[i, 3] <- paste0("(", round(exp(ci[2,1]),2), " - ", round(exp(ci[2,2]),2), ")")
    glmUnivarClinical[i, 4] <- round(pval,3)
}

print(glmUnivarClinical)

with(dataTrain, table(lives, nash))
with(dataTrain, table(lives, hiv))
with(dataTrain, table(lives, endemic))
with(dataTrain, table(lives, hepBeAg))

predClinical <- varsClinical[!varsClinical %in% c("nash", "hiv", "endemic","hepBeAg")]


# Fit the full model
dataClinical <- dataTrain %>%
    dplyr::select(one_of(predClinical, "lives"))

fullModel.Clinical <- glm(lives ~ .,
                          dataClinical,
                          family = binomial)

# Model via stepwise regression
modelClinicalA <- stepAIC(fullModel.Clinical,
                          direction = "both")

summary(modelClinicalA)

# Model via best subsets
dataClinical.Xy <- dataClinical %>%
    mutate_all(as.character) %>%
    mutate_all(as.numeric) %>%
    rename(y = lives) %>%
    dplyr::select(-c(mcv, protein_total, hepCvAb))  # max p = 15

modelClinical.bestglm <- bestglm(dataClinical.Xy,
                                 family = binomial,
                                 TopModels = 1,
                                 IC = "BIC",
                                 method = "exhaustive")

modelClinicalB <- modelClinical.bestglm$BestModel
summary(modelClinicalB)

modelClinicalB <- glm(lives ~ age + diabetes + hypertension + symptoms,
                      dataClinical,
                      family = "binomial")

# Random forest for feature selection
set.seed(seed)
rForestClinical <- randomForest(lives ~ ., dataClinical, importance = TRUE)
importance(rForestClinical)
varImpPlot(rForestClinical)

# Model via automated model selection and model-averaging
modelClinical.glmulti <- glmulti(lives ~ creatinine + hepBcAb + protein_total + leuk +
                                     diabetes + obesity + mcv + symptoms +
                                     hypertension + alcohol + hemochromatosis + renalInsufficiency,
                             data = dataClinical,
                             level = 1,
                             crit = "bic",
                             confsetsize = 1,
                             fitfunction = "glm",
                             family = binomial)
print(modelClinical.glmulti)

modelClinicalC <- glm(lives ~ diabetes + symptoms + hypertension + creatinine,
                      dataClinical,
                      family = "binomial")

summary(modelClinicalC)

# Multicollinearity
vif(modelClinicalA)
vif(modelClinicalB)
vif(modelClinicalC)

# Autocorrelation
durbinWatsonTest(modelClinicalA)
durbinWatsonTest(modelClinicalB)
durbinWatsonTest(modelClinicalC)

# Check goodness of fit with Hosmer lemeshow
gof(modelClinicalA, plotROC = FALSE)$gof[1]
gof(modelClinicalB, plotROC = FALSE)$gof[1]
gof(modelClinicalC, plotROC = FALSE)$gof[1]

# Review influential observations with Cook's Distance
cooksD.ClinicalA <- cooks.distance(modelClinicalA)
plot(cooksD.ClinicalA, pch = "*", cex = 1.5)
abline(h = 4*mean(cooksD.ClinicalA, na.rm = TRUE), col = "red")
influentialRows.ClinicalA <- which(cooksD.ClinicalA > 4*mean(cooksD.ClinicalA, na.rm = TRUE))
dataClinical %>%  # print influential observations
    dplyr::select(one_of(names(modelClinicalA$model)[-1])) %>%
    slice(influentialRows.ClinicalA)

cooksD.ClinicalB <- cooks.distance(modelClinicalB)
plot(cooksD.ClinicalB, pch = "*", cex = 1.5)
abline(h = 4*mean(cooksD.ClinicalB, na.rm = TRUE), col = "red")
influentialRows.ClinicalB <- which(cooksD.ClinicalB > 4*mean(cooksD.ClinicalB, na.rm = TRUE))
dataClinical %>%  # print influential observations
    dplyr::select(one_of(names(modelClinicalB$model)[-1])) %>%
    slice(influentialRows.ClinicalB)

cooksD.ClinicalC <- cooks.distance(modelClinicalC)
plot(cooksD.ClinicalC, pch = "*", cex = 1.5)
abline(h = 4*mean(cooksD.ClinicalC, na.rm = TRUE), col = "red")
influentialRows.ClinicalC <- which(cooksD.ClinicalC > 4*mean(cooksD.ClinicalC, na.rm = TRUE))
dataClinical %>%  # print influential observations
    dplyr::select(one_of(names(modelClinicalC$model)[-1])) %>%
    slice(influentialRows.ClinicalC)

# Accuracy
probTrainClinicalA <- predict(modelClinicalA, type = "response")
modelClinicalA.cut <- optimalCutoff(dataClinical$lives, probTrainClinicalA)
1 - misClassError(dataClinical$lives, probTrainClinicalA, threshold = modelClinicalA.cut)

probTrainClinicalB <- predict(modelClinicalB, type = "response")
modelClinicalB.cut <- optimalCutoff(dataClinical$lives, probTrainClinicalB)
1 - misClassError(dataClinical$lives, probTrainClinicalB, threshold = modelClinicalB.cut)

probTrainClinicalC <- predict(modelClinicalC, type = "response")
modelClinicalC.cut <- optimalCutoff(dataClinical$lives, probTrainClinicalC)
1 - misClassError(dataClinical$lives, probTrainClinicalC, threshold = modelClinicalC.cut)

# AUC
gof(modelClinicalA)$auc
gof(modelClinicalB)$auc
gof(modelClinicalC)$auc

# PPV/NPV
precision(dataClinical$lives, probTrainClinicalA, threshold = modelClinicalA.cut)
npv(dataClinical$lives, probTrainClinicalA, threshold = modelClinicalA.cut)

precision(dataClinical$lives, probTrainClinicalB, threshold = modelClinicalB.cut)
npv(dataClinical$lives, probTrainClinicalB, threshold = modelClinicalB.cut)

precision(dataClinical$lives, probTrainClinicalC, threshold = modelClinicalC.cut)
npv(dataClinical$lives, probTrainClinicalC, threshold = modelClinicalC.cut)

# Best Clinical model
modelClinical <- modelClinicalB



# BUILD THE TUMOR/LIVER MODEL ----

# Correlation matrix for BCLC variables
corrplot(cor(dataTrain %>%
                 dplyr::select(one_of(varsBCLC)) %>%
                 dplyr::select_if(is.numeric),
             use = "pairwise.complete.obs"),
         method = "number",
         type = "lower", order = "FPC", tl.col = "black", tl.srt = 45)

# Univariate logistic regressions
glmUnivarBCLC <- cbind(varsBCLC,
                       setNames(data.frame(matrix(nrow = length(varsBCLC),
                                                  ncol = 3)),
                                c("OR_", "CI_", "pval")))

for (i in seq_along(varsBCLC)) {
    fit <- glm(as.formula(paste0("lives ~ ", varsBCLC[i])),
               data = dataTrain,
               family = binomial)
    beta <- coef(fit)
    ci <- confint(fit)
    pval <- coef(summary(fit))[2, 4]
    glmUnivarBCLC[i, 2] <- round(exp(beta[2]), 2)
    glmUnivarBCLC[i, 3] <- paste0("(", round(exp(ci[2,1]),2), " - ", round(exp(ci[2,2]),2), ")")
    glmUnivarBCLC[i, 4] <- round(pval,3)
}

print(glmUnivarBCLC)

with(dataTrain, table(lives, encephalopathy))

predBCLC <- varsBCLC[!varsBCLC %in% c("encephalopathy")]

# Fit the full model
dataBCLC <- dataTrain %>%
    dplyr::select(one_of(predBCLC, "lives"))

fullModel.BCLC <- glm(lives ~ ., dataBCLC,
                      family = binomial)
summary(fullModel.BCLC)

# Model via stepwise regression
modelBCLCA <- stepAIC(fullModel.BCLC,
                      direction = "both")

summary(modelBCLCA)

# Model via best subsets
dataBCLC.Xy <- dataBCLC %>%
    mutate_all(as.character) %>%
    mutate_all(as.numeric) %>%
    dplyr::select(-c(performance, alp)) %>%  # max p = 15
    rename(y = lives)

modelBCLC.bestglm <- bestglm(dataBCLC.Xy,
                             family = binomial, TopModels = 1,
                             IC = "BIC", method = "exhaustive")

modelBCLCB <- modelBCLC.bestglm$BestModel

summary(modelBCLCB)

modelBCLCB <- glm(lives ~ nodule_dim + bilirubin_total + hb + ascites,
                 dataBCLC,
                 family = "binomial")

# Random forest method
set.seed(seed)
rForestBCLC <- randomForest(lives ~ ., dataBCLC, importance = TRUE)
importance(rForestBCLC)
varImpPlot(rForestBCLC)

# Model via automated model selection and model-averaging
modelBCLC.glmulti <- glmulti(lives ~ hb + performance + albumin + alp +
                                 nodule_dim + metastasis + plt + ast +
                                 ascites + bilirubin_total + alt + pvt,
                             data = dataBCLC,
                             level = 1,
                             crit = "bic",
                             confsetsize = 1,
                             fitfunction = "glm",
                             family = binomial)
print(modelBCLC.glmulti)

modelBCLCC <- glm(lives ~ metastasis + hb + bilirubin_total,
                  data = dataBCLC,
                  family = binomial)
summary(modelBCLCC)


# Multicollinearity
vif(modelBCLCA)
vif(modelBCLCB)
vif(modelBCLCC)

# Autocorrelation
durbinWatsonTest(modelBCLCA)
durbinWatsonTest(modelBCLCB)
durbinWatsonTest(modelBCLCC)

# Check goodness of fit with Hosmer lemeshow
gof(modelBCLCA, plotROC = FALSE)$gof[1]
gof(modelBCLCB, plotROC = FALSE)$gof[1]
gof(modelBCLCC, plotROC = FALSE)$gof[1]

# Review influential observations with Cook's Distance
cooksD.BCLCA <- cooks.distance(modelBCLCA)
plot(cooksD.BCLCA, pch = "*", cex = 1.5)
abline(h = 4*mean(cooksD.BCLCA, na.rm = TRUE), col = "red")
influentialRows.BCLCA <- which(cooksD.BCLCA > 4*mean(cooksD.BCLCA, na.rm = TRUE))
dataBCLC %>%  # print influential observations
    dplyr::select(one_of(names(modelBCLCA$model)[-1])) %>%
    slice(influentialRows.BCLCA)

cooksD.BCLCB <- cooks.distance(modelBCLCB)
plot(cooksD.BCLCB, pch = "*", cex = 1.5)
abline(h = 4*mean(cooksD.BCLCB, na.rm = TRUE), col = "red")
influentialRows.BCLCB <- which(cooksD.BCLCB > 4*mean(cooksD.BCLCB, na.rm = TRUE))
dataBCLC %>%  # print influential observations
    dplyr::select(one_of(names(modelBCLCB$model)[-1])) %>%
    slice(influentialRows.BCLCB)

cooksD.BCLCC <- cooks.distance(modelBCLCC)
plot(cooksD.BCLCC, pch = "*", cex = 1.5)
abline(h = 4*mean(cooksD.BCLCC, na.rm = TRUE), col = "red")
influentialRows.BCLCC <- which(cooksD.BCLCC > 4*mean(cooksD.BCLCC, na.rm = TRUE))
dataBCLC %>%  # print influential observations
    dplyr::select(one_of(names(modelBCLCC$model)[-1])) %>%
    slice(influentialRows.BCLCC)

# Accuracy
probTrainBCLCA <- predict(modelBCLCA, type = "response")
modelBCLCA.cut <- optimalCutoff(dataBCLC$lives, probTrainBCLCA)
1 - misClassError(dataBCLC$lives, probTrainBCLCA, threshold = modelBCLCA.cut)

probTrainBCLCB <- predict(modelBCLCB, type = "response")
modelBCLCB.cut <- optimalCutoff(dataBCLC$lives, probTrainBCLCB)
1 - misClassError(dataBCLC$lives, probTrainBCLCB, threshold = modelBCLCB.cut)

probTrainBCLCC <- predict(modelBCLCC, type = "response")
modelBCLCC.cut <- optimalCutoff(dataBCLC$lives, probTrainBCLCC)
1 - misClassError(dataBCLC$lives, probTrainBCLCC, threshold = modelBCLCC.cut)

# AUC
gof(modelBCLCA)$auc
gof(modelBCLCB)$auc
gof(modelBCLCC)$auc

# PPV/NPV
precision(dataBCLC$lives, probTrainBCLCA, threshold = modelBCLCA.cut)
npv(dataBCLC$lives, probTrainBCLCA, threshold = modelBCLCA.cut)

precision(dataBCLC$lives, probTrainBCLCB, threshold = modelBCLCB.cut)
npv(dataBCLC$lives, probTrainBCLCB, threshold = modelBCLCB.cut)

precision(dataBCLC$lives, probTrainBCLCC, threshold = modelBCLCC.cut)
npv(dataBCLC$lives, probTrainBCLCC, threshold = modelBCLCC.cut)

# Best BCLC model
modelBCLC <- modelBCLCB



# BUILD THE BIOMARKER MODEL ----

# Correlation matrix for clinical variables
with(dataTrain, cor(afp, ggt, use = "pairwise.complete.obs"))


# Univariate logistic regressions

glmUnivarMarker <- cbind(varsMarker,
                         setNames(data.frame(matrix(nrow = length(varsMarker),
                                                    ncol = 3)),
                                  c("OR_", "CI_", "pval")))

for (i in seq_along(varsMarker)) {
    fit <- glm(as.formula(paste0("lives ~ ", varsMarker[i])),
               data = dataTrain,
               family = binomial)
    beta <- coef(fit)
    ci <- confint(fit)
    pval <- coef(summary(fit))[2, 4]
    glmUnivarMarker[i, 2] <- round(exp(beta[2]), 2)
    glmUnivarMarker[i, 3] <- paste0("(", round(exp(ci[2,1]),2), " - ", round(exp(ci[2,2]),2), ")")
    glmUnivarMarker[i, 4] <- round(pval,3)
}

print(glmUnivarMarker)

predMarker <- varsMarker



# Fit models
dataMarker <- dataTrain %>%
    dplyr::select(one_of(predMarker, "lives"))

modelMarkerA <- glm(lives ~ afp,
                    data = dataMarker,
                    family = binomial)
modelMarkerB <- glm(lives ~ ggt,
                    data = dataMarker,
                    family = binomial)
modelMarkerC <- glm(lives ~ afp + ggt,
                    data = dataMarker,
                    family = binomial)

# Multicollinearity
vif(modelMarkerC)

# Autocorrelation
durbinWatsonTest(modelMarkerA)
durbinWatsonTest(modelMarkerB)
durbinWatsonTest(modelMarkerC)

# Check goodness of fit with Hosmer lemeshow
gof(modelMarkerA, plotROC = FALSE)$gof[1]
gof(modelMarkerB, plotROC = FALSE)$gof[1]
gof(modelMarkerC, plotROC = FALSE)$gof[1]

# Review influential observations with Cook's Distance
cooksD.MarkerA <- cooks.distance(modelMarkerA)
plot(cooksD.MarkerA, pch = "*", cex = 1.5)
abline(h = 4*mean(cooksD.MarkerA, na.rm = TRUE), col = "red")
influentialRows.MarkerA <- which(cooksD.MarkerA > 4*mean(cooksD.MarkerA, na.rm = TRUE))
dataMarker %>%  # print influential observations
    dplyr::select(one_of(names(modelMarkerA$model)[-1])) %>%
    slice(influentialRows.MarkerA)

cooksD.MarkerB <- cooks.distance(modelMarkerB)
plot(cooksD.MarkerB, pch = "*", cex = 1.5)
abline(h = 4*mean(cooksD.MarkerB, na.rm = TRUE), col = "red")
influentialRows.MarkerB <- which(cooksD.MarkerB > 4*mean(cooksD.MarkerB, na.rm = TRUE))
dataMarker %>%  # print influential observations
    dplyr::select(one_of(names(modelMarkerB$model)[-1])) %>%
    slice(influentialRows.MarkerB)

cooksD.MarkerC <- cooks.distance(modelMarkerC)
plot(cooksD.MarkerC, pch = "*", cex = 1.5)
abline(h = 4*mean(cooksD.MarkerC, na.rm = TRUE), col = "red")
influentialRows.MarkerC <- which(cooksD.MarkerB > 4*mean(cooksD.MarkerB, na.rm = TRUE))
dataMarker %>%  # print influential observations
    dplyr::select(one_of(names(modelMarkerC$model)[-1])) %>%
    slice(influentialRows.MarkerC)

# Accuracy
probTrainMarkerA <- predict(modelMarkerA, type = "response")
modelMarkerA.cut <- optimalCutoff(dataMarker$lives, probTrainMarkerA)
1 - misClassError(dataMarker$lives, probTrainMarkerA, threshold = modelMarkerA.cut)

probTrainMarkerB <- predict(modelMarkerB, type = "response")
modelMarkerB.cut <- optimalCutoff(dataMarker$lives, probTrainMarkerB)
1 - misClassError(dataMarker$lives, probTrainMarkerB, threshold = modelMarkerB.cut)

probTrainMarkerC <- predict(modelMarkerC, type = "response")
modelMarkerC.cut <- optimalCutoff(dataMarker$lives, probTrainMarkerC)
1 - misClassError(dataMarker$lives, probTrainMarkerC, threshold = modelMarkerC.cut)

# AUC
gof(modelMarkerA)$auc
gof(modelMarkerB)$auc
gof(modelMarkerC)$auc

# PPV/NPV
precision(dataMarker$lives, probTrainMarkerA, threshold = modelMarkerA.cut)
npv(dataMarker$lives, probTrainMarkerA, threshold = modelMarkerA.cut)

precision(dataMarker$lives, probTrainMarkerB, threshold = modelMarkerB.cut)
npv(dataMarker$lives, probTrainMarkerB, threshold = modelMarkerB.cut)

precision(dataMarker$lives, probTrainMarkerC, threshold = modelMarkerC.cut)
npv(dataMarker$lives, probTrainMarkerC, threshold = modelMarkerC.cut)

# Best Marker model
modelMarker <- modelMarkerA



# EVALUATE ADDITIVE VALUE OF TUMOR/LIVER & MOLECULAR BIOMARKER MODELS ----

# Recall individual models
modelClinical$formula
modelBCLC$formula
modelMarker$formula

# Build added models
probTrainClinical <- predict(modelClinical, type = "response")

modelCB <- glm(lives ~ age + diabetes + hypertension + symptoms + creatinine + leuk +
                   nodule_dim + hb + ascites + bilirubin_total,
               data = dataTrain,
               family = binomial)
probTrainCB <- predict(modelCB, type = "response")

summary(modelCB)

modelCBM <- glm(lives ~ age + diabetes + hypertension + symptoms + creatinine + leuk +
                    nodule_dim + hb + ascites + bilirubin_total +
                    afp,
                data = dataTrain,
                family = binomial)
probTrainCBM <- predict(modelCBM, type = "response")

summary(modelCBM)

# Likelihood ratio test
anova(modelClinical, modelCB, test ="Chisq")
anova(modelCB, modelCBM, test ="Chisq")

# Compare ROC curves
rocClinical <- roc(dataClinical$lives, probTrainClinical)
rocCB <- roc(dataTrain$lives, probTrainCB)
rocCBM <- roc(dataTrain$lives, probTrainCBM)

gof(modelClinical, plotROC = FALSE)$auc
gof(modelCB, plotROC = FALSE)$auc
gof(modelCBM, plotROC = FALSE)$auc

ggroc(data = list(C = rocClinical,
                  CB = rocCB,
                  CBM = rocCBM),
      size = 1,
      aes = c("color")) +
    theme_light() +
    ggthemes::scale_color_tableau(name = "Model")

roc.test(rocClinical, rocCB)
roc.test(rocCB, rocCBM)
roc.test(rocClinical, rocCBM)

# Accuracy
modelClinical.cut <- optimalCutoff(dataTrain$lives, probTrainClinical)
1 - misClassError(dataTrain$lives, probTrainClinical, threshold = modelClinical.cut)

modelCB.cut <- optimalCutoff(dataTrain$lives, probTrainCB)
1 - misClassError(dataTrain$lives, probTrainCB, threshold = modelCB.cut)

modelCBM.cut <- optimalCutoff(dataTrain$lives, probTrainCBM)
1 - misClassError(dataTrain$lives, probTrainCBM, threshold = modelCBM.cut)

# PPV/NPV
precision(dataTrain$lives, probTrainClinical, threshold = modelClinical.cut)
npv(dataTrain$lives, probTrainClinical, threshold = modelClinical.cut)

precision(dataTrain$lives, probTrainCB, threshold = modelCB.cut)
npv(dataTrain$lives, probTrainCB, threshold = modelCB.cut)

precision(dataTrain$lives, probTrainCBM, threshold = modelCBM.cut)
npv(dataTrain$lives, probTrainCBM, threshold = modelCBM.cut)



# INTERNALLY VALIDATE ON TEST DATA ----

# Make model predictions on test data
probTestClinical <- predict(modelClinical, type = "response", newdata = dataTest)
probTestCB <- predict(modelCB, type = "response", newdata = dataTest)
probTestCBM <- predict(modelCBM, type = "response", newdata = dataTest)

# Compare ROC curves
rocTestClinical <- roc(dataTest$lives, probTestClinical)
rocTestCB <- roc(dataTest$lives, probTestCB)
rocTestCBM <- roc(dataTest$lives, probTestCBM)

rocTestClinical$auc
rocTestCB$auc
rocTestCBM$auc

ggroc(data = list(C = rocTestClinical,
                  CB = rocTestCB,
                  CBM = rocTestCBM),
      size = 1,
      aes = c("color", "linetype")) +
    theme_light() +
    scale_linetype(name = "Model") +
    ggthemes::scale_color_tableau(name = "Model")

roc.test(rocTestClinical, rocTestCB)
roc.test(rocTestCB, rocTestCBM)
roc.test(rocTestClinical, rocTestCBM)

# Accuracy
modelClinical.testCut <- optimalCutoff(dataTest$lives, probTestClinical)
1 - misClassError(dataTest$lives, probTestClinical, threshold = modelClinical.testCut)

modelCB.testCut <- optimalCutoff(dataTest$lives, probTestCB)
1 - misClassError(dataTest$lives, probTestCB, threshold = modelCB.testCut)

modelCBM.testCut <- optimalCutoff(dataTest$lives, probTestCBM)
1 - misClassError(dataTest$lives, probTestCBM, threshold = modelCBM.testCut)

# PPV/NPV
precision(dataTest$lives, probTestClinical, threshold = modelClinical.testCut)
npv(dataTest$lives, probTestClinical, threshold = modelClinical.testCut)

precision(dataTest$lives, probTestCB, threshold = modelCB.testCut)
npv(dataTest$lives, probTestCB, threshold = modelCB.testCut)

precision(dataTest$lives, probTestCBM, threshold = modelCBM.testCut)
npv(dataTest$lives, probTestCBM, threshold = modelCBM.testCut)


# K-FOLD CROSS VALIDATION ----

ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     savePredictions = TRUE)

modelCV.C <- train(lives ~ age + diabetes + hypertension + symptoms + creatinine + leuk,
                   data = data.imp,
                   method = "glm",
                   family = "binomial",
                   trControl = ctrl,
                   tuneLength = 5)
modelCV.C$results

predictionsCV.C <- predict(modelCV.C, data.imp)
confusionC <- caret::confusionMatrix(predictionsCV.C, data.imp$lives, positive = "1")
print(confusionC)

modelCV.CB <- train(lives ~ age + diabetes + hypertension + symptoms + creatinine + leuk +
                        nodule_dim + hb + ascites + bilirubin_total,
                    data = data.imp,
                    method = "glm",
                    family = "binomial",
                    trControl = ctrl,
                    tuneLength = 5)

modelCV.CB$results

predictionsCV.CB <- predict(modelCV.CB, data.imp)
confusionCB <- caret::confusionMatrix(predictionsCV.CB, data.imp$lives, positive = "1")
print(confusionCB)

modelCV.CBM <- train(lives ~ age + diabetes + hypertension + symptoms + creatinine + leuk +
                         nodule_dim + hb + ascites + bilirubin_total + afp,
                     data = data.imp,
                     method = "glm",
                     family = "binomial",
                     trControl = ctrl,
                     tuneLength = 5)

modelCV.CBM$results

predictionsCV.CBM <- predict(modelCV.CBM, data.imp)
confusionCBM <- caret::confusionMatrix(predictionsCV.CBM, data.imp$lives, positive = "1")
print(confusionCBM)


modelCB.imp <- glm(lives ~ age + diabetes + hypertension + symptoms + creatinine +
                       leuk + nodule_dim + hb + ascites + bilirubin_total,
                   data.imp,
                   family = "binomial")
summary(modelCB.imp)

exp(cbind(coef(modelCB.imp), confint(modelCB.imp)))



# NOMOGRAM ----

# Nomogram function
## Create simple Fagan nomograms as ggplot objects
##   Based on Perl web-implementation (https://araw.mede.uic.edu/cgi-bin/testcalc.pl)
##   Authors: AM. Chekroud* & A. Schwartz (* adam dot chekroud at yale . edu)
##   December 2016

nomogrammer <- function(Prevalence,
                        Sens = NULL,
                        Spec = NULL,
                        Plr = NULL,
                        Nlr = NULL,
                        Detail = FALSE,
                        NullLine = FALSE,
                        LabelSize = (14/5),
                        Verbose = FALSE){

    ## Function inputs:
    # Prevalence (prior probability) as a number between 0 and 1
    # Either
    # Sens & Spec
    # model sensitivity and specificity as a number between 0 and 1
    # Or
    # Likelihood ratios
    # Positive and Negative LRs (numeric)

    ## Function options:
    # Detail: If true, will overlay key statistics onto the plot
    # NullLine: If true, will add a line from prior prob through LR = 1
    # LabelSize: Tweak this number to change the label sizes
    # Verbose: Print out relevant metrics in the console

    ## Function returns:
    # ggplot object



    ######################################
    ########## Libraries & Functions #####
    ######################################

    ## Libraries
    require(ggplot2)
    require(scales)






    ## Helper functions
    ##   (defined inside nomogrammer, so remain local only & wont clutter user env)
    odds         <- function(p){
        # Function converts probability into odds
        o <- p/(1-p)
        return(o)
    }

    logodds      <- function(p){
        # Function returns logodds for a probability
        lo <- log10(p/(1-p))
        return(lo)
    }

    logodds_to_p <- function(lo){
        # Function goes from logodds back to a probability
        o <- 10^lo
        p <- o/(1+o)
        return(p)
    }

    p2percent <- function(p){
        # Function turns numeric probability into string percentage
        # e.g. 0.6346111 -> 63.5%
        scales::percent(signif(p, digits = 3))}


    ######################################
    ########## Calculations     ##########
    ######################################

    ## Checking inputs

    ## Prevalence
    # needs to exist
    if(missing(Prevalence)){
        stop("Prevalence is missing")
    }
    # needs to be numeric
    if(!is.numeric(Prevalence)){stop("Prevalence should be numeric")}
    # needs to be a prob not a percent
    if((Prevalence > 1) | (Prevalence <= 0)){stop("Prevalence should be a probability (did you give a %?)")}

    # Did user give sens & spec?
    if(missing(Sens) | missing(Spec)){
        sensspec <- FALSE
    } else{ sensspec <- TRUE}
    # if yes, make sure they are numbers
    if(sensspec == TRUE){
        if(!is.numeric(Sens)){stop("Sensitivity should be numeric")}
        if(!is.numeric(Spec)){stop("Specificity should be numeric")}
        # numbers that are probabilities not percentages
        if((Sens > 1) | (Sens <= 0)){stop("Sensitivity should be a probability (did you give a %?)")}
        if((Spec > 1) | (Spec <= 0)){stop("Specificity should be a probability (did you give a %?)")}
    }


    # Did user give PLR & NLR?
    if(missing(Plr) | missing(Nlr)){
        plrnlr <- FALSE
    } else{plrnlr <- TRUE}
    # if yes, make sure they are numbers
    if(plrnlr == TRUE){
        if(!is.numeric(Plr)){stop("PLR should be numeric")}
        if(!is.numeric(Nlr)){stop("NLR should be numeric")}
        # numbers that vaguely make sense
        if(Plr < 1){stop("PLR shouldn't be less than 1")}
        if(Nlr < 0){stop("NLR shouldn't be below zero")}
        if(Nlr > 1){stop("NLR shouldn't be more than 1")}
    }

    # Did they give a valid sensspec and plrnlr? If yes, ignore the LRs and tell them
    if((sensspec == TRUE) && (plrnlr == TRUE) ){
        warning("You provided sens/spec as well as likelihood ratios-- I ignored the LRs!")
    }


    ## If sens/spec provided, we calculate posterior probabilities & odds using sens & spec
    ##  otherwise, if plr and nlr provided, we calculate posteriors using them
    ##  if neither exist, then return an error
    if(sensspec == TRUE){
        prior_prob  <- Prevalence
        prior_odds  <- odds(prior_prob)
        sensitivity <- Sens
        specificity <- Spec
        PLR <- sensitivity/(1-specificity)
        NLR <- (1-sensitivity)/specificity
        post_odds_pos  <- prior_odds * PLR
        post_odds_neg  <- prior_odds * NLR
        post_prob_pos  <- post_odds_pos/(1+post_odds_pos)
        post_prob_neg  <- post_odds_neg/(1+post_odds_neg)
    } else if(plrnlr == TRUE){
        prior_prob  <- Prevalence
        prior_odds  <- odds(prior_prob)
        PLR <- Plr
        NLR <- Nlr
        sensitivity <- (PLR*(1-NLR))/(PLR-NLR)    ## TODO: check Adam's math!
        specificity <- (1-PLR)/(NLR-PLR)          ## TODO: check Adam's math!
        post_odds_pos  <- prior_odds * PLR
        post_odds_neg  <- prior_odds * NLR
        post_prob_pos  <- post_odds_pos/(1+post_odds_pos)
        post_prob_neg  <- post_odds_neg/(1+post_odds_neg)
    } else{
        stop("Couldn't find sens & spec, or positive & negative likelihood ratios")
    }



    ######################################
    ########## Plotting (prep)  ##########
    ######################################


    ## Set common theme preferences up front
    theme_set(theme_bw() +
                  theme(axis.text.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.title.x = element_blank(),
                        axis.title.y = element_text(angle = 0),
                        axis.title.y.right = element_text(angle = 0),
                        axis.line = element_blank(),
                        panel.grid = element_blank(),
                        legend.position = "none"
                  )
    )

    ## Setting up the points of interest along the y-axes

    # Select probabilities of interest (nb as percentages)
    ticks_prob    <- c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 30,
                       40, 50, 60, 70, 80, 90, 95, 99)
    # Convert % to odds
    ticks_odds    <- odds(ticks_prob/100)
    # Convert % to logodds
    ticks_logodds <- logodds(ticks_prob/100)

    # Select the likelihood ratios of interest (for the middle y-axis)
    ticks_lrs     <- sort(c(10^(-3:3), 2*(10^(-3:2)), 5*(10^(-3:2))))
    # Log10 them since plot is in logodds space
    ticks_log_lrs <- log10(ticks_lrs)




    ## Fixing particular x-coordinates
    left     <- 0
    right    <- 1
    middle   <- 0.5
    midright <- 0.75

    ## Lay out the four key plot points
    ##  (the start and finish of the positive and negative lines)

    # Initially these are expressed as probabilities
    df <- data.frame(x=c(left, right, left, right),
                     y=c(prior_prob, post_prob_pos, prior_prob, post_prob_neg),
                     line = c("pos", "pos", "neg", "neg"))

    adj_min      <- range(ticks_logodds)[1]
    adj_max      <- range(ticks_logodds)[2]
    adj_diff     <- adj_max - adj_min
    scale_factor <- abs(adj_min) - adj_diff/2
    #df$lo_y <- ifelse(df$x==left,(10/adj_diff)*logodds(1-df$y)-1,logodds(df$y))

    # Convert probabilities to logodds for plotting
    df$lo_y  <- ifelse(df$x==left,logodds(1-df$y)-scale_factor,logodds(df$y))
    # zero         <- data.frame(x = c(left,right),
    #                            y = c(0,0),
    #                            line = c('pos','pos'),
    #                            lo_y = c(-scale_factor,0))





    rescale   <- range(ticks_logodds) + abs(adj_min) - adj_diff/2
    rescale_x_breaks  <- ticks_logodds + abs(adj_min) - adj_diff/2



    ######################################
    ########## Plot             ##########
    ######################################


    p <- ggplot(df) +
        geom_line(aes(x = x, y = lo_y, color = line), size = 1) +
        geom_vline(xintercept = middle) +
        annotate(geom = "text",
                 x = rep(middle+.075, length(ticks_log_lrs)),
                 y = (ticks_log_lrs-scale_factor)/2,
                 label = ticks_lrs,
                 size = rel(LabelSize)) +
        annotate(geom="point",
                 x = rep(middle, length(ticks_log_lrs)),
                 y = (ticks_log_lrs-scale_factor)/2,
                 size = 1) +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0),
                           limits = rescale,
                           breaks = -rescale_x_breaks,
                           labels = ticks_prob,
                           name = "prior \n prob.",
                           sec.axis = sec_axis(trans = ~.,
                                               name = "posterior \n prob.",
                                               labels = ticks_prob,
                                               breaks = ticks_logodds))

    ## Optional overlay text: prevalence, PLR/NLR, and posterior probabilities
    detailedAnnotation <- paste(
        paste0("prevalence = ", p2percent(prior_prob)),
        paste("PLR =", signif(PLR, 3),", NLR =", signif(NLR, 3)),
        paste("post. pos =", p2percent(post_prob_pos),
              ", neg =", p2percent(post_prob_neg)),
        sep = "\n")


    ## Optional amendments to the plot

    ## Do we add the null line i.e. LR = 1, illustrating an uninformative model
    if(NullLine == TRUE){
        ## If yes, first calculate the start and end points
        uninformative <- data.frame(
            x = c(left,right),
            lo_y = c( (logodds(1-prior_prob) - scale_factor) , logodds(prior_prob))
        )

        p <- p + geom_line(aes(x = x, y = lo_y), data = uninformative,
                           color = "gray",
                           lty = 2,
                           inherit.aes = FALSE)
    }


    ## Do we add the detailed stats to the top right?
    if(Detail == TRUE){
        p <- p + annotate(geom = "text",
                          x = midright,
                          y = 2,
                          label = detailedAnnotation,
                          size = rel(LabelSize))
    }

    if(Verbose == TRUE){
        writeLines(
            text = c(
                paste0("prevalence = ", p2percent(prior_prob)),
                paste("PLR =", signif(PLR, 3)),
                paste("NLR =", signif(NLR, 3)),
                paste("posterior probability (positive) =", p2percent(post_prob_pos)),
                paste("posterior probability (negative) =", p2percent(post_prob_neg)),
                paste("sensitivity =", p2percent(sensitivity)),
                paste("specificity =", p2percent(specificity))
                # sep = "\n"
            )
        )
    }


    return(p)

}


nomogrammer(Prevalence = prop.table(table(data.imp$lives))[2],
            Sens = confusionC$byClass[1], Spec = confusionC$byClass[2])

nomogrammer(Prevalence = prop.table(table(data.imp$lives))[2],
            Sens = confusionCB$byClass[1], Spec = confusionCB$byClass[2])

nomogrammer(Prevalence = prop.table(table(data.imp$lives))[2],
            Sens = confusionCBM$byClass[1], Spec = confusionCBM$byClass[2])


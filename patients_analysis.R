# INITIALIZATION =====

### Load libraries
library(dplyr)
library(tidyr)
library(janitor)
library(lubridate)
library(ggplot2)
library(SmartEDA)
library(ggstatsplot)
library(GGally)
library(viridis)
library(car)


### Load data
wd <- getwd()
patients_original <- read.csv(paste(wd, "/Patients_HCPAnywhere.csv", sep = ""), header = TRUE, sep = ",")
diagnostics <- read.csv(paste(wd, "/HCPAnywhere_Diagnostics.csv", sep = ""), header = TRUE, sep = ",")
plomb <- read.csv(paste(wd, "/HCPAnywhere_Plomb_Data.csv", sep = ""), header = TRUE, sep = ",")

# NETTOYAGE DONNEES =====
# Elimination de variables protégées/sensibles sur les patients
patients_original<-patients_original %>%
  select (-housenum, -streetname, -rownames)

diagnostics<-diagnostics %>%
  mutate(temps_suivi_jours = as.integer(difftime(last_activity, inclusion, units = "days"))) %>%
  mutate(year_last_activity = year(last_activity)) %>%
  select(-inclusion, -last_activity, -rownames) 

# Identifier et Eliminer les duplicats
check_duplicated_idmr <- function(df) {
  duplicates <- df %>%
    group_by(idmr) %>%
    filter(n() > 1)
  return(duplicates)
}

check_duplicated_idmr(patients_original)
check_duplicated_idmr(diagnostics)
check_duplicated_idmr(plomb)

# Supprimmer les lignes avec identifiants identiques, car on ne peut pas déterminer qui est qui pour faire correctement les jointures entre tables

diagnostic_filtered <- diagnostics %>%
  filter (!(idmr %in% c("3433643965", "3335633931")))

# Verifier valeurs manquantes 
verif_missing<- function(df1, df2){
  # Faire une jointure complète sur la colonne 'idmr'
  full_merge <- full_join(df1, df2, by = "idmr")
  # Compter les valeurs manquantes pour chaque variable
  na_counts <- full_merge %>%
    summarise(across(everything(), ~ sum(is.na(.))))
  return(na_counts) 
}

verif_missing(patients_original, diagnostic_filtered) # les 26 valeurs manquantes sont la table diagnostic

# Merge patients et diagnostics data 

patients_full<-patients_original %>%
  right_join(diagnostic_filtered, by = "idmr") # ne joint que les lignes pour lesquelles on a pas les données de diagnostic 

str(patients_full)

## EDA rapide =====
# Vérification de la distribution des variables et abscence d'outliers avant la transformation des variables

#ExpReport(patients_full, op_file="EDA/_patients_original_smartEDA.html")

### TRANSFORMATION DE VARIABLES =====

## Age (discretization) ----
# Calcul de l'age lors de la derniere activité a partir de l'année de naissance et découpage en classes
patients<- patients_full %>%
  mutate (age=(year_last_activity-birthyear))

# Distribution variable age
plt_age<- ggplot(patients, aes(x = age)) +
  geom_histogram(binwidth = 2, fill = "blue", color = "black") +
  labs(title = "Histogramme 'Age'", x = "Age", y = "Fréquence") +
  geom_vline(xintercept = c(30, 50, 70, 90), linetype = "dotted", color = "red", size=2)

plt_age

# Decoupage en classes
patients <- patients %>%
  mutate(age_class = cut(age, breaks = c(0, 30, 50, 70, 90, Inf), labels = c("0-30", "31-50", "51-70","71-90", "91+"))) 


# Vérification qu'aucune classe ne comporte moins de 10 individus
table_age <- patients %>%
  count(age_class)
table_age

## Temps de suivi (discretization) ----

# Distribution de la variable 'temps_suivi'
plt_temps_suivi <- ggplot(patients, aes(x = temps_suivi_jours)) +
  geom_histogram(binwidth = 10, fill = "blue", color = "black") +
  labs(title = "Histogramme 'Temps de Suivi'", x = "Temps de Suivi (jours)", y = "Fréquence")+
  geom_vline(xintercept = c(0, 365, 730), linetype = "dotted", color = "red", size=2)

plt_temps_suivi

# Découpage en classes
patients <- patients %>%
  mutate(temps_suivi_class = cut(temps_suivi_jours, breaks = c(0, 365, 730, Inf), labels = c("0-1an", "1-2ans", "+2ans")))

# Vérification qu'aucune classe ne comporte pas moins de 10 individus
table_temps_suivi <- patients %>%
  count(temps_suivi_class) 

table_temps_suivi

## Alive ----
table_alive <- patients %>%
  count(alive) 
table_alive
# Variable problématique car la n=9 pour la catégorie décédé. Elle ne pourra pas être utilisée pour les analyses. Elle n'apporte pas beaucoup à l'analyse de toute façon. 

## DIAGNOSTIC (recodage)--- 

patients %>%
  distinct(diagnostic)

# Recodage 
patients <- patients %>%
  mutate(diagnostic = case_when(
    diagnostic == 304 ~ 'EB_simpl',
    diagnostic == 158684 ~ 'EBS_atresie',
    diagnostic == 257 ~ 'EBS_distrophie',
    diagnostic == 305 ~ 'EB_jonct',
    diagnostic == 79403 ~ 'ABJ_atresie',
    diagnostic == 251393 ~ 'EBJ_local'
    )) %>%
  mutate(diagnostic_general = case_when(
    diagnostic %in% c('EB_simpl', 'EBS_atresie', 'EBS_distrophie') ~ 'EBS',
    diagnostic %in% c('EB_jonct', 'ABJ_atresie', 'EBJ_local') ~ 'EBJ'))

## Level (recodage) ----

patients %>%
  count(level) 

# Recodage des valeurs de level
patients <- patients %>%  
  mutate(level = case_when(
    level == "UNK" ~ "Inconnu",
    level == "PRO" ~ "Probable",
    level == "CON" ~ "Confirmé",
    level == "ONG" ~ "EnCours",
    level == '' ~ "Non renseigné"
  )) 

## Idmr (recodage) ----
patients <- patients %>%
  mutate(idmr = as.character(idmr))

## EDA sur variables transformées ---
#ExpReport(patients, op_file="EDA/_patients_transform_smartEDA.html")

## Suppression des variables transformées/inutiles ----
patients_clean <- patients %>%
  select (-c('birthyear', "year_last_activity", 'alive'))


## ANALYSE UNIVARIEE ====

## EDA automatique sur variables transformées  ----
#ExpReport(patients_clean, op_file="EDA/patients_clean_smartEDA.html")

## ANALYSE BIVARIEE =====
 patients_noID<- patients_clean %>%
  select (-idmr)

## Exploratoire - Plot de toutes les paires de variables ----
ggpairs(patients_noID)

## Test Chi2 pour identifier correlations entre variables qualitatives
CatVar <- c('age_class', 'temps_suivi_class', 'diagnostic', 'diagnostic_general', 'level')

chi_squared_test <- function(var1, var2, data) {
  table <- table(data[[var1]], data[[var2]])
  test <- chisq.test(table)
  return(test$p.value)
}

chi2_all<- function(df, CatVar) {
  # Tableau pour conserver les p-values
  p_values <- matrix(NA, nrow = length(CatVar), ncol = length(CatVar))
  rownames(p_values) <- CatVar
  colnames(p_values) <- CatVar
  
  # Loop through all pairs of variables and perform the Chi-squared test
  for (i in 1:length(CatVar)) {
    for (j in i:length(CatVar)) {
      if (i != j) {
        p_values[i, j] <- chi_squared_test(CatVar[i], CatVar[j], df)
      }
    }
  }
  return(p_values)
}

chi2_all(patients_noID, CatVar) # message d'erreur car certaines categories n<5

## Plots individuels intéressants ---

stacked_barplot<- function(df, varX, varY, palette){
  StkPlot <- ggbarstats(
    data = df,
    x = !! sym(varX),
    y = !! sym(varY),
    package = "wesanderson",
    palette = palette)
  return(StkPlot)
}

boxplot <- function (df, continuousVar, CatVar) {
  ggbetweenstats(
    data  = df,
    x     = !! sym(CatVar),
    y     = !! sym(continuousVar),
    title = paste("Distribution de", continuousVar, "en fonction de", CatVar, sep=' ')
  )
}
    

patients_noID <- patients_noID %>% 
  mutate(level = factor(level, levels = c("Non renseigné", "Inconnu", "EnCours", "Probable", "Confirmé"), ordered = TRUE))

# Summary
summary(patients_noID)

# PLot types et sous-types de diagnostic
stacked_barplot(patients_noID, 'diagnostic', 'diagnostic_general', 'AsteroidCity2')

# Plot diagnostic et niveau
stacked_barplot(patients_noID, 'level', 'diagnostic_general', 'Zissou1')
stacked_barplot(patients_noID, 'level', 'diagnostic', 'Zissou1')

# Plot niveau de diagnostic et temps de suivi
stacked_barplot(patients_noID, 'level', 'temps_suivi_class', 'Zissou1')
boxplot(patients_noID, 'temps_suivi_jours', 'level')

# Plot repartition age et diagnostic
boxplot(patients_noID, 'age', 'diagnostic_general')
boxplot(patients_noID, 'age', 'diagnostic')

stacked_barplot(patients_noID, 'diagnostic', 'age_class', 'AsteroidCity2')

# Stacked barplot non proportionnel age et diagnostic
# 
# ggplot(patients_noID, aes(fill=diagnostic, y=value, x=diagnostic_general)) + 
#   geom_bar(position="stack", stat="identity") +
#   scale_fill_viridis(discrete = T) +
#   ggtitle(".....") +
#   theme_minimal() +
#   xlab("")

### ANALYSE NIVEAUX DE PLOMB =======

## Jointure et transformation ------

plomb <- plomb %>%
  mutate(idmr = as.character(idmr))

verif_missing(patients_clean, plomb)

# Il faut effectuer un inner join car il y a des valeurs manquantes de part et d'autre

plomb_full<- patients_clean %>%
  inner_join(plomb, by='idmr')

# Résumé des données
summary(plomb_full)

# Transformer en format long pour colonnes W0 à W6
plomb_long <- plomb_full %>%
  pivot_longer(cols = starts_with("W"), 
               names_to = "semaine", 
               values_to = "niveau_plomb")

## Verifier distribution normale de la variable niveau_plomb ------

test_normality_niveauPlomb <- function(df) {
  # Histogramme
  hist<- ggplot(df, aes(x = niveau_plomb)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "blue", color = "black") +
    geom_density(color = "red") +
    labs(title = "Histogramme", x = "Variable Numérique", y = "Densité")
  print(hist)
  # Q-Q Plot
  qqplot<- ggplot(df, aes(sample = niveau_plomb)) +
    stat_qq() +
    stat_qq_line() +
    labs(title = "Q-Q Plot", x = "Théorique", y = "Observé")
  print(qqplot)
  # Test de Shapiro-Wilk
  shapiro_test <- shapiro.test(df$niveau_plomb)
  print(shapiro_test)
}

test_normality_niveauPlomb(plomb_long)

# Identifier outlier - verifier démarche
max_outlier <- plomb_long %>%
  filter(niveau_plomb == max(niveau_plomb)) %>%
  select(idmr, semaine, niveau_plomb)

print(max_outlier)

# Retirer individu outlier? on ne peux pas eliminer des données sans bonne raison
# plomb_long_filtered <- plomb_long %>%
#   filter(idmr != '3033313335')
# 
# test_normality_niveauPlomb(plomb_long_filtered)

## Vérification homogénéité de variance -----

leveneTest(niveau_plomb ~ semaine, data = plomb_long)
leveneTest(niveau_plomb ~ diagnostic_general, data = plomb_long)
leveneTest(niveau_plomb ~ diagnostic, data = plomb_long)

 # --> prvilégier tests statistiques non paramétriques

## Evolution du niveau de plomb pour l'ensemble des patients - Echantillons appariés ----

ggwithinstats(
  data    = plomb_long,
  x       = semaine,
  y       = niveau_plomb,
  title   = "Evolution des niveaux de plomb de tous les patients",
  ylab            = "Niveau de Plomb (ug/dL)",
  xlab            = "Semaine de traitement",
)

### Test statistique non-paramétrique 
kruskal.test(niveau_plomb ~ semaine, data = plomb_long)

## Evolution du niveau de plomb pour patients selon diagnostic général ----

grouped_ggwithinstats(
  data            = plomb_long,
  x               = semaine,
  y               = niveau_plomb,
  type            = "np",
  ylab            = "Niveau de Plomb (ug/dL)",
  xlab            = "Semaine de traitement",
  grouping.var    = diagnostic_general
)



## Evolution du niveau de plomb pour patients selon diagnostic détaillé ----

grouped_ggwithinstats(
  data            = plomb_long,
  x               = semaine,
  y               = niveau_plomb,
  type            = "np",
  ylab            = "Niveau de Plomb (ug/dL)",
  xlab            = "Semaine de traitement",
  grouping.var    = diagnostic
)

## Evolution du niveau de plomb pour patients selon classe age ----

# pas de correlation entre age et niveau de plomb détecté en général
ggscatterstats(
  data  = plomb_long,
  x     = age,
  y     = niveau_plomb,
  xlab  = "age",
  ylab  = "niveau de plomb",
  title = "Age x niveau plomb"
)

grouped_ggwithinstats(
  data            = plomb_long,
  x               = semaine,
  y               = niveau_plomb,
  type            = "np",
  ylab            = "Niveau de Plomb (ug/dL)",
  xlab            = "Semaine de traitement",
  grouping.var    = age_class
)


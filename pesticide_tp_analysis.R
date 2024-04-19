
## Data Wrangling
# Set working directory
setwd("C:/Users/Melina Padron/OneDrive/Spring 2024 Research")


# Load packages
library(tidyverse)
library(tidymodels)
library(ggbiplot)
library(kernlab)
library(gridExtra)
library(grid)
library(knitr)
library(ggrepel)
library(cowplot)
library(FSA)

# Read data
parent_table <- read_delim("Research Data/Table3a_RSQA_DetectedOrganics_Parentonly.txt")
degradate_table <- read_delim("Research Data/Table3b_RSQA_Detectedorganics_Degradatesonly.txt")
rsqa_flow <- read_csv("Research Data/RSQA_Flow_Data.csv")
watershed <- read.csv("Research Data/RSQA_WatershedData.csv")

# Add column that says whether it is a parent or pesticide degradate
pesticide <- full_join(parent_table, degradate_table)
pesticide <- pesticide %>%
  rename("Pesticide Type" = "Pesticide degradate") %>%
  mutate("Pesticide Type" = ifelse(`Pesticide Type` == "Pesticide", "Parent", "Degradate"))

# Join flow data and pesticide data by station number and date
pesticide_and_flow <- inner_join(rsqa_flow, pesticide, by = c("TSITE_NO" = "USGS station number", "Date" = "Sample date (YYYYMMDD)"))

# Join pesticide and flow data with watershed data by station number
pfw <- inner_join(pesticide_and_flow, watershed, by = c("TSITE_NO" = "TSITE_NO"))

# Select detected pesticides
original <- pfw %>%
  filter(Remark == "DET",
         !is.na(BaseFlow.Frac)) %>%
  mutate(flow_type = ifelse(BaseFlow.Frac <= 0.5, "Runoff", "Groundwater"))

# Count unique parent pesticides
unique_parent_pesticides <- original %>%
  filter(`Pesticide Type` == "Parent") %>%
  distinct(`Chemical name`) %>%
  count()

# Count unique degradates
unique_degradates <- original %>%
  filter(`Pesticide Type` == "Degradate") %>%
  distinct(`Chemical name`) %>%
  count()

# Print the counts
print(unique_parent_pesticides)
print(unique_degradates)

## Preliminary Statistical Analysis
# Wilcoxon test for pesticide concentrations
wilcoxon_ptype <- wilcox.test(`Result (ug/L)` ~ `Pesticide Type`, data = original)

# Extract relevant information from the test result
results <- data.frame(
  Statistic = wilcoxon_ptype$statistic,
  P_value = wilcoxon_ptype$p.value,
  Alternative = wilcoxon_ptype$alternative,
  stringsAsFactors = FALSE
)

# Create the table using knitr::kable()
knitr::kable(results, 
             format = "simple")

# Prepare data for Kruskal-Wallis test
original$"Pesticide Type:flow_type" <- with(original, paste(`Pesticide Type`, flow_type, sep = ":"))

# Kruskal-Wallis test
kruskal <- kruskal.test(`Result (ug/L)` ~ `Pesticide Type:flow_type`, data = original)

# Extract relevant information from the test result
results_kruskal <- data.frame(
  Statistic = kruskal$statistic,
  P_value = sprintf("%.2e", kruskal$p.value),
  Df = kruskal$parameter
  
)

# Create the table using knitr::kable()
knitr::kable(results_kruskal, 
             format = "simple")

# Dunn Test
dunnTest(`Result (ug/L)` ~ `Pesticide Type:flow_type`, 
         data = original,
         method = "bonferroni",
         list=TRUE)

# Create a data frame to store the provided data
results_table <- data.frame(
  Comparison = c("Degradate:Groundwater - Degradate:Runoff",
                 "Degradate:Groundwater - Parent:Groundwater",
                 "Degradate:Runoff - Parent:Groundwater",
                 "Degradate:Groundwater - Parent:Runoff",
                 "Degradate:Runoff - Parent:Runoff",
                 "Parent:Groundwater - Parent:Runoff"),
  Z = c(4.447179, 6.527394, 1.144261, 1.567790, -3.061596, -4.904178),
  P_unadj = c(8.700525e-06, 6.692390e-11, 2.525155e-01, 1.169301e-01, 2.201604e-03, 9.381954e-07),
  P_adj = c(5.220315e-05, 4.015434e-10, 1.000000e+00, 7.015804e-01, 1.320962e-02, 5.629173e-06)
)

# Display the table using knitr::kable()
knitr::kable(results_table, 
             format = "simple")

## Exploratory Data Analysis (EDA)
# Select relevant variables and only look at Pesticide Degradates
pd_data <- original %>%
  filter(`Pesticide Type` == "Degradate") %>%
  select(`Result (ug/L)`, Depth.M, Q.M3.S, BaseFlow.Frac, StreamDensity,
         SiteID.Flow, Date, `Chemical name`, flow_type)

# Scatter plots for predictor variables vs PD concentrations
pd_data %>%
  select(`Result (ug/L)`, Depth.M, Q.M3.S, BaseFlow.Frac, StreamDensity) %>%
  pivot_longer(-`Result (ug/L)`, names_to = "name", values_to = "metric") %>%
  ggplot(aes(x = metric, y = `Result (ug/L)`)) +
  geom_point(alpha = 0.10, na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, na.rm = TRUE) +
  facet_wrap(vars(name), scales = "free_x") +
  labs(x = "Predictor Values",
       y = "Pesticide Concentrations (ug/L)") +
  theme_classic() +
  theme(strip.background = element_blank())

# Visualize outliers with box plot
pd_data %>%
  ggplot(aes(x = "", y = `Result (ug/L)`)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(x = "", y = "Pesticide Concentrations (ug/L) (Log 10 Scale)") +
  theme_classic()
    

# Remove outliers
pd_data_clean <- pd_data %>%
  filter(`Result (ug/L)` >= quantile(`Result (ug/L)`, 0.25) - 1.5 * IQR(`Result (ug/L)`) &
           `Result (ug/L)` <= quantile(`Result (ug/L)`, 0.75) + 1.5 * IQR(`Result (ug/L)`))

# Create scatterplots that compare all predictor variables vs PD concentrations
# with clean dataset
pd_data_clean %>%
  select(`Result (ug/L)`, Depth.M, Q.M3.S, BaseFlow.Frac, StreamDensity) %>%
  pivot_longer(-`Result (ug/L)`, names_to = "name", values_to = "metric") %>%
  ggplot(aes(x = metric, y = `Result (ug/L)`)) +
  geom_point(alpha = 0.10, na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, na.rm = TRUE) +
  facet_wrap(vars(name), scales = "free_x") +
  labs(x = "Predictor Values",
       y = "Pesticide Concentrations (ug/L)") +
  theme_classic() +
  theme(strip.background = element_blank())

## Principal Component Analysis (PCA)
# Create metadata for pca plots
# Identify the relevant rows
pd_data <- pd_data_clean %>% na.omit()
relevant_siteID_flow <- unique(pd_data$SiteID.Flow)

# Filter metadata to retain only the relevant rows
metadata <- pd_data_clean %>%
  filter(SiteID.Flow %in% relevant_siteID_flow) %>%
  na.omit()

# Join relevant metadata with pd_data
metadata <- semi_join(metadata, pd_data, by = c("SiteID.Flow", "Date"))

# Select predictors and outcome variables
numeric_clean <- pd_data %>%
  select(Depth.M, Q.M3.S, BaseFlow.Frac, StreamDensity) %>%
  na.omit()

# Perform PCA
pca_result <- prcomp(numeric_clean, center = TRUE, scale. = TRUE)

# Explore PCA results
head(pca_result)
summary(pca_result)

# Visualize PCA results
# Scree plot
pca_result %>%
  # Extract eigenvalues
  tidy(matrix = "eigenvalues") %>%
  mutate(cumulative_percent = cumsum(percent)) %>%
  ggplot(aes(PC)) + 
  geom_col(aes(y = percent), fill = "skyblue", alpha = 0.7) +
  geom_line(aes(y = cumulative_percent), color = "red", size = 1) +
  geom_point(aes(y = cumulative_percent), color = "red", size = 2) +
  scale_x_continuous(
    limits = c(0, 5),
    breaks = 1:4
  ) +
  scale_y_continuous(
    name = "Variance explained %",
    label = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(y = "Variance explained %",
       x = "Principal Component") +
  theme_classic()

# Biplot
arrow_style <- arrow(angle = 20,
                     length = grid::unit(8, "pt"),
                     ends = "first", type = "closed")
  
pca_result %>%
  tidy(matrix = "rotation") %>%
  pivot_wider(
    names_from = "PC", values_from = "value",
    names_prefix = "PC"
  ) %>%
  ggplot(aes(PC1, PC2)) +
  geom_segment(
    xend = 0, yend = 0,
    arrow = arrow(
      angle = 20, length = grid::unit(8, "pt"),
      ends = "first", type = "closed"
    )
  ) +
  geom_text_repel(aes(label = column), size = 3) +
  xlim(-1, 1) + ylim(-1, 1) +
  coord_fixed() +
  labs(
    x = "PC1",
    y = "PC2"
  ) +
  theme_classic()

# Loading plot with numeric variables
loadings_tp <- pca_result %>%
  # Extract rotation matrix
  tidy(matrix = "rotation") %>%
  pivot_wider(
    names_from = "PC", values_from = "value",
    names_prefix = "PC"
  ) %>%
  select(' ' = 'column', 'PC1', 'PC2')

# Display table and adjust settings
knitr::kable(loadings_tp, format = "simple")

# Biplot - Grouping by baseflow
ggbiplot(pca_result, choices = c(1, 2), var.axes = FALSE,
         groups = paste(metadata$flow_type), ellipse = TRUE) + 
  theme_classic() +
  xlab("PC1") +
  ylab("PC2") +
  labs(color = "Flow Condition")

## Prepare data for Predictive Modeling
# Extract the first 2 PCs
pca_vars <- as.data.frame(pca_result$x[, 1:2])

# Combine the first 2 PCs with the outcome variable
pca_data_final <- cbind(pca_vars, pd_data$`Result (ug/L)`, pd_data$`Chemical name`)

# Rename outcome variable
pca_data_final <- pca_data_final %>%
  rename(outcome =  "pd_data$`Result (ug/L)`")

# Make code reproducible
set.seed(123)

# Split the data into training and test sets
dat_split <- initial_split(pca_data_final, prop = 0.8)  # 80% training, 20% testing
dat_train <- training(dat_split)
dat_test <- testing(dat_split)

## Predictive Models
## Multiple Linear Regression

# Define the recipe
rec_lr <- dat_train %>%
  recipe(outcome ~ PC1 + PC2)

# Specify the linear regression model
model_lr <- linear_reg() %>%
  set_engine("lm") %>%
  set_mode("regression")

# Create the workflow
wf_lr <- workflow() %>%
  add_recipe(rec_lr) %>%
  add_model(model_lr)

# Fit the model
lr_results <- fit(wf_lr, data = dat_train)

# Use cross-validation
# Make code reproducible
set.seed(123)

# Check performance using cross-validation
folds <- vfold_cv(dat_train, v = 10)
lr_results <- fit_resamples(wf_lr, resamples = folds)
lr_metrics <- lr_results %>%
  collect_metrics()
lr_metrics

## k-NN Model
# Define the recipe
rec_knn <- dat_train %>%
  recipe(outcome ~ PC1 + PC2)

# Specify the k-NN model
model_knn <- nearest_neighbor(neighbors = 10) %>%
  set_engine("kknn") %>%
  set_mode("regression")

# Create the workflow
wf_knn <- workflow() %>%
  add_model(model_knn) %>%
  add_recipe(rec_knn)

# Fit the model
knn_results <- fit(wf_knn, data = dat_train)

# Check performance using cross-validation
knn_results <- fit_resamples(wf_knn, resamples = folds)
knn_metrics <- knn_results %>%
  collect_metrics()
knn_metrics

# Tune k-NN 
# Set up knn model
model_knn <- nearest_neighbor(neighbors = tune("k")) %>%
  set_engine("kknn") %>%
  set_mode("regression")

# Create the workflow
wf_knn <- workflow() %>%
  add_model(model_knn) %>%
  add_recipe(rec_knn)

# Tune for the optimal number of neighbors
res_knn <- tune_grid(wf_knn, resamples = folds,
                     grid = tibble(k = c(3, 5, 10, 15, 20, 25, 30, 35)))

# Find best number of neighbors using RMSE as the main metric
res_knn %>%
  select_best(metric = "rmse")

# Now, use tuned k-NN parameters
# Define the recipe
rec_knn <- dat_train %>%
  recipe(outcome ~ PC1 + PC2)

# Specify the k-NN model
model_knn <- nearest_neighbor(neighbors = 35) %>%
  set_engine("kknn") %>%
  set_mode("regression")

# Create the workflow
wf_knn <- workflow() %>%
  add_model(model_knn) %>%
  add_recipe(rec_knn)

## Check performance using cross-validation using same folds
knn_results <- fit_resamples(wf_knn, resamples = folds)
knn_metrics <-knn_results %>%
  collect_metrics()
knn_metrics

## Random Forest Model
# Define the recipe
rec_rf <- dat_train %>%
  recipe(outcome ~ PC1 + PC2)

# Specify the random forest model
model_rf <- rand_forest(mtry = 2,
                        min_n = 1) %>%
  set_engine("ranger") %>%
  set_mode("regression")

# Create the workflow
wf_rf <- workflow() %>%
  add_recipe(rec_rf) %>%
  add_model(model_rf)

# Check performance using cross-validation using same folds
rf_results <- fit_resamples(wf_rf, resamples = folds)
rf_metrics <- rf_results %>%
  collect_metrics()

# Try a grid of tuning parameters
model_rf <- rand_forest(mtry = tune("mtry"),
                        min_n = tune("min_n")) %>%
  set_engine("ranger") %>%
  set_mode("regression")

# Create the workflow
wf_rf <- workflow() %>%
  add_recipe(rec_rf) %>%
  add_model(model_rf)

# Fit model over grid of tuning parameters
res_rf <- tune_grid(wf_rf, resamples = folds,
                    grid = expand.grid(mtry = c(1, 2),
                                       min_n = c(3, 5, 7)))

# Find best parameters using RMSE as the main metric
res_rf %>%
  select_best(metric = "rmse")

# Now, use best parameters
rec_rf <- dat_train %>%
  recipe(outcome ~ PC1 + PC2)

# Specify the random forest model
model_rf <- rand_forest(mtry = 1,
                        min_n = 7) %>%
  set_engine("ranger") %>%
  set_mode("regression")

# Create the workflow
wf_rf <- workflow() %>%
  add_recipe(rec_rf) %>%
  add_model(model_rf)

# Check performance using cross-validation
rf_results <- fit_resamples(wf_rf, resamples = folds)
rf_metrics <- rf_results %>%
  collect_metrics()
rf_metrics

## SVM (Support Vector Model)
# Define the recipe
rec_svm <- dat_train %>%
  recipe(outcome ~ PC1 + PC2)

# Specify the SVM model
model_svm <- svm_rbf() %>%
  set_engine("kernlab") %>%
  set_mode("regression")

# Create the workflow
wf_svm <- workflow() %>%
  add_recipe(rec_svm) %>%
  add_model(model_svm)

# Check performance using cross-validation using same folds\
wf_results <- fit_resamples(wf_svm, resamples = folds)

wf_metrics <-wf_results %>%
  collect_metrics()
wf_metrics

## Overview of all models
# Combine metrics into a summary table
all_metrics <- bind_rows(
  lr_metrics %>% mutate(model = "Linear Regression"),
  knn_metrics %>% mutate(model = "KNN"),
  rf_metrics %>% mutate(model = "Random Forest"),
  wf_metrics %>% mutate(model = "SVM"),
)
# View summary table
all_metrics %>%
  pivot_wider(names_from = model, values_from = mean) %>%
  select(.metric,`Linear Regression`, `KNN`, `Random Forest`, `SVM`)

# Manually clean up table from above
summary_table <- data.frame(Metric = c("RMSE", "rsq"),
                            Linear_Regression = c("26.6", "0.0431"),
                            KNN = c("26.6","0.0516"),
                            Random_Forest = c("28", "0.0354"),
                            SVM = c('28.7','0.0351'))

# Display table and adjust settings
knitr::kable(summary_table, format = "simple")

## Evaluate models
# Make linear regression model the "final model"
lr_final <- last_fit(wf_lr, split = dat_split)

# Summary statistics for MLR
lr_metrics <- lr_final %>%
  collect_metrics()

# Make k-NN model the "final model"
knn_final <- last_fit(wf_knn, split = dat_split)

# Summary statistics for k-NN
knn_metrics <- knn_final %>%
  collect_metrics() 

# Make Random Forest model the "final model"
rf_final <- last_fit(wf_rf, split = dat_split)

# Summary statistics for Random Forest
rf_metrics <- rf_final %>%
  collect_metrics()

# Make SVM model the "final model"
svm_final <- last_fit(wf_svm, split = dat_split)

# Summary statistics for SVM
svm_metrics <- svm_final %>%
  collect_metrics()

## Visualize results
# Combine the plots
plot_grid(
  lr_final %>%
    collect_predictions() %>%
    ggplot(aes(.pred, outcome)) +
    geom_point(alpha = 0.10, na.rm = TRUE) +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    labs(title = "Linear Regression ",
         x = " ",
         y = "Observed Concentrations (ug/L)") +
    theme_classic(),
  
  knn_final %>%
    collect_predictions() %>%
    ggplot(aes(.pred, outcome)) +
    geom_point(alpha = 0.10, na.rm = TRUE) +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    labs(title = "k-NN ",
         x = " ",
         y = " ") +
    theme_classic(),
  
  rf_final %>%
    collect_predictions() %>%
    ggplot(aes(.pred, outcome)) +
    geom_point(alpha = 0.10, na.rm = TRUE) +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    labs(title = "Random Forest",
         x = "Predicted Concentrations (ug/L)",
         y = "Observed Concentrations (ug/L)") +
    theme_classic(),
  
  svm_final %>%
    collect_predictions() %>%
    ggplot(aes(.pred, outcome)) +
    geom_point(alpha = 0.10, na.rm = TRUE) +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    labs(title = "SVM",
         x = "Predicted Concentrations (ug/L)",
         y = " ") +
    theme_classic(),
  
  ncol = 2)


## Model Performance 
# Add predictions from final model to test data
model_test <- knn_final %>%
  extract_fit_parsnip() %>%
  augment(new_data = dat_test)

## Identify specific pesticides
# Find pesticides with the min differences between observed and predicted values
min_pest <- model_test %>%
  mutate(Difference = abs(.resid)) %>%
  arrange(Difference) %>%
  rename(`Chemical Name` = "pd_data$`Chemical name`",
         Observed = outcome,
         Predicted = .pred) %>%
  select(`Chemical Name`, `Predicted`, `Observed`, `Difference`) %>%
  distinct(`Chemical Name`,.keep_all = TRUE) %>%
  slice(1:20)

# Display table
knitr::kable(min_pest, format = "simple")

# Find pesticides with the max differences between observed and predicted values
max_pest <- model_test %>%
  mutate(Difference = abs(.resid)) %>%
  arrange(-Difference) %>%
  rename(`Chemical Name` = "pd_data$`Chemical name`",
         Observed = outcome,
         Predicted = .pred) %>%
  select(`Chemical Name`, `Predicted`, `Observed`, `Difference`) %>%
  distinct(`Chemical Name`,.keep_all = TRUE) %>%
  slice(1:20)

# Display table
knitr::kable(max_pest, format = "simple")


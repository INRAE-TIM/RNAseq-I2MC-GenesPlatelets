# Loading assets ###############################################################
source("./R/Data_load.R")


# Data visualization ###########################################################
# plot(sqrt(pData$TP),
#      Counts[Counts$gene_name == "Cyp4a14", -1],
#      pch = 20,
#      col = as.factor(pData$diet))
# legend("topright",
#        as.character(unique(pData$diet)),
#        col = 1:2,
#        pch = 20)
#
# LM.all |>
#   group_by(gene_name, diet) |>
#   slice(1) |>
#   ggplot() +
#   geom_histogram(aes(x = r.sqr)) +
#   facet_wrap(~ diet)




# Diet independent genes #######################################################
## Tables ####
LM.mean.perDiet <- fit |>
  pivot_longer(!gene_name, names_to = "samples", values_to = "fitted") |>
  separate(samples, into = c("diet", "samples"), sep = "\\.") |>
  separate(samples,
           into = c("treatment", "rep_name"),
           sep = "_") |>
  left_join(pData) |>
  group_by(gene_name, diet) |>
  nest() |>
  dplyr::mutate(models = map(data, ~ lm(fitted ~ TP, data = .))) |>
  ungroup()


LM.mean.ALL <- fit |>
  pivot_longer(!gene_name, names_to = "samples", values_to = "fitted") |>
  separate(samples, into = c("diet", "samples"), sep = "\\.") |>
  separate(samples,
           into = c("treatment", "rep_name"),
           sep = "_") |>
  left_join(pData) |>
  dplyr::mutate(diet = "ALL") |>
  group_by(gene_name, diet) |>
  nest() |>
  mutate(models = map(data, ~ lm(fitted ~ TP, data = .))) |>
  ungroup()

# Best model chosed from the maximum r²
LM.all <- rbind.data.frame(LM.mean.perDiet, LM.mean.ALL) |>
  group_by(gene_name, diet) |>
  mutate(r.sqr = map(models, ~ summary(.)$r.squared),
         r.sqr = as.numeric(r.sqr)) |>
  mutate(sumUps = map(models, ~ tidy(.))) |>
  unnest(sumUps) |>
  group_by(term) |>
  mutate(adj.p = p.adjust(p.value, method = "BH")) |>
  group_by(gene_name) |>
  mutate(bestR = ifelse(r.sqr == max(r.sqr), "best", NA)) |>
  ungroup()


diet.indep.genes <- LM.all |>
  filter(bestR == "best") |>
  # if bestR = ALL ==> diet does not influence gene expression (models 1, 2 or 3)
  filter(diet == "ALL") |>
  select(gene_name) |>
  unlist() |>
  unique()


diet.indep.table <- LM.all |>
  filter(gene_name %in% diet.indep.genes &
           term == "TP" & diet == "ALL") |>
  mutate(modelTP = ifelse(adj.p > 0.01, 1, ifelse(estimate < 0, 2, 3)))



## Models ####
### Model 1 - flat
unaffected.genes.model <- diet.indep.table |>
  filter(modelTP == 1) |>
  arrange(desc(estimate)) |>
  left_join(Counts) |>
  column_to_rownames("gene_name") |>
  select(contains("ND") | contains("HFC"))


### Model 2 - both down (nonexistant)
m2 <- diet.indep.table |>
  filter(modelTP == 2) |>
  arrange(desc(estimate)) |>
  left_join(Norm.Counts) |>
  column_to_rownames("gene_name") |>
  select(contains("ND") | contains("HFC"))


### Model 3 - both up (nonexistant)
m3 <- diet.indep.table |>
  filter(modelTP == 3) |>
  arrange(desc(estimate)) |>
  left_join(Norm.Counts) |>
  column_to_rownames("gene_name") |>
  select(contains("ND") | contains("HFC"))




# Diet dependent genes #########################################################
## Tables ####
diet.dep.table <- LM.all |>
  select(!bestR) |>
  filter(!gene_name %in% diet.indep.genes &
           term == "TP" & diet != "ALL") |>
  group_by(gene_name, diet) |>
  mutate(slope = ifelse(
    r.sqr > 0.8 & adj.p < 0.01 & estimate > 0,
    "up",
    ifelse(r.sqr > 0.8 &
             adj.p < 0.01 & estimate < 0, "down", "flat")
  )) |>
  ungroup()


two.flat <- diet.dep.table |>
  select(gene_name, diet, slope) |>
  pivot_wider(names_from = diet, values_from = slope) |>
  group_by(gene_name) |>
  filter(ND == "flat" & HFC == "flat") |>
  pivot_longer(names_to = "diet", cols = c(HFC, ND)) |>
  left_join(y = LM.all, by = c("gene_name", "diet")) |>
  select(gene_name, diet, term, estimate) |>
  filter(term == "(Intercept)") |>
  pivot_wider(names_from = diet, values_from = estimate) |>
  mutate(slopes = "flat") |>
  ungroup()


two.down <- diet.dep.table |>
  select(gene_name, diet, slope) |>
  pivot_wider(names_from = diet, values_from = slope) |>
  group_by(gene_name) |>
  filter(ND == "down" & HFC == "down") |>
  pivot_longer(names_to = "diet", cols = c(HFC, ND)) |>
  left_join(y = LM.all, by = c("gene_name", "diet")) |>
  select(gene_name, diet, term, estimate) |>
  filter(term == "(Intercept)") |>
  pivot_wider(names_from = diet, values_from = estimate) |>
  mutate(slopes = "down") |>
  ungroup()


two.up <- diet.dep.table |>
  select(gene_name, diet, slope) |>
  pivot_wider(names_from = diet, values_from = slope) |>
  group_by(gene_name) |>
  filter(ND == "up" & HFC == "up") |>
  pivot_longer(names_to = "diet", cols = c(HFC, ND)) |>
  left_join(y = LM.all, by = c("gene_name", "diet")) |>
  select(gene_name, diet, term, estimate) |>
  filter(term == "(Intercept)") |>
  pivot_wider(names_from = diet, values_from = estimate) |>
  mutate(slopes = "up") |>
  ungroup()



## Models ####
### Symmetric slopes models ####
symmetrics.models <- two.flat |>
  # m04 - both flat - ND higher
  filter(ND > HFC) |>
  select(gene_name, ND) |>
  rename(estimate = ND) |>
  mutate(
    diet = "ND",
    slope = "flat",
    model_name = "4",
    .after = gene_name
  ) |>
  # m05 - both flat - HFC higher
  bind_rows(
    two.flat |>
      filter(ND < HFC) |>
      select(gene_name, ND) |>
      rename(estimate = ND) |>
      mutate(
        diet = "HFC",
        slope = "flat",
        model_name = "5",
        .after = gene_name
      )
  ) |>
  # m06 - both down - ND higher
  bind_rows(
    two.down |>
      filter(ND > HFC) |>
      select(gene_name, ND) |>
      rename(estimate = ND) |>
      mutate(
        diet = "ND",
        slope = "down",
        model_name = "6",
        .after = gene_name
      )
  ) |>
  # m07 - both down - HFC higher
  bind_rows(
    two.down |>
      filter(ND < HFC) |>
      select(gene_name, ND) |>
      rename(estimate = ND) |>
      mutate(
        diet = "HFC",
        slope = "down",
        model_name = "7",
        .after = gene_name
      )
  ) |>
  # m08 - both up - ND higher
  bind_rows(
    two.up |>
      filter(ND > HFC) |>
      select(gene_name, ND) |>
      rename(estimate = ND) |>
      mutate(
        diet = "ND",
        slope = "up",
        model_name = "8",
        .after = gene_name
      )
  ) |>
  # m09 - both up - HFC higher
  bind_rows(
    two.up |>
      filter(ND < HFC) |>
      select(gene_name, ND) |>
      rename(estimate = ND) |>
      mutate(
        diet = "HFC",
        slope = "up",
        model_name = "9",
        .after = gene_name
      )
  )


### Asymmetric slopes models ####
asymmetrics.models <- diet.dep.table |>
  group_by(gene_name) |>
  # m10 - ND sensible down
  filter((diet == "ND" &
            slope == "down") |
           (diet == "HFC" & slope == "flat")) |>
  filter(n() == 2) |>
  select(gene_name, diet, slope, estimate) |>
  mutate(model_name = "10", .after = slope) |>
  # m11 - ND sensible up
  bind_rows(
    diet.dep.table |>
      group_by(gene_name) |>
      filter((diet == "ND" &
                slope == "up") |
               (diet == "HFC" & slope == "flat")) |>
      filter(n() == 2) |>
      select(gene_name, diet, slope, estimate) |>
      mutate(model_name = "11", .after = slope) |>
      ungroup()
  ) |>
  # m12 - HFC sensible down
  bind_rows(
    diet.dep.table |>
      group_by(gene_name) |>
      filter((diet == "ND" & slope == "flat") |
               (diet == "HFC" &
                  slope == "down")) |>
      filter(n() == 2) |>
      select(gene_name, diet, slope, estimate) |>
      mutate(model_name = "12", .after = slope) |>
      ungroup()
  ) |>
  # m13 - HFC sensible up
  bind_rows(
    diet.dep.table |>
      group_by(gene_name) |>
      filter((diet == "ND" & slope == "flat") |
               (diet == "HFC" &
                  slope == "up")) |>
      filter(n() == 2) |>
      select(gene_name, diet, slope, estimate) |>
      mutate(model_name = "13", .after = slope) |>
      ungroup()
  ) |>
  # m14 - ND down + HFC up
  bind_rows(
    diet.dep.table |>
      group_by(gene_name) |>
      filter((diet == "ND" & slope == "down") |
               (diet == "HFC" &
                  slope == "up")) |>
      filter(n() == 2) |>
      select(gene_name, diet, slope, estimate) |>
      mutate(model_name = "14", .after = slope) |>
      ungroup()
  ) |>
  # m15 - ND up + HFC down
  bind_rows(
    diet.dep.table |>
      group_by(gene_name) |>
      filter((diet == "ND" & slope == "up") |
               (diet == "HFC" &
                  slope == "down")) |>
      filter(n() == 2) |>
      select(gene_name, diet, slope, estimate) |>
      mutate(model_name = "15", .after = slope) |>
      ungroup()
  ) |>
  ungroup() |>
  filter(!slope == "flat")


### Genes number per model of both diet dependent tables ####
bind_rows(symmetrics.models, asymmetrics.models) |>
  group_by(model_name) |>
  summarize(size = n())




# Heatmap ######################################################################
HM.model <- function(model, model_no) {
  gene_order <- model |>
    arrange(estimate) |>
    filter(model$model_name == model_no)
  
  
  Counts <- column_to_rownames(Counts, "gene_name")
  
  data.HM <- cbind.data.frame(Counts[gene_order$gene_name, row.names(pData[grepl("ND", pData$diet), ])
                                     [order(pData[grepl("ND", pData$diet), "TP"])]], Counts[gene_order$gene_name, row.names(pData[grepl("HFC", pData$diet), ])
                                                                                            [order(pData[grepl("HFC", pData$diet), "TP"])]])
  
  pheatmap(
    data.HM,
    scale = "row",
    color = colorRampPalette(c("blue", "yellow"))(1000),
    cluster_rows = F,
    cluster_cols = F,
    show_rownames = F,
    main = paste("Model", model_no, sep = " "),
    legend = F,
    show_colnames = F,
    gaps_col = 15
  )
}



## Profile of HM per model ####
profiles <- Counts |>
  pivot_longer(!gene_name) |>
  left_join(rownames_to_column(pData, var = "name")) |>
  group_by(gene_name) |>
  mutate(z.score = as.numeric(scale(value))) |>
  group_by(gene_name, diet, treatment) |>
  summarise(z.score = mean(z.score), TP = mean(TP)) |>
  ungroup() |>
  left_join(
    symmetrics.models |>
      bind_rows(asymmetrics.models) |>
      select(gene_name, model_name) |>
      group_by(gene_name) |>
      summarize(model_name = mean(as.numeric(model_name)))
  ) |>
  group_by(model_name) |>
  filter(n() > 10)


profiles.plot <- profiles |>
  ggplot(aes(x = sqrt(TP), y = z.score, colour = diet)) +
  stat_summary(
    fun.data = mean_ci,
    geom = "errorbar",
    size = .5,
    na.rm = TRUE
  ) +
  geom_smooth(method = "lm") +
  facet_wrap(~ model_name, ncol = 6)




# Graphs files generation (PDF) #########################################################################
repertory <- "~/Severin_GenesPlatelets/Analysis/Outputs/Heatmaps.pdf"
symmetricals <- c(4, 5, 7, 8, 9) # because the model 6 doesn't contains any gene, it needs to be counted out for the for loop to run

pdf(file = repertory)

profiles.plot

for (i in symmetricals) {
  HM.model(symmetrics.models, model_no = i)
}

for (i in 10:15) {
  HM.model(asymmetrics.models, model_no = i)
}

dev.off()




# Tables files generation (XLSX) ##############################################################
export <- createWorkbook(title = "Tables of genes by Models")

addWorksheet(export, sheetName = "Diet Independent genes")
writeData(export, sheet = "Diet Independent genes", x = diet.indep.table)

addWorksheet(export, sheetName = "Diet Dependent genes")
writeData(export, sheet = "Diet Dependent genes", x = diet.dep.table)

addWorksheet(export, sheetName = "Symmetrical Slopes models")
writeData(export, sheet = "Symmetrical Slopes models", x = symmetrics.models)

addWorksheet(export, sheetName = "Asymmetrical Slopes models")
writeData(export, sheet = "Asymmetrical Slopes models", x = asymmetrics.models)

saveWorkbook(wb = export, file = "./Outputs/Tables_of_genes_by_Models.xlsx")

# Open your file in xlsx format, don't forget to save it !
openXL(export)

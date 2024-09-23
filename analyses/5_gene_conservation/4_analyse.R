library(ggplot2)
pdf("output/plots/5_mean_conservation_scores.pdf", height=14)

for(seqtype in c("proteins", "cds")) {
  mean_conservation_scores <- data.frame(gene = character(),
                                         mean_score = numeric(),
                                         stringsAsFactors = FALSE)
  
  # Loop through each file
  for (file in list.files(paste0("output/", seqtype, "/al2co_scores"))) {
    tryCatch({
      data <- read.csv(file.path(paste0("output/", seqtype, "/al2co_scores"), file), header = FALSE, sep = "")
      mean_score <- mean(data$V3, na.rm = TRUE)
      
      mean_conservation_scores <- rbind(mean_conservation_scores,
                                        data.frame(gene = gsub("\\.\\w+$", "", file),
                                                   mean_score = mean_score,
                                                   stringsAsFactors = FALSE))
      
    }, error=function(cond) {
      print(paste0("Error with file ", file))
    })
  }
  
  hs_internal_names = read.csv("../1_core_response/input/hs_core_genes_internal_names.csv")
  mean_conservation_scores = merge(mean_conservation_scores, hs_internal_names, by.x = "gene", by.y = "target", all.x = TRUE)
  
  core_genes = mean_conservation_scores[!is.na(mean_conservation_scores$internal.name),]
  baseline_genes = mean_conservation_scores[is.na(mean_conservation_scores$internal.name),]
  
  print(paste0(seqtype, " conservation (HS core genes vs random baseline)"))
  print(t.test(core_genes$mean_score, baseline_genes$mean_score))
  print(paste0(seqtype, " median conservation score"))
  print(paste0("core genes: ", median(core_genes$mean_score)))
  print(paste0("random baseline: ", median(baseline_genes$mean_score)))
  
  # Plot mean_conservation_scores
  plot <- ggplot(mean_conservation_scores, aes(x = jitter(as.numeric(factor(gene))), y = mean_score)) +
      geom_point(aes(color = ifelse(!is.na(internal.name), "core gene", "random baseline")), show.legend = FALSE) +
      geom_text(aes(label = internal.name), color = "red", vjust = -0.5) +
      labs(x = "Gene", y = "Mean Conservation Score") +
      ggtitle(seqtype) +
      theme_minimal()
  print(plot)

}
dev.off()

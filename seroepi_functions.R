require(ggplot2)
require(tidyverse)
require(dplyr)
require(patchwork)
require(ggridges)

parseModelledEstimates <- function(global_post, region_post, fixOnames=FALSE, median=FALSE) {
  
  # get global posterior
  bayes_global_post <- read_csv(global_post)
  
  # get region posterior
  bayes_region_post <- read_csv(region_post)
  
  # update O names if required
  if (fixOnames) {
    bayes_global_post <- bayes_global_post %>%
     # mutate(locus=replace(locus, locus=="O1.v1", "O1v1")) %>%
     # mutate(locus=replace(locus, locus=="O1.v2", "O1v2")) %>%
     # mutate(locus=replace(locus, locus=="O1.v3", "O1v3")) %>%
      mutate(locus=replace(locus, locus=="O2a.v1", "O2a")) %>%
     # mutate(locus=replace(locus, locus=="O2a.v3", "O2av3")) %>%
      mutate(locus=replace(locus, locus=="O2afg.v2", "O2afg"))
    
    bayes_region_post <- bayes_region_post %>%
     # mutate(locus=replace(locus, locus=="O1.v1", "O1v1")) %>%
     # mutate(locus=replace(locus, locus=="O1.v2", "O1v2")) %>%
     # mutate(locus=replace(locus, locus=="O1.v3", "O1v3")) %>%
      mutate(locus=replace(locus, locus=="O2a.v1", "O2a")) %>%
     # mutate(locus=replace(locus, locus=="O2a.v3", "O2av3")) %>%
      mutate(locus=replace(locus, locus=="O2afg.v2", "O2afg"))
  }
  
  # get global medians, means and 95% CI; rank by median value
  bayes_global_est <- bayes_global_post %>% group_by(locus) %>% 
    summarise(median=median(prob), mean=mean(prob), 
              lower=quantile(prob, 0.025), upper=quantile(prob, 0.975)) %>%
    arrange(-mean) %>%
    mutate(rank=row_number())
  
  if (median) {
    bayes_global_est <- bayes_global_est %>%
      arrange(-median) %>%
      mutate(rank=row_number())
  }
  
  # get ranks based on adjusted global estimates
  locus_rank <- bayes_global_est %>% 
    select(locus, rank)
  
  bayes_region_est <- bayes_region_post %>% group_by(locus, subgroup) %>%
    summarise(median=median(prob), mean=mean(prob),
              lower=quantile(prob,0.025), upper=quantile(prob, 0.975))
  
  return(list(locus_rank=locus_rank, 
              global_est=bayes_global_est, 
              global_post=bayes_global_post,
              region_est=bayes_region_est, 
              region_post=bayes_region_post
  ))
}

locus_barplot <- function(estimates, ranks, maxRank=30, lines_every=10,
                          col="skyblue", error_bars=TRUE,
                          xtitle="Global prevalence (%)",
                          title="a) Global prevalence",
                          median=FALSE) {
  
  if (median) { estimates <- estimates %>% mutate(mean=median) }
  
  plot <- estimates %>% 
    filter(rank<=maxRank) %>%
    ggplot(aes(x=mean*100, y=locus)) +
    geom_col(fill=col) + 
    scale_y_discrete(limits=rev(ranks$locus[1:maxRank]))+ 
    theme_bw() + 
    labs(y=NULL, x=xtitle, title=title)
  
  if (error_bars) {
    plot <- plot + geom_errorbar(aes(xmin = lower*100, xmax=upper*100))
  }
  
  # add horizontal separator lines
  if (!is.null(lines_every)) {
    if (lines_every < maxRank) {
      lines_at <- seq((lines_every+0.5), maxRank, by=lines_every)
      for (y in lines_at) {
        plot <- plot + geom_hline(yintercept=y, col="darkgrey")
      }
    }
  }
  
  print(plot)
  
  return(plot)
}


locus_ridgesplot <- function(posterior, ranks, maxRank=30, lines_every=10,
                          col="#d17187", scale=3, axis_font_size=9,
                          rel_min_height = 0.01, linewidth=0.2,
                          xlim=c(0,10), xbreaks=seq(0,10,2),
                          xtitle="Global prevalence (%)",
                          title="a) Global prevalence") {

  plot <- posterior %>% 
    left_join(ranks) %>%
    filter(rank<=maxRank) %>%
    ggplot(aes(x=prob*100, y=locus)) +
    geom_density_ridges(scale=scale, 
                        fill=col, 
                        rel_min_height=rel_min_height,
                        linewidth=linewidth) + 
    scale_y_discrete(limits=rev(ranks$locus[1:maxRank])) + 
    scale_x_continuous(limits=xlim, breaks=xbreaks) + 
    theme_bw() + 
    theme(axis.text.x = element_text(size=axis_font_size), 
          axis.text.y = element_text(size=axis_font_size)) + 
    labs(y=NULL, x=xtitle, title=title)
  
  # add horizontal separator lines
  if (!is.null(lines_every)) {
    if (lines_every < maxRank) {
      lines_at <- seq((lines_every+1), maxRank, by=lines_every)
      for (y in lines_at) {
        plot <- plot + geom_hline(yintercept=y, col="darkgrey")
      }
    }
  }
  
  print(plot)
  
  return(plot)
}

comparative_locus_ridgesplot <- function(posterior1, posterior2, ranks,
                             type1="cluster-adjusted", type2="raw counts",
                             maxRank=30, lines_every=10,
                             col1="#d17187", col2="navy", axis_font_size=9,
                             scale=5, alpha=0.5, rel_min_height = 0.01,
                             xlim=c(0,10), xbreaks=seq(0,10,2),
                             xtitle="Global prevalence estimate (%)",
                             title="Global prevalence", legend_title="Prevalence") {
  
  plot <- posterior1 %>% mutate(type=type1) %>%
    bind_rows(posterior2 %>% mutate(type=type2)) %>%
    left_join(ranks) %>%
    filter(rank<=maxRank) %>%
    ggplot(aes(x=prob*100, y=locus)) +
    geom_density_ridges(aes(fill=type), scale=scale, alpha=alpha, col=NA,
                        rel_min_height = rel_min_height) + 
    scale_y_discrete(limits=rev(ranks$locus[1:maxRank])) + 
    scale_x_continuous(limits=xlim, breaks=xbreaks) + 
    scale_fill_manual(values=c(col1, col2)) +
    theme_bw() + 
    theme(axis.text.x = element_text(size=axis_font_size), 
          axis.text.y = element_text(size=axis_font_size)) + 
    labs(y=NULL, x=xtitle, title=title, fill=legend_title)

  # add horizontal separator lines
  if (!is.null(lines_every)) {
    if (lines_every < maxRank) {
      lines_at <- seq((lines_every+1), maxRank, by=lines_every)
      for (y in lines_at) {
        plot <- plot + geom_hline(yintercept=y, col="darkgrey")
      }
    }
  }

  print(plot)
  
  return(plot)
}


region_heatmap <- function(estimates, global, ranks,
                             maxRank=30, min_col="white", max_col="#d17187", 
                             rel_min_height = 0.01, linewidth=0.2,
                             label_size=2.5, axis_font_size=9, xaxis_angle=45,
                             region_order=c("Global", "Western Africa", "Southern Africa", "Eastern Africa", "Southern Asia"),
                             xtitle="Region",
                             title="b) Regional prevalence",
                             legend_title="Mean\nprevalence (%)", median=FALSE,
                            includeGlobal=TRUE) {
  
  if (is.null(region_order)) {region_order <- unique(estimates$subgroup)}
  
  if (includeGlobal) {estimates <- bind_rows(estimates, global %>% mutate(subgroup="Global"))}
  else {region_order <- region_order[-1]}
  
  if (median) {estimates <- estimates %>% mutate(mean=median)}
  
  plot <- estimates %>%
    ggplot(aes(x=subgroup, y=locus, fill=mean*100)) + 
    geom_tile() + 
    scale_y_discrete(limits=rev(ranks$locus[1:maxRank])) +
    scale_x_discrete(limits=region_order) +
    scale_fill_gradient(low=min_col, high=max_col) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle=xaxis_angle, hjust=1, size=axis_font_size), 
          axis.text.y = element_text(size=axis_font_size), legend.position="right") + 
    labs(x=xtitle, y="", fill=legend_title, title=title) + 
    geom_text(aes(x=subgroup, y=locus, label=round(mean*100,1)), size=label_size)
  
  print(plot)
  
  return(plot)
}


# get cumulative coverage estimates by adding posteriors
getCumulativeCoverage <- function(ranks, posterior, maxRank=30) {
  
  coverage <- tibble(
    rank = integer(0),
    mean = double(0),
    median = double(0),
    lower = double(0),
    upper = double(0)
  )
  
  for (i in 1:maxRank) {
    k <- ranks %>% filter(rank <=i) %>% pull(locus)
    draws <- posterior %>% filter(locus %in% k) %>% group_by(draw_id) %>%
      summarise(cumsum=sum(prob))
    coverage <- tibble(rank = i, 
                       mean = mean(draws$cumsum), 
                       median = median(draws$cumsum), 
                       lower = quantile(draws$cumsum, 0.025),
                       upper = quantile(draws$cumsum, 0.975),
                       locus = k[i]
    ) %>%
      bind_rows(coverage)
  }
  
  coverage <- coverage %>% arrange(rank)
  
  return(coverage)
}

getCumulativeCoverageGlobalRegional <- function(bayes, ranks, maxRank=30) {
  
  # global coverage
  coverage <- getCumulativeCoverage(ranks %>% ungroup(), bayes$global_post %>% ungroup(), maxRank=maxRank) %>% 
    mutate(subgroup="Global")
  
  for (region in unique(bayes$region_post$subgroup)) {
    regional <- getCumulativeCoverage(ranks %>% ungroup(), bayes$region_post %>% filter(subgroup==region) %>% ungroup(), maxRank=maxRank) %>% 
                mutate(subgroup=region)
    coverage <- bind_rows(coverage, regional)
  }
  
  return(coverage)
}


# input = raw data object with K_locus, Cluster, and Site columns
# output = per-site raw and cluster-adjusted counts and proportions
#          (suitable for input to function coverageAnalysis)
prepCoverageData <- function(data) {
  dat2 <- data %>%
    raw_adj_prop(grouping_vars = c("K_locus"), summarise_by = "K_locus", 
                 adj_vars = c("Cluster", "Site")) %>%
    mutate(raw_sum=nrow(data)) %>%
    mutate(adj_sum=length(unique(data$Cluster))) %>%
    mutate(raw_prop = raw_count/raw_sum) %>%
    mutate(adj_prop = adj_count/adj_sum) %>%
    mutate(rawRank=row_number(desc(raw_count))) %>%
    mutate(adjRank=row_number(desc(adj_count))) %>%
    mutate(raw_lower = raw_prop - 1.96*sqrt(raw_prop*(1-raw_prop)/nrow(data))) %>%
    mutate(raw_upper = raw_prop + 1.96*sqrt(raw_prop*(1-raw_prop)/nrow(data))) %>%
    mutate(adj_lower = adj_prop - 
             1.96*sqrt(adj_prop*(1-adj_prop)/length(unique(data$Cluster)))) %>%
    mutate(adj_upper = adj_prop + 
             1.96*sqrt(adj_prop*(1-adj_prop)/length(unique(data$Cluster)))) %>%
    rename(locus=K_locus) 
  
  return(dat2)
}

prepCoverageDataO <- function(data) {
  dat2 <- data %>%
    raw_adj_prop(grouping_vars = c("O_genotype"), summarise_by = "O_genotype", 
                 adj_vars = c("Cluster", "Site")) %>%
    mutate(raw_sum=nrow(data)) %>%
    mutate(adj_sum=length(unique(data$Cluster))) %>%
    mutate(raw_prop = raw_count/raw_sum) %>%
    mutate(adj_prop = adj_count/adj_sum) %>%
    mutate(rawRank=row_number(desc(raw_count))) %>%
    mutate(adjRank=row_number(desc(adj_count))) %>%
    mutate(raw_lower = raw_prop - 1.96*sqrt(raw_prop*(1-raw_prop)/nrow(data))) %>%
    mutate(raw_upper = raw_prop + 1.96*sqrt(raw_prop*(1-raw_prop)/nrow(data))) %>%
    mutate(adj_lower = adj_prop - 
             1.96*sqrt(adj_prop*(1-adj_prop)/length(unique(data$Cluster)))) %>%
    mutate(adj_upper = adj_prop + 
             1.96*sqrt(adj_prop*(1-adj_prop)/length(unique(data$Cluster)))) %>%
    rename(locus=O_genotype) 
  
  return(dat2)
}

## function to filter data
filterData <- function(data, studies=studies, years=years) {
  
  data_NNS <- data %>% 
    filter(is.na(Year) | Year %in% years) %>%
    filter(species =="Klebsiella pneumoniae") %>%
    filter(Study %in% studies) %>%
    filter(neonatal == "neonate")
  
  data_NNS_qc <- data_NNS %>%
    filter(total_size >= 5E6 & total_size <= 6.2E6) %>%
    filter(contig_count <= 1000) %>%
    mutate(O_type = str_replace(O_type, "unknown \\(", "")) %>%
    mutate(O_type = str_replace(O_type, "\\)", "")) %>%
    mutate(O_locus_variant = case_when(
      O_type %in% c("O1a", "O1ab") ~ paste0("O1 (", O_locus, ")"),
      O_type %in% c("O2a", "O2afg") ~ paste0(O_type, " (", O_locus, ")"),
      TRUE ~ O_locus
    ))
  
  # filter to sites with N≥10 included genomes
  sites_n10 <- data_NNS_qc %>% 
    group_by(Study, Country, Site) %>% 
    summarise(n=n()) %>% 
    filter(n>=10)
  
  data_NNS_sites10 <- data_NNS_qc %>% filter(Site %in% sites_n10$Site)
  
  return(list(filtered=data_NNS_sites10, qc=data_NNS_qc, sites_n10=sites_n10, data_NNS=data_NNS))
}

coverageAnalysis <- function(global_prevalence, ranks=bayes_global_K_adjRanks, maxRank=30, labelK=T) {
  
  cov_table <- global_prevalence %>% left_join(ranks) %>% 
    arrange(rank) %>% 
    filter(rank <= maxRank) %>%
    mutate(adj_cov = cumsum(adj_prop)) %>%
    mutate(raw_cov = cumsum(raw_prop)) %>%
    mutate(raw_lower=cumsum(raw_lower)) %>%
    mutate(raw_upper=cumsum(raw_upper)) %>%
    mutate(adj_lower=cumsum(adj_lower)) %>%
    mutate(adj_upper=cumsum(adj_upper)) %>% 
    mutate(adj_upper = if_else(adj_upper>1, 1, adj_upper)) %>%
    mutate(raw_upper = if_else(raw_upper>1, 1, raw_upper))
  
  if (!(ranks$locus[maxRank] %in% cov_table$locus)) {
    last_row <- cov_table[nrow(cov_table),]
    last_row <- last_row %>% mutate(locus=ranks$locus[maxRank], rank=maxRank)
    cov_table <- bind_rows(cov_table, last_row)
  }
  
  cov_plot <- cov_table %>% 
    select(-c(adj_prop, raw_prop)) %>% 
    pivot_longer(cols=c(adj_cov, raw_cov, raw_lower, raw_upper, adj_lower, adj_upper), 
                 names_sep="_", names_to=c("type", "est")) %>%
    select(locus, rank, type, est, value) %>%
    pivot_wider(id_cols=c(locus, rank, type), names_from=est) %>%
    mutate(type=fct_relevel(type, c("raw", "adj"))) %>%
    ggplot(aes(x=rank, y=cov*100, group=type)) +
    geom_line(aes(linetype=type, col=type)) + 
    geom_ribbon(aes(ymin=lower*100, ymax=upper*100, fill=type), alpha=0.1, linetype=0) +
    theme_bw() +
    ylim(0,100) + 
    labs(x=NULL, y="Cumulative coverage (%)", fill="", colour="", linetype="") +
    theme(axis.text.x = element_text(angle=90))
  
  if (labelK) {
    #Klist <- cov_table$locus
    #names(Klist) <- cov_table$rank
    #Klist <- Klist[as.character(1:maxRank)]
    #names(Klist) <- 1:maxRank
    #cov_plot <- cov_plot + scale_x_continuous(breaks=seq(1,nrow(cov_table)), labels=cov_table$locus)
    cov_plot <- cov_plot + scale_x_continuous(breaks=seq(1,maxRank), labels=ranks$locus[1:maxRank])
  }
  
  cov_table_long <- cov_table %>% 
    select(locus, rank, adj_cov, raw_cov) %>% 
    pivot_longer(cols=c(adj_cov, raw_cov), names_sep="_", names_to=c("type", "est")) 
  
  return(list(cov_table=cov_table, cov_table_long=cov_table_long, cov_plot=cov_plot))
}

plotCoverage <- function(coverage, ranks, maxRank=20, y_title="Cumulative coverage (%)", alpha=0.08, xintercept=c(5,10,15,20), cols=region_cols, col_title="", linetype_title="") {
  
  plot <- coverage %>% 
    mutate(upper=if_else(upper>1,1,upper)) %>% # max value for plotting interval using geom_ribbon is 1
    filter(rank <= maxRank) %>%
    ggplot(aes(x=rank, y=mean*100, group=subgroup)) +
    geom_hline(yintercept=70, linetype=2) + 
    geom_vline(xintercept=xintercept, col="grey") +
    geom_line(aes(col=subgroup), lwd=1) + 
    geom_ribbon(aes(ymin=lower*100, ymax=upper*100, fill=subgroup), alpha=alpha, linetype=0) +
    theme_bw() +
    ylim(0,100) + 
    labs(x=NULL, y=y_title, fill=col_title, colour=col_title, linetype="") +
    theme(axis.text.x = element_text(angle=90), panel.grid.minor.x = element_blank()) +
    scale_x_continuous(breaks=seq(1,maxRank), labels=ranks$locus[1:maxRank]) + 
    scale_color_manual(values=cols) + 
    scale_fill_manual(values=cols)
  
  print(plot)
  
  return(plot)
}

timelineBarplot <- function(data, region, xlimits=c(2012,2024), xbreaks=seq(2013,2023,2), ylim=c(0,250),
                            legend_title=NULL, xlab="Year", ylab="", legend.position="right",
                            cols=unname(alphabet()), legend_size=unit(6,"pt"), axis_text=8, title_size=10) {
  plot <- data %>% 
    filter(Region==region) %>%
    mutate(label=paste0(Country, ", ",Study)) %>%
    ggplot(aes(x=Year, fill=label)) + 
    geom_bar() + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = axis_text, 
                                   colour = "black"), 
        axis.title = element_text(size = axis_text, colour = "black"),
        axis.text.y = element_text(size = axis_text, colour = "black"),
        plot.title = element_text(hjust=0, size=title_size),
        legend.position = legend.position, 
        legend.text = element_text(size=7), legend.key.size=legend_size) +
    labs(x=xlab, y=ylab, fill=legend_title, title=region) + 
    scale_x_continuous(limits=xlimits, breaks=xbreaks) + 
    ylim(ylim) + 
    scale_fill_manual(values=cols)

    print(plot)
  
  return(plot)
}
library(shiny)
library(shinydashboard)
library(tidyverse)
library(foreach)
library(glue)
library(DT)
library(recipes)
library(heatmaply)
library(RColorBrewer)
library(plotly)
library(embed)

#### UI

ui <- dashboardPage(
  dashboardHeader(title = "TAC-seq expression app"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Input", tabName = "input"),
      menuItem("Normalization", tabName = "normalization"),
      menuItem("Visualization", tabName = "visualization"),
      menuItem("Scoring", tabName = "scoring"),
      tags$a(
        href = "https://github.com/seqinfo/tac-seq-expression-app", "GitHub",
        align = "center", style = "
        position:absolute;
        bottom:0;
        width:100%;
        weight:40px;   /* Height of the footer */
        color: white;
        padding: 10px;
        background-color: black;
        z-index: 1000;"
      )
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(
        tabName = "input",
        h1("Input"),
        p("Use TAC-seq data analysis output file as an input."),
        fileInput("counts", label = "Choose input file(s):", multiple = TRUE,
          accept = "text"
        ),
        selectInput(
          "target_list", label = "Choose target list or file:",
          choices = list(
            "Choose one" = "",
            "READY 61 targets" = "data/targets/READY_61targets.tsv",
            "READY 72 targets" = "data/targets/READY_72targets.tsv"
          )
        ),
        fileInput("target_file", label = "", accept = "text"),
        h3("Counts"),
        plotOutput("count_plot"),
        dataTableOutput("counts"),
        h3("Targets"),
        dataTableOutput("targets")
      ),
      tabItem(
        tabName = "normalization",
        h1("Normalization"),
        p("Molecule counts are normalized by geometric mean of housekeeping genes."),
        h3("Housekeeping genes"),
        plotOutput("housekeepers"),
        h3("Normalized counts"),
        p("Samples with geometric mean of zero are removed."),
        dataTableOutput("norm_biomarkers"),
        uiOutput("download")
      ),
      tabItem(
        tabName = "visualization",
        h1("Visualization"),
        p("Visualizing the normalized molecule counts of targeted biomarkers."),
        selectInput(
          "control_list", label = "Choose control list or file (optional):",
          choices = list(
            "Choose one" = "",
            "READY controls" = "data/controls/READY_controls.tsv",
            "READY controls (small set)" = "data/controls/READY_controls_small.tsv",
            "READY HRT controls (400k reads)" = "data/controls/READY_HRT_controls_400k_reads.tsv",
            "READY HRT controls (all reads)" = "data/controls/READY_HRT_controls_all_reads.tsv",
            "READY HRT controls (read count)" = "data/controls/READY_HRT_controls_read_count.tsv",
            "READY PCOS controls (400k reads)" = "data/controls/READY_PCOS_controls_400k_reads.tsv",
            "READY PCOS controls (filtered outliers)" = "data/controls/READY_PCOS_controls_400k_filtered.tsv"
          )
        ),
        fileInput("control_file", label = "", accept = "text"),
        h3("Controls"),
        dataTableOutput("controls"),
        h3("PCA"),
        plotlyOutput("pca", width = "auto", height = "auto"),
        h3("UMAP"),
        plotlyOutput("umap", width = "auto", height = "auto"),
        h3("Heatmap"),
        plotlyOutput("heatmap", width = "auto", height = "auto")
      ),
      tabItem(
        tabName = "scoring",
        h1("Scoring"),
        p("Score and classify the input samples depending on the reference groups."),
        selectInput(
          "scoring_ref", label = "Choose reference for scoring:",
          choices = list(
            "Choose one" = "",
            "READY PCOS controls (400k reads)" = "data/controls/READY_PCOS_controls_400k_reads.tsv",
            "READY PCOS controls (filtered outliers)" = "data/controls/READY_PCOS_controls_400k_filtered.tsv"
          )
        ),
        selectInput(
          "gene_subset", label = "Choose the gene subset to apply:",
          choices = list(
            "Choose one" = "",
            "LH+10 vs. LH+7 differential genes" = "data/targets/sig_set72.tsv",
            "All READY 72 targets" = "data/targets/READY_72targets.tsv",
            "All READY 61 targets" = "data/targets/READY_61targets.tsv"
          )
        ),
        h3("Scoring results"),
        dataTableOutput("scores_table"),
        uiOutput("download_score"),
        h3("PCA of filtered set"),
        plotlyOutput("pca_filt", width = "auto", height = "auto"),
      )
    )
  )
)

#### Server

server <- function(input, output) {

# counts ------------------------------------------------------------------

  counts <- reactive({
    count_file <- input$counts

    req(count_file, targets())

    n_targets = nrow(targets())

    validate(
      need(
        try(
          counts <- count_file$datapath %>%
            set_names(nm = count_file$name) %>%
            map_dfr(read_tsv, .id = "file") %>%
            filter(!str_detect(sample, "Undetermined"),  # remove undetermined samples
                   locus != "unmatched") %>%  # remove unmatched loci
            right_join(targets(), by = c("locus" = "target")) %>%
            complete(nesting(file, sample), nesting(locus, type)) %>%
            group_by(sample) %>%
            # filter(n() == n_targets,  # remove duplicated samples
            #        !any(is.na(molecule_count))) %>%  # remove samples with missing targets
            mutate(hk_geo_mean = exp(mean(log(molecule_count[type == "housekeeper"]))),  # geometric mean of housekeeping genes
                   norm_molecule_count = molecule_count / hk_geo_mean) %>%
            ungroup()
        ),
        "Incorrect input file(s). Please choose correct TAC-seq count file(s) with columns 'sample', 'locus', and 'molecule_count'."
      )
    )

    validate(
      need(
        counts %>%
          count(sample) %>%
          filter(n != n_targets) %>%
          nrow() == 0,
        counts %>%
          count(sample) %>%
          filter(n != n_targets) %>%
          transmute(duplicated = str_c(sample, " is duplicated")) %>%
          pull()
      ),
      need(
        !anyNA(counts$molecule_count),
        counts %>%
          filter(is.na(molecule_count), !is.na(sample)) %>%
          transmute(missing = str_c(sample, " is missing a ", locus, " locus")) %>%
          pull()
      )
    )

    counts
  })

  output$counts <- renderDataTable(counts())


# targets -----------------------------------------------------------------

  targets <- reactive({
    target_file <- input$target_file
    target_list <- input$target_list

    if (isTruthy(target_file)) {
      validate(
        need(
          try(
            targets <- read_tsv(target_file$datapath) %>%
              select(target, type)
          ),
          "Incorrect target file. Please choose correct target file with columns 'target', and 'type'."
        )
      )
    } else if (isTruthy(target_list)) {
      targets <- read_tsv(target_list)
    } else {
      req(FALSE)
    }

    targets
  })

  output$targets <- renderDataTable(
    targets()
  )


# plot counts -------------------------------------------------------------

  output$count_plot <- renderPlot(
    counts() %>%
      ggplot(aes(sample, molecule_count)) +
      geom_boxplot() +
      geom_point(aes(color = locus), show.legend = FALSE) +
      scale_y_log10() +
      facet_wrap(vars(file), scales = "free") +
      labs(x = NULL, y = "molecule count", color = NULL) +
      theme(axis.text.x = element_text(vjust = 0.5, angle = 90))
  )


# plot housekeepers -------------------------------------------------------

  output$housekeepers <- renderPlot(
    counts() %>%
      filter(type == "housekeeper") %>%
      ggplot(aes(sample, molecule_count)) +
      geom_point(aes(color = locus)) +
      geom_errorbar(aes(y = hk_geo_mean, ymin = hk_geo_mean, ymax = hk_geo_mean),
                    size = 1) +
      facet_wrap(vars(file), scales = "free") +
      labs(x = NULL, y = "molecule count") +
      theme(axis.text.x = element_text(vjust = 0.5, angle = 90))
  )


# biomarker counts --------------------------------------------------------

  norm_biomarkers <- reactive({
    counts() %>%
      filter(type == "biomarker") %>%
      pivot_wider(id_cols = c(file, sample), names_from = locus,
                  values_from = norm_molecule_count) %>%
      filter(across(where(is.numeric), is.finite))
  })

  output$norm_biomarkers <- renderDataTable(
    norm_biomarkers(),
    options = list(scrollX = TRUE)
  )


# download normalized counts --------------------------------------------------

  output$download <- renderUI({
    req(norm_biomarkers())
    downloadButton("file", "Download normalized counts")
  })

  output$file <- downloadHandler(
    filename = "TAC-seq_normalized_counts.tsv",
    content = function(file) {
      norm_biomarkers() %>%
        write_tsv(file)
    }
  )


# gene subset list --------------------------------------------------------

  gene_subset <- reactive({
    req(input$gene_subset, norm_biomarkers())

    biomarkers <-
      norm_biomarkers() %>%
      select(where(is.numeric)) %>%
      colnames

    read_tsv(input$gene_subset) %>%
      filter(target %in% biomarkers)
  })


# scoring control set -----------------------------------------------------

  scoring_ref <- reactive({
    req(counts(), input$scoring_ref)
    read_tsv(input$scoring_ref)
  })


# sample gene subsetter ---------------------------------------------------

  filt_biomarkers <- reactive({
    req(norm_biomarkers(), gene_subset())

    filter_list <-
      gene_subset() %>% 
      filter(target %in% colnames(norm_biomarkers())) %>%
      .$target
     
    norm_biomarkers() %>%
      select(c("sample", filter_list))
  })


# scored counts -----------------------------------------------------------

  mah.dist <- function(x, centr, covar)
    sweep(as.matrix(x), 2L, centr, "-") %>% 
      { setNames(rowSums(. %*% solve(covar) * .), rownames(x)) }
  
  mah_dist_p <- function(sample_dat, centroids, covs, dgf)
    foreach(centr = centroids$group, .combine = "full_join") %do% {
      sample_dat %>%
        mutate("mah.dist_{centr}" := 
                 mah.dist(sample_dat %>% select(matches("PC[0-9]+$")),
                          centroids[centroids$group == centr %>% as.character, -1] %>% 
                            as.numeric, 
                          # check if covars are per group or just general
                          if (is.list(covs)) covs[[centr]] else covs)) %>%
        mutate("p.chisq_{centr}" :=
                 1 - pchisq(get(glue("mah.dist_{centr}")), dgf))
    }


  centrs <- function(comp_proj, 
                     inner_quant = 0, # don't filter IQR values
                     label_predict = "")
    comp_proj %>%  # calculate the centroids per reference groups >(label_one, label_two)
      filter(group != label_predict) %>%
      select(c(group, where(is.numeric))) %>%
      group_by(group) %>%
      # centroids from the inner quantile values of the dataset
      summarize_all( ~ mean(.x[.x >= quantile(.x, inner_quant) & 
                            .x <= quantile(.x, 1 - inner_quant)]))


  clust_dist <- function(comp_proj, 
                         label_predict,
                         n_dim = 3, 
                         inner_quant = 0) {
    comp_proj <-
      comp_proj %>% 
      # select the n_dim principal components
      select(-num_range("PC", (n_dim + 1):999))
  
    covs <- 
      comp_proj %>% # covariance matrix over the whole reference PCA
        filter(!(group %in% c(label_predict))) %>%
        select(where(is.numeric)) %>%
        cov

    comp_proj %>%
      filter(group %in% c(label_predict)) %>%
      mah_dist_p(centrs(comp_proj, label_predict = label_predict), 
                 covs, n_dim)
  }


  project_pca <- function(ref_set, 
                          sample_set, 
                          ref_groups_list = "all") {
    ref_set <-
      ref_set %>%
      filter(group %in%  # filter out groups if needed
             if (ref_groups_list != "all") c(ref_groups_list) 
             else ref_set$group %>% unique)

    ref_set %>%
      recipe() %>%
      step_normalize(all_numeric()) %>%
      step_pca(all_numeric()) %>%
      prep(strings_as_factors = F) %>% 
      bake(new_data = 
           sample_set %>%
             mutate(group = NA) %>%
             bind_rows(ref_set)) %>%
      mutate(group = group %>% as.character) %>%
      mutate(group = ifelse(is.na(group), "Input", group))
  }


  # TODO: put counts filtering to another object, 
  # as it is actually not needed in scoring functions
  scores <- reactive({

    # prepare reference set 
    sig_genes <- gene_subset()
    ref_set <- 
      scoring_ref() %>%
      select(c("sample", "group", 
               sig_genes$target %>% as.character)) %>%
      mutate(group = recode(group, 
                            "prolif" = "pre", 
                            "PE" = "pre", 
                            "LH+2" = "pre", 
                            "pre-receptive" = "pre",
                            "LH+7" = "rec",
                            "receptive" = "rec",
                            "receptive HRT" = "rec",
                            "post-receptive" = "post",
                            "LH+10" = "post"))

    # generate the projected PCA objects for all groups, pre and post comparisons
    proj_pca_all <-
      project_pca(ref_set,
                  filt_biomarkers())
    proj_pca_pre <-
      project_pca(ref_set,
                  filt_biomarkers(),
                  ref_groups_list = c("pre", "rec"))
    proj_pca_post <- 
      project_pca(ref_set,
                  filt_biomarkers(),
                  ref_groups_list = c("post", "rec"))

    # calculate the mah. distances, p-values and normalize
    dists_all <-
      clust_dist(proj_pca_all, "Input") %>%
      mutate(norm_pre =  `p.chisq_pre`  / (`p.chisq_rec` + `p.chisq_pre`)) %>%
      mutate(norm_post = `p.chisq_post` / (`p.chisq_rec` + `p.chisq_post`)) %>%
      mutate(norm_rec =  `p.chisq_rec`  / (`p.chisq_rec` + `p.chisq_pre`))
    dists_pre <-
      clust_dist(proj_pca_pre, "Input") %>%
      mutate(norm_pre =  `p.chisq_pre`  / (`p.chisq_rec` + `p.chisq_pre`)) %>%
      mutate(norm_rec =  `p.chisq_rec`  / (`p.chisq_rec` + `p.chisq_pre`))
    dists_post <-
      clust_dist(proj_pca_post, "Input") %>%
      mutate(norm_post = `p.chisq_post` / (`p.chisq_rec` + `p.chisq_post`)) %>%
      mutate(norm_rec =  `p.chisq_rec`  / (`p.chisq_rec` + `p.chisq_post`))

    sigif <- 0.05

    dists_all %>%
      transmute(
        `Sample Name` = sample,
        `Pre-receptive (%)` = 
          ifelse(norm_pre >= norm_post & `p.chisq_pre` > sigif & !is.na(norm_pre),
                 round( dists_pre$norm_pre * 100, 2), NA_real_),
        `Receptive (%)` = 
          case_when(`p.chisq_pre` < sigif & `p.chisq_rec` < sigif & `p.chisq_post` < sigif ~ NA_real_,
                    is.na(norm_pre) & is.na(norm_rec) & is.na(norm_post) ~ NA_real_,
                    norm_pre >= norm_post ~ round( dists_pre$norm_rec * 100, 2),
                    T ~ round( dists_post$norm_rec * 100, 2)),
        `Post-receptive (%)` = 
          ifelse(norm_pre < norm_post & `p.chisq_post` > sigif & !is.na(norm_post),
                 round( dists_post$norm_post * 100, 2), NA_real_)
        ) %>%
      mutate(
        `Predicted Group` =
          case_when(
            is.na(`Receptive (%)`) ~ "No result",
            !is.na(`Pre-receptive (%)`) & `Receptive (%)` < 25 ~ "Pre-receptive",
            !is.na(`Pre-receptive (%)`) & `Receptive (%)` < 50 ~ "Early-receptive",
            !is.na(`Post-receptive (%)`) & `Receptive (%)` < 25 ~ "Post-receptive",
            !is.na(`Post-receptive (%)`) & `Receptive (%)` < 50 ~ "Late-receptive",
            T ~ "Receptive"
          ),

        `Receptometer score` = case_when(
          is.na(`Receptive (%)`) ~ NA_real_,
          is.na(`Post-receptive (%)`) ~ 
            round( 100 - (`Pre-receptive (%)` / (`Pre-receptive (%)` * 2.4 + `Receptive (%)`)) * 100),
          is.na(`Pre-receptive (%)`) ~ 
            round( 200 - (`Post-receptive (%)` / (`Post-receptive (%)` * 2.4 + `Receptive (%)`)) * 100))
        )
  })

  output$scores_table <- renderDataTable(
    scores(),
    options = list(scrollX = TRUE)
  )


# download scores --------------------------------------------------

  output$download_score <- renderUI({
    req(scores())
    downloadButton("score_file", "Download scored results")
  })

  output$score_file <- downloadHandler(
    filename = "TAC-seq_scores.tsv",
    content = function(file) {
      scores() %>%
        write_tsv(file)
    }
  )


# controls ----------------------------------------------------------------

  controls <- reactive({
    control_file <- input$control_file
    control_list <- input$control_list

    req(targets())

    biomarkers <- targets() %>%
      filter(type == "biomarker") %>%
      pull(target)

    if (isTruthy(control_file)) {
      validate(
        need(
          try(
            controls <- read_tsv(control_file$datapath) %>%
              select(sample, group, !!biomarkers)
          ),
          "Incorrect control file. Please choose correct control file with columns 'sample', 'group', and column for each 'target'."
        )
      )
    } else if (isTruthy(control_list)) {
      validate(
        need(
          try(
            controls <- read_tsv(control_list) %>%
              select(sample, group, !!biomarkers) %>%
              mutate(group = factor(group, c("pre-receptive",
                                             "early-receptive",
                                             "receptive",
                                             "receptive HRT",
                                             "late-receptive",
                                             "post-receptive",
                                             "menstrual blood",
                                             "polyp",
                                             "PE",
                                             "LH+2",
                                             "LH+7",
                                             "LH+10")))
          ),
          "Missing target(s) in controls. Please choose correct controls or targets."
        )
      )
    } else {
      return(NULL)
    }

    controls
  })

  output$controls <- renderDataTable(controls(), options = list(scrollX = TRUE))


# train and test data -----------------------------------------------------

  train_data <- reactive(req(controls()))
  test_data <- reactive({
    req(controls())

    test_data <- bind_rows(norm_biomarkers(), controls())

    validate(
      need(
          test_data %>%
            count(sample) %>%
            filter(n > 1) %>%
            nrow() == 0,
          test_data %>%
            count(sample) %>%
            filter(n > 1) %>%
            transmute(duplicated = str_c(sample, " name collides with a control sample")) %>%
            pull()
      )
    )

    test_data
  })


# heatmap -----------------------------------------------------------------

  output$heatmap <- renderPlotly(
    test_data() %>%
      na_if(0) %>%
      mutate(across(where(is.numeric), log)) %>%
      select(-file) %>%
      column_to_rownames("sample") %>%
      heatmaply(colors = rev(brewer.pal(n = 7, name = "RdYlBu")),
                scale = "column", show_dendrogram = c(TRUE, FALSE),
                hide_colorbar = TRUE)
  )


# Filtered PCA ---------------------------------------------------------------------

  output$pca_filt <- renderPlotly({
    sig_genes <- gene_subset()
    ref_set <- 
      scoring_ref() %>%
      select(c("sample", "group", 
               sig_genes$target %>% as.character)) %>%
      mutate(group = recode(group, 
                            "prolif" = "pre", 
                            "PE" = "pre", 
                            "LH+2" = "pre", 
                            "pre-receptive" = "pre",
                            "LH+7" = "rec",
                            "receptive" = "rec",
                            "receptive HRT" = "rec",
                            "post-receptive" = "post",
                            "LH+10" = "post"))

    ref_pca <-
      project_pca(ref_set,
                  filt_biomarkers())

    plot_ly() %>%
      add_trace(data = ref_pca %>%
                  filter(group != "Input"),
                x = ~PC1, y = ~PC2, z = ~PC3,
                color = ~group, text = ~sample,
                marker = list(size = 5), 
                type = "scatter3d", mode = "markers") %>%
      add_trace(data = ref_pca %>%
                  filter(group == "Input"),
                x = ~PC1, y = ~PC2, z = ~PC3,
                color = ~group, text = ~sample,
                marker = list(size = 3, symbol = "x"), 
                type = "scatter3d", mode = "markers")
  })


# PCA ---------------------------------------------------------------------

  output$pca <- renderPlotly({
    ref_pca <-
      train_data() %>%
      recipe() %>%
      step_normalize(all_numeric()) %>%
      step_pca(all_numeric(), num_comp = 3) %>%
      prep(strings_as_factors = FALSE) %>% 
      bake(new_data = test_data()) %>%
      mutate(group = group %>% as.character) %>%
      mutate(group = ifelse(is.na(group), "Input", group))

    plot_ly() %>%
      add_trace(data = ref_pca %>%
                  filter(group != "Input"),
                x = ~PC1, y = ~PC2, z = ~PC3,
                color = ~group, text = ~sample,
                marker = list(size = 5), 
                type = "scatter3d", mode = "markers") %>%
      add_trace(data = ref_pca %>%
                  filter(group == "Input"),
                x = ~PC1, y = ~PC2, z = ~PC3,
                color = ~group, text = ~sample,
                marker = list(size = 3, symbol = "x"), 
                type = "scatter3d", mode = "markers")
  })


# UMAP --------------------------------------------------------------------

  output$umap <- renderPlotly({
    ref_umap <-
      train_data() %>%
      recipe() %>%
      step_normalize(all_numeric()) %>%
      step_string2factor(group) %>%
      step_umap(all_numeric(), outcome = vars(group), seed = c(1, 1)) %>%
      prep(strings_as_factors = FALSE) %>%
      bake(new_data = test_data()) %>%
      mutate(group = group %>% as.character) %>%
      mutate(group = ifelse(is.na(group), "Input", group))

    plot_ly() %>%
      add_trace(data = ref_umap %>% # reference samples
                  filter(group != "Input"),
                x = ~umap_1, y = ~umap_2,
                color = ~group, text = ~sample,
                marker = list(size = 7), 
                type = "scatter", mode = "markers") %>%
      add_trace(data = ref_umap %>% # input samples
                  filter(group == "Input"),
                x = ~umap_1, y = ~umap_2,
                color = ~group, text = ~sample,
                marker = list(size = 7, symbol = "x"), 
                type = "scatter", mode = "markers")
  })
}


shinyApp(ui = ui, server = server)

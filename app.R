library(shiny)
library(shinydashboard)
library(tidyverse)
library(DT)
library(recipes)
library(heatmaply)
library(RColorBrewer)
library(plotly)
library(embed)


ui <- dashboardPage(
  dashboardHeader(title = "TAC-seq expression app"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Input", tabName = "input"),
      menuItem("Normalization", tabName = "normalization"),
      menuItem("Visualization", tabName = "visualization"),
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
            "READY PCOS controls (400k reads)" = "data/controls/READY_PCOS_controls_400k_reads.tsv"
          )
        ),
        fileInput("control_file", label = "", accept = "text"),
        h3("Controls"),
        dataTableOutput("controls"),
        h3("Heatmap"),
        plotlyOutput("heatmap", width = "auto", height = 900),
        h3("PCA"),
        plotlyOutput("pca", width = "auto", height = 900),
        h3("UMAP"),
        plotlyOutput("umap", width = "auto", height = 900)
      )
    )
  )
)

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


# download ----------------------------------------------------------------

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


# PCA ---------------------------------------------------------------------

  output$pca <- renderPlotly({
    train_data() %>%
      recipe() %>%
      step_normalize(all_numeric()) %>%
      step_pca(all_numeric(), num_comp = 2) %>%
      prep(strings_as_factors = FALSE) %>%
      bake(new_data = test_data()) %>%
      ggplot(aes(PC1, PC2, color = group, shape = group, sample = sample)) +
      geom_point() +
      scale_shape_manual(values = 1:7)
    ggplotly()
  })


# UMAP --------------------------------------------------------------------

  output$umap <- renderPlotly({
    train_data() %>%
      recipe() %>%
      step_normalize(all_numeric()) %>%
      step_string2factor(group) %>%
      step_umap(all_numeric(), outcome = vars(group), seed = c(1, 1)) %>%
      prep(strings_as_factors = FALSE) %>%
      bake(new_data = test_data()) %>%
      ggplot(aes(umap_1, umap_2, color = group, shape = group, sample = sample)) +
      geom_point() +
      scale_shape_manual(values = 1:7)
    ggplotly()
  })
}


shinyApp(ui = ui, server = server)

# app.R
# Oncology Treatment Access Navigator - NSCLC (Phase 1 prototype)

library(shiny)
library(bslib)
library(DT)

# ---- 1. Load helper functions and data
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)

rules <- readr::read_csv("data/therapy_rules.csv",
                         col_types = readr::cols(.default = "c", min_value = "d"),
                         show_col_types = FALSE)
catalog       <- readr::read_csv("data/biomarker_catalog.csv", show_col_types = FALSE)
payer_rules   <- readr::read_csv("data/payer_rules.csv", show_col_types = FALSE)
demo_patients <- readr::read_csv("data/demo_patients.csv",
                                 col_types = readr::cols(.default = "c", age = "i", pdl1_tps = "d"),
                                 show_col_types = FALSE)
demo_maf <- read_maf("data/demo_cohort.maf")

dna_rna_markers <- catalog[catalog$biomarker_id != "PDL1", ]
marker_choices  <- setNames(dna_rna_markers$biomarker_id, dna_rna_markers$display_name)
stage_choices   <- c("Advanced / metastatic (IIIB-IV)", "Early stage, resected (IB-IIIA)")
payer_choices   <- c("Original Medicare", "Medicare Advantage", "Medicaid", "Commercial")

# ---- 2. User interface
ui <- page_sidebar(
  title = "Oncology Treatment Access Navigator - NSCLC (Phase 1 prototype)",
  theme = bs_theme(bootswatch = "flatly"),
  fillable = FALSE,
  
  sidebar = sidebar(
    width = 340,
    radioButtons("source", "Patient data",
                 choices = c("Demo patient" = "demo", "Upload a MAF file" = "upload")),
    conditionalPanel(
      "input.source == 'demo'",
      selectInput("demo_id", "Demo patient", choices = demo_patients$patient_id),
      textOutput("demo_story")
    ),
    conditionalPanel(
      "input.source == 'upload'",
      fileInput("maf_file", "NGS results file (.maf, .tsv or .txt)",
                accept = c(".maf", ".tsv", ".txt")),
      checkboxGroupInput("tested_upload", "Biomarkers covered by this test",
                         choices = marker_choices, selected = marker_choices)
    ),
    hr(),
    selectInput("histology", "Histology", choices = c("Nonsquamous", "Squamous")),
    selectInput("stage", "Stage", choices = stage_choices),
    radioButtons("line", "Line of therapy (advanced disease)",
                 choices = c("First-line", "Subsequent"), inline = TRUE),
    numericInput("pdl1", "PD-L1 TPS (%) - leave blank if not tested",
                 value = NA, min = 0, max = 100),
    selectInput("payer", "Payer", choices = payer_choices)
  ),
  
  navset_tab(
    nav_panel(
      "Biomarkers & Gaps",
      br(),
      uiOutput("gap_alert"),
      h5("Biomarker results"),
      DTOutput("status_table", fill = FALSE),
      br(),
      h5("Variants in the NGS report"),
      DTOutput("variant_table", fill = FALSE)
    ),
    nav_panel(
      "Therapy Options",
      br(),
      p("FDA-approved options that match the biomarker results. Chemotherapy and",
        "chemo-immunotherapy regimens that need no biomarker are not listed in Phase 1."),
      DTOutput("therapy_table", fill = FALSE)
    ),
    nav_panel(
      "Coverage & Prior Auth",
      br(),
      card(card_header("Biomarker testing coverage"), textOutput("test_note")),
      h5("Billing and typical prior-auth pattern"),
      DTOutput("coverage_dt", fill = FALSE),
      br(),
      h5("Prior-auth document checklist"),
      uiOutput("checklist")
    ),
    nav_panel(
      "About",
      br(),
      h4("About this prototype"),
      p("Portfolio project by DeAnna Genthner (UMGC, Bioinformatics). Uses synthetic patients only."),
      p("Therapy rules are hand-built from FDA prescribing information (DailyMed / openFDA).",
        "Each rule must be checked against the current label; see the 'Rule verified' column."),
      p("NCCN Guidelines are copyrighted and are not reproduced here. Clinicians should consult",
        "the current NCCN Guidelines for NSCLC directly."),
      p(strong("Not a medical device. For education and decision support only."),
        "All treatment decisions belong to the treating clinician.")
    )
  )
)

# ---- 3. Server (the logic that reacts to inputs)
server <- function(input, output, session) {
  
  current_patient <- reactive({
    demo_patients[demo_patients$patient_id == input$demo_id, ]
  })
  
  output$demo_story <- renderText(current_patient()$story)
  
  # When a demo patient is chosen, fill the sidebar with their details
  observeEvent(list(input$source, input$demo_id), {
    req(input$source == "demo")
    p <- current_patient()
    updateSelectInput(session, "histology", selected = p$histology)
    updateSelectInput(session, "stage", selected = p$stage_group)
    if (p$line %in% c("First-line", "Subsequent")) {
      updateRadioButtons(session, "line", selected = p$line)
    }
    pdl1_value <- p$pdl1_tps
    if (is.na(pdl1_value)) {
      pdl1_value <- ""
    }
    updateNumericInput(session, "pdl1", value = pdl1_value)
    updateSelectInput(session, "payer", selected = p$payer)
  })
  
  # Variants from the demo cohort or from an uploaded file
  variants <- reactive({
    if (input$source == "demo") {
      maf <- demo_maf[demo_maf$Tumor_Sample_Barcode == input$demo_id, ]
    } else {
      req(input$maf_file)
      maf <- try(read_maf(input$maf_file$datapath), silent = TRUE)
      validate(need(!inherits(maf, "try-error"),
                    "Could not read that file. It needs MAF columns: Hugo_Symbol, Variant_Classification, HGVSp_Short."))
    }
    annotate_variants(maf)
  })
  
  found_ids <- reactive({
    ids <- variants()$biomarker_id
    unique(ids[!is.na(ids)])
  })
  
  tested_ids <- reactive({
    if (input$source == "demo") {
      base <- expand_tested(current_patient()$biomarkers_tested, catalog)
    } else {
      base <- input$tested_upload
    }
    ids <- union(base, found_ids())
    if (!is.na(input$pdl1)) {
      ids <- union(ids, "PDL1")
    }
    ids
  })
  
  stage_setting <- reactive(setting_for_stage(input$stage))
  status  <- reactive(build_biomarker_status(catalog, found_ids(), tested_ids(), input$pdl1))
  gaps    <- reactive(find_testing_gaps(rules, catalog, stage_setting(), input$histology, tested_ids()))
  matches <- reactive(match_therapies(rules, status(), stage_setting(), input$histology, input$line))
  
  # ---- Tab 1: Biomarkers & Gaps
  output$gap_alert <- renderUI({
    g <- gaps()
    if (nrow(g) == 0) {
      div(class = "alert alert-success",
          "No testing gaps: every biomarker with an FDA-approved option in this setting has a result.")
    } else {
      div(class = "alert alert-warning",
          strong("Testing gaps: "),
          "no result on file for ", paste(g$display_name, collapse = ", "),
          ". Treatment selection may be premature.")
    }
  })
  
  output$status_table <- renderDT({
    s <- status()
    shown <- data.frame(Biomarker = s$display_name,
                        Result = s$result,
                        Value = ifelse(is.na(s$value), "", paste0(s$value, "%")))
    dt <- datatable(shown, rownames = FALSE, options = list(dom = "t", pageLength = 20))
    formatStyle(dt, "Result",
                backgroundColor = styleEqual(c("Positive", "Negative", "Not tested"),
                                             c("#d4edda", "#f8f9fa", "#fff3cd")))
  })
  
  output$variant_table <- renderDT({
    v <- variants()
    v$biomarker_id <- ifelse(is.na(v$biomarker_id), "Not actionable in this rule set", v$biomarker_id)
    datatable(v, rownames = FALSE, options = list(dom = "t"))
  })
  
  # ---- Tab 2: Therapy Options
  output$therapy_table <- renderDT({
    m <- format_matches(matches(), catalog)
    validate(need(nrow(m) > 0,
                  "No biomarker-matched FDA-approved options for this patient in the current rule set."))
    datatable(m, rownames = FALSE, options = list(pageLength = 15))
  })
  
  # ---- Tab 3: Coverage & Prior Auth
  output$test_note <- renderText(biomarker_test_note(input$payer))
  
  output$coverage_dt <- renderDT({
    ct <- coverage_table(matches(), payer_rules, input$payer)
    validate(need(nrow(ct) > 0, "No matched therapies to check yet."))
    datatable(ct, rownames = FALSE, options = list(pageLength = 15))
  })
  
  output$checklist <- renderUI({
    items <- pa_checklist(matches(), input$line, stage_setting())
    tags$ul(lapply(items, tags$li))
  })
}

shinyApp(ui, server)

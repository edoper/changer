# changer_app — UI

library(shiny)
library(shinythemes)
library(shinyWidgets)
library(DT)
library(shinycssloaders)  # (7) spinners
library("shinyalert")

app_title <- "CHANGER v1"

shinyUI(
  navbarPage(
    theme = shinytheme("flatly"),
    title = app_title,
    id = "tabs",
    selected = "home",
    header = tags$head(
      tags$style(HTML(".card-soft {background:#fff; border:1px solid #e5e7eb; border-radius:10px;box-shadow:0 2px 8px rgba(0,0,0,0.06);}")),
      tags$style(HTML("
      .st-hero { 
      position: relative; 
      min-height: calc(100vh - 120px); /* da altura para que quepa el mapa */
      }
      .st-hero::after {
      content: '';
      position: absolute; inset: 0;
      background-image: url('chile.png');
      background-repeat: no-repeat;
      /* centra verticalmente y encaja por ALTURA */
      background-position: right 2% center;
      background-size: auto 92%;
      opacity: .15; 
      filter: saturate(115%) contrast(110%);
      pointer-events: none;
      }
      /* Ajustes móviles/tablet: reduce tamaño para no tapar la búsqueda */
      @media (max-width: 992px){
      .st-hero { min-height: 70vh; }
      .st-hero::after {
      background-position: center bottom 8%;
      background-size: auto 65%;
      }
      }
      "))
    ),
    
    # [ST] Search Tab ----------------------------------------------------------
    tabPanel(
      title = "Home", value = "home",
      div(class = "st-hero",
          fluidRow(
            column(
              width = 8, offset = 2, align = "center",
              br(), br(), h3(app_title),
              p("CHilean Aggregated National GEnomics Resource."),
              br(),
              selectizeInput("gene_search", NULL, choices = character(0), selected = NULL,
                             options = list(placeholder = "Search gene symbol (e.g., AACS)", create = FALSE),
                             width = "60%"),
              actionButton("go_search", "Search", class = "btn btn-success"),
              br(), br(),
              p("CHANGER is a valuable resource for genetic research and a reference for local variant interpretation, by providing a more relevant genetic landscape of the Chilean general population.")
            )
          )
      )
    ),
    
    # [RT] Results Tab (PA, PB, PC, PD) ---------------------------------------
    tabPanel(
      title = "Results",
      value  = "results",
      fluidRow(
        # ROW 1: PC (2/12), PB (10/12)
        column(
          width = 2,
          wellPanel(
            class = "card-soft",
            h3("Gene summary"),
            DTOutput("gene_summary_card")
          )
        ),
        column(
          width = 10,
          wellPanel(
            class = "card-soft",
            uiOutput("pb_title"),
            withSpinner(plotOutput("pb_lollipop", height = "360px"))
          )
        )
      ),
      fluidRow(
        # ROW 2: PA (2/12), PD (10/12)
        column(
          width = 2,
          wellPanel(
            class = "card-soft",
            h3("Filters"),
            checkboxGroupInput(
              inputId = "filters",
              label   = "Show only:",
              inline  = FALSE,
              choices = c(
                "Loss of function (HC)" = "f_lofhc",
                "Missense damaging"    = "f_dmg",
                "Not in gnomAD"        = "f_nognomad",
                "Common (AF ≥ 1%)"     = "f_common"
              ),
              selected = NULL
            ),
            uiOutput("consequence_picker")
          )
        ),
        column(
          width = 10,
          wellPanel(
            class = "card-soft",
            h3("Variants Table View"),
            withSpinner(DTOutput("variants_table"))
          )
        )
      ),
      fluidRow(
        column(
          width = 12, align = "right", style = "margin-bottom:16px;",
          downloadButton("dl_gene_tsv", "Download TSV", class = "btn btn-default", style = "margin-right:8px;"),
          actionButton("new_search", "New Search", class = "btn btn-success"),
          br()
        )
      )
    ),
    # [AB] About Tab -----------------------------------------------------------
    tabPanel(
      title = "About",
      value  = "about",
      div(class = "st-hero",
          ##########
          fluidRow(
            column(width = 10, offset = 1, align = "center",
                   wellPanel("Citation and Contact",
                             style = "background-color: #b2baba;
                                        color: white;
                                        border-top-color: #b2baba;
                                        border-left-color: #b2baba;
                                        border-right-color: #b2baba;
                                        box-shadow: 3px 3px 3px #d8d8d8;
                                        margin-bottom: 0px; padding:5px"), 
                   wellPanel(p("Evelin González, Camilo Villamán, Boris Rebolledo-Jaramillo, Carlos Hernandez, Dominga Berrios, Gabriela Moreno, Cecilia Poli, Juan Francisco Calderón, Paula Muñoz-Venturelli, Mario I. Fernández, Juan Alberto Lecaros, Ricardo Armisen, Gabriela M. Repetto, Eduardo Pérez-Palma",
                               a("A novel genetic reference from 902 unrelated Chilean exomes to enhance South American medical genetics and research.", target="_blank", href="https://genome.cshlp.org/content/30/1/62"),
                               "MedRxiv, Dic 2025."),
                             br(), 
                             #p(icon("envelope", lib = "glyphicon"),"eduardoperez@udd.cl"),
                             style = "background-color: #ffffff;
                                         border-color: #b2baba;
                                         box-shadow: 3px 3px 3px #d8d8d8;
                                         margin-top: 0px")
            ) # Column
          ),#Fluid row
          fluidRow(
            column(width = 10, offset = 1, 
                   wellPanel("Abstract", align = "center",
                             style = "background-color: #2b3e51;
                                         color: white;
                                         border-top-color: #2b3e51;
                                         border-left-color: #2b3e51;
                                         border-right-color: #2b3e51;
                                         box-shadow: 3px 3px 3px #d8d8d8;
                                         margin-bottom: 0px; padding:5px"), 
                   wellPanel(br(),p("Studying underrepresented populations, such as those in South America, is essential for a comprehensive understanding of global human genetics. Historically, genetic research has focused predominantly on populations of European descent, resulting in disparities in knowledge and applicability of genetic findings across different ethnic groups. By focusing on the Chilean population, our work aims to bridge these gaps, providing more accurate reference data for genetic research and enhancing diagnostic and therapeutic strategies for underrepresented communities. This is particularly important for rare disorders and cancer research, where precise variant interpretation can significantly impact patient care and outcomes."),
                             p(""),
                             p("Here we present CHANGER (Chilean Aggregated National Genomics Resource), a comprehensive aggregation of 902 unrelated Chilean exomes designed to create a detailed reference of genetic variation in Chile. By incorporating data from multiple cohorts, we identified a total of 774,110 unique genetic variants, with an average of 42,601 aggregated variants per individual. We found a total of 132,363 novel variants, of which 31,470 were common within our cohort (frequency > 1 %). We provide variant annotation and direct comparison with other publicly available general population references."),
                             p(""),
                             p("CHANGER will become a valuable resource for genetic research and a reference for variant interpretation, by providing a more relevant genetic landscape of the Chilean general population. Our results are available at: changer.cl."),
                             p(""),
                             style = "background-color: #ffffff;
                                         border-color: #2b3e51;
                                         box-shadow: 3px 3px 3px #d8d8d8;
                                         margin-top: 0px")
            ) # column
          )
      )
    ),
    # Footer 
    fluidRow(column(width = 12, align = "center",
                    wellPanel(style = "background:#ffffff; padding: 20px 20px; border: 1px solid white; border-top: 1px solid grey",
                              p("Disclaimer: No Medical Advice. The information and Content provided are not substitutes for medical or professional
                              care, and you should not use the information in place of a visit, call, consultation or the advice of your
                              physician or other healthcare provider. You are liable or responsible for any advice, course of treatment, 
                              diagnosis or any other information, services or product obtained through this site.")
                    ),
                    fluidRow(column(width = 4, align = "center",
                                    p("Powered by Shiny")
                    ),
                    column(width = 4, align = "center",
                           p("Copyright ", icon("copyright"), "2025", "[Page under construction]")
                    ),
                    column(width = 4, align = "center",
                           p("Contact | ", 
                             a(icon("envelope"), target="_blank", href="mailto:eduardoperez@udd.cl"),
                             " | ",
                             a(icon("twitter"), target="_blank", href="https://x.com/ICIMUDD"),
                             " | ",
                             a(icon("github"), target="_blank", href="https://github.com/edoper/changer"))
                    )
                    )#Fluid row
    )#Column
    )#Fluid row 
    ###########footer a mano
  )
)

#' start the candiSNP app in a browser
#'
#' @examples
#' \dontrun{ #starts shiny app
#'  app()
#' }
#' @export
app <- function(){
  message("Starting candiSNP app in browser. When done use CTRL-C to quit the process and get the console back.")

  ui <- shiny::fluidPage(

    # Application title
    shiny::titlePanel("CandiSNP"),

    # Sidebar with a slider input for number of bins
    shiny::sidebarLayout(
      shiny::sidebarPanel(

        shiny::selectizeInput('genome', "Select genome", choices = get_installed_genomes() ),

        shiny::fileInput("upload", "Upload a file", accept = "csv"),

        shiny::sliderInput("allele_range",
                    "Allele Frequency Range:",
                    min = 0,
                    max = 1,
                    value = c(0,1) ),

        colourpicker::colourInput("non_syn_cod_ctga", "Non-Synonymous Coding CT/GA", "red", allowTransparent = TRUE),#CDS
        colourpicker::colourInput("non_syn_cod", "Non-Synonymous Coding", "red", allowTransparent = TRUE),
        colourpicker::colourInput("annot", "Annotated Region", "steelblue", allowTransparent = TRUE),
        colourpicker::colourInput("non_annot", "Non-Annotated Region", "steelblue", allowTransparent = TRUE),
        colourpicker::colourInput("cds", "CDS", "steelblue", allowTransparent = TRUE),
        shiny::radioButtons("hide_centromere", "Centromere SNPS:", c("Visible", "Hidden") ),
        shiny::downloadButton("download", label="Download annotations")
      ),


      # Show a plot of the generated distribution
      shiny::mainPanel(
        plotly::plotlyOutput("candiplot",width = "800px", height = "1500px")
      )
    )
  )


  server <- function(input, output) {

    output$files <- shiny::renderTable(input$upload)

    data <- shiny::reactive({
      shiny::req(input$upload)
    })

    snp_eff <- shiny::reactive({
      file <- data()$datapath
      shiny::withProgress(message = "Running SNPEff", value = 0, {
        shiny::incProgress(0.3, detail = "finding effects...")
        do_snp_eff(file, input$genome)
      })
    })




    output$candiplot <- plotly::renderPlotly({

      p <- do_centromere_filter(snp_eff(),  input$hide_centromere, input$genome) |>
        ggplot2::ggplot() +
        ggplot2::aes(.data$Pos, .data$Allele_Freq, label=.data$info) +
        ggplot2::geom_jitter(ggplot2::aes(colour=.data$type)) +
        ggplot2::facet_grid(.data$Chr ~ . , scales = "free_x", space="free") +
        ggplot2::ylim(input$allele_range) +
        ggplot2::scale_colour_manual(values = c(input$non_syn_cod_ctga, input$non_syn_cod, input$annot, input$non_annot, input$cds), drop=FALSE) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none") +
        ggplot2::theme(panel.spacing = ggplot2::unit(2, "lines"))
      plotly::ggplotly(p)
    })

    output$download <- shiny::downloadHandler(
      filename = function() {
        paste0(input$dataset, "candisnp_out.csv")
      },
      content = function(file) {
        readr::write_csv(snp_eff(), file)
      }
    )


  }

  # Run the application
  shiny::shinyApp(ui = ui, server = server)
}

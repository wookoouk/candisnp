
#' install the default genome annotations from SNPEff
#'
#' installs the following genome annotations from SNPEff:
#'      "athalianaTair9",
#'      "athalianaTair10",
#'      "gmax1.09v8",
#'      "grape",
#'      "maizeZmB73",
#'      "rice7"
#'
#' these being the ones originally available in candiSNP
#' @examples
#' \dontrun{ #installs 300Mb database
#'  install_default_genomes()
#'  }
#' @returns No value
#' @export
install_default_genomes <- function() {
  message("Installing default genomes. This can take a few minutes.")
  genomes = c(
    "athalianaTair9",
    "athalianaTair10",
    "gmax1.09v8",
    "grape",
    "maizeZmB73",
    "rice7"
  )
  download_annotation(genomes = genomes)
}

#' download annotations from SNPEff
#'
#' @param genomes character vector of annotation keys (genome names) as described in snpeff.config
#' @returns no value
#' @examples
#' \dontrun{ ##installs very large database
#' download_annotation('GRCh37.75')
#' }
#' @export
download_annotation <- function(genomes = NULL) {


  snp_eff <- system.file(file.path("snpEff", "snpEff.jar"), package="candiSNP")

  for (g in genomes){
    # States ERROR if no connection is made at first attempt, retries and is ok if no FATAL ERROR
    print(glue::glue("java -jar {snp_eff} download {g}"))
    system(glue::glue("java -jar {snp_eff} download {g}"))
    unlink(paste0("snpEff_v3_6_", g, ".zip"))

  }
}

#' runs SNPeff with a candiSNP input file and a preinstalled genome
#'
#' @param file path to .csv file (Chr, Pos, Ref, Alt, Allele_freq)
#' @param genome genome name (annotation key) for annotation to use
#' @export
#' @examples
#' \dontrun{ #cant work without installing a large db
#'  do_snp_eff("athal_sample.csv", "athalianaTair10")
#'  }
#' @returns dataframe
do_snp_eff <- function(file, genome) {

  snp_eff <- system.file(file.path( "snpEff", "snpEff.jar"), package="candiSNP")
  conf <- system.file(file.path( "snpEff", "snpEff.config"), package="candiSNP")
  tmpfile <- tempfile()

  raw <- readr::read_csv(file, show_col_types = FALSE)

  raw |>
    dplyr::select(.data$Chr, .data$Pos, .data$Ref, .data$Alt) |>
    dplyr::arrange(.data$Chr, .data$Pos, .data$Ref, .data$Alt) |>
    readr::write_tsv(file = tmpfile)


  message(glue::glue("Running SNPEff."))
  #cmd <- glue::glue("java -jar -Xmx4g {snp_eff} -c {conf} -i txt -o txt -noLog  -noStats -canon -snp -no-downstream -no-upstream -no-utr {genome} {tmpfile}")
  args <- c("-jar", "-Xmx2g", snp_eff, "-c", conf, "-i", "txt", "-o", "txt", "-noLog", "-noStats", "-canon", "-snp", "-no-downstream", "-no-upstream", "-no-utr", genome, tmpfile)
  #d <- system(cmd, ignore.stderr = TRUE, intern = TRUE)
  d <- system2("java", args, stdout=TRUE, stderr=FALSE)
  message(d[1])
  readr::read_delim(I(d), comment = "#", col_names = FALSE,show_col_types = FALSE) |>
    dplyr::select(.data$X1, .data$X2, .data$X3, .data$X4, .data$X11, .data$X16, .data$X17) |>
    dplyr::rename(Chr=.data$X1, Pos=.data$X2, Ref=.data$X3, Alt=.data$X4, gene=.data$X11, effect=.data$X16, nucs=.data$X17) |>
    dplyr::mutate(info = paste(.data$gene, .data$nucs, .data$effect)) |>
    dplyr::left_join(raw,by = c("Chr", "Pos", "Ref", "Alt") ) |>
    get_snp_type()
}


#' determines the candiSNP snp type for each SNP from SNPEff data
#' adds a new column 'type'
#'
#' @param df dataframe of parsed SNPEff data
#' @returns dataframe
get_snp_type <- function(df){

  df |> dplyr::mutate(
    is_synonymous = dplyr::case_when(
      effect == "NON_SYNONYMOUS_CODING" ~ FALSE,
      effect == "INTERGENIC" ~ FALSE,
      effect == "INTRON" ~ FALSE,
      TRUE ~ TRUE
    ),
    is_CTGA = dplyr::case_when(
      Ref == "C" & Alt == "T" ~ TRUE,
      Ref == "G" & Alt == "A" ~ TRUE,
      TRUE ~ FALSE
    ),
    in_cds = dplyr::case_when(
      effect == "INTERGENIC" ~ FALSE,
      effect == "INTRON" ~ FALSE,
      TRUE ~ TRUE
    ),

    type = dplyr::case_when(
        !is_synonymous & is_CTGA & in_cds ~ "NON_SYNONYMOUS_CODING_CT_GA",
        !is_synonymous & ! is_CTGA & in_cds ~ "NON_SYNONYMOUS_CODING",
        !is_synonymous & ! is_CTGA & in_cds ~ "ANNOTATED_REGION",
        in_cds ~ "CDS",
        !in_cds ~ "NON_ANNOTATED_REGION",
        TRUE ~ "NON_ANNOTATED_REGION"
      )
  ) |>
    dplyr::mutate(
      type = factor(.data$type, ordered=TRUE, levels = c("NON_SYNONYMOUS_CODING_CT_GA",
                                                   "NON_SYNONYMOUS_CODING",
                                                   "ANNOTATED_REGION",
                                                   "NON_ANNOTATED_REGION",
                                                   "CDS")
                    )
    )

}

#' list installed genomes
#'
#' @export
#' @examples
#'  get_installed_genomes()
#' @returns character vector of genomes available
get_installed_genomes <- function() {
  #c("athal10", "athal9")
  g <- list.dirs(system.file(file.path("snpEff", "data"), package="candiSNP"), full.names = FALSE )
  g <- g[g != ""]
  if (length(g) == 0){
    return(c("no genomes installed! Please see documentation before proceeding"))
  }
  g

}

#' removes points within +/- 0.5 Mb of the centromere in Arabidopsis plots in the app
#'
#' @param df parsed SNPEff data
#' @param centromere_state "Hidden" or "Visible" as to whether centromere SNPs should be shown
#' @param genome SNPEff genome in use.
#' @returns dataframe
#' @importFrom rlang .data
do_centromere_filter <- function(df, centromere_state, genome) {

  if (centromere_state == "Visible" | ! genome %in% c("athalianaTair9", "athalianaTair10", "athaliana130")){
    return(df)
  }
  else {
    df |>
      dplyr::mutate(
        is_centromere = dplyr::case_when(
          Chr == 1 & (Pos < 15086545 + 500000 & Pos > 15086545 - 500000) ~ TRUE,
          Chr == 2 & (Pos < 3608429 + 500000 & Pos > 3608429 - 500000) ~ TRUE,
          Chr == 3 & (Pos < 14209452 + 500000 & Pos > 14209452 - 500000) ~ TRUE,
          Chr == 4 & (Pos < 3956521 + 500000 & Pos > 3956521 - 500000) ~ TRUE,
          Chr == 5 & (Pos < 11725524 + 500000 & Pos > 11725524 - 500000) ~ TRUE,
          TRUE ~ FALSE
        )
      ) |>
      dplyr::filter(! .data$is_centromere)

  }
}


help_text <- function() {
  shiny::tags$div(
    shiny::HTML('<hr>
        <div id="help">

        <h3>Help</h3>
        <h4>1. Pre-processing your SNP data</h4>

        <p><span style="color:#34BEDA">Candi</span><span style="color:#FFB416">SNP</span> works on pre-processed,
      filtered, high-confidence SNPs that you provide, there are lots of ways of preparing SNP data. Here are some
      hints on how to go about this:</p>
        <ul>
        <li>The simplest way is with the spreadsheet or VCF file that is provided by your sequence service, if you
      employed one to do SNP calls for you. These files can be edited in Excel or another spreadsheet program
      to include the columns with headers "Chr", "Pos", "Ref", "Alt" and "Allele_Freq", as appropriate. They
      can be in any order, any other columns present are just ignored. Chromosome names are important. See the
      "Selecting a genome" section below for exact details of the names you must use for the. Save the file as a
      "comma-separated values file" for direct use in <span style="color:#34BEDA">Candi</span><span
      style="color:#FFB416">SNP</span></li>
        <li>If you have an alignment file, such as a BAM or SAM file and you need to determine the SNPs from this
      yourself, you can use many published workflows such as those at <a
      href="http://usegalaxy.org">Galaxy</a> or tools in <a href="https://www.bioconductor.org">Bioconductor</a>.
      </li>
        </ul>


        <h4>2. Selecting a genome</h4>

        <p>The genome annotations for TAIR10, TAIR9, Rice genome v7, Tomato genome v2.40, Glycine max genome 1.09v8,
      Grape genome v1 and Maize B73 v5b are made available in <span style="color:#34BEDA">Candi</span><span
      style="color:#FFB416">SNP</span> when you use the <span style="font: monotype;">install_default_genomes() function</span> . You can install other genomes if you wish (see the package documentation for that). In any case it is essential that the chromosome names in column "Chr" match
      that
      in the database. Here are the chromosome names for each default genome build:</p>

        <table border="1">
        <tr>
        <th>Genome</th>
        <th>Chromosome names</th>
        </tr>
        <tr>
        <td>TAIR10</td>
        <td>1, 2, 3, 4, 5, M, C</td>
        </tr>
        <tr>
        <td>TAIR9</td>
        <td>1, 2, 3, 4, 5, M, C</td>
        </tr>
        <tr>
        <td>Rice genome v7</td>
        <td>1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, Un, Sy</td>
        </tr>
        <tr>
        <td>Tomato genome v2.40</td>
        <td>SL2.40ch00, SL2.40ch01, SL2.40ch02, SL2.40ch03, SL2.40ch04, SL2.40ch05, SL2.40ch06, SL2.40ch07,
      SL2.40ch08, SL2.40ch09, SL2.40ch10, SL2.40ch11, SL2.40ch12
      </td>
        </tr>
        <tr>
        <td>Glycine max genome 1.09v8</td>
        <td>Gm01, Gm02, Gm03, Gm04, Gm05, Gm06, Gm07, Gm08, Gm09, Gm10, Gm11, Gm12, Gm13, Gm14, Gm15, Gm16,
      Gm17, Gm18, Gm19, Gm20 and 21 to 2287
      </td>
        </tr>
        <tr>
        <td>Grape genome v1</td>
        <td>1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 1_random, 10_random, 11_random,
      12_random, 13_random, 16_random, 17_random, 18_random, 3_random, 4_random, 5_random, 7_random,
      9_random
      </td>
        </tr>
        <tr>
        <td>Maize B73 v5b</td>
        <td>1, 4, 2, 3, 5, 7, 8, 6, 9, 10, UNKNOWN, Mt, Pt</td>
        </tr>
        </table>

        <h4>3. Upload</h4>

        <p>Once your file is prepared select the genome you wish to use in the sidebar of the "Main" panel, use the "Upload a file" button to select your SNP file. Upload will begin and your data will be analysed. When it is done, an interactive plot will be
      displayed on the screen.</p>

        <h4>4. Control Panel</h4>

        <p>The control panel on the left allows you to customise and filter the plot</p>
        <ul>
            <li>"Allele Frequency Range" - set the upper and lower limit of the Allele Frequency in which to show
                SNPs. Will dynamically change the scale of the plot
            </li>
            <li>"SNP types" - These four coloured boxes allow you to set the colours of the spots representing SNPs of different classes</li>

            <li>"Centromere  SNPs" - for the genome assemblies with known centromere positions, setting this to "Hidden"
                will hide SNPs within 0.5 Mb +/- the centromeres
            </li>
            <li>"Download Annotationss" - will send <span style="color:#34BEDA">Candi</span><span
                    style="color:#FFB416">SNP</span> annotated table to your computer
            </li>
            <li>Other plot parameters can be changed using the plot tool button that appears above the plot area. You can save a copy of the plot from this panel.</li>
        </ul>
        <h4>5. The Plot</h4>

        <p>You will see a track for every contig or chromosome of your chosen assembly. </p>

        <ul>
            <li>The spots represent SNPs, these are plotted on the X-axis according to position on the chromosome and on
                the Y-axis according to allele frequency
            </li>
            <li>Hovering over spots with the mouse pointer brings up a blue box containing summary information for that
                SNP
            </li>
            <li>Spots are coloured according to the colour scheme set in the control panel</li>
        </ul>

        <h4>6. Plotting Tips</h4>
        <ul>
            <li>I just see a splodge of spots...</li>
            <p>If you have too many spots to see the distinctions, first set the upper and lower limit of the Allele
                Frequency slider, this should reduce some. Also note that you can set the transparency of the spots in
                the plot - use the colour selector and press the "more" button, the slider that appears lets you set the
                transparency of the current colour. Lastly, try reducing the number of SNPs in your file - there is a
                fair chance that you have lots of false positive SNPs if you have too many to plot clearly.
            </p>
        </ul>'
    ))
}


intro_text <- function() {
  shiny::tags$div(
    shiny::HTML(
      '<h1><span><span style="color:#34BEDA">Candi</span><span style="color:#FFB416">SNP</span> <br />identifying candidate SNPs in genomes</h1>
          <a href="http://www.eric-carle.com/"><img id="logo_pic" src="https://raw.githubusercontent.com/danmaclean/candisnp/awesome/img/logo.png" title="&lsaquo;Miam miam!&rsaquo; dit la chenille qui fait des trous" style="height: 200px; width: auto; box-shadow: 10px 10px 5px #888888; border-radius: 20px; float: right;"></img></a>

          <div id="intro" class="div-text"><p><span style="color:#34BEDA">Candi</span><span style="color:#FFB416">SNP</span> classifies, annotates and visualises SNPs on genomes. Provide it with a list of SNP positions and allele frequencies and <span style="color:#34BEDA">Candi</span><span style="color:#FFB416">SNP</span> will return the type of each SNP and an interactive  visualisation that you can explore to identify potential causative mutations. </p>
          <p id="cite_us">
          If you use <span style="color:#34BEDA">Candi</span><span style="color:#FFB416">SNP</span> please cite: Etherington, Monaghan <em>et al</em> "Mapping mutations in plant genomes with the user-friendly web application CandiSNP." <a href="http://www.plantmethods.com/content/10/1/306/abstract">Plant Methods 2014, 10:306. doi:10.1186/s13007-014-0041-7</a>.
        </p>

          </div>

          <hr width=1000 align=left>
          <div id="version_text">version 0.4.0 - "le ecureuil emotif"</div>
        '


    )

  )
}

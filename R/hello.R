
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
#'  install_default_genomes()
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

#' download annoations from SNPEff
#'
#' @param genomes character vector of annotation keys (genome names) as described in snpeff.config
#' @returns no value
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
  args <- c("-jar", "-Xmx4g", snp_eff, "-c", conf, "-i", "txt", "-o", "txt", "-noLog", "-noStats", "-canon", "-snp", "-no-downstream", "-no-upstream", "-no-utr", genome, tmpfile)
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

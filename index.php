<html>
<!DOCTYPE html>
<head>
    <meta charset="utf-8">
    <script src="//code.jquery.com/jquery-1.10.2.js"></script>
    <script src="//code.jquery.com/ui/1.11.2/jquery-ui.js"></script>


    <link rel="stylesheet" href="css/bootstrap/css/bootstrap.min.css"/>
    <script type="text/javascript" src="css/bootstrap/js/bootstrap.min.js"></script>

    <link rel='stylesheet' href='css/Flat-UI/dist/css/flat-ui.css'/>

    <link rel="stylesheet" href="//code.jquery.com/ui/1.11.2/themes/smoothness/jquery-ui.css">

    <script src="js/d3/d3.js"></script>
    <script type="text/javascript" src="js/d3-tip.js"></script>

    <link rel='stylesheet' href='css/spectrum.css'/>

    <script src="js/dropzone.js"></script>
    <link rel='stylesheet' href='css/dropzone.css'/>

    <script src="js/genomics_lib.js"></script>
    <script src="js/scatter.js"></script>
    <script src='js/spectrum.js'></script>
    <script src='js/spin.min.js'></script>
    <script src='js/Blob.js'></script>
    <script src='js/FileSaver.min.js'></script>

    <link rel='stylesheet' href='css/custom.css'/>


    <title>CandiSNP</title>

</head>
<body>
<div id="wrapper">
    <div id="header">
        <h1><span id="c">C</span>andi<span style="color:#FFB416">SNP:</span> <br/>identifying candidate SNPs in genomes
        </h1>
        <a href="http://www.eric-carle.com/home.html"><img id="logo_pic" src="img/logo.png"
                                                           title="&lsaquo;Miam miam!&rsaquo; dit la chenille qui fait des trous"></img></a>

        <div id="intro" class="div-text"><p><span style="color:#34BEDA">Candi</span><span
                    style="color:#FFB416">SNP</span> classifies, annotates and visualises SNPs on genomes. Provide it
                with a list of SNP positions and allele frequencies and <span style="color:#34BEDA">Candi</span><span
                    style="color:#FFB416">SNP</span> will return the type of each SNP and an interactive visualisation
                that you can explore to identify potential causative mutations. </p>

            <p id="cite_us">
                If you use <span style="color:#34BEDA">Candi</span><span style="color:#FFB416">SNP</span> please cite:
                Etherington, Monaghan <em>et al</em> "Mapping mutations in plant genomes with the user-friendly web
                application CandiSNP." <a href="http://www.plantmethods.com/content/10/1/306/abstract">Plant Methods
                    2014, 10:306. doi:10.1186/s13007-014-0041-7</a>.
            </p>

        </div>

        <hr width=1000 align=left>

        <div id="version_text">version 0.3.0 - "la chenille charmant"</div>

    </div>


    <?php if (isset($_GET['session']) && isset($_GET['species'])) : ?>

        <script>
            var sessionID = "<?php echo $_GET['session']; ?>";
            var species = "<?php echo $_GET['species']; ?>";
            $.getJSON('http://candisnp.tsl.ac.uk:8080/' + sessionID, function (json) {
                var outputDiv = $("#output");
                outputDiv.css("display", "block");
//                pageData = json;
                pageData = JSON.parse(JSON.parse(json));
                draw(pageData, species);
                $("html, body").delay(100).animate({scrollTop: outputDiv.offset().top}, 2000);
            });
        </script>

    <?php else : ?>
        <div id="input">
            <div class="blue_box">
                <h2 class="input_h2">1. Select the organism &amp; genome version</h2>

                <p>CandiSNP will classify your SNPs from the genome annotations you choose</p>
                <!-- The species selection -->
                <div class="select2-container select select-primary">
                    <select class="select2-choice" id="species_select"
                            title="Select the genome on which SNPs were called">
                        <option value="athalianaTair10">Arabidopsis TAIR 10</option>
                        <option value="athalianaTair9">Arabidopsis TAIR 9</option>
                        <option value="rice7">Rice genome v7</option>
                        <option value="SL2.40">Tomato genome v2.40</option>
                        <option value="gmax1.09v8">Glycine max genome 1.09v8</option>
                        <option value="grape">Grape genome v1</option>
                        <option value="maizeZmB73">Maize B73 v5b</option>
                    </select>

                </div>
            </div>

            <div class="blue_box">
                <h2 class="input_h2">2. Click to select or drop a data file to begin</h2>
                <a href="data/sample_data.csv">Download some sample data to examine and use</a>

                <div class="dropzone" id="my-second-awesome-dropzone">
                </div>

            </div>
        </div>
        <!--input-->
    <?php endif; ?>


    <div id="output">
        <div class="grey_box">
            <h2 class="input_h2">Control Panel</h2>

            <div id="spot_colours">
                <ul>Choose spot colours</ul>
                <!-- The colour pickers -->
                <li><input type='text' id="non_synonymous_in_coding_CT_GA_picker"/> Non-synonymous coding CT/GA</li>
                <li><input type='text' id="non_synonymous_coding_region_picker"/> Non-synonymous coding</li>
                <li><input type='text' id="annotated_region_picker"/> Annotated region</li>
                <li><input type='text' id="non_annotated_region_picker"/> Non-annotated region</li>
                </ul>
            </div>
            <div id="slider_readout">
                <p>Set allele frequency range</p>

                <div id="slider-range"></div>
                <input type="text" id="amount" readonly>

            </div>

            <div id="hide_centromeres_container"></div>
            <div id="download_buttons">
                <div class="btn btn-block btn-lg btn-primary" id="save_as_svg">Download SVG</div>
                <div class="btn btn-block btn-lg btn-primary" id="save_as_png">Download PNG</div>
                <div class="btn btn-block btn-lg btn-primary" id="save_as_table">Download Annotated SNPs</div>
            </div>
        </div>
    </div>
    <!--output-->

    <div id="results">
    </div>

    <hr>
    <div id="help">

        <h3>Help</h3>
        <h4>1. Pre-processing your SNP data</h4>

        <p><span style="color:#34BEDA">Candi</span><span style="color:#FFB416">SNP</span> works on pre-processed,
            filtered, high-confidence SNPs that you provide, there are lots of ways of preparing SNP data. Here are some
            hints on how to go about this:
        <ul>
            <li>The simplest way is with the spreadsheet or VCF file that is provided by your sequence service, if you
                employed one to do SNP calls for you. These files can be edited in Excel or another spreadsheet program
                to include the columns with headers 'Chr', 'Pos', 'Ref', 'Alt' and 'Allele_Freq', as appropriate. They
                can be in any order, any other columns present are just ignored. Chromosome names are important. See the
                'Selecting a genome' section below for exact details of the names you must use. Save the file as a
                'comma-separated values file' for direct use in <span style="color:#34BEDA">Candi</span><span
                    style="color:#FFB416">SNP</span></li>
            <li>If you have an alignment file, such as a BAM or SAM file and you need to determine the SNPs from this
                yourself, you can use many published workflows such as those at <a
                    href="http://usegalaxy.org">Galaxy</a>. We used the following command-line script in our work <a
                    href="https://github.com/danmaclean/candisnp/blob/master/pileup_to_snps.rb">linked here</a>.
            </li>
        </ul>
        </p>

        <h4>2. Selecting a genome</h4>

        <p>The genome annotations for TAIR10, TAIR9, Rice genome v7, Tomato genome v2.40, Glycine max genome 1.09v8,
            Grape genome v1 and Maize B73 v5b are currently available in <span style="color:#34BEDA">Candi</span><span
                style="color:#FFB416">SNP</span>. It is essential that the chromosome names in column 'Chr' match that
            in the database. Here are the chromosome names for each genome build:</p>

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

        <p>Once your file is prepared click or drop your file into the 'Click to select or drop a data file to begin'
            area. Upload will begin and your data will be analysed. When it is done, an interactive plot will be
            displayed on the screen.</p>

        <h4>4. Control Panel</h4>

        <p>When the plot appears, a control panel to determine it's data and appearance is also displayed, it has four
            sections</p>
        <ul>
            <li>'Choose Spot Colours' - set the colours of the spots representing SNPs of different classes</li>
            <li>'Set Allele Frequency Range' - set the upper and lower limit of the Allele Frequency in which to show
                SNPs
            </li>
            <li>'Hide Centromere Region SNPs' - for the genome assemblies with known centromere positions, clicking this
                button will hide SNPs within 0.5 Mb +/- the centromeres
            </li>
            <li>'Download SVG/PNG/Annotated SNPs' - will send <span style="color:#34BEDA">Candi</span><span
                    style="color:#FFB416">SNP</span> images or an annotated table to your computer
            </li>
        </ul>
        <h4>5. The Plot</h4>

        <p>You will see a plot like this for every contig or chromosome of your chosen assembly. </p>
        <img src="img/example_plot.png"></img>

        <p>
        <ul>
            <li>The spots represent SNPs, these are plotted on the X-axis according to position on the chromosome and on
                the Y-axis according to allele frequency
            </li>
            <li>Hovering over spots with the mouse pointer brings up a blue box containing summary information for that
                SNP
            </li>
            <li>Spots are coloured according to the colour scheme set in the control panel</li>
        </ul>
        </p>
        <h4>6. Plotting Tips</h4>
        <ul>
            <li>I just see a splodge of spots...</li>
            <p>If you have too many spots to see the distinctions, first set the upper and lower limit of the Allele
                Frequency slider, this should reduce some. Also note that you can set the transparency of the spots in
                the plot - use the colour selector and press the 'more' button, the slider that appears lets you set the
                transparency of the current colour. Lastly, try reducing the number of SNPs in your file - there is a
                fair chance that you have lots of false positive SNPs if you have too many to plot clearly.
            </p>
    </div>


    <script>
        $("#non_synonymous_in_coding_CT_GA_picker").spectrum(colour_settings("#ff0000", "NON_SYNONYMOUS_CODING_CT_GA"));

        $("#non_synonymous_coding_region_picker").spectrum(colour_settings("#ff0000", "NON_SYNONYMOUS_CODING"));

        $("#annotated_region_picker").spectrum(colour_settings("#aaaaaa", "ANNOTATED_REGION"));

        $("#non_annotated_region_picker").spectrum(colour_settings("#aaaaaa", "NON_ANNOTATED_REGION"));
    </script>


    <!-- Hidden <FORM> to submit the SVG data to the server, which will convert it to SVG/PDF/PNG downloadable file.
         The form is populated and submitted by the save_png() JavaScript below. -->
    <form id="svgform" method="post" action="cgi/download.cgi">
        <input type="hidden" id="output_format" name="output_format" value="">
        <input type="hidden" id="data" name="data" value="">
    </form>


    <script>
        $("#save_as_svg").click(function () {
            save_svg();
        });
    </script>

    <script>
        $('#save_as_png').click(function () {
            save_png();
        });
    </script>

    <script>
        $('#save_as_table').click(function () {
            save_table();
        });
    </script>


    <script>
        var ohcanada = new Audio('easter_eggs/oh_canada.mp3');
        $("#c").mouseover(function () {
            ohcanada.play();
        }).mouseout(function () {
            ohcanada.pause();
        })
    </script>


</div>
</body>
</html>

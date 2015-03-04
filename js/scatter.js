//get the chromosome list
function get_chromosomes(data) {
    var hash = {};
    for (var i in data) {
        hash[data[i].chromosome] = 1;
    }
    return Object.keys(hash).sort()
}

//get the data for a given chromosome
function select_for(chr, data) {
    var result = [];
    for (var i in data) {
        if (data[i].chromosome == chr) {
            result.push(data[i]);
        }
    }
    return result;
}

//reformat big numbers to SI units
function bases_to_unit_bases(bases) {
    var units = ['Kbp', 'Mbp', 'Gbp'];
    var u = -1;
    do {
        bases /= 1000;
        ++u;
    } while (bases >= 1000);
    return bases.toFixed(1) + ' ' + units[u];
}

//get the svg as a string

function get_svg_string(plots) {
    // Get the d3js SVG element
    var xml = '<svg xmlns="http://www.w3.org/2000/svg" id="candisnp_output" width="1000" height="' + (plots.length * 300) + '" >'; //+
    for (var i = 0; i < plots.length; i++) {
        xml = xml + (new XMLSerializer).serializeToString(plots[i]);
    }
    xml = xml + '</svg>';
    return xml;
}


//send the svg to a file
function save_svg() {
    // Get the d3js SVG element
    var body = document.body;
    var plots = body.getElementsByClassName("plot_element");
    var xml = get_svg_string(plots);
    var blob = new Blob([xml], {type: "image/svg+xml"});
    saveAs(blob, "CandiSNP_Output.svg");
    return true;
}

//send the svg to the external CGI script and convert to PNG
function save_png() {
    var body = document.body;
    var plots = body.getElementsByClassName("plot_element");
    var form = document.getElementById("svgform");
    form['output_format'].value = 'png';
    form['data'].value = get_svg_string(plots);
    form.submit();
}


//save the json data as a csv file
function save_table() {
    console.log(pageData[0]);
    var string = ' "Chr", "Pos", "Alt", "Ref", "Allele_Freq", "Is_CTGA", "Is_Synonymous", "In_CDS", "Change", "Effect"\n';
    //console.log(string);
    for (i in pageData) {
        var snp = pageData[i];
        string = string + snp_to_string(snp);
    }
    var blob = new Blob([string], {type: "text/csv"});
    saveAs(blob, "CandiSNP_Output.csv");
    return true;
}

function snp_to_string(snp) {
    var string = "";
    string = string + [snp.chromosome, snp.position, snp.alternate_base, snp.reference_base, snp.allele_freq, snp.is_ctga, snp.is_synonymous, snp.in_cds, snp.change, snp.effect].join(",") + "\n";
    return string;
}


function ratio_current_chromosome_to_longest(species, chr) {

    console.log('here', species, chr);


    var lengths = genome_lengths(species);
    var current = lengths[chr];
    var keys = Object.keys(lengths);
    var values = keys.map(function (v) {
        return lengths[v];
    });
    var longest = values[0];
    return current / longest;
}

//splits data objext into single sets per chromosome
//adds species specific buttons
function draw(data) {

    console.log('going to try and draw', data);

    var chromosomes = get_chromosomes(data);
    var species = $('#species_select').val();

    if (has_centromere_positions(species)) {
        add_centromere_select(species);
    }

    var svg = d3.select("#results")
        .append('svg')
        .attr('width', 1000)
        .attr('height', (300 * chromosomes.length));

    var count = 0;
    for (i in chromosomes) {
        var chr = chromosomes[i];
        var centromere_positions = get_centromere_positions(species, chr);
        var single_chr_data = select_for(chr, data);
        draw_single(species, chr, centromere_positions, single_chr_data, svg, count);
        count = count + 1;
    }

}

// actually draws the panels in the plots
function draw_single(species, chr, centromere_positions, data, svg, count) {

    var plot_id = "chr_" + chr;
    var ratio = ratio_current_chromosome_to_longest(species, chr);
    var margin = 50, width = (ratio * 1000), height = 300;
    //var x_extent = d3.extent(data, function(d){return d.position});
    var x_extent = [1, genome_lengths(species)[chr]];
    var y_extent = d3.extent(data, function (d) {
        return d.allele_freq
    });


    var x_scale = d3.scale.linear()
        .range([margin, width - margin])
        .domain(x_extent);

    var y_scale = d3.scale.linear()
        .range([height - margin, margin])
        .domain(y_extent);


    var x_axis = d3.svg.axis().scale(x_scale).tickFormat(function (d) {
        return bases_to_unit_bases(d)
    }).ticks(7);
    var y_axis = d3.svg.axis().scale(y_scale).orient("left").ticks(5);

    var tip = d3.tip()
        .attr('class', 'd3-tip')
        .offset([-10, 0])
        .html(function (d) {
            return "Position: " + d.position + " <br />  \
	  Allele Frequency: " + d.allele_freq + "<br /> \
	  Locus: " + d.gene + " <br /> \
	  Change: " + d.change + "<br /> \
	  Reference base: " + d.reference_base + "<br /> \
	  Alternate base: " + d.alternate_base + "<br />";
        })


    var g = svg.append("g");

    g.attr('class', 'plot_element');
    g.attr("transform", "translate(0," + ( count * 300) + ")");

    g.attr("id", plot_id)
        .attr('width', width)
        .attr('height', height)
        .selectAll("circle")
        .data(data)
        .enter()
        .append("circle");

    g.selectAll("circle")
        .attr("cx", function (d) {
            return x_scale(d.position)
        })
        .attr("cy", function (d) {
            return y_scale(d.allele_freq)
        })
        .on('mouseover', tip.show)
        .on('mouseout', tip.hide);

    g.selectAll("circle")
        .attr("r", 0);

    //add class to circles
    g.selectAll("circle")
        .attr("class", function (d) {
            return get_snp_type(d)
        })
        .style("fill", function (d) {
            var snp_type = get_snp_type(d);
            return default_colour(snp_type);
        });

    //add within centromere information
    if (centromere_positions) {
        g.selectAll("circle")
            .classed('in_centromere', function (d) {
                return in_centromere(d, centromere_positions);
            })
    }


    g.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + (height - margin) + " )")
        .call(x_axis);

    g.append("text")
        .attr("class", "x label")
        .attr("x", 50)
        .attr("y", height - 6)
        .text("Chromosome/contig: " + chr);

    g.append("g")
        .attr("class", "y axis")
        .attr("transform", "translate(" + margin + ",0 )")
        .call(y_axis);

    g.append("text")
        .attr("class", "y label")
        .attr("text-anchor", "end")
        .attr("y", 35)
        .attr("x", "110")
        //.attr("transform", "rotate(-90)")
        .text("Allele Frequency");

    d3.selectAll('.axis path')
        .style("fill", "none")
        .style("stroke", "#000")
        .style("shape-rendering", "crispEdges");

    d3.selectAll('.axis line')
        .style("fill", "none")
        .style("stroke", "#000")
        .style("shape-rendering", "crispEdges");


    g.selectAll("circle")
        .transition()
        .delay(function (d, i) {
            return i / data.length * 1000;
        })
        .attr("r", 5);

    g.call(tip);

}

//add the centromere select button
function add_centromere_select(species) {
    if ($('#hide_centromeres').length == 0) {
        $('#hide_centromeres_container').append('<form>Hide Centromere Region SNPs <input type="checkbox" id="hide_centromeres" name="hide_centromeres" value="true"></form>');
        add_centromere_listener(species);
    }
}

//function that creates button to obscure centromeres when the genome selected allows it
function add_centromere_listener(species) {


    // now add the listener for the select box that does the hiding
    $('#hide_centromeres').change(function () {
        if (this.checked) {
            d3.selectAll(".in_centromere")
                .transition()
                .duration(250)
                .style("opacity", 0);

        }
        else {
            d3.selectAll("circle")
                //.filter(function(d,i){
                //  return (d.position >= centromere_range[0] && d.position <= centromere_range[1]);
                //})
                .transition()
                .duration(250)
                .style("opacity", 1);
        }
    });
}

//function that controls the slider behaviour for the allele freq filter
$(function () {
    $("#slider-range").slider({
        range: true,
        min: 1.0,
        max: 100.0,
        orientation: "horizontal",
        values: [75, 100],
        slide: function (event, ui) {
            var max = (ui.values[1] / 100);
            var min = (ui.values[0] / 100)
            $("#amount").val(min + " - " + max);
            d3.selectAll("circle")
                .filter(
                function (d, i) {
                    return !(d.allele_freq >= min && d.allele_freq <= max )

                })
                .transition()
                .duration(250)
                .style("opacity", 0);

            d3.selectAll("circle")
                .filter(
                function (d, i) {
                    return (d.allele_freq >= min && d.allele_freq <= max )

                })
                .transition()
                .duration(250)
                .style("opacity", 1);


        }
    });
    $("#amount").val("0.75 - 1");
});

//decide the snp type
function get_snp_type(d) {

    if (d.is_synonymous == "FALSE" && d.is_ctga == "TRUE" && d.in_cds == "TRUE") {
        return "NON_SYNONYMOUS_CODING_CT_GA";
    }
    else if (d.is_synonymous == 'FALSE' && d.is_ctga == "FALSE" && d.in_cds == "TRUE") {
        return "NON_SYNONYMOUS_CODING";
    }
    else if (d.is_synonymous == 'TRUE' && d.is_ctga == "FALSE" && d.in_cds == "TRUE") {
        return "ANNOTATED_REGION";
    }
    else if (d.in_cds == "FALSE") {
        return "NON_ANNOTATED_REGION";

    }
    else {
        return "NON_ANNOTATED_REGION";
    }
}

//format the popout string
function format_popup(d) {
    // Example object in d  "chromosome":"1","position":825457,"reference_base":"G","alternate_base":"A","allele_freq":0.764705882,"in_cds":"TRUE","is_synonymous":"TRUE","is_ctga":"TRUE","change":"R/R","gene":"AT1G03360","summary":"SYNONYMOUS_CODING"

    return d.gene + " " + d.position + "<tspan> " + d.reference_base + " -> " + d.alternate_base + " : " + d.change + "</tspan>";
}

function default_colour(snp_type) {

    if (snp_type == "NON_SYNONYMOUS_CODING_CT_GA") {
        return "#ff0000";
    }
    else if (snp_type == "NON_SYNONYMOUS_CODING") {
        return "#ff0000";
    }
    else if (snp_type == "ANNOTATED_REGION") {
        return "#aaaaaa";
    }
    else {
        return "#aaaaaa";
    }

}

function set_spot_colour(colour, snp_type) {
    d3.selectAll("." + snp_type)
        .transition()
        .duration(200)
        .style("fill", colour)
}


function colour_settings(default_colour, snp_type) {
    return {
        color: default_colour,
        showAlpha: true,
        showPalette: true,
        showPaletteOnly: true,
        togglePaletteOnly: true,
        togglePaletteMoreText: 'more',
        togglePaletteLessText: 'less',
        clickoutFiresChange: true,
        allowEmpty: true,
        change: function (colour) {
            if (colour == null) {
                set_spot_colour('rgba(0,0,0,0)', snp_type);
            }
            else {
                set_spot_colour(colour.toRgbString(), snp_type);
            }
        },
        palette: [
            ["#FF0000", "#AAAAAA", "#000000", "#FFFFFF"],
            ["#225EA8", "#41B6C4", "#A1DAB4", "#FFFFCC"],
            ["#D7191C", "#FDAE61", "#2C7BB6", "#ABD9E9"],
            ["#E31A1C", "#FD8D3C", "#FECC5C", "#FFFFB2"],
            ["#CB181D", "#FB6A4A", "#FCAE91", "#FEE5D9"],


        ]
    };
}

function all_in_array(needed, to_check) {

    var total = 0;
    for (var i = 0; i < needed.length; i++) {
        if (one_in_array(needed[i], to_check)) {
            total = total + 1;
        }
    }

    return total == needed.length;

}


function one_in_array(item, list) {

    for (var i = 0; i < list.length; i++) {
        if (list[i].toUpperCase() === item.toUpperCase()) {
            return true;
        }
    }
    return false;
}

function chomp(raw_text) {
    return raw_text.replace(/(\n|\r)+$/, '');
}

function headers_checker(file_contents, done) {
    // do file munging
    var lines = file_contents.split("\n");
    var headers = chomp(lines[0]).split(",");
    if (all_in_array(["Chr", "Pos", "Alt", "Ref", "Allele_Freq"], headers)) {
        return done();
    }
    else {
        return done("File doesn't have proper headings. Please check the help and try again");
    }
}


function file_ok(file, done) {
    if (window.File && window.FileReader && window.FileList && window.Blob) {
        //file is a proper html5 File object
        var reader = new FileReader();
        reader.onload = function (e) {
            var file_contents = reader.result;
            headers_checker(file_contents, done);
        };
        reader.readAsText(file);
    } else {
        done('The File APIs are not fully supported by your browser. Please try a different browser');
    }


}

function add_species_to_form() {
    this.on("sending", function (file, xhr, formData) {
        formData.append("species", $('#species_select').val())
    })
}

var pageData = null;

Dropzone.options.mySecondAwesomeDropzone = {
    paramName: "file",  // The name that will be used to transfer the file
    maxFilesize: 200, // MB
    url: "cgi/snpeff.cgi",
    addRemoveLinks: true,
    acceptedFiles: "text/csv",
    uploadprogress: true,
    maxFiles: 1,
    // params: {"organism": $('#species_select').val(), "file": "filename" },
    init: add_species_to_form,
    accept: function (file, done) {
        file_ok(file, done);
        var opts = spinner_opts();
        var target = document.getElementById('my-second-awesome-dropzone');
        var spinner = new Spinner(opts).spin(target);

    },
    success: function (file, response_json) {
        var out = $("#output");
        //console.log(response_json);
        $('.spinner').remove();
        out.css("display", "block");
        //d3.json(response_json, draw);
        var data = JSON.parse(response_json);
        pageData = data;
        draw(data);
        $("html, body").delay(100).animate({scrollTop: out.offset().top}, 2000);
    },
    canceled: function () {
        $('.spinner').remove();
    },
    reset: function () {
        $("svg").remove();
        $("#output").css("display", "none");
    }
};

function spinner_opts() {
    return {
        lines: 13, // The number of lines to draw
        length: 40, // The length of each line
        width: 23, // The line thickness
        radius: 56, // The radius of the inner circle
        corners: 1, // Corner roundness (0..1)
        rotate: 8, // The rotation offset
        direction: 1, // 1: clockwise, -1: counterclockwise
        color: '#000', // #rgb or #rrggbb or array of colors
        speed: 0.5, // Rounds per second
        trail: 22, // Afterglow percentage
        shadow: true, // Whether to render a shadow
        hwaccel: false, // Whether to use hardware acceleration
        className: 'spinner', // The CSS class to assign to the spinner
        zIndex: 2e9, // The z-index (defaults to 2000000000)
        top: '50%', // Top position relative to parent
        left: '50%' // Left position relative to parent
    };
}

function in_centromere(d, centromere_range) {
    return !!(d.position >= centromere_range[0] && d.position <= centromere_range[1]);

}
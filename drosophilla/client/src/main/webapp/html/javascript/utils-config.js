// var data_url = "https://ondex.rothamsted.ac.uk/drosophilla_data/";
var data_url = "http://babvs43.rothamsted.ac.uk:8080/drosophilla_data/"; // change from 8080 to port 80
// var data_url = "http://qtlnetminer-test.rothamsted.ac.uk/drosophilla_data/"; // now, using HTTP instead of HTTPS
var applet_url = "https://ondex.rothamsted.ac.uk/OndexWebBeta/";
var species = 'Drosophilla melanogaster genome';
var reference_genome = true;
var multiorganisms = false;
// adjusts scaling in the GViewer in case the first chromosome is not the longest.
var longest_chr= 35000000/*32079331*/;

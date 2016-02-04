var GENEMAP = GENEMAP || {};

GENEMAP.XmlReader = function() {
  if (!(this instanceof arguments.callee)) {
    return new arguments.callee();
  }

  // read a chromosome band object from an XML element
  var _readBand = function(elt) {
    return {
      index: elt.getAttribute("index"),
      start: elt.getElementsByTagName("start")[0].childNodes[0].nodeValue,
      end: elt.getElementsByTagName("end")[0].childNodes[0].nodeValue,
      color: elt.getElementsByTagName("color")[0].childNodes[0].nodeValue,
    };
  }

  // read a chromosome JS object from an XML element
  var _readChromosome = function(elt) {
    var chromosome = {
      index:elt.getAttribute("index"),
      length:elt.getAttribute("length"),
      number:elt.getAttribute("number"),
      bands:[]
    }

    var bandElements = elt.getElementsByTagName("band");

    for (var j = 0; j < bandElements.length; j++) {
      var band = _readBand(bandElements[j])
      chromosome.bands.push(band);
    }

    return chromosome;
  }

  // reads the genome data from a basemap XML document
  var _readBasemapXML = function(xml) {
    var genome = {};
    genome.chromosomes = [];

    var chromosomeElements = xml.getElementsByTagName("chromosome");
    for (var i = 0; i < chromosomeElements.length; i++){
      var chromosome = _readChromosome(chromosomeElements[i]);
      genome.chromosomes.push(chromosome);
    }

    return genome;
  }

  return {

    readBasemapXML: function(path, callbackFn) {
      d3.xml(path, function(error, xml) {
        var data = _readBasemapXML(xml)
        callbackFn(data);
      });
    }
  }
}
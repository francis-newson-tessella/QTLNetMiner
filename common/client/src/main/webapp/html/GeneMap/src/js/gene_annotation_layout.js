var GENEMAP = GENEMAP || {};

//Produce layout for gene annotations
//GENEMAP.GeneClusterer is used to cluster genes if necessary
//Labella is used to generate layout of nodes

GENEMAP.GeneAnnotationLayout = function (userConfig) {

  var defaultConfig = {
    longestChromosome: 100,
    layout: {
      width: 10, //not used
      height: 100,
      x: 0, //not used
      y: 0, //not used
    },
    autoLabels: true,
    manualLabels: true,
    annotationMarkerSize: 5,
    annotationLabelSize: 5,
    doCluster: true,
    nClusters: 6,
    maxAnnotationLayers: 3,
    scale: 1,
  };

  var config = _.merge({}, defaultConfig, userConfig);

  var buildYScale = function () {
    return d3.scale.linear()
      .range([0, config.layout.height])
      .domain([0, config.longestChromosome]);
  };

  var shouldRecluster = function (nodes) {
    if (!config.doCluster) {
      return false;
    }
    var layers = nodes.map(function (d) {
      return d.getLayerIndex();
    });
    var maxLayer = Math.max.apply(null, layers);
    return ( maxLayer > config.maxAnnotationLayers );
  }

  var generateLayoutParameters = function (chromosome) {

    //How much space do we have?
    var availableWidth = config.layout.width;
    var availableHeight = config.layout.height;
    var labelLength = 20;
    var minDisplayedFontSize = 12;
    var maxDisplayedFontSize = 14;
    var fontCoordRatio = 3.5;

    //How many rows of labels do we show?

    //Can we show 1 row?
    var par1 = {};
    par1.spaceForLabel = availableWidth - ( 0.1 + 0.1 / config.scale)
    par1.setFontSize = par1.spaceForLabel / labelLength * fontCoordRatio;
    par1.possible = par1.setFontSize * config.scale > minDisplayedFontSize;
    par1.nodeSpacing = par1.setFontSize;
    par1.density = 1.0;

    //Can we show 2 rows?
    var par2 = {};
    par2.spaceForLabel = availableWidth / 3;
    par2.setFontSize = par2.spaceForLabel / labelLength * fontCoordRatio;
    par2.possible = par2.setFontSize * config.scale > minDisplayedFontSize;

    par2.setFontSize = Math.min(par2.setFontSize, maxDisplayedFontSize / config.scale);
    par2.nodeSpacing = par2.setFontSize;
    par2.spaceForLabel = 1.3 * labelLength * par2.setFontSize / fontCoordRatio;
    par2.density = 0.9;

    //var nRows = par1.possible + par2.possible;
    var nRows = par1.possible + par2.possible;
    var par = (nRows == 2) ? _.clone(par2) : _.clone(par1);
    par.setFontSize = Math.min(par.setFontSize, maxDisplayedFontSize / config.scale);
    par.layerGap = Math.min(5 * par.setFontSize, availableWidth / 3);
    par.nodeSpacing = Math.max(2, par.setFontSize);
    par.lineSpacing = 1;
    par.scale = config.scale;
    par.availableHeight = availableHeight;
    if (nRows == 0) {
      par.nLabels = 0;
    } else if (nRows == 1) {
      par.nLabels = 0.4 * availableHeight / (par.nodeSpacing + par.lineSpacing);
    } else if (nRows == 2) {
      par.nLabels = 0.6 * availableHeight / (par.nodeSpacing + par.lineSpacing);
    }

    if (false && chromosome.number == "2B") {
      log.info("par1", par1);
      log.info("par2", par2);
      log.info("par", par);
    }
    return par;
  };

  var setMerge = function(autoLabels, manualLabels){
    var nodeSet = new Set(manualLabels);
    autoLabels.forEach( function(gene){ nodeSet.add(gene)});
    var nodeSource = Array.from(nodeSet);
    return nodeSource;
  }

  var generateNodes = function( force,y, par, nodeSource ){
    if( nodeSource.length < 1){
      return []
    };

    var nodes = nodeSource.map(function (d) {
      return new labella.Node(y(d.midpoint), par.setFontSize, d);
    } );

    try {
      force.nodes(nodes).compute();
      }

    catch (e) {
      if (e instanceof RangeError) {
        log.error(e.message);
        log.info(par);
        return null;
      }
      else {
        throw e;
      }
    }

    var layers = nodes.map(function (d) {  return d.getLayerIndex(); } );
    var maxLayer = Math.max.apply(null, layers);
    if (maxLayer < 2) {
      return nodes;
    } else {
      return null;
    }
  }

  //Use labella to generate layout nodes for each gene
  //or cluster of genes
  var generateChromosomeLayout = function(chromosome){

    var par = generateLayoutParameters(chromosome);

    var y = buildYScale();

    forceConfig = {
      nodeSpacing: par.nodeSpacing,
      lineSpacing: par.lineSpacing,
      algorithm: 'overlap',
      minPos: 0,
      maxPos: par.availableHeight,
      density: par.density,
    };

    var force = new labella.Force( forceConfig );

    //Decide which labels to display
    //Start with no labels displayed
    allGenes = chromosome.annotations.allGenes;
    allGenes.forEach( function(gene){ gene.displayed = false}) ;

    //Include all genes set to visible
    var manualLabels = config.manualLabels
      ? allGenes.filter( function(gene){ return gene.visible;})
      : [];

    //Automatically show some additional labels
    var autoLabels =  config.autoLabels
    ? allGenes.slice(0, par.nLabels).filter( function(gene){return !gene.hidden;})
    : [];

    var nodeSource = setMerge(autoLabels, manualLabels);

    if ( nodeSource.length < 1 || par.nLabels < 1){
      return [];
    }

    var nodes = generateNodes( force, y, par, nodeSource);

    if ( ! nodes) {
      nAutoLabels = Math.max(par.nLabels - manualLabels.length, 0 );

      autoLabels =  config.autoLabels
        ? allGenes.slice(0, nAutoLabels).filter( function(gene){return !gene.hidden;})
        : [];

      nodeSource = setMerge(autoLabels, manualLabels);
      nodes = generateNodes( force, y, par, nodeSource);
    }

    if ( ! nodes ){
      var geneClusterer = GENEMAP.GeneClusterer().nClusters(par.nLabels)
      try {
      nodeSource = geneClusterer.createClustersFromGenes(manualLabels);
        }
      catch (e) {
        log.info( manualLabels);
        log.info( par.nLabels);
        throw e;
      }
      nodes = generateNodes( force, y, par, nodeSource);
    }

    if ( !nodes){
      log.error( 'Null nodes');
    }

    nodeSource.forEach( function(gene){
      gene.displayed = true
      gene.fontSize = par.setFontSize;
    });

    //Compute paths
    renderConfig  = {
      direction: 'right',
      layerGap: par.layerGap,
      nodeHeight: par.spaceForLabel,
    };

    var renderer = new labella.Renderer(renderConfig);
    renderer.layout(nodes);

    nodes.forEach( function(node){
     node.data.path = renderer.generatePath(node);
    });


    if ( false &&  chromosome.number == "2B") {
      log.info( nodes );
    }


    return nodes;

  }

  //Produce list of clusters (which could be single genes)
  //for a given chromosome
  var generateChromosomeClusters = function(chromosome){
    var geneClusterer = GENEMAP.GeneClusterer();
    //Run clustering algorithm so we can use the clusters later when drawing
    genes = chromosome.annotations.genes;
    geneClusters = geneClusterer.createClustersFromGenes(genes);
    return geneClusters;
  }

  my = {};

  my.layoutChromosome = function(chromosome){
    chromosome.layout.annotationNodes = (
    generateChromosomeLayout(chromosome) || chromosome.layout.annotationNodes);
  }

  my.computeChromosomeClusters = function(chromosome){
    chromosome.layout.annotationClusters = generateChromosomeClusters(chromosome);
    chromosome.layout.annotationDisplayClusters = chromosome.layout.annotationClusters.slice();
  };

  my.expandAllChromosomeClusters = function(chromosome) {
    chromosome.layout.annotationDisplayClusters = chromosome.annotations.genes;
  };

  my.collapseAllChromosomeClusters = function(chromosome) {
    chromosome.layout.annotationDisplayClusters = chromosome.layout.annotationClusters.slice();
  };

  my.expandAChromosomeCluster= function( chromosome, cluster) {
    chromosome.layout.annotationDisplayClusters = chromosome.layout.annotationClusters.slice();

    //add each gene as it's own cluster
    cluster.genesList.forEach( function(gene){
      chromosome.layout.annotationDisplayClusters.push(gene) ;
    } );

    //delete the original cluster
    var clusterIndex = chromosome.layout.annotationDisplayClusters.indexOf(cluster);
    chromosome.layout.annotationDisplayClusters.splice(clusterIndex, 1);
  };

  return my;
}

var GENEMAP = GENEMAP || {};

//Produce layout for gene annotations
//GENEMAP.GeneClusterer is used to cluster genes if necessary
//Labella is used to generate layout of nodes

GENEMAP.GeneBandsLayout = function (userConfig) {

  var defaultConfig = {
    longestChromosome: 100,
    layout: {
      width: 10, //not used
      height: 100,
      x: 0, //not used
      y: 0, //not used
    },
    doCluster : true,
    nClusters: 6,
    scale: 1,
  };

  var config = _.merge({}, defaultConfig, userConfig);
  var y;

  var buildYScale = function () {
    return d3.scale.linear()
      .range([0, config.layout.height])
      .domain([0, config.longestChromosome]);
  };

  var shouldRecluster = function(nodes) {
    return config.doCluster;
  }


  var createNode = function(cluster){

    result = _.pick(cluster, 'start', 'end', 'midpoint', 'color');
    result.data = cluster;

    if (cluster.type == "gene") {
      result.color = cluster.color;
    }
    else if (cluster.type == "geneslist"){
      result.color = '#0000FF';
    }
    else{
      log.error( "unregconized cluster type");
      log.info( cluster);
    }
      return result;
  }

  var generateChromosomeLayout = function(chromosome){
    y = buildYScale();

    //Start by constructing nodes directly from genes
    var nodeSource = chromosome.layout.geneBandDisplayClusters;
    var nodes = nodeSource.map( createNode );
    return nodes;

  }

  var mergeGenes = function( genes){
    var id = genes.reduce(function (sum, current) {
      return sum + current.id.toString();
    }, "");

    var end =  genes.reduce( function(max,current){
      return Math.max(max, current.end);
    }, 0);
    var start = genes.reduce( function(min,current){
      return Math.min(min, current.start);
    }, Infinity);

    var midpoint = (start + end)  /2;

    var resultCluster = {
      type: "geneslist",
      id: id,
      genesList: genes,
      start: start,
      midpoint: midpoint,
      end: end,
    };
    return resultCluster;
  }

  var mergeIdenticalPositions = function (genes){
    var geneClusters = [];

    var iGene = 0;
    while (iGene < genes.length) {
      iDiff = iGene;
      while (iDiff < genes.length &&
      (genes[iGene].midpoint == genes[iDiff].midpoint)) {
        iDiff++;
      }
      nMatching = iDiff - iGene;

      if (nMatching == 1) {
        geneClusters.push(genes[iGene]);
        iGene++;
      }
      else {
        var genesList = genes.slice(iGene, iDiff);
        var mergedCluster = mergeGenes(genesList);
        geneClusters.push(mergedCluster);
        iGene = iDiff;
      }
    }
    return geneClusters;
  }

  var mergeNearbyClusters = function(rawClusters, threshold){
    var clusters = rawClusters.slice();


    while( true){

      if (clusters.length < 2){
        break;
      }

      var minDistance = Infinity;
      var iMin = 0;
      var jMin = 0;

      //Find min clusters
      for ( var iClus = 0 ; iClus < clusters.length ; iClus++ ){
        for( var jClus = 0 ; jClus < iClus ; jClus++ ){
          var iCluster = clusters[iClus];
          var jCluster = clusters[jClus];
          var distance = Math.abs(iCluster.midpoint - jCluster.midpoint);

          if (distance < minDistance){
            iMin = iClus;
            jMin = jClus;
            minDistance = distance;
          }
        }
      }

      if (minDistance > threshold){
        break;
      }

      //Merge clusters
      var iCluster = clusters[iMin];
      var jCluster = clusters[jMin];
      var iGenes = ( ( iCluster.type  == "geneslist" ) ? iCluster.genesList.slice() : [iCluster]);
      var jGenes = ( ( jCluster.type  == "geneslist" ) ? jCluster.genesList.slice() : [jCluster]);
      var newGenes = iGenes.concat(jGenes);
      var newCluster = mergeGenes(newGenes);

      clusters.splice(iMin, 1);
      clusters.splice(jMin, 1);
      clusters.push(newCluster);

    }

    return clusters;
  }

//Produce list of clusters (which could be single genes)
//for a given chromosome
  var generateChromosomeClusters = function(chromosome) {
    log.info( "generateChromosomeClusters");

    var genes = chromosome.annotations.allGenes.slice();
    genes.sort(function (lhs, rhs) {
        return lhs.midpoint - rhs.midpoint
      });

    geneClusters = mergeIdenticalPositions(genes);

    minClusterSeparation = config.longestChromosome / config.layout.height ;
    geneClusters = mergeNearbyClusters( geneClusters, minClusterSeparation);

    geneClusters.sort(function (lhs, rhs) {
      return lhs.midpoint < rhs.midpoint
    });

    return geneClusters;
  }


  my = {};

  my.layoutChromosome = function(chromosome){
    chromosome.layout.geneBandNodes = generateChromosomeLayout(chromosome)
  }

  my.computeChromosomeClusters = function(chromosome){
    ly = chromosome.layout;
    ly.geneBandClusters = generateChromosomeClusters(chromosome);
    ly.geneBandDisplayClusters = ly.geneBandClusters.slice();
  };

  my.expandAllChromosomeClusters = function(chromosome) {
    ly = chromosome.layout;
    ly.geneBandDisplayClusters = chromosome.annotations.allGenes;
  };

  my.collapseAllChromosomeClusters = function(chromosome) {
    ly = chromosome.layout;
    ly.geneBandDisplayClusters = ly.geneBandClusters.slice();
  };

  my.expandAChromosomeCluster= function( chromosome, cluster) {
    ly = chromosome.layout;
    ly.geneBandDisplayClusters = ly.geneBandClusters.slice();

    //add each gene as it's own cluster
    cluster.genesList.forEach( function(gene){
      ly.geneBandDisplayClusters.push(gene) ;
    } );

    //delete the original cluster
    var clusterIndex = ly.geneBandDisplayClusters.indexOf(cluster);
    ly.geneBandDisplayClusters.splice(clusterIndex, 1);
  };

  return my;
}

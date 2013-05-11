var aai = require('./adaptive_ai').GA;

var genome = new aai.Genome();

genome.setChromosomeCount(4);		// 4 states

for (var i=0; i<genome.getChromosomeCount(); i++) {
	var chrom = genome.getChromosome(i);
	chrom.setGeneCount(4);			// 4 transition states

	for (var j=0; j<chrom.getGeneCount(); j++) {
		var gene = chrom.getGene(j);
		gene.setLength(8);			// 7 sensor values (N-1)
		gene.setMutationChance(0.5);
	}
}

for (var i=0; i<8; i++) {
	console.log('Generation ' + (i+1) + ':');
	printGenome(genome);
	genome.mutate();
}

function printGenome(genome) {
	for (var c=0; c<genome.getChromosomeCount(); c++) {
		var chrom = genome.getChromosome(c);

		for (var g=0; g<chrom.getGeneCount(); g++) {
			var gene = chrom.getGene(g);

			var os = 'c'+c+' '+'g'+g+': ';

			for (var i=0, l=gene._sequence.length; i<l; i++) {
				os += gene._sequence[i].toFixed(4) + ' ';
			}

			console.log(os);
		}
	}
	console.log('\n');
}

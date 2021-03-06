var aai = require('./adaptive_ai').GA;
var org = require('./adaptive_ai').Org;

var nStates = 3;

var orgs = []
for (var i=0;i<10;i++) {
	var organism  = new org.Organism();
	orgs.push(organism);

	organism.setStateCount(nStates);
	organism.setStateName(0, '.');
	organism.setStateName(1, '-');
	organism.setStateName(2, 'X');
	organism.setSensorCount(1);

	for (j=0; j<nStates; j++)
		for (k=0; k<nStates; k++)
			organism.setTransition(j, k, .2, [0]);
}

var gen = 0;
for (var i=0; i<100000; i++) {
	update(i, false);
}

for (var i=0; i<10000; i++) {
	update(i, true);
}

function update(i, outFlag) {
	var mod = i%30;

	if (i%100 === 99) {
		gen += 1;

		for (var j=0;j<10;j++) {
			var org = orgs[j]
			org.mutate();
		}
	}

	var ostring = 'gen['+(gen+1)+']: ';

	for (var j=0;j<10;j++) {
		var org = orgs[j]

		if (mod === 29) {
			org.setSensorValue(0,100);
		} else if (mod === 0) {
			org.setSensorValue(0,0);
		}

		var idx = org.getCurrentState();
	
		org.updateState();

		if (outFlag === true)
			ostring += organism.getStateName(idx) + '\t';
	}

	if (outFlag === true)
		console.log(ostring);
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

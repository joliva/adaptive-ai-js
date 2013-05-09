/*****************************************************************************
			Copyright (c) 2013 by John Oliva - All Rights Reserved
*****************************************************************************
	File:				adaptive_ai.js
	Purpose:			Implementation of Adaptive AI FSM as a Common.JS module.
*****************************************************************************/

var ADAPTAI_DEFAULTCHANCE = 0.001;
var ADAPTAI_DEFAULTRATE = 0.1;

Random = function() {
	return Math.random();
}

// --------------------------------------------------------------------------
//	Gene
// --------------------------------------------------------------------------
Gene = function (gene) {
	if (gene === undefined) {
		this._sequence = [];		// array of float
		this._mutationChance = ADAPTAI_DEFAULTCHANCE;
		this._mutationRate = ADAPTAI_DEFAULTRATE;
		this._sequenceLength = 0;
	} else {
		this.copy(gene);
	}
}

Gene.prototype.getSequence = function() {
	return this._sequence;
}

Gene.prototype.setElement = function(idx, element) {
	if (idx >= 0 && idx < this._sequenceLength) {
		this._sequence[idx] = element;
		return true;
	}

	return false;
}

Gene.prototype.getElement = function(idx) {
	if (idx >= 0 && idx < this._sequenceLength) {
		return this._sequence[idx];
	}

	return 0.0;
}

Gene.prototype.setLength = function(length) {
	this._sequenceLength = length;

	this._sequence=[];

	for (var i=0; i<length; i++) {
		this._sequence[i] = 0.0;
	}

	return true;
}

Gene.prototype.getLength = function() {
	return this._sequenceLength;
}

Gene.prototype.setMutationChance = function(chance) {
	// Crop to [0, 1]
	if (chance < 0.0) {
		chance = 0.0;
	} else if (chance > 1.0) {
		chance = 1.0;
	}

	this._mutationChance = chance;

	return true;
}

Gene.prototype.getMutationChance = function() {
	return this._mutationChance;
}

Gene.prototype.setMutationRate = function(mutationRate) {
	this._mutationRate = mutationRate;
	return true;
}

Gene.prototype.getMutationRate = function() {
	return this._mutationRate;
}

Gene.prototype.mutate = function() {
	if (this._sequenceLength <= 0)
		return false;

	for (var i=0, l=this._sequenceLength; i<l; i++) {
		if (Random() <= this._mutationChance)
			this._sequence[i] = this._sequence[i] + (2.0 * Random() - 1.0) * this._mutationRate;
	}

	return true;
}

Gene.prototype.mutateMutationFactors = function(chance, rate) {
	if (Random() <= chance) {
		this._mutationChance += (Random() * 2.0 - 1.0) * rate;
		this._mutationRate += (Random() * 2.0 - 1.0) * rate;
	}

	return true;
}

Gene.prototype.copy = function(gene) {
	this._sequenceLength = gene._sequenceLength;
	this._mutationChance = gene._mutationChance;
	this._mutationRate = gene._mutationRate;
	this._sequence = gene._sequence.slice();

	return this;
}

Gene.prototype.add = function(gene) {
	var newGene = new Gene();

	if (this._sequenceLength != gene._sequenceLength)
		return newGene;

	newGene.setLength(this._sequenceLength);

	// Arithmetic average of elements:
	for (var i=0, l=this._sequenceLength; i<l; i++) {
		newGene._sequence[i] = (this._sequence[i] + gene._sequence[i]) / 2.0;
	}

	return newGene;
}

Gene.prototype.save = function(fstream) {
}

Gene.prototype.load = function(fstream) {
}

// --------------------------------------------------------------------------
//	Chromosome
// --------------------------------------------------------------------------
Chromosome = function (chromosome) {
	if (chromosome === undefined) {
		this._geneList = [];
		this._geneCount = 0;
		this._crossover = true;
		this._crossoverMutationChance = ADAPTAI_DEFAULTCHANCE;
	} else {
		this.copy(chromosome);
	}
}

Chromosome.prototype.getGeneList = function() {
	return this._geneList;
}

Chromosome.prototype.setGene = function(idx, gene) {
	if (idx < 0 || idx >= this._geneCount) {
		return false;
	}

	this._geneList[idx] = new Gene(gene);

	return true;
}

Chromosome.prototype.getGene = function(idx) {
	if (idx < 0 || idx >= this._geneCount) {
		throw new RangeError('Invalid index into geneList');
	}

	return this._geneList[idx];
}

Chromosome.prototype.setGeneCount = function (length) {
	if (length < 0)
		return false;

	this._geneList = [];
	this._geneCount = length;

	for (var i=0; i<length; i++) {
		this._geneList.push(new Gene());
	}

	return true;
}

Chromosome.prototype.getGeneCount = function() {
	return this._geneCount;
}

Chromosome.prototype.copy = function(chrom) {
	this._geneCount = chrom._geneCount;
	this._crossover = chrom._crossover;
	this._crossoverMutationChance = chrom._crossoverMutationChance;

	// Copy gene info:
	for (var i=0; i<this._geneCount; i++) {
		this._geneList.push(new Gene(chrom._geneList[i]));
	} 

	return this;
}

Chromosome.prototype.add = function(chrom) {
	var Temp = new Chromosome();

	if (this._geneCount != chrom._geneCount)
		return Temp;

	Temp.setGeneCount (this._geneCount);

	if (this._crossover) {
		// 50% chance of inheriting crossover trait & mutation rate from either parent:
		if (Random () < 0.5) {
			Temp._crossover = this._crossover;
			Temp._crossoverMutationChance = this._crossoverMutationChance;
		} else {
			Temp._crossover = chrom._crossover;
			Temp._crossoverMutationChance = chrom._crossoverMutationChance;
		}

		for (var i=0; i<this._geneCount; i++) {
			// 50% chance of inheriting gene from either parent:
			if (Random () < 0.5) {
				Temp._geneList[i] = new Gene(this._geneList[i]);
			} else {
				Temp._geneList[i] = new Gene(chrom._geneList[i]);
			}
		} 
	} else {
		// 50% chance of inheriting crossover trait from either parent:
		if (Random () < 0.5) {
			Temp._crossover = this._crossover;
		} else {
			Temp._crossover = chrom._crossover;
		}

		// Mutation chance is numerical average of parents:
		Temp._crossoverMutationChance = (this._crossoverMutationChance + chrom._crossoverMutationChance) / 2.0;

		// Genes are numerical average of parents:
		for (var i=0; i < this._geneCount; i++) {
			Temp._geneList[i] = this._geneList[i].add(chrom._geneList[i]);
		} 
	}

	return Temp;
}

Chromosome.prototype.setCrossoverState = function(state) {
	// state is boolean
   this._crossover = state;

   return true;
}

Chromosome.prototype.getCrossoverState = function() {
   return this._crossover;
}

Chromosome.prototype.setCrossoverMutationChance = function(chance) {
   // Clamp mutation chance:
   if (chance < 0.0) {
      chance = 0.0;
	} else if (chance > 1.0) {
      chance = 1.0;
	}

   this._crossoverMutationChance = chance;

   return true;
}

Chromosome.prototype.getCrossoverMutationChance = function() {
	return this._crossoverMutationChance;
}

Chromosome.prototype.mutateChromosome = function() {
	if (Random () <= this._crossoverMutationChance) {
		if (this._crossover) {
			this._crossover = false;
		} else {
			this.crossover = true;
		}
	}

	return true;
}

Chromosome.prototype.mutateGenes = function() {
	for (var i=0; i<this._geneCount; i++)
		this._geneList[i].mutate();

	return true;
}

Chromosome.prototype.mutate = function() {
	return (this.mutateChromosome() && this.mutateGenes());
}

Chromosome.prototype.mutateMutationFactors = function(chance, rate) {
	if (Random () <= chance) {
		this.setCrossoverMutationChance(this._crossoverMutationChance + (Random() * 2.0 - 1.0) * rate);
	}

	for (var i=0; i<this._geneCount; i++) {
		this._geneList[i].mutateMutationFactors(chance, rate);
	}

	return true;
}

Chromosome.prototype.save = function(fstream) {
}

Chromosome.prototype.load = function(fstream) {
}

// --------------------------------------------------------------------------
//	Genome
// --------------------------------------------------------------------------

Genome = function (genome) {
	if (genome === undefined) {
		this._chromosomeList = NULL;
		this._chromosomeCount = 0;
	} else {
		this.copy(genome);
	}
}

Genome.prototype.setChromosome = function(idx, chrom) {
	if (idx < 0 || idx >= this._chromosomeCount)
		return false;

	this._chromosomeList[idx] = new Chromosome(chrom);

	return true;
}

Genome.prototype.getChromosome = function(idx) {
	if (idx < 0 || idx >= this._chromosomeCount) {
		throw new RangeError('Invalid index into chromosomeList');
	}

	return this._chromosomeList[idx];
}

Genome.prototype.setChromosomeCount = function(count) {
	if (count < 0)
		return false;

	this._chromosomeList = [];
	this._chromosomeCount = count;

	for (var i=0; i<count; i++) {
		this._chromosomeList.push(new Chromosome());
	}

	return true;	
}

Genome.prototype.getChromosomeCount = function() {
	return this._chromosomeCount;
}

// --------------------------------------------------------------------------
exports.GA={};
exports.GA.Gene = Gene;
exports.GA.Chromosome = Chromosome;
exports.GA.Genome = Genome;


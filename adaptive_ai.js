/*****************************************************************************
			Copyright (c) 2013 by John Oliva - All Rights Reserved
*****************************************************************************
	File:				adaptive_ai.js
	Purpose:			Implementation of Adaptive AI FSM as a Common.JS module.
*****************************************************************************/

// constants
var ADAPTAI_DEFAULT_CHANCE = 0.001;
var ADAPTAI_DEFAULT_RATE = 0.1;

Random = function() {
	return Math.random();
}

GA = {
	defaultChance: ADAPTAI_DEFAULT_CHANCE,
	defaultRate: ADAPTAI_DEFAULT_RATE
}

// --------------------------------------------------------------------------
//	Genetic Algorithm: Gene
// --------------------------------------------------------------------------
Gene = function (gene) {
	if (gene === undefined) {
		this._sequence = [];		// array of float
		this._mutationChance = GA.defaultChance;
		this._mutationRate = GA.defaultRate;
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
//	Genetic Algorithm: Chromosome
// --------------------------------------------------------------------------
Chromosome = function (chromosome) {
	if (chromosome === undefined) {
		this._geneList = [];
		this._geneCount = 0;
		this._crossover = true;
		this._crossoverMutationChance = ADAPTAI_DEFAULT_CHANCE;
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

	this._geneList = [];

	// Copy gene info:
	for (var i=0, l=this._geneCount; i<l; i++) {
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

		for (var i=0, l=this._geneCount; i<l; i++) {
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
	for (var i=0, l=this._geneCount; i<l; i++)
		this._geneList[i].mutate();

	return true;
}

Chromosome.prototype.mutate = function() {
	return (this.mutateChromosome() && this.mutateGenes());
}

Chromosome.prototype.setMutationFactors = function(chance, rate) {
	for (var i=0, l=this._geneCount; i<l; i++) {
		var gene = this._geneList[i];
		gene.setMutationChance(chance);
		gene.setMutationRate(rate);
	}

	return true;
}

Chromosome.prototype.mutateMutationFactors = function(chance, rate) {
	if (Random () <= chance) {
		this.setCrossoverMutationChance(this._crossoverMutationChance + (Random() * 2.0 - 1.0) * rate);
	}

	for (var i=0, l=this._geneCount; i<l; i++) {
		this._geneList[i].mutateMutationFactors(chance, rate);
	}

	return true;
}

Chromosome.prototype.save = function(fstream) {
}

Chromosome.prototype.load = function(fstream) {
}

// --------------------------------------------------------------------------
//	Genetic Algorithm: Genome
// --------------------------------------------------------------------------

Genome = function (genome) {
	if (genome === undefined) {
		this._chromosomeList = [];
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

Genome.prototype.copy = function(genome) {
	this.setChromosomeCount(genome._chromosomeCount);

	this._chromosomeList = [];

	for (var i=0, l=this._chromosomeCount; i<l; i++) {
		this._chromosomeList.push(new Chromosome(genome._chromosomeList[i]));
	}

	return this;
}

Genome.prototype.add = function(genome) {
	var Temp = new Genome();

	if (this._chromosomeCount !== genome._chromosomeCount) {
		throw new RangeError('Attempt to add genomes with differing chromosome counts.');
	}

	Temp.setChromosomeCount(this._chromosomeCount);

	for (var i=0, l=this._chromosomeCount; i<l; i++) {
		Temp._chromosomeList[i] = this._chromosomeList[i].add(genome._chromosomeList[i]);
	}

	return Temp;
}

Genome.prototype.mutate = function() {
	for (var i=0, l=this._chromosomeCount; i<l; i++) {
		this._chromosomeList[i].mutate();
	}

	return true;
}

Genome.prototype.setMutationFactors = function(chance, rate) {
	for (var i=0, l=this._chromosomeCount; i<l; i++) {
		this._chromosomeList[i].setMutationFactors(chance, rate);
	}

	GA.defaultChance = chance;
	GA.defaultRate = rate;

	return true;
}

Genome.prototype.mutateMutationFactors = function(chance, rate) {
	for (var i=0, l=this._chromosomeCount; i<l; i++) {
		this._chromosomeList[i].mutateMutationFactors(chance, rate);
	}

	return true;
}

Genome.prototype.save = function(fstream) {
}

Genome.prototype.load = function(fstream) {
}

// --------------------------------------------------------------------------
//	Organism
// --------------------------------------------------------------------------

Organism = function(organism) {
	if (organism === undefined) {
		this._states = [];
		this._sensors = [];
		this._stateCount = 0;
		this._sensorCount = 0;
		this._currentState = 0;
		this._orgGenome = new Genome();
	} else {
		this.copy(organism);
	}
}

Organism.prototype.copy = function(organism) {
	this._states = [];
	this._sensors = [];
	this._stateCount  = organism._stateCount;
	this._sensorCount = organism._sensorCount;

	States  = new State  [StateCount];
	Sensors = new Sensor [SensorCount];

	for (var i=0, l=this._stateCount; i<l; i++) {
		this._states.push(new Organism.State(organism._states[i]));
	}

	for (var i=0, l=this._sensorCount; i<l; i++) {
		this._sensors.push(new Organism.Sensor(organism._sensors[i]));
	}

	this._orgGenome.copy(organism._orgGenome);

	return this;
}

Organism.prototype.add = function(organism) {
	temp = new Organism();

	if (this._stateCount != organism._stateCount || this._sensorCount != organism._sensorCount) {
		throw new RangeError('Incompatible # of states/sensors for copying Organsim');
	}

	temp.copy(this);

	temp._orgGenome = this._orgGenome.add(organism._orgGenome);

	// Mutate the offspring's genome:
	temp.mutate ();

	return temp;
}

Organism.prototype.getStateName = function(idx) {
	if (idx < 0 || idx >= this._stateCount)
		return '';

	return this._states[idx]._name;
}

Organism.prototype.setStateName = function(idx, name) {
	if (idx < 0 || idx >= this._stateCount)
		return false;

	this._states[idx].setName(name);

	return true;
}

Organism.prototype.setStateCount = function(count) {
	if (count < 0)
		return false;

	// One chromosome for each state:
	this._stateCount = count;

	this._states=[];

	for (var i=0, l=this._stateCount; i<l; i++) {
		this._states.push(new Organism.State());
	}

	if (this._orgGenome.setChromosomeCount(this._stateCount) === false) {
		return false;
	}

	// StateCount genes for each state:
	for (var i=0, l=this._stateCount; i<l; i++) {
		var chrom = this._orgGenome.getChromosome(i);
		chrom.setGeneCount(this._stateCount);
	}

	// Now set sensor count since it depends on state count:
	return this.setSensorCount(this._sensorCount);
}

Organism.prototype.getStateCount = function() {
	return this._stateCount;
}

Organism.prototype.getStateIndex = function(name) {
	for (var i=0, l=this.stateCount; i<l; i++) {
		if (name === this._states[i]._name) {
			return i;
		}
	}

	return null;
}

Organism.prototype.getSensorValue = function(arg) {
	// inspect argument type to infer argument semantics
	var argtype = typeof arg;

	if (argtype === 'string') {
		// arg is inferred to be 'name'
		var name = arg;

		for (var i=0, l=this._sensorCount; i<l; i++) {
			if (name === this._sensors[i]._name)
				return this._sensors[i]._value;
		}
	} else if (argtype === 'number') {
		// arg is inferred to be 'index'
		var idx = arg;

		if (idx < 0 || idx >= this._sensorCount)
			return 0.0;

		return this._sensors[idx]._value;
	}

	return 0.0;
}

Organism.prototype.setSensorValue = function(arg, value) {
	// inspect argument type to infer argument semantics
	var argtype = typeof arg;

	if (argtype === 'string') {
		// arg is inferred to be 'name'
		var name = arg;

		for (var i=0, l=this._sensorCount; i<l; i++) {
			if (name === this._sensors[i]._name) {
				this._sensors[i].setValue(value);
	
				return true;
			}
		}

		return false;
	} else if (argtype === 'number') {
		// arg is inferred to be 'index'
		var idx = arg;

		if (idx < 0 || idx >= this._sensorCount)
			return false;

		this._sensors[idx].setValue(value);

		return true;
	}

	return false;
}

Organism.prototype.getSensorName = function(idx) {
	if (idx < 0 || idx >= this._sensorCount)
		return null;

	return this._sensors[idx]._name;
}

Organism.prototype.setSensorName = function(idx, name) {
	if (idx < 0 || idx >= this._sensorCount)
		return false;

	this._sensors[idx].setName(name);

	return true;
}

Organism.prototype.setSensorCount = function(count) {
	if (count < 0) {
		return false;
	}
	
	this._sensorCount = count;
	this._sensors= [];

	for (var i=0, l=this._sensorCount; i<l; i++) {
		this._sensors.push(new Organism.Sensor());
	}

	for (var i=0, l=this._stateCount; i<l; i++) {
		var chrom = this._orgGenome.getChromosome(i);

		for (var j=0, ll=this._stateCount; j<ll; j++) {
			var gene = chrom.getGene(j);

			// Allocate enough room for each coefficient (base + sensor coeff's):
			gene.setLength(1 + this._sensorCount);
		}
	}

	return true;
}

Organism.prototype.getSensorCount = function() {
	return this._sensorCount;
}

Organism.prototype.getSensorIndex = function(name) {
	for (var i=0, l=this._sensorCount; i<l; i++) {
		if (name === this._sensors[i]._name) {
			return i;
		}
	}

	return null;
}

Organism.prototype.setTransition  = function(idx1, idx2, baseChance, sensorCoeff) {
	if (idx1 < 0 || idx1 >= this._stateCount || idx2 < 0 || idx2 >= this._stateCount)
		return false;

	var chrom = this._orgGenome.getChromosome(idx1);
	var gene = chrom.getGene(idx2);

	gene.setElement(0, baseChance);

	for (var i=0, l=this._sensorCount; i<l; i++)
		gene.setElement(1 + i, sensorCoeff[i]);

	return true;
}

Organism.prototype.getCurrentState = function() {
	return this._currentState;
}

Organism.prototype.setCurrentState = function(arg) {
	// inspect argument type to infer argument semantics
	var argtype = typeof arg;

	if (argtype === 'string') {
		// arg is inferred to be 'name'
		var name = arg;

		for (var i=0, l=this._stateCount; i<l; i++) {
			if (this._states[i]._name === name) {
				this._currentState = i;

				return true;
			}
		}

		return false;
	} else if (argtype === 'number') {
		// arg is inferred to be 'index'
		var idx = arg;

		if (idx < 0 || idx >= this._stateCount)
			return false;

		this._currentState = idx;

		return true;
	}

	return false;
}

Organism.prototype.mutate = function() {
	return this._orgGenome.mutate();
}

Organism.prototype.setMutationFactors = function(chance, rate) {
	return this._orgGenome.setMutationFactors(chance, rate);
}

Organism.prototype.mutateMutationFactors = function(chance, rate) {
	return this._orgGenome.mutateMutationFactors(chance, rate);
}

Organism.prototype.updateState = function() {
	var prob = [];

	for (var i=0, l=this._stateCount; i<l; i++) {
		prob.push(0.0);	
	}

	var totalProb = 0.0;

	var chrom = this._orgGenome.getChromosome(this._currentState);

	for (var i=0, l=this._stateCount; i<l; i++) {
		var gene = chrom.getGene(i);

		// Base chance:
		prob[i] = gene.getElement(0);

		totalProb += prob[i];

		// Sensor coefficients:
		for (var j=1, ll=this._sensorCount; j<=ll; j++) {
			var p = gene.getElement(j) * this._sensors[j - 1]._value;
			
			prob[i] += p;
			totalProb += p;
		}
	}

	// Normalize probability:
	for (var i=0, l=this._stateCount; i<l; i++) {
		prob[i] /= totalProb;
	}

	// cdf - cumulative distribution function
	for (var i=this._stateCount - 1; i >= 0; i--) {
		for (var j=0; j<i; j++)
			prob[i] += prob[j];
	}

	var choice = Random();

	var nextState = this._currentState;

	for (var i=0, l=this._stateCount; i<l; i++) {
		if (choice <= prob[i]) {
			nextState = i;
			break;
		}
	}

	this._currentState = nextState;

	return true;
}

Organism.prototype.save = function(fstream) {
}

Organism.prototype.load = function(fstream) {
}

// --------------------------------------------------------------------------
//	Organism: State
// --------------------------------------------------------------------------

Organism.State = function(state) {
	if (state === undefined) {
		this._name = 'NOT_SET';
	} else {
		this.copy(state);
	}
}

Organism.State.prototype.setName = function(newName) {
	this._name = newName;
	return true;
}

Organism.State.prototype.copy = function(state) {
	this._name = state.name;
	return this;
}

Organism.State.prototype.save = function(fstream) {
}

Organism.State.prototype.load = function(fstream) {
}

// --------------------------------------------------------------------------
//	Organism: Sensor
// --------------------------------------------------------------------------

Organism.Sensor = function(sensor) {
	if (sensor === undefined) {
		this._value = 0.0;
		this._name = 'NOT_SET';
	} else {
		this.copy(sensor);
	}
}

Organism.Sensor.prototype.setName = function(newName) {
	this._name = newName;
	return true;
}

Organism.Sensor.prototype.getName = function() {
	return this._name
}

Organism.Sensor.prototype.setValue = function(value) {
	this._value = value;

	return true;
}

Organism.Sensor.prototype.getValue = function() {
	return this._value;
}

Organism.Sensor.prototype.copy = function(sensor) {
	this._name = sensor._name;
	this._value = sensor._value;

	return this;
}

Organism.Sensor.prototype.save = function(fstream) {
}

Organism.Sensor.prototype.load = function(fstream) {
}

// --------------------------------------------------------------------------

exports.GA = GA;
exports.GA.Gene = Gene;
exports.GA.Chromosome = Chromosome;
exports.GA.Genome = Genome;

exports.Org={};
exports.Org.Organism = Organism;
exports.Org.Organism.State = Organism.State;
exports.Org.Organism.Sensor = Organism.Sensor;


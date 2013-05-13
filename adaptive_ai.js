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
//	Genetic Algorithm: Gene
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
//	Genetic Algorithm: Chromosome
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

	this._geneList = [];

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

	for (var i=0; i<this._chromosomeCount; i++) {
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

	for (var i=0; i<this._chromosomeCount; i++) {
		Temp._chromosomeList[i] = this._chromosomeList[i].add(genome._chromosomeList[i]);
	}

	return Temp;
}

Genome.prototype.mutate = function() {
	for (var i=0; i<this._chromosomeCount; i++) {
		this._chromosomeList[i].mutate();
	}

	return true;
}

Genome.prototype.mutateMutationFactors = function(chance, rate) {
	for (var i=0; i<this._chromosomeCount; i++) {
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
		this._states.push(new Organism.State(organism._states[i]);
	}

	for (var i=0, l=this._sensorCount; i<l; i++) {
		this._sensors.push(new Organism.Sensor(organism._sensors[i]);
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

	this._states[idx]._setName(name);

	return true;
}

Organism.prototype.setStateCount = function(count) {
	if (count < 0)
		return false;

	// One chromosome for each state:
	this._stateCount = count;

	this._states=[];

	for (var i=0, l=this._stateCount; i<l; i++) {
		this.states.push(new Organism.State());
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

/*
int Organism::GetStateCount () const {
	return StateCount;
}

int Organism::GetStateIndex (std::string Name) const {
	for (int i = 0; i < StateCount; i++) {
		if (Name == States [i].Name) {
			return i;
		}
	}

	return -1;
}

float Organism::GetSensorValue (std::string Name) const {
	for (int i = 0; i < SensorCount; i++) {
		if (Name == Sensors [i].Name)
			return Sensors [i].Value;
	}

	return 0.0F;
}

float Organism::GetSensorValue (int Index) const {
	if (Index < 0 || Index >= SensorCount)
		return 0.0F;

	return Sensors [Index].Value;
}

bool Organism::SetSensorValue (std::string Name, float Value) {
	for (int i = 0; i < SensorCount; i++) {
		if (Name == Sensors [i].Name) {
			Sensors [i].SetValue (Value);

			return true;
		}
	}

	return false;
}

bool Organism::SetSensorValue (int Index, float Value) {
	if (Index < 0 || Index >= SensorCount)
		return false;

	Sensors [Index].SetValue (Value);

	return true;
}

std::string Organism::GetSensorName (int Index) const {
	if (Index < 0 || Index >= SensorCount)
		return NULL;

	return Sensors [Index].Name;
}

bool Organism::GetSensorName (int Index, std::string* Name) const {
	if (Index < 0 || Index >= SensorCount)
		return false;

	*Name = Sensors [Index].Name;

	return true;
}

bool Organism::SetSensorName (int Index, std::string Name) {
	if (Index < 0 || Index >= SensorCount)
		return false;

	Sensors [Index].SetName (Name);

	return true;
}

bool Organism::SetSensorCount (int Count) {

	if (Count < 0) {
		return false;
	}
	
	SensorCount = Count;

	if (Count > 0) {
		if (Sensors == NULL) {
			Sensors = new Sensor [SensorCount];
		}
	}

	for (int i = 0; i < StateCount; i++) {
		Chromosome &Chrom = OrgGenome.GetChromosome (i);

		for (int j = 0; j < StateCount; j++) {
			Gene &G = Chrom.GetGene (j);

			// Allocate enough room for each coefficient
			// (base + sensor coeff's):
			G.SetLength (1 + SensorCount);
		}
	}

	return true;
}

int Organism::GetSensorCount () const {
	return SensorCount;
}

int Organism::GetSensorIndex (std::string Name) const {
	for (int i = 0; i < SensorCount; i++) {
		if (Name == Sensors [i].Name) {
			return i;
		}
	}

	return -1;
}

bool Organism::SetTransition (int Index1, int Index2, float BaseChance, const float *SensorCoeff) {
	if (Index1 < 0 || Index1 >= StateCount || Index2 < 0 || Index2 >= StateCount)
		return false;

	Chromosome &Chrom = OrgGenome.GetChromosome (Index1);

	Gene &G = Chrom.GetGene (Index2);

	G.SetElement (0, BaseChance);

	for (int i = 0; i < SensorCount; i++)
		G.SetElement (1 + i, SensorCoeff [i]);

	return true;
}

int Organism::GetCurrentState () const {
	return CurrentState;
}

bool Organism::GetCurrentState (std::string* Name) const {
	*Name = States [CurrentState].Name;
	return true;
}

bool Organism::SetCurrentState (std::string Name) {
	for (int i = 0; i < StateCount; i++) {
		if (States [i].Name == Name) {
			CurrentState = i;

			return true;
		}
	}

	return false;
}

bool Organism::SetCurrentState (int Index) {
	if (Index < 0 || Index >= StateCount)
		return false;

	CurrentState = Index;

	return true;
}

bool Organism::Mutate () {
	return OrgGenome.Mutate ();
}

bool Organism::MutateMutationFactors (float Chance, float Rate) {
	return OrgGenome.MutateMutationFactors (Chance, Rate);
}

bool Organism::UpdateState () {
	float *Prob = new float [StateCount];

	float TotalProb = 0.0F;

	int i;

	Chromosome &Chrom = OrgGenome.GetChromosome (CurrentState);

	for (i = 0; i < StateCount; i++) {
		Gene &G = Chrom.GetGene (i);

		// Base chance:
		Prob [i]  = G.GetElement (0);

		TotalProb  += Prob [i];

		// Sensor coefficients:
		for (int j = 1; j <= SensorCount; j++) {
			float p = G.GetElement (j) * Sensors [j - 1].Value;
			
			Prob [i] += p;

			TotalProb  += p;
		}
	}

	// Normalize probability:
	for (i = 0; i < StateCount; i++)
		Prob [i] /= TotalProb;

	// cdf - cumulative distribution function
	for (i = StateCount - 1; i >= 0; i--) {
		for (int j = 0; j < i; j++)
			Prob [i] += Prob [j];
	}

	float Choice = Random ();

	int NextState = CurrentState;

	for (i = 0; i < StateCount; i++) {
		if (Choice <= Prob [i]) {
			NextState = i;

			break;
		}
	}

	CurrentState = NextState;

	delete [] Prob;

	return true;
}

bool Organism::Save (std::fstream &File) const {
	// File format:
	//			 StateCount			sizeof (int)
	//			 SensorCount		  sizeof (int)
	//			 CurrentState		 sizeof (int)
	//			 States				 sizeof (State) * StateCount
	//			 Sensors				sizeof (Sensor) * SensorCount
	//			 OrgGenome			 varies
	File.write ((const char *) &StateCount,	sizeof (int));
	File.write ((const char *) &SensorCount,  sizeof (int));
	File.write ((const char *) &CurrentState, sizeof (int));
	
	// Save states:
	int i;
	for (i = 0; i < StateCount; i++) {
		if (!States [i].Save (File))
			return false;
	}

	// Save sensors:
	for (i = 0; i < SensorCount; i++) {
		if (!Sensors [i].Save (File))
			return false;
	}

	// Save the genome:
	if (!OrgGenome.Save (File))
		return false;

	if (!File.good ())
		return false;

	return true;
}

bool Organism::Load (std::fstream &File) {
	File.read ((char *) &StateCount,	sizeof (int));
	File.read ((char *) &SensorCount,  sizeof (int));
	File.read ((char *) &CurrentState, sizeof (int));
	
	// Allocate memory for states:
	delete [] States;
	States = new State [StateCount];

	// Load states:
	int i;
	for (i = 0; i < StateCount; i++) {
		if (!States [i].Load (File))
			return false;
	}

	// Allocate memory for the sensors:
	delete [] Sensors;
	Sensors = new Sensor [SensorCount];

	// Load sensors:
	for (i = 0; i < SensorCount; i++) {
		if (!Sensors [i].Load (File))
			return false;
	}

	// Load the genome:
	if (!OrgGenome.Load (File))
		return false;

	if (!File.good ())
		return false;

	return true;
}
*/

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
exports.GA={};
exports.GA.Gene = Gene;
exports.GA.Chromosome = Chromosome;
exports.GA.Genome = Genome;
exports.Org={};
exports.Org.Organism = Organism;
exports.Org.Organism.State = Organism.State;
exports.Org.Organism.Sensor = Organism.Sensor;


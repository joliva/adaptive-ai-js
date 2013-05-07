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

Gene = function (gene) {
		gene = gene || {};
		this._sequence = gene._sequence || null,
		this._mutationChance = gene._mutationChance || ADAPTAI_DEFAULTCHANCE,
		this._mutationRate = gene._mutationRate || ADAPTAI_DEFAULTRATE,
		this._sequenceLength = gene._sequenceLength || 0;
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

Gene.prototype.set = function(gene) {
	setLength(gene._sequenceLength);

	this._mutationChance = gene._mutationChance;
	this._mutationRate = gene._mutationRate;

	for (var i=0, l=this._sequenceLength; i<l; i++)
		this._sequence[i] = gene._sequence[i];

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

exports.Gene = Gene

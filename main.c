文件名：Classifier.h

#ifndef CLASSIFIER_H
#define CLASSIFIER_H

#include <iostream>
#include <stdio.h>

#define SELF                                    0
#define NONSELF                                 1
#define MASKVALUE                               2

// detector class
class Detector
{
        public:
                Detector(const unsigned int length);
                Detector::~Detector(void);

                unsigned int length;
                unsigned int *value;
                double threshold;
                unsigned int type;

                void save(FILE *outputStream);
                void show(void) { save(stdout); };
};

#endif

/**********************************************/
//文件名：Classifier.cpp

#include "Classifier.h"

// detector class public methods

Detector::Detector(const unsigned int length)
{
        this->length = length;
        threshold = 0.0;
        value = new unsigned int [length];
        type = 0;
}

Detector::~Detector(void)
{
        delete [] value;
}

void Detector::save(FILE *outputStream)
{
        register unsigned int i;

        fprintf(outputStream, \
                "%-3d %-.10f %-1d\n", \
                length, \
                threshold, \
                type \
                );

        for(i = 0; i < length; i++)
                        fprintf(outputStream, "%-1d ", value[i]);
        fprintf(outputStream, "\n");

        fflush(outputStream);
}

/**********************************************/
//文件名：EvolutionaryAlgorithm.h

#ifndef EVOLUTIONARYALGORITHM_H
#define EVOLUTIONARYALGORITHM_H

#include "Classifier.h"
#include <stdio.h>

// genome class
class Genome
{
        public:
                Genome(const unsigned int length);
                ~Genome(void);

                unsigned int size;
                unsigned int *locus;
                unsigned int type;
                double mutationProbability;
                double crossoverProbability;
                double fitness, scaledFitness;
                unsigned int thresholdLength, patternLength;
                double generalityBias;
                double typeBias;

                void copyGenome(Genome *genome);
                void uniformCrossover(Genome *genome1, Genome *genome2);
                void mutateBinary(void);
                void rand()omiseBinary(void);
                void setDetector(Detector *detector);
                void save(FILE *outputStream);
                void show(void) { save(stdout); };
};

// species class
class Species
{
        public:
                Species(const unsigned int speciesSize, const unsigned int genomeLength);
                ~Species(void);

                unsigned int speciesSize;
                Genome **genome;
                unsigned int fittestIndividual;
                double speciesScaledFitnessSum;
                double meanSpeciesFitness;

                Genome *FPSelection(void);
                void rand()omise(void);
                void copySpecies(Species *species);
                void save(FILE *outputStream);
                void show(void) { save(stdout); };
};

#endif

/**********************************************************/
//文件名：EvolutionaryAlgorithm.cpp

#include "EvolutionaryAlgorithm.h"
#include <stdlib.h>

// genome class public methods
Genome::Genome(const unsigned int length)
{
        thresholdLength = 8;
        patternLength = length;
        size = thresholdLength + 2 * patternLength;
        locus = new unsigned int [size];
        mutationProbability = 2.0 / double(size);
        crossoverProbability = 0.6;
        generalityBias = typeBias = 0.5;
        fitness = 0.0;
        type = 0;
}

Genome::~Genome(void)
{
        delete [] locus;
}

void Genome::copyGenome(Genome *genome)
{
        register unsigned int i = size;
        register unsigned int *from = genome->locus;
        register unsigned int *to = locus;

        while(i--)
                        to[i] = from[i];
        mutationProbability = genome->mutationProbability;
        crossoverProbability = genome->crossoverProbability;
        fitness = genome->fitness;
        scaledFitness = genome->scaledFitness;
        generalityBias = genome->generalityBias;
        size = genome->size;
        patternLength = genome->patternLength;
        thresholdLength = genome->thresholdLength;
        type = genome->type;
}

void Genome::uniformCrossover(Genome *genome1, Genome *genome2)
{
        register unsigned int i = size;
        register unsigned int *from1 = genome1->locus;
        register unsigned int *from2 = genome2->locus;
        register unsigned int *to = locus;
        register double cp = crossoverProbability;

        while(i--)
        {
                if(drand()48() < cp)
                        to[i] = from1[i];
                else
                        to[i] = from2[i];
        }
        if(drand48() < cp)
                type = genome1->type;
        else
                type = genome2->type;
}

void Genome::mutateBinary(void)
{
        register unsigned int i = size;
        register unsigned int *loci = locus;
        register double mp = mutationProbability;

        while(i--)
                if(drand48() < mp)
                        loci[i] = 1 - loci[i];
        if(drand48() < mp)
                        type = 1 - type;
}

void Genome::randomiseBinary(void)
{
        register unsigned int index, i;

        index = 0;

        i = thresholdLength;
        while(i--)
                locus[index++] = int((double(rand()) * 2.0) / double(RAND_MAX + 1.0));
        i = patternLength;
        while(i--)
                locus[index++] = int((double(rand()) * 2.0) / double(RAND_MAX + 1.0));
        i = patternLength;
        while(i--)
                if(drand48() < generalityBias)
                        locus[index++] = 0;
                else
                        locus[index++] = 1;
        if(drand48() < typeBias)
                type = SELF;
        else
                type = NONSELF;
}

void Genome::save(FILE *outputStream)
{
        register unsigned int i;
        Detector *detector = new Detector(patternLength);

        fprintf(outputStream, \
                        "%-3d %-3d %-3d %-1d %-10f %-10f %-10f %-10f %-10f %-10f\n", \
                        size, \
                        thresholdLength, \
                        patternLength, \
                        type, \
                        fitness, \
                        scaledFitness, \
                        mutationProbability, \
                        crossoverProbability, \
                        generalityBias, \
                        typeBias \
                        );

        for(i = 0; i < size; i++)
                fprintf(outputStream, "%-2d ", locus[i]);
        fprintf(outputStream, "\n");

        setDetector(detector);
        detector->save(outputStream);

        delete detector;

        fflush(outputStream);
}

void Genome::setDetector(Detector *detector)
{
        register unsigned int i, loci = 0, sum, lastLoci;

        // set activation threshold
        // gray coding for threshold gene
        sum = lastLoci = locus[loci++];
        while(loci < thresholdLength)
        {
                sum = (sum << 1) | (lastLoci ^ locus[loci]);
                lastLoci = locus[loci++];
        }
        detector->threshold = double(sum) / 255.0; // !!!!!!!!!!!!!!!!!!!!
        for(i = 0; i < patternLength; i++)
                detector->value[i] = locus[loci++];
        for(i = 0; i < patternLength; i++)
                if(!locus[loci++])
                        detector->value[i] = MASKVALUE;
        detector->type = type;
        detector->length = patternLength;
}

// species class public methods
Species::Species(const unsigned int speciesSize, \
                const unsigned int genomeLength)
{
        register unsigned int i = speciesSize;

        this->speciesSize = speciesSize;
        fittestIndividual = 0;
        speciesScaledFitnessSum = meanSpeciesFitness = 0.0;
        genome = new Genome * [speciesSize];
        while(i--)
                genome[i] = new Genome(genomeLength);
}

Species::~Species(void)
{
        register unsigned int i = speciesSize;

        while(i--)
                delete genome[i];
        delete genome;
}

Genome *Species::FPSelection(void)
{
        register unsigned int i = 0;
        register double dtmp1, dtmp2;

        dtmp1 = drand48() * speciesScaledFitnessSum;
        dtmp2 = 0.0;
        while((i < speciesSize) && ((dtmp2 = dtmp2 + genome[i]->scaledFitness) \
                        < dtmp1))
                i++;
        return((i < speciesSize) ? genome[i] : genome[i - 1]);
}

void Species::randomise(void)
{
        register unsigned int i = speciesSize;

        while(i--)
                genome[i]->randomiseBinary();
}

void Species::save(FILE *outputStream)
{
        fprintf(outputStream, \
                        "%-4d %-4d %-5.10f %-.10f\n", \
                        speciesSize, \
                        fittestIndividual, \
                        speciesScaledFitnessSum, \
                        meanSpeciesFitness \
                        );

        genome[fittestIndividual]->save(outputStream);

        fflush(outputStream);
}

void Species::copySpecies(Species *species)
{
        register unsigned int i = species->speciesSize;

        speciesSize = i;
        while(i--)
                genome[i]->copyGenome(species->genome[i]);
        fittestIndividual = species->fittestIndividual;
        speciesScaledFitnessSum = species->speciesScaledFitnessSum;
        meanSpeciesFitness = species->meanSpeciesFitness;
}

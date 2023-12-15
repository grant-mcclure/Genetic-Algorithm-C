// Include everything necessary here
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

extern double Objective_function(int NUM_VARIABLES, double x[NUM_VARIABLES]);  //Import OF

double generate_random(double min, double max)
{
    if (max < min)
    {
        printf("Max cannot be less than the minimum value\n");
        exit(1);
    }

    double random_number = (double)rand() / RAND_MAX;    // Generate a random value between 0 and 1
    random_number = (random_number * (max - min)) + min; //scale to in rangee


    return random_number;
}

unsigned int generate_int()
{
    unsigned int random_int = rand();

    if (random_int > UINT_MAX || random_int < 0) //double check correct size
    {
        printf("Not an unsigned integer\n"); //checl error
        exit(1);
    }
    return random_int;
}

void generate_population(int POPULATION_SIZE, int NUM_VARIABLES, double population[POPULATION_SIZE][NUM_VARIABLES], double Lbound[NUM_VARIABLES], double Ubound[NUM_VARIABLES])
{
    for (int i = 0; i < POPULATION_SIZE; i++)
    {
        for (int j = 0; j < NUM_VARIABLES; j++)
        {
            population[i][j] = generate_random(Lbound[j], Ubound[j]); //random generic population within bounds
        }
    }
}
//call objecive fucntion
void compute_objective_function(int POPULATION_SIZE, int NUM_VARIABLES, double population[POPULATION_SIZE][NUM_VARIABLES], double fitness[POPULATION_SIZE])
{
    for (int i = 0; i < POPULATION_SIZE; i++)
    {
        // Call the Objective_function from OF.c with the appropriate arguments
        fitness[i] = Objective_function(NUM_VARIABLES, population[i]); //computes fitness and stores in array
    }

   
}

void crossover(int POPULATION_SIZE, int NUM_VARIABLES, double fitness[POPULATION_SIZE], double new_population[POPULATION_SIZE][NUM_VARIABLES], double population[POPULATION_SIZE][NUM_VARIABLES], double crossover_rate)
{
   int selected_indices[POPULATION_SIZE]; //create an array of selected for indices
    
    // Select random indices for crossover
    for (int i = 0; i < POPULATION_SIZE; i++) {
        selected_indices[i] = rand() % POPULATION_SIZE;
    }

    // Perform crossover
    for (int i = 0; i < POPULATION_SIZE; i += 2) {
        if ((double)rand() / RAND_MAX < crossover_rate) {
            int parent1_index = selected_indices[i]; //create parent1 index
            int parent2_index = selected_indices[i + 1]; //create oarent 2 index
            int crosspoint = rand() % (NUM_VARIABLES - 1) + 1; // Ensure crosspoint is  not at end of array

            // create new individuals using one point crossover
            for (int j = 0; j < crosspoint; j++) {
                new_population[i][j] = population[parent1_index][j]; //mix genes
                new_population[i + 1][j] = population[parent2_index][j]; //mix genes
            }
            for (int j = crosspoint; j < NUM_VARIABLES; j++) {  
                new_population[i][j] = population[parent2_index][j]; //mix other gene
                new_population[i + 1][j] = population[parent1_index][j];//miz other genes after crosspoint
            }
        } else {
            // If crossover doesn't happen, copy the parents to the new population
            for (int j = 0; j < NUM_VARIABLES; j++) {
                new_population[i][j] = population[selected_indices[i]][j];
                new_population[i + 1][j] = population[selected_indices[i + 1]][j];
            }
        }
    }
}

void mutate(int POPULATION_SIZE, int NUM_VARIABLES, double new_population[POPULATION_SIZE][NUM_VARIABLES], double population[POPULATION_SIZE][NUM_VARIABLES], double Lbound[NUM_VARIABLES], double Ubound[NUM_VARIABLES], double mutate_rate)
{


    //this is directly translated from the python mutate code provided
    
    //calcluate number of genes in popluation
    int total_genes = POPULATION_SIZE * NUM_VARIABLES;

    // calculate mutate rate
    int genes_to_mutate = (int)(total_genes * mutate_rate);

    // Perform mutations
    for (int i = 0; i < genes_to_mutate; i++)
    {
        //find index to mutate
        int gene_index = generate_int() % total_genes;

        // Calculate the row and column of the gene to mutate
        int row = gene_index / NUM_VARIABLES;
        int col = gene_index % NUM_VARIABLES;

        // Mutate the gene within the bounds
        new_population[row][col] = generate_random(Lbound[col], Ubound[col]);
    }
}

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include "functions.h"

// function prototypes
void printProgressBar(int iteration, int total);
double calculateDistance(double x, double y);
double adaptiveCrossoverRate(int generation, int MAX_GENERATIONS, double startRate, double endRate);
double adaptiveMutationRate(int generation, int MAX_GENERATIONS, double startRate, double endRate);

int main(int argc, char *argv[])
{
    /*all storage variabels*/

    srand((unsigned int)time(NULL)); // seed rng
    int best_index = -1;
    double best_fitness = __DBL_MAX__;
    int closest_solution_index = -1;
    double global_best_fitness = __DBL_MAX__;
    int no_improvement_counter = 0;
    double previous_best_fitness = __DBL_MAX__;

    // <YOUR CODE: Handle the possible errors in input data given by the user and say how to execute the code>
    if (argc != 6) // check correct number of arguments
    {
        printf("Usage: %s POPULATION_SIZE MAX_GENERATIONS crossover_rate mutate_rate stop_criteria\n", argv[0]);
        printf("%d", argc);
        return 1;
    }

    // <YOUR CODE: Assign all inputs given by the user argv[i]> like

    // cconvert argymnets to coorenct date types
    unsigned int POPULATION_SIZE = atoi(argv[1]);
    unsigned int MAX_GENERATIONS = atoi(argv[2]);
    double crossover_rate = atof(argv[3]);
    double mutate_rate = atof(argv[4]);
    double stop_criteria = atof(argv[5]);

    const int NO_IMPROVEMENT_THRESHOLD = MAX_GENERATIONS/3 ; // set an arbirtay convergnece criteria

    // possive erors
    if (mutate_rate >= 1 || crossover_rate >= 1)
    {
        printf("incorrect mutate or crossover rate\n");
        exit(1);
    }
    if (POPULATION_SIZE > __UINT64_MAX__ || MAX_GENERATIONS > __UINT64_MAX__)
    {
        printf("population size or max generations too big too big\n");
        exit(1);
    }
    if(POPULATION_SIZE < 0 || MAX_GENERATIONS <0 || crossover_rate <0 || mutate_rate < 0 || stop_criteria< 0){
        printf("no negative input arguments\n");
        exit(1);
    }

    // POPULATION_SIZE, MAX_GENERATIONS, crossover_rate, mutate_rate, stop_criteria

    // ###################################################################################
    // you dont need to change anything here
    // the number of variables
    int NUM_VARIABLES = 2;
    // the lower bounds of variables
    double Lbound[] = {-5.0, -5.0};
    // the upper bounds of variable
    double Ubound[] = {5.0, 5.0};
    // ###################################################################################

    // <YOUR CODE: Here make all the initial print outs>

    printf("genetic algorithm initiated\n");
    printf("====================================\n");
    printf("number of varibles: %d\n", NUM_VARIABLES);
    for (int i = 0; i < 1; i++)
    {
        printf("lower bound: %.10lf , %.10lf\n", Lbound[0], Lbound[1]);
    }

    for (int i = 0; i < 1; i++)
    {
        printf("Upper bound: %.10lf , %.10lf\n", Ubound[0], Ubound[1]);
    }

    printf("\n");

    printf("POPULATION_SIZE : %d\n", POPULATION_SIZE);
    printf("MAX GENERATIONS : %d\n", MAX_GENERATIONS);
    printf("crossover rate: %lf\n", crossover_rate);
    printf("mutate_rate : %lf\n", mutate_rate);
    printf("stop criteria : %.20lf\n", stop_criteria);

    printf("results\n");
    printf("====================================\n");

    clock_t start_time, end_time;
    double cpu_time_used;
    start_time = clock();

    // <YOUR CODE: Declare all the arrays you need here>
    double population[POPULATION_SIZE][NUM_VARIABLES];
    double fitness[POPULATION_SIZE];
    double new_population[POPULATION_SIZE][NUM_VARIABLES];
    double best_fitness_values[MAX_GENERATIONS];
    double final_population[POPULATION_SIZE][NUM_VARIABLES];

    // <YOUR CODE: Call generate_population function to initialize the "population"> like:
    generate_population(POPULATION_SIZE, NUM_VARIABLES, population, Lbound, Ubound);

    
    double fitness_probs[POPULATION_SIZE];
    

    /*mutation rates, are adaptive thorughout depending on where solution is foib*/
    double initial_mutate_rate = mutate_rate; // these can chamge depending on scenatio
    double final_mutate_rate = 0.001;         //

    double initial_crossover_rate = crossover_rate;
    double final_crossover_rate = 0.6;

    /*issue loop only running once put in break points*/
    for (int generation = 0; generation < MAX_GENERATIONS; generation++)
    {

        //   <YOUR CODE: Compute the fitness values using objective function for
        //   each row in "population" (each set of variables)> like:
        //   compute_objective_function(POPULATION_SIZE, NUM_VARIABLES, population, fitness);

        // call functios

        compute_objective_function(POPULATION_SIZE, NUM_VARIABLES, population, fitness);
        mutate_rate = adaptiveMutationRate(generation, MAX_GENERATIONS, initial_mutate_rate, final_mutate_rate);
        crossover_rate = adaptiveCrossoverRate(generation, MAX_GENERATIONS, initial_crossover_rate, final_crossover_rate);

        /* Elitism - this section takes the genes o the fittest individual and pastes into the next generation */

        int elite_index = 0;                // find index of most elite individual
        double elite_fitness = __DBL_MAX__; // set to super high indivual - will change

        for (int i = 0; i < POPULATION_SIZE; i++)
        {
            if (fitness[i] < elite_fitness)
            {
                elite_fitness = fitness[i]; // find indi wih best fitness
                elite_index = i;
            }
        }

        // keep best indiviaial - store
        double elite_individual[NUM_VARIABLES];
        for (int i = 0; i < NUM_VARIABLES; i++)
        {
            elite_individual[i] = population[elite_index][i];
        }

        // <YOUR CODE: Here implement the logic of finding best solution with minimum fitness value
        // and the stopping criteria>

        /* code to calculate fintess probabilties affects how fitness is ppassed on*/

        for (int i = 0; i < POPULATION_SIZE; i++)
        {

            if (fitness[i] < best_fitness)
            {
                best_fitness = fitness[i]; // find best i=fintess
                best_index = i;
            }
        }
        double sum_fitness_probs = 0;
        for (int i = 0; i < POPULATION_SIZE; i++)
        {
            fitness_probs[i] = 1 / (fitness[i] + 1e-5); // find i=fitness probs
            sum_fitness_probs += fitness_probs[i];
        }
        for (int i = 0; i < POPULATION_SIZE; i++)
        {
            fitness_probs[i] /= sum_fitness_probs;
        }

        // Selection of parents based on fitness probabilities and creation of a new population
        for (int i = 0; i < POPULATION_SIZE; i++)
        {

            double random_value = generate_random(0, 1); // Generate a random value between 0 and 1
            double cumulative_prob = 0.0;

            int selected_index = -1; // Initialize with an invalid value

            for (int j = 0; j < POPULATION_SIZE; j++)
            {
                cumulative_prob += fitness_probs[j];
                if (random_value <= cumulative_prob)
                {
                    selected_index = j;
                    break; // if seected index found
                }
            }

            if (selected_index != -1)
            {
                // Copy the selected individual to the new population
                for (int k = 0; k < NUM_VARIABLES; k++)
                {
                    new_population[i][k] = population[selected_index][k];
                }
            }
        }

        // <YOUR CODE: Here call the crossover function>
        crossover(POPULATION_SIZE, NUM_VARIABLES, fitness, new_population, population, crossover_rate);

        // <YOUR CODE: Here call the mutation function>
        mutate(POPULATION_SIZE, NUM_VARIABLES, new_population, population, Lbound, Ubound, mutate_rate);
        // recompite from loop
        for (int i = 0; i < POPULATION_SIZE; i++)
        {
            for (int j = 0; j < NUM_VARIABLES; j++)
            {
                population[i][j] = new_population[i][j]; // copy new pop into pop (maybe unexxasy)
            }
        }

        for (int i = 0; i < NUM_VARIABLES; i++)
        {
            population[0][i] = elite_individual[i]; // replace first individial with most elite
        }

        previous_best_fitness = best_fitness;

        /*stop criteria*/

        double fitness_change = fabs(previous_best_fitness - best_fitness);
        if (best_fitness < global_best_fitness)
        {
            // There is an improvement
            global_best_fitness = best_fitness;
            previous_best_fitness = best_fitness; // Update previous best to current best
            no_improvement_counter = 0;           // Reset the counter since we found a better fitness
            // printf("found better solution\n");
        }
        else if (fitness_change < stop_criteria)
        {
            // No significant improvement
            no_improvement_counter++; // Increment the no improvement counter
            // printf("no improvment\n");
        }

        // Check if the no improvement threshold has been reached
        if (no_improvement_counter >= NO_IMPROVEMENT_THRESHOLD)
        {
            printf("Stopping: No signifcant improvment in best fitness for %d generations.\n", NO_IMPROVEMENT_THRESHOLD);
            break; // Exit the generation loop
        }

        // for plotting (done  by chatgpt)
        double generation_best_fitness = __DBL_MAX__;
        for (int i = 0; i < POPULATION_SIZE; i++)
        {
            if (fitness[i] < generation_best_fitness)
            {
                generation_best_fitness = fitness[i];
            }
        }
        if (generation == MAX_GENERATIONS - 1)
        {
            for (int i = 0; i < POPULATION_SIZE; i++)
            {
                for (int j = 0; j < NUM_VARIABLES; j++)
                {
                    final_population[i][j] = population[i][j];
                }
            }
        }

        // Update the array with the best fitness value for this generation
        best_fitness_values[generation] = generation_best_fitness;
        FILE *best_fitness_file = fopen("best_fitness_values.txt", "w");
        if (best_fitness_file)
        {
            for (int generation = 0; generation < MAX_GENERATIONS; generation++)
            {
                fprintf(best_fitness_file, "%f\n", best_fitness_values[generation]);
            }
            fclose(best_fitness_file);
        }
        // Open a file to save the final population
        FILE *final_population_file = fopen("final_population.txt", "w");
        if (final_population_file)
        {
            for (int i = 0; i < POPULATION_SIZE; i++)
            {
                for (int j = 0; j < NUM_VARIABLES; j++)
                {
                    fprintf(final_population_file, "%f ", population[i][j]);
                }
                fprintf(final_population_file, "\n");
            }
            fclose(final_population_file);
        }

        // Find the solution closest to (0, 0) using distance formual
        // Variables to hold the minimum distance and its index
        // Variables to hold the minimum values and their indices
        double min_x = fabs(final_population[0][0]);
        double min_y = fabs(final_population[0][1]);
        int min_index_x = 0;
        int min_index_y = 0;

        // Loop through the population to find the smallest absolute values
        for (int i = 0; i < POPULATION_SIZE; i++)
        {
            if (fabs(final_population[i][0]) < min_x)
            {
                min_x = fabs(final_population[i][0]);
                min_index_x = i;
            }
            if (fabs(final_population[i][1]) < min_y)
            {
                min_y = fabs(final_population[i][1]);
                min_index_y = i;
            }
        }

        printProgressBar(generation, MAX_GENERATIONS); // prints progress bar
    }



   

    // Find the solution closest to (0, 0) using distance formula
    // This should be after the final population has been completely assigned
    double closest_solution_distance = __DBL_MAX__;

    for (int i = 0; i < POPULATION_SIZE; i++)
    {
        double x = population[i][0];
        double y = population[i][1];
        double distance = calculateDistance(x, y);

        if (distance < closest_solution_distance)
        {
            closest_solution_distance = distance;
            closest_solution_index = i;
        }
    }

    

    // <YOUR CODE: Jump to this part of code if the stopping criteria is met before MAX_GENERATIONS is met>

    // ###################################################################################
    // You dont need to change anything here
    // Here we print the CPU time taken for your code

    printf("genetic algorithm complete\n");
    printf("====================================\n");
    printf("number of varibles: %d\n", NUM_VARIABLES);
    for (int i = 0; i < 1; i++)
    {
        printf("lower bound: %.10lf , %.10lf\n", Lbound[0], Lbound[1]);
    }

    for (int i = 0; i < 1; i++)
    {
        printf("Upper bound: %.10lf ,%.10lf\n", Ubound[0], Ubound[1]);
    }

    printf("\n");

    printf("POPULATION_SIZE : %d\n", POPULATION_SIZE);
    printf("MAX GENERATIONS : %d\n", MAX_GENERATIONS);
    printf("crossover rate: %lf\n", crossover_rate);
    printf("mutate_rate : %lf\n", mutate_rate);
    printf("stop criteria : %lf\n", stop_criteria);

    end_time = clock();
    cpu_time_used = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    printf("CPU time: %f seconds\n", cpu_time_used);
    // ###################################################################################

    // <Here print out the best solution and objective function value for the best solution like the format>

    printf("Best fitness value: %.10lf\n", best_fitness);
    if (closest_solution_index != -1)
    {
        printf("Best soltuon): (%.10lf, %.10lf)",
               population[closest_solution_index][0],
               population[closest_solution_index][1]);
    }

    printf("\n");

    

    return 0;
}

// chat gpt did this
void printProgressBar(int iteration, int total)
{
    int progressBarWidth = 50;
    float progress = (float)iteration / total;
    int numBarChars = (int)(progress * progressBarWidth);

    printf("Progress: [");
    for (int i = 0; i < numBarChars; i++)
    {
        putchar('=');
    }
    for (int i = numBarChars; i < progressBarWidth; i++)
    {
        putchar(' ');
    }
    printf("] %.1f%%\r", progress * 100);
    fflush(stdout);
}

double calculateDistance(double x, double y)
{
    return sqrt(x * x + y * y); // disrnace formula
}
double adaptiveMutationRate(int generation, int MAX_GENERATIONS, double startRate, double endRate)
{
    // Linearly decrease the mutation rate from startRate to endRate over the generations
    double rate = startRate + (endRate - startRate) * ((double)generation / (double)MAX_GENERATIONS);
    return rate < endRate ? endRate : rate; // Ensure that the rate never goes below the endRate
}

double adaptiveCrossoverRate(int generation, int MAX_GENERATIONS, double startRate, double endRate)
{
    // Linearly decrease the crossover rate from startRate to endRate over the generations
    double rate = startRate + (endRate - startRate) * ((double)generation / (double)MAX_GENERATIONS);
    return rate < endRate ? endRate : rate; // Ensure that the rate never goes below the endRate
}

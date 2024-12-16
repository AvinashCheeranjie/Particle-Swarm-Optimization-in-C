#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "utility.h"
#include "OF_lib.h"

// Helper function to generate random numbers in a range
double random_double(double min, double max) {
    return min + (max - min) * ((double)rand() / RAND_MAX);
}

// Calculates inertia weight as a decreasing linear function of the iteration number
double calc_inertia(int iter, int MAX_ITERATIONS, double w_min, double w_max) {
    int dec_stage = 0.75*MAX_ITERATIONS;
    if (iter <= dec_stage)
        return w_min + (w_max - w_min) * (dec_stage - iter) / dec_stage;
    else
        return w_min;
}

// Box Muller transform to generate Gaussian white noise
double gauss(double mean, double std_dev) {
    double u1 = random_double(0.0, 1.0);
    double u2 = random_double(0.0, 1.0);
    double z0 = sqrt(-2.0*log(u1))*cos(2*M_PI*u2);
    return z0 * std_dev + mean;
}

double pso(ObjectiveFunction objective_function, int NUM_VARIABLES, Bound *bounds, int NUM_PARTICLES, int MAX_ITERATIONS, double *best_position, char *objective_function_name) {
    // Algorithm Parameters
    double w = 0.729; // Inertia weight for higher dimensions
    double w_min = 0.4; 
    double w_max = 0.9;
    double c1 = 1.49; // Cognitive coefficient
    double c2 = 1.50; // Social coefficient
    double c3 = 0.2; // Adaptive coefficient
    double epsilon = 1e-17; // Tolerance 
    double sig_max =1.0;
    double sig_min = 0.1;

    // Particle Parameters
    double **x = (double **)malloc(NUM_PARTICLES * sizeof(double *));                // Positions
    double **v = (double **)malloc(NUM_PARTICLES * sizeof(double *));                // Velocities
    double **p = (double **)malloc(NUM_PARTICLES * sizeof(double *));                // Personal best positions
    double *personal_best_scores = (double *)malloc(NUM_PARTICLES * sizeof(double)); // f_p
    double *g = (double *)malloc(NUM_VARIABLES * sizeof(double));                    // Global best positions

    // Swarm Initialization
    double global_best_score = INFINITY; 

    for (int i = 0; i < NUM_PARTICLES; i++) {
        x[i] = (double *)malloc(NUM_VARIABLES * sizeof(double));
        v[i] = (double *)malloc(NUM_VARIABLES * sizeof(double));
        p[i] = (double *)malloc(NUM_VARIABLES * sizeof(double));

        for (int j = 0; j < NUM_VARIABLES; j++) {
            x[i][j] = random_double(bounds[j].lowerBound, bounds[j].upperBound);
            double v_max = 0.5*fabs(bounds[j].upperBound - bounds[j].lowerBound);
            v[i][j] = random_double(-v_max, v_max); 
            p[i][j] = x[i][j];
        }
        personal_best_scores[i] = objective_function(NUM_VARIABLES, x[i]);
        if (personal_best_scores[i] < global_best_score) {
            global_best_score = personal_best_scores[i];
            memcpy(g, p[i], sizeof(double) * NUM_VARIABLES);
        }
    }

    // PSO loop
    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
        if(iter%1000 == 0 && MAX_ITERATIONS > 1000) printf("Iteration: %d, Best Score: %0.16lf\n", iter, global_best_score);
        if(iter%10 == 0 && MAX_ITERATIONS < 1000) printf("Iteration: %d, Best Score: %0.16lf\n", iter, global_best_score);

        double Pc = c3*(1+cos(iter*M_PI/MAX_ITERATIONS)); // Probability function for Gaussian pertubation

        // Dynamically adjust inertia weight for lower dimensions (smaller search space)
        if(NUM_VARIABLES <= 50 || strcmp(objective_function_name, "rastrigin")==0) {
            w = calc_inertia(iter, MAX_ITERATIONS, w_min, w_max);
        } 
            
        // Tolerance check
        if(fabs(global_best_score) < epsilon){
            printf("Optimal Fitness Found. Number of Iterations Executed: %d\n\n", iter);
            break;
        }

        for (int i = 0; i < NUM_PARTICLES; i++) {            
            for (int j = 0; j < NUM_VARIABLES; j++) {
                double r1 = random_double(0.0, 1.0);
                double r2 = random_double(0.0, 1.0);
                
                // Randomly apply Gaussian pertubation closer to convergence
                if(iter > 0.7*MAX_ITERATIONS && random_double(0.0,1.0) <= Pc && NUM_VARIABLES <= 50){
                    double r3 = random_double(0.0, 1.0);
                    double sigma = sig_max-(sig_max-sig_min)*iter/MAX_ITERATIONS;
                    x[i][j] *= r3*gauss(0, sigma*sigma);   
                    x[i][j] += v[i][j];   
                }
                
                if(iter > 0.5*MAX_ITERATIONS && (random_double(0.0,1.0) <= Pc) && (NUM_VARIABLES > 50)){
                    double r3 = random_double(0.0, 1.0 && NUM_VARIABLES > 50);
                    double sigma = sig_max-(sig_max-sig_min)*iter/MAX_ITERATIONS;
                    x[i][j] *= r3*gauss(0, sigma*sigma);   
                    x[i][j] += v[i][j];   
                } 

                // Update velocity & position
                v[i][j] = w*v[i][j] + c1*r1*(p[i][j] - x[i][j]) + c2*r2*(g[j] - x[i][j]); 
                x[i][j] += v[i][j]; 

                // Clip position within bounds
                if (x[i][j] < bounds[j].lowerBound){ 
                    x[i][j] = bounds[j].lowerBound;
                    v[i][j] = 0.0;
                }
                if (x[i][j] > bounds[j].upperBound){
                    x[i][j] = bounds[j].upperBound;
                    v[i][j] = 0.0;
                }
            }

            double score = objective_function(NUM_VARIABLES, x[i]);

            if (score < personal_best_scores[i]) {
                personal_best_scores[i] = score;
                memcpy(p[i], x[i], sizeof(double) * NUM_VARIABLES);
            }

            if (score < global_best_score) {
                global_best_score = score;
                memcpy(g, x[i], sizeof(double) * NUM_VARIABLES);
            }  
                  
        }
        
    }

    memcpy(best_position, g, sizeof(double) * NUM_VARIABLES);
    
    // Free all allocated memory
    for (int i = 0; i < NUM_PARTICLES; i++) {
        free(x[i]);
        free(v[i]);
        free(p[i]);
    }
    free(x);
    free(v);
    free(p);
    free(personal_best_scores);
    free(g);

    return global_best_score;
}


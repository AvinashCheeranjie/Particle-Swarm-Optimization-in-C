# Particle Swarm Optimization (PSO) in C

[Particle Swarm Optimization](https://en.wikipedia.org/wiki/Particle_swarm_optimization) (PSO) is a population-based optimization algorithm inspired by the social behavior of birds flocking or fish schooling. Developed by Kennedy and Eberhart in 1995, PSO mimics how individuals in a group adapt by sharing information to reach an optimal solution in a multi-dimensional search space. In PSO, a ”swarm” of particles represents potential solutions to the optimization problem. Each particle adjusts its position based on its own best experience and the best experience of the swarm, gradually converging to an optimal or near optimal solution. This algorithm is widely used for solving complex optimization problems due to its simplicity, efficiency, and ability to avoid local optima in many cases. In this project, I developed a very basic PSO algorithm and tested it on the following benchmark objective functions:

1. [Griewank](https://www.sfu.ca/%7Essurjano/griewank.html)
2. [Levy](https://www.sfu.ca/%7Essurjano/levy.html)
3. [Rastrigin](https://www.sfu.ca/%7Essurjano/rastr.html)
4. [Rosenbrock](https://www.sfu.ca/%7Essurjano/rosen.html)
5. [Schwefel](https://www.sfu.ca/%7Essurjano/schwef.html)
6. [Dixon-Price](https://www.sfu.ca/%7Essurjano/dixonpr.html)
7. [Michalewicz](https://www.sfu.ca/%7Essurjano/michal.html) with $m=10$
8. [Styblinski-Tang](https://www.sfu.ca/%7Essurjano/stybtang.html)

## Features
- Adaptive inertia weight to encourage exploitation closer to optimal convergence
- Random application of Gaussian perturbation to prevent premature convergence to local minima 


## Usage 
- Clone the repository
- Run `make pso` in the terminal to compile and build the program
- Execute the program with a command line of the following format:
  ```bash
  ./pso <objective function> <number of dimensions> <lower bound> <upper bound> <number of particles> <number of iterations>
  ```
  For example:
  ```bash
  ./pso griewank 10 -600 600 450 100000
  ```


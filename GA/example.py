import pygad
import numpy

def fitness_func(solution, solution_idx):
    output = numpy.sum(solution*function_inputs)
    fitness = 1.0 / numpy.abs(output - desired_output)
    return fitness

function_inputs = [4,-2,3.5,5,-11,-4.7]
desired_output = 44

fitness_function = fitness_func
num_genes = len(function_inputs)


ga_instance = pygad.GA(num_generations=50,
                       num_parents_mating=4,
                       fitness_func=fitness_function,
                       sol_per_pop=8,
                       num_genes=num_genes,
                       init_range_low=-2,
                       init_range_high=5,
                       parent_selection_type="sss",
                       keep_parents=1,
                       crossover_type="single_point",
                       mutation_type="random",
                       mutation_percent_genes=10)

ga_instance.run()

solution, solution_fitness, solution_idx = ga_instance.best_solution()
print("\nParameters of the best solution : {solution}".format(solution=solution))
print("\nFitness value of the best solution = {solution_fitness}".format(solution_fitness=solution_fitness))

prediction = numpy.sum(numpy.array(function_inputs)*solution)
print("\nPredicted output based on the best solution : {prediction}\n".format(prediction=prediction))

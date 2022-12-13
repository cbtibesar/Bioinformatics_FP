"""
This file provides the different functions for plotting the distances by generation for the different evolutionary
model simulations
"""

import numpy as np
from matplotlib import pyplot as plt
import math

"""
Function for plotting a single scatter plot of the genetic distance by generation
"""
def plot_data(distance_by_generation, graph_title, color, generations_per_item):
    ## create a numpy array from the genetic distance by generation values and the number of generations

    y = np.array(distance_by_generation)
    x = np.array(range(len(distance_by_generation)))

    ## set the plot title and labels
    plt.title(graph_title)
    plt.xlabel("Generation (x 10^" + str(math.log(generations_per_item, 10)) + ")")
    plt.ylabel("Genetic Distance")

    ## plot the graph and show
    plt.scatter(x, y, s=5, color=color)
    plt.show()


"""
Function for plotting an overlay of all scatter plots of the genetic distance by generation from the different models
"""
def plot_all_data(distance_by_generation_JC, distance_by_generation_K2P, distance_by_generation_HKY85, distance_by_generation_GTR, graph_title, generations_per_item):
    ## create numpy arrays for graphing like when plotting a single graph
    y_JC = np.array(distance_by_generation_JC)
    x_JC = np.array(range(len(distance_by_generation_JC)))

    y_K2P = np.array(distance_by_generation_K2P)
    x_K2P = np.array(range(len(distance_by_generation_K2P)))

    y_HKY85 = np.array(distance_by_generation_HKY85)
    x_HKY85 = np.array(range(len(distance_by_generation_HKY85)))

    y_GTR = np.array(distance_by_generation_GTR)
    x_GTR = np.array(range(len(distance_by_generation_GTR)))

    ## set the title and labels
    plt.title(graph_title)
    plt.xlabel("Generation (x 10^" + str(math.log(generations_per_item, 10)) + ")")
    plt.ylabel("Genetic Distance")

    ## plot each graph over each other with different colors and labels
    plt.scatter(x_JC, y_JC, s=5, color="blue", label="JC")
    plt.scatter(x_K2P, y_K2P, s=5, color="green", label="K2P")
    plt.scatter(x_HKY85, y_HKY85, s=5, color="red", label="HKY85")
    plt.scatter(x_GTR, y_GTR, s=5, color="orange", label="GTR")
    ## show the legend in the lower right corner of the graph
    plt.legend(loc="lower right")
    ## show the graph
    plt.show()

"""
This file provides the different functions for plotting the distances by generation for the different evolutionary
model simulations
"""

import numpy as np
from matplotlib import pyplot as plt

def plot_data(distance_by_generation, graph_title):
    y = np.array(distance_by_generation)
    x = np.array(range(len(distance_by_generation)))
    plt.title(graph_title)
    plt.xlabel("Generation")
    plt.ylabel("Distance")
    plt.scatter(x, y)
    plt.show()


"""
Initialize the problem.
"""

include("Problem.jl")

# -------   Problem Parameters   -------
T = 0.1

problem = Problem.Gaussian(T)

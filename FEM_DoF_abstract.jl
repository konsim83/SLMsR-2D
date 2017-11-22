"""

Define an abstrcat type Dof that serves as a parent type for type
safety. This can be improved. Note that the syntax changes to Julia
v0.6.

"""
abstract Dof

abstract Dof_square <: Dof

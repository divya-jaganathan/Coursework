#!/usr/bin/env python
# -*- coding: utf-8 -*-
# CFT - 25 Feb, 2019: Solved Elastostatics for Beam with a hole

import fenics as fn
import mshr as ms
import matplotlib.pyplot as plt


# Material properties

lmbda=10.0
mu=20.0

def eps(u):
	return fn.sym(fn.grad(u))

def sigma(u):
	return (lmbda*fn.tr(eps(u))*fn.Identity(2)+2*mu*eps(u))

# Geometry Specification/Mesh Generation

L=10.0
H=2.0

bar=ms.Rectangle(fn.Point(0,0),fn.Point(L,H))
hole=ms.Circle(fn.Point(L/2,H/2),H/3)
domain=bar-hole
mesh=ms.generate_mesh(domain,30)

VFS=fn.VectorFunctionSpace(mesh,'P',2)

tol=10**-14

# Boundary conditions

def left(x, on_boundary):
	return on_boundary and x[0]<tol

bc=fn.DirichletBC(VFS, fn.Constant((0,0)),left)
fext=fn.Expression(('0','-0.1'),degree=1)

# Weak formulation and solving it

u=fn.TrialFunction(VFS)
v=fn.TestFunction(VFS)

lhs=fn.inner(eps(v),sigma(u))*fn.dx
rhs=fn.inner(v,fext)*fn.dx

u=fn.Function(VFS)
fn.solve(lhs==rhs, u,bc)

# Plot and save result

fn.plot(u,mode='displacement')
plt.savefig('Exercise 3: Beam_with_hole')


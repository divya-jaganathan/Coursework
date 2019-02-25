#!/usr/bin/env python
# -*- coding: utf-8 -*-
# CFT - 21 Feb, 2019: Solved Elastostatics for Beam

import fenics as fn
import mshr as ms
import matplotlib.pyplot as plt

lmbda=1.0
mu=2.0

def eps(u):
	return fn.sym(fn.grad(u))

def sigma(u):
	return (lmbda*fn.tr(eps(u))*fn.Identity(2)+2*mu*eps(u))

L=10.0
H=2.0

domain=ms.Rectangle(fn.Point(0,0),fn.Point(L,H))
mesh=ms.generate_mesh(domain,10)

VFS=fn.VectorFunctionSpace(mesh,'P',2)

tol=10**-14

# Boundary conditions

def left(x, on_boundary):
	return on_boundary and x[0]<tol

bc=fn.DirichletBC(VFS, fn.Constant((0,0)),left)
fext=fn.Expression(('0','-1.0'),degree=1)

u=fn.TrialFunction(VFS)
v=fn.TestFunction(VFS)

lhs=fn.inner(eps(v),sigma(u))*fn.dx
rhs=fn.inner(v,fext)*fn.dx

u=fn.Function(VFS)
fn.solve(lhs==rhs, u,bc)

fn.plot(u)
plt.savefig('Exercise 2: Beam')


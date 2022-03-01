# Molecular Dynamics C++ code

   A LJ fluid solved with Langevin dynamics using a velocity Verlet integrator.  
   Interactions are handled via a head-and-list neighbour list.  

###  USAGE:

Compile using ```$ make```. This will create the file "md".  

A file called ```read.in``` must be present with parameters:  
```bash
   [numberParticles]
   [Lx] [Ly] [Lz]
   [cutOff]
```
A file called pos.init must be present containing [numberParticle] positions  

```bash
   [X] [Y] [Z]
   ...
```
Writes positions to pos.out and velocities to vel.out   

/* Raul P. Pelaez 2021-2022
   A LJ fluid solved with Langevin dynamics using a velocity Verlet integrator.
   Interactions are handled via a head-and-list netihgbour list.

   USAGE:
   A file called read.in must be present with parameters:
   [numberParticles]
   [Lx] [Ly] [Lz]
   [cutOff]
   
   A file called pos.init must be present containing [numberParticle] positions
   [X] [Y] [Z]
   ...

   Writes positions to pos.out and velocities to vel.out   
 */

#include"utils/vector.h"
#include"utils/helper.h"
#include"utils/Timer.cpp"

#include<iostream>
#include<numeric>
#include<vector>
#include<string>
#include<fstream>
#include<random>
#include<algorithm>

using namespace utils;

struct Parameters{
  real3 boxSize;
  int numberParticles;
  real cutOff = 2.5;
  real temperature = 0.1;
  real friction = 1.0;
  real dt = 0.01;

  std::string outFile = "pos.out";
  std::string outFileVel = "vel.out";
  
  int numberSteps = 5000;
  int printSteps =  100;
  int relaxSteps = 0;
};

struct HeadAndList{
  
  std::vector<int> head, list;
  int3 cells;
  real3 boxSize;

  void initialize(int numberParticles, real3 i_boxSize, real cutOff){
    this->boxSize = i_boxSize;
    this->cells = make_int3(boxSize/cutOff);
    cells.x = std::max(cells.x, 3);
    cells.y = std::max(cells.y, 3);
    cells.z = std::max(cells.z, 3);
    head.resize(cells.x*cells.y*cells.z);
    list.resize(numberParticles+1);
    std::fill(head.begin(), head.end(), 0);
 } 

  void updateList(std::vector<real3> &pos, real3 i_boxSize, real cutOff){
    initialize(pos.size(), i_boxSize, cutOff);
    for(int i = 0; i<pos.size(); i++){
      int icell = getCellIndex(getCell(pos[i]));
      list[i+1] = head[icell];
      head[icell] = i+1;
    }
  }
  
  int3 getCell(real3 pos){
    auto posInBox = apply_mic(pos, this->boxSize);
    return make_int3((posInBox/this->boxSize + 0.5)*make_real3(this->cells));
  }

  int getCellIndex(int3 ic){
    return ic.x + (ic.y + ic.z*cells.y)*cells.x;
  }

  int getNeighbourCellIndex(int ic, int3 cell_i){
    const int3 cell_j = pbc_cell(cell_i + indexToCell3D(ic, cells), cells);
    const int jcell = getCellIndex(cell_j);
    return jcell;
  }
  
};
 
struct Particles{
  std::vector<real3> pos, forces, velocities;
};

Parameters readParameters(std::string fileName){
  std::ifstream in(fileName);
  Parameters par;
  in>>par.numberParticles;
  in>>par.boxSize.x>>par.boxSize.y>>par.boxSize.z;
  in>>par.cutOff;
  return par;
}

Particles initializeParticles(Parameters par, std::string fileName){
  Particles particles;
  particles.pos.resize(par.numberParticles);
  particles.forces.resize(par.numberParticles, real3());
  particles.velocities.resize(par.numberParticles, real3());
  std::ifstream in(fileName);
  for(int i=0 ; i< par.numberParticles; i++){
    in>>particles.pos[i].x>>particles.pos[i].y>>particles.pos[i].z;
  }
  return particles;
}

real3 computeLJForce(real3 pi, real3 pj, const Parameters &par){
  real3 rij = apply_mic(pj - pi, par.boxSize);
  real r2 = dot(rij, rij);
  if(r2>0 and r2<par.cutOff*par.cutOff){
    real invr2 = real(1.0)/r2;
    real invr6 = invr2*invr2*invr2;  
    real fmoddivr = (real(-48.0)*invr6 + real(24.0))*invr6*invr2;  
    return fmoddivr*rij;
  }
  else return real3();
}

void computeShortRangedForces(Particles &particles, const Parameters &par){
  static HeadAndList nl;
  nl.updateList(particles.pos, par.boxSize, par.cutOff);
  constexpr int numberNeighbourCells = 27;
  for(int i = 0; i<par.numberParticles; i++){
    real3 fi = real3();
    const real3 pi = particles.pos[i];
    const int3 cell_i = nl.getCell(particles.pos[i]);
    for(int ic = 0; ic<numberNeighbourCells; ic++){
      int jcell = nl.getNeighbourCellIndex(ic, cell_i);
      int j = nl.head[jcell];
      while(j!=0){
	if(i != j-1){
	  fi += computeLJForce(pi, particles.pos[j-1], par);
	}
	j = nl.list[j];
      }
    }
    particles.forces[i] += fi;
  }
}

void computeForces(Particles &particles, const Parameters &par){
  std::fill(particles.forces.begin(), particles.forces.end(), real3());
  computeShortRangedForces(particles, par);
}

void firstVerletUpdate(Particles & particles, const Parameters &par, std::mt19937 &gen){
  std::normal_distribution<real> dis(0.0, 1.0);
  real noiseAmplitude = sqrt(par.temperature*par.friction*par.dt);  
  for(int i = 0; i< par.numberParticles; i++){
    const real3 fi = particles.forces[i];
    real3 vi = particles.velocities[i];
    const real3 noise = make_real3(dis(gen), dis(gen), dis(gen))*noiseAmplitude;
    vi += (fi - par.friction*vi)*par.dt*0.5 + noise;
    particles.pos[i] += vi*par.dt*0.5;
    particles.velocities[i] = vi;    
  }
}

void secondVerletUpdate(Particles & particles, const Parameters &par, std::mt19937 &gen){
  std::normal_distribution<real> dis(0.0, 1.0);
  real noiseAmplitude = sqrt(par.temperature*par.friction*par.dt);
  for(int i = 0; i< par.numberParticles; i++){
    const real3 fi = particles.forces[i];
    real3 vi = particles.velocities[i];
    const real3 noise = make_real3(dis(gen), dis(gen), dis(gen))*noiseAmplitude;    
    vi = (vi + fi*par.dt*0.5 + noise)/(1+par.friction*par.dt*0.5);
    particles.velocities[i] = vi;    
  }
}

real computeKineticEnergy(std::vector<real3> &vel){
  //This runs in parallel if execution is uncommented
  real K = 0.5*std::transform_reduce(//std::execution::par_unseq,
				     vel.begin(), vel.end(), 0.0,
				     std::plus<real>(),
				     [](real3 v){return dot(v,v);});  
  return K;
}

void writeSimulationState(Particles &particles, Parameters par){
  static std::ofstream out(par.outFile);
  static std::ofstream vout(par.outFileVel);
  out<<"#\n";
  vout<<"#\n";
  real kineticEnergy = computeKineticEnergy(particles.velocities);
  std::cout<<"Instantaneous temperature: "<<kineticEnergy/(par.numberParticles*1.5)<<std::endl;
  for(int i = 0; i<par.numberParticles; i++){
    auto p = particles.pos[i];
    out<<p.x<<" "<<p.y<<" "<<p.z<<"\n";
    auto v = particles.velocities[i];
    vout<<v.x<<" "<<v.y<<" "<<v.z<<"\n";
  }
  out<<std::flush;
  vout<<std::flush;
}

void doStep(Particles &particles, Parameters par){
  static std::mt19937 rng(std::random_device{}());
  firstVerletUpdate(particles, par, rng);
  computeForces(particles, par);
  secondVerletUpdate(particles, par, rng);
}

int main(int argc, char* argv[]){
  auto par = readParameters("read.in");
  auto particles = initializeParticles(par, "pos.init");
  computeForces(particles, par);
  Timer tim;
  tim.tic();
  for(int i = 0; i<par.numberSteps + par.relaxSteps; i++){
    doStep(particles, par);
    if(i >= par.relaxSteps and par.printSteps and i%par.printSteps == 0){
      writeSimulationState(particles, par);
    }
  }
  auto elapsed = tim.toc();
  std::cerr<<"Time per step: "<<elapsed/par.numberSteps<<std::endl;
  return 0;
}

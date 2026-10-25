#include "kernel.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

/*******time parameteres******/
//dt set by CFL: dt <= 0.4*h/cvel, h_water~0.0736, cvel=80 -> dt<=0.00037
const int TS = 10000;  //number of time steps
const double dt = 0.00005; //0.0003
const int frame_every = 1000; //write a frame (NNN.dat) every this many steps

/*****cvel warm-up ramp*******/
//a freshly-placed (random) particle field has large local density variance;
//starting stiff (target cvel) reacts to that variance with huge pressure spikes.
//Ramping cvel up from a soft start gives the field time to settle first.
const double cvel_ramp_start = 10.0;
const int cvel_ramp_steps = 3000;

/*****water fluid paramenters******/
const double rhowater = 998.29; // kg/m^3
const double Volwater= 0.5; //m^3
//Monaghan stiffness rule: cvel >= 10*Vmax; Vmax~6.8 m/s from the drop impact -> cvel>=68
const double cvelwater = 80.0; //speed of sound of the fluid
const double Gammawater = 7.0;  //adiabatic constant
const double muwater = 3.5; //viscosity
const double buoyancywater = 0.0;
const double sigmawater = 0.0728; //surface tension constant
const double thresholdwater = 6.768;//was 7.065 for poly6; see the note on kernel_vol_const
//buoyancy must stay 0: the explicit (rho-rhozero)*(-g) term double-counts the
//buoyancy that the pressure gradient already produces from the rhozero
//difference, and because it is divided by density it makes the *effective*
//gravity -9.82*(2 - rhozero/rho). Any particle below half rest density then
//accelerates upward, which is a runaway (dilute -> push up -> more dilute).

/*****petroleum (crude oil) fluid parameters******/
//medium crude, ~31 API at 20°C. Replaces the old light-oil phase, which ran at an
//artificial 400 kg/m^3; at 870 the density contrast with water is only ~13%, so
//the phases separate much more slowly than before.
const double rhopetroleum = 870.0;   // kg/m^3
const double Volpetroleum = 0.5;     // m^3
const double cvelpetroleum = 80.0;   // speed of sound of the fluid, matched to cvelwater for consistency
const double Gammapetroleum = 7.0;   // adiabatic constant
//same scale as muwater (3.5 <-> water's 1 mPa·s): medium crude is ~10 mPa·s
const double mupetroleum = 35.0;     // viscosity
const double buoyancypetroleum = 0.0;  //see buoyancywater
const double sigmapetroleum = 0.027; // surface tension constant, N/m (crude: ~0.025-0.030)
const double thresholdpetroleum = 5.748;//was 6.0 for poly6; see the note on kernel_vol_const

/*****petroleum-water interface parameters******/
//repulsive force between unlike phases, keeps them from interpenetrating/mixing;
//real crude-water interfacial tension is ~0.02-0.03 N/m
const double sigma_interphase = 0.05;  // try 0.02 N/m

/*****mucus fluid paramenters******/
const double rhomucus = 1500.00; // kg/m^3
const double Volmucus = 0.15; //m^3 
const double mumucus = 36.0; //viscosity
const double buoyancymucus = 0.0;  
const double sigmamucus = 6.0; //surface tension constant
const double thresholdmucus = 5.0;

/******Smoothing length / particle spacing ratio*******/
//h = (kernel_vol_const*315*Vol/(64*pi*N))^(1/3), i.e. h/s = (c*315/(64*pi))^(1/3)
//with s the rest particle spacing. 5.1063 gives h/s = 2.0, about 33 neighbours.
//h is the SUPPORT RADIUS, so in the cubic spline's usual 2H convention this is
//H/s = 1.0 -- comfortably below the H/s ~ 1.5 where the spline starts to suffer
//the pairing instability.
//At this ratio the cubic-spline density sum on a rest-spaced lattice returns
//0.99997*rhozero, i.e. the reference state sits essentially exactly at rho0 and
//the EOS restores symmetrically in both directions. (poly6 returned 1.010 here
//and 0.963 at the original 3.05761318, where the 3.7% deficit drove the rest
//pressure negative and it was then clamped to zero, leaving the fluid with no
//restoring stiffness at its own equilibrium.)
//The surface-tension thresholds are tied to this choice: the colour-field
//normal at a flat free surface is 0.958x its poly6 value, so thresholdwater and
//thresholdpetroleum were rescaled by that factor to keep surface detection identical.
const double kernel_vol_const = 5.1063;

/******Number of particles *********/
const int N1 = 5500;  //number of particles
const int N2 = 5500;

/******Water drop released above the equilibrium pool*******/
//Shape of the drop: "sphere", "cube", or "none". Use "none" for the settling
//pass that generates 000.data: the file reader expects exactly N1+N2 rows, so
//the run that writes it must not carry any drop particles.
//Otherwise the particles sit on a
//simple-cubic lattice at the water particle's own rest spacing
//s = (Volwater/N1)^(1/3), so the drop starts at the pool's density rather than
//artificially compressed or dilute, and the offsets are exactly symmetric under
//reflection through the drop centre in all three axes (no random jitter).
const char drop_shape[] = "none";
//Drop radius in lattice cells. The kernel support is 2 cells (h/s = 2), so a
//drop of radius 1 or 2 is essentially all surface: its density sum never
//reaches rhozero, the EOS pressure clamps to 0, and it disperses in flight
//instead of holding together as a drop.
//Resulting Ndrop -- sphere: 7, 33, 123, 257 for R = 1, 2, 3, 4
//                   cube (2R+1)^3: 27, 125, 343, 729
const int drop_radius_cells = 3;
const double dropX = 0.0, dropY = 0.0, dropZ = 3.0;

struct Offset { double x, y, z; };

//Lattice offsets for the drop, in units of the water particle spacing. Ndrop is
//taken from this same function, so the particle count can never disagree with
//the placement loop.
vector<Offset> drop_offsets(){
  vector<Offset> off;
  if( strcmp(drop_shape, "none") == 0 ) return off;//no drop: Ndrop == 0
  const int R = drop_radius_cells;
  const bool sphere = ( strcmp(drop_shape, "sphere") == 0 );
  for(int ix=-R; ix<=R; ix++)
    for(int iy=-R; iy<=R; iy++)
      for(int iz=-R; iz<=R; iz++){
        if( sphere && (ix*ix + iy*iy + iz*iz > R*R) ) continue;
        Offset o;
        o.x = ix;  o.y = iy;  o.z = iz;
        off.push_back(o);
      }
  return off;
}

const int Ndrop = (int)drop_offsets().size();

const int N = N1+N2+Ndrop;
const double VolT = Volwater + Volpetroleum;
//smoothing length of n particles filling volume vol (see the note on kernel_vol_const)
inline double smoothing_length(double vol, int n){
  return cbrt(315*vol*kernel_vol_const/(64*M_PI*n));
}
const double h_water_ref = smoothing_length(Volwater, N1); //water smoothing length, reused as a lengthscale outside the Particle class
const double s_water = cbrt(Volwater/N1); //rest spacing of the water particles
const double box_dim = cbrt(VolT/4.0); //box half-width, sized so that (2*box_dim)*(2*box_dim)*box_dim == VolT

/******Boundary repulsion parameters******/
//Smooth push-back as a particle nears a wall. Range is HALF a smoothing length,
//i.e. about one particle spacing: at a full h (two spacings) the wall holds the
//near-wall layers ~1.5 spacings clear of the floor, which opens a gap and lifts
//the whole column. D must still be large enough to carry the hydrostatic load
//of the column: rho*g*depth spread over one particle's footprint s^2 needs
//~140 m/s^2 for a 0.63 m pool, and D = 500 puts the bottom layer's standoff at
//about half a spacing. Fast particles that beat it are caught by the
//(dissipative) clamp in useboundaries().
const double boundary_r0 = 0.5*h_water_ref;
const double boundary_D  = 500.0;
double wall_repulsion(double d, double r0, double D); //defined near useboundaries()

/******Periodic boundary in x******/
//When true, the x walls are removed: a particle leaving through x = x0+Lx
//re-enters at x0 and vice versa, and every pair separation uses the minimum
//image in x. Requires Lx >= 3*h so the neighbour grid has at least 3 distinct
//cells across x (otherwise the 3x3x3 stencil would visit a cell twice).
const bool periodic_x = true;
double box_x0 = 0.0, box_Lx = 0.0; //set by boundaries::set_boundaries

//map x back into [box_x0, box_x0+box_Lx)
inline double wrap_x(double x){
  if(!periodic_x) return x;
  x -= box_Lx*floor( (x-box_x0)/box_Lx );
  if(x >= box_x0+box_Lx) x = box_x0;//floor rounding can land exactly on the upper edge
  return x;
}

//minimum-image x separation: the nearest periodic copy of the pair
inline double min_image_x(double dx){
  if(!periodic_x) return dx;
  return dx - box_Lx*round(dx/box_Lx);
}

#include "rocks.h" //the rocks: needs the wall and periodic definitions above

//caps a single surface-tension contribution's acceleration magnitude: the CSF
//normal estimate is noisy at a sparsely-sampled free surface, and the kernel
//Laplacian is largest at small separations, so an occasional pair can otherwise
//"pop" a particle with an outsized, spurious kick
const double surface_tension_amax = 50.0; //m/s^2

/******Relaxation damping******/
//Artificial velocity damping, a = -friction*v, in 1/s. This is a RELAXATION
//knob for settling the field to equilibrium, not a physical fluid property --
//turn it down or to 0 once you want real dynamics, since at 20/s it suppresses
//genuine fluid motion along with the numerical noise.
//The kernel inconsistency that used to leak energy in is gone: density and the
//pressure gradient now both come from kernel_cubic, so the pressure force is
//the exact gradient of the discrete thermal energy and the scheme conserves
//energy up to integrator error. h is uniform across both phases, so no grad-h
//(Omega) correction terms are needed for that to hold.
//One non-conservative term remains -- clamping negative pressure to 0 rectifies
//the EOS, pushing particles apart under compression without pulling them back
//under expansion. Prefer Monaghan artificial stress if that shows up.
//With the leak fixed, friction should no longer be needed to hold the field
//steady; drop it toward 0 and let av_alpha do the work.
const double friction = 0.0;

/******Monaghan artificial viscosity******/
//Pi_ij, the standard WCSPH dissipation term. It enters the momentum equation in
//exactly the same place as the pressure term, so it is simply added into
//(prhoi+prhoj) below; Pi_ij is symmetric in i,j, so the pair force stays
//antisymmetric and momentum is conserved exactly.
//
//The difference from the global `friction` above is what makes it useful here:
//Pi_ij is nonzero ONLY for pairs that are APPROACHING each other. It therefore
//dissipates relative motion -- acoustic ringing, particle disorder, the noise
//that the non-conservative Mueller force model keeps feeding in -- while a
//uniformly translating blob (a falling drop, a sloshing pool) has no
//approaching pairs and is left essentially untouched. `friction` cannot make
//that distinction: it damps the drop and the noise at the same rate.
//
//  alpha : linear term, the one that does the damping. Free-surface WCSPH
//          typically runs alpha ~ 0.01-0.1. alpha = 1.0 is the classic shock
//          value and would be heavily over-damped here. Set 0 to disable.
//  beta  : quadratic term, suppresses particle interpenetration at high Mach.
//          The flow here is Mach ~ 6.5/80 = 0.08, so 0 is fine; use beta = 2*alpha
//          if fast particles start passing through each other.
//  eps   : softens the denominator so mu stays finite as r -> 0.
const double av_alpha = 0.05;
const double av_beta  = 0.0;
const double av_eps   = 0.01;

class boundaries;

class Particle {
private:
  double h, mass, rhozero, density, mu, sigma, threshold, buoyancy, cvel, Gamma; //fluid properties
  unsigned int phase;
  double normal[3], pressure; //calculated properties
  vector<int> neighbourlist;
public:
  double r[3];
  double v[3];
  double a[3];
  void setparams(int phases);
  void setphase(unsigned int ph);
  void setdensity(double rho){density = rho;};
  void ramp_cvel(double frac);
  void reinit();
  double kin_energy();
  friend void press_calc(Particle *body);
  friend void color_calc(Particle *body);
  friend void forces(Particle *body, boundaries& limit);
  friend void kinetic_E(Particle *body, int count);
};

void Particle::setparams(int phases){
  if( phases == 1){
    phase = 1;
    h = smoothing_length(VolT, N);
    rhozero = rhowater;
    mass = rhozero * VolT/N;
    mu = muwater;
    sigma = sigmawater;
    threshold = thresholdwater;
    buoyancy = buoyancywater;
    cvel = cvelwater;
    Gamma = Gammawater;
  }
  else{
    if( (double)rand()/(RAND_MAX) <= 0.5 )
      setphase(1);
    else
      setphase(2);
  }
}

void Particle::setphase(unsigned int ph){
  if(ph == 1){
    phase = 1;
    h = smoothing_length(Volwater, N1);
    rhozero = rhowater;
    mass = rhozero * Volwater/N1;
    mu = muwater;
    sigma = sigmawater;
    threshold = thresholdwater;
    buoyancy = buoyancywater;
    cvel = cvelwater;
    Gamma = Gammawater;
  }
  else{
    phase = 2;
    h = smoothing_length(Volpetroleum, N2);
    rhozero = rhopetroleum;
    mass = rhozero * Volpetroleum/N2;
    mu = mupetroleum;
    sigma = sigmapetroleum;
    threshold = thresholdpetroleum;
    buoyancy = buoyancypetroleum;
    cvel = cvelpetroleum;
    Gamma = Gammapetroleum;
  }
}

void Particle::ramp_cvel(double frac){
  double target = (phase == 1) ? cvelwater : cvelpetroleum;
  cvel = cvel_ramp_start + (target - cvel_ramp_start)*frac;
}

void Particle::reinit(){
  pressure = 0.0;
  fill( normal, normal+3, 0.0 );
  fill( a, a+3, 0.0 );
  neighbourlist.resize(0);
}

double Particle::kin_energy(){
  //mass-weighted, so the sum is a real energy in J and phases with different
  //particle masses are compared on the same footing
  return 0.5*mass*( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
}

class boundaries{
private:
  string figure;
  double boundaryX[2], boundaryY[2], boundaryZ[2], rad;
public:
  void set_boundaries(string fig, double dim, double zheight);
  friend void useboundaries(Particle& body, boundaries& limit);
  friend void forces(Particle *body, boundaries& limit);
  friend int main();
};

void boundaries::set_boundaries(string fig, double dim, double zheight){
  figure = fig;
  if(fig == "cube"){
    boundaryX[0] = -dim;
    boundaryX[1] = dim;
    boundaryY[0] = -dim;
    boundaryY[1] = dim;
    boundaryZ[0] = 0.0;
    boundaryZ[1] = zheight;
    box_x0 = boundaryX[0];
    box_Lx = boundaryX[1] - boundaryX[0];
  }
  if(fig == "sphere"){
    rad = dim;
  }
}

void initial_cond(string fig, Particle *body, int phases){
  
  if(fig == "cube"){
    //dim must match the box set by set_boundaries so particles are placed
    //inside the actual domain instead of outside its walls
    const double dim = box_dim;
    //Lattice, not uniform random. Random placement puts many pairs at
    //near-zero separation, where the pressure gradient is largest:
    //two nearly-coincident particles see a huge acceleration and the run starts by exploding
    //and never recovers. The lattice starts every particle at exactly rest
    //density; the small jitter breaks the lattice symmetry so it does not
    //behave like a crystal. The top layer is partial, which is the free surface.
    double spacing = cbrt(VolT/N);
    int nx = (int)(2.0*dim/spacing);
    int ny = nx;
    double jit = 0.05*spacing;
    for(int i=0; i<N; i++ ){
      body[i].setdensity(0.0);
      body[i].setparams(phases);
      body[i].reinit();
      int ix = i % nx;
      int iy = (i/nx) % ny;
      int iz = i/(nx*ny);
      body[i].r[0] = -dim + (ix+0.5)*spacing + jit*(2*((double)rand()/(RAND_MAX)) - 1.0);
      body[i].r[1] = -dim + (iy+0.5)*spacing + jit*(2*((double)rand()/(RAND_MAX)) - 1.0);
      body[i].r[2] =        (iz+0.5)*spacing + jit*(2*((double)rand()/(RAND_MAX)) - 1.0);
      fill( body[i].v, body[i].v+3, 0.0 );
    }
  }
  
  if(fig == "layered"){
    //petroleum starts at the bottom, water on top: the opposite of the buoyant
    //equilibrium, so the swap is visible as the simulation runs.
    //dim must match the box set by set_boundaries so the placement volume
    //(and therefore the SPH density estimate) is consistent with Volwater/Volpetroleum
    const double dim = box_dim;
    double zmid = dim/2.0;
    for(int i=0; i<N; i++ ){
      body[i].setdensity(0.0);
      body[i].reinit();
      fill( body[i].v, body[i].v+3, 0.0 );
      body[i].r[0] = (2*((double)rand()/(RAND_MAX)) - 1.0)*dim;//(-dim,dim)
      body[i].r[1] = (2*((double)rand()/(RAND_MAX)) - 1.0)*dim;//(-dim,dim)
      if(i < N1){
        body[i].setphase(1);//water: upper band
        body[i].r[2] = zmid + ((double)rand()/(RAND_MAX))*(dim-zmid);
      }
      else{
        body[i].setphase(2);//petroleum: lower band
        body[i].r[2] = ((double)rand()/(RAND_MAX))*zmid;
      }
    }
  }

  if(fig == "file"){
    //Parse line by line with the commas turned into whitespace. Reading the
    //separator with `file>>comma` into a string only works when the commas are
    //themselves surrounded by spaces: against "x,y,z" the string extraction
    //runs to the next whitespace and swallows the entire rest of the line, so
    //each particle ends up holding one field from seven consecutive lines and
    //phase_int parses "-0" as 0 -> setphase(0) -> every particle becomes petroleum.
    //Going through a stringstream makes the reader independent of the spacing.
    ifstream file;
    file.open("000.data");
    if(!file.is_open()){
      cerr<<"initial_cond: could not open 000.data"<<endl;
      exit(1);
    }
    int Nfile = N1+N2;//particles stored in the equilibrium file
    string line;
    int nread = 0;
    while(nread < Nfile && getline(file, line)){
      for(size_t c=0; c<line.size(); c++) if(line[c] == ',') line[c] = ' ';
      istringstream ss(line);
      double x,y,z,vx,vy,vz;
      int phase_int;
      if( !(ss>>x>>y>>z>>vx>>vy>>vz>>phase_int) ) continue;//skip blank/short lines
      if(phase_int != 1 && phase_int != 2){
        cerr<<"initial_cond: 000.data line "<<nread+1<<": phase is "<<phase_int
            <<", expected 1 or 2"<<endl;
        exit(1);
      }
      body[nread].setdensity(0.0);
      body[nread].reinit();
      body[nread].r[0]=x;  body[nread].r[1]=y;  body[nread].r[2]=z;
      body[nread].v[0]=vx; body[nread].v[1]=vy; body[nread].v[2]=vz;
      body[nread].setphase((unsigned int)phase_int);
      nread++;
    }
    file.close();
    if(nread != Nfile){
      cerr<<"initial_cond: 000.data supplied "<<nread<<" particles, expected "
          <<Nfile<<" (columns must be x,y,z,vx,vy,vz,phase)"<<endl;
      exit(1);
    }

    //Water drop: a symmetric cluster released above the pool, built from the
    //same offsets that determined Ndrop and scaled by the water particle's rest
    //spacing so the drop begins at the pool's own density. Deliberately not
    //jittered -- the symmetry is the point, so the splash stays axisymmetric.
    vector<Offset> off = drop_offsets();
    for(size_t k=0; k<off.size(); k++){
      int i = Nfile + (int)k;
      body[i].setdensity(0.0);
      body[i].reinit();
      fill( body[i].v, body[i].v+3, 0.0 );
      body[i].setphase(1);//water
      body[i].r[0] = dropX + off[k].x*s_water;
      body[i].r[1] = dropY + off[k].y*s_water;
      body[i].r[2] = dropZ + off[k].z*s_water;
    }
  }
}

void press_calc(Particle *body){
  //to make the code faster I calculate the density and the pressure in the same loop
#pragma omp parallel for
  for(int i=0; i<N; i++){
    body[i].reinit();
    body[i].density = body[i].mass*kernel_cubic( body[i].h, 0.0 );
  }

  //uniform grid neighbour search: cell size = largest smoothing length, so any
  //pair within a particle's own h is guaranteed to lie within the surrounding
  //3x3x3 block of cells. Turns the neighbour search from O(N^2) into O(N).
  double xmin=body[0].r[0], xmax=body[0].r[0];
  double ymin=body[0].r[1], ymax=body[0].r[1];
  double zmin=body[0].r[2], zmax=body[0].r[2];
  double hmax=body[0].h;
  for(int k=1; k<N; k++){
    if(body[k].r[0]<xmin) xmin=body[k].r[0];
    if(body[k].r[0]>xmax) xmax=body[k].r[0];
    if(body[k].r[1]<ymin) ymin=body[k].r[1];
    if(body[k].r[1]>ymax) ymax=body[k].r[1];
    if(body[k].r[2]<zmin) zmin=body[k].r[2];
    if(body[k].r[2]>zmax) zmax=body[k].r[2];
    if(body[k].h>hmax) hmax=body[k].h;
  }
  //A non-finite position makes the (int) casts below produce arbitrary indices
  //and cells[...] then reads out of bounds, which shows up as a bare segfault
  //with no clue where it came from. Fail loudly instead.
  if( !(xmin<=xmax) || !(ymin<=ymax) || !(zmin<=zmax) || !(hmax>0.0) ){
    cerr<<"press_calc: particle positions are not finite (the run has diverged)."
        <<" x["<<xmin<<","<<xmax<<"] y["<<ymin<<","<<ymax<<"] z["<<zmin<<","<<zmax
        <<"] hmax="<<hmax<<endl;
    exit(1);
  }
  double cellsize = hmax;
  int nx = (int)((xmax-xmin)/cellsize) + 2;
  double cellsize_x = cellsize;
  if(periodic_x){
    //the grid spans exactly one period, so cell 0 and cell nx-1 are neighbours.
    //cellsize_x = Lx/nx >= hmax keeps the 3-cell stencil wide enough.
    xmin = box_x0;
    nx = (int)(box_Lx/hmax);
    if(nx < 3){
      cerr<<"press_calc: periodic x needs box_Lx >= 3*h (Lx="<<box_Lx<<", h="<<hmax<<")"<<endl;
      exit(1);
    }
    cellsize_x = box_Lx/nx;
  }
  int ny = (int)((ymax-ymin)/cellsize) + 2;
  int nz = (int)((zmax-zmin)/cellsize) + 2;
  vector< vector<int> > cells(nx*ny*nz);
  for(int k=0; k<N; k++){
    int cx = (int)((body[k].r[0]-xmin)/cellsize_x);
    int cy = (int)((body[k].r[1]-ymin)/cellsize);
    int cz = (int)((body[k].r[2]-zmin)/cellsize);
    cx = clamp(cx, 0, nx-1);//rounding at the extremes must not index out of bounds
    cy = clamp(cy, 0, ny-1);
    cz = clamp(cz, 0, nz-1);
    cells[cx + nx*(cy + ny*cz)].push_back(k);
  }

#pragma omp parallel for
  for(int i=0; i<N; i++) {
    double dr0, dr1, dr2, normdr;
    double rhoi = 0.0;//accumulate locally, flush once: see the note in forces()
    int cx = (int)((body[i].r[0]-xmin)/cellsize_x);
    int cy = (int)((body[i].r[1]-ymin)/cellsize);
    int cz = (int)((body[i].r[2]-zmin)/cellsize);
    cx = clamp(cx, 0, nx-1);//must match the binning clamp above
    cy = clamp(cy, 0, ny-1);
    cz = clamp(cz, 0, nz-1);
    for(int dx=-1; dx<=1; dx++){
      int nxc = cx+dx;
      if(periodic_x) nxc = (nxc + nx) % nx;//wrap across the periodic seam
      else if(nxc<0 || nxc>=nx) continue;
      for(int dy=-1; dy<=1; dy++){
        int nyc = cy+dy;
        if(nyc<0 || nyc>=ny) continue;
        for(int dz=-1; dz<=1; dz++){
          int nzc = cz+dz;
          if(nzc<0 || nzc>=nz) continue;
          vector<int>& cell = cells[nxc + nx*(nyc + ny*nzc)];
          for(size_t idx=0; idx<cell.size(); idx++){
            int j = cell[idx];
            if(j<=i) continue;//keep the same i<j pairing as the original O(N^2) loop
            dr0 = min_image_x( body[i].r[0] - body[j].r[0] );
            dr1 = body[i].r[1] - body[j].r[1];
            dr2 = body[i].r[2] - body[j].r[2];
            normdr = sqrt( dr0*dr0 + dr1*dr1 + dr2*dr2 );
            if(normdr <= body[i].h){
              body[i].neighbourlist.push_back(j);
              if(body[i].phase == body[j].phase){
                rhoi += body[j].mass*kernel_cubic( body[i].h, normdr );
                #pragma omp atomic
                body[j].density += body[i].mass*kernel_cubic( body[j].h, normdr );
              }
            }
          }
        }
      }
    }
    #pragma omp atomic
    body[i].density += rhoi;
  }
  //Shepard renormalization: the same-phase-only density sum above
  //under-counts near the interface, since part of a particle's support
  //radius is legitimately occupied by the other phase but contributes
  //nothing to its same-phase density. gamma_i approximates the fraction of
  //the kernel's normalized volume integral actually sampled by same-phase
  //neighbours (1.0 for a fully-sampled bulk particle, <1.0 near an
  //interface); dividing the raw density by it corrects the under-count
  //instead of just clamping the resulting negative pressure afterward.
  vector<double> gamma(N);
#pragma omp parallel for
  for(int i=0; i<N; i++){
    gamma[i] = body[i].mass/body[i].density*kernel_cubic( body[i].h, 0.0 );
  }
#pragma omp parallel for
  for(int i=0; i<N; i++){
    double dr0, dr1, dr2, normdr;
    double gi = 0.0;//accumulate locally, flush once: see the note in forces()
    for(size_t k=0; k<body[i].neighbourlist.size(); k++){
      int neigh = body[i].neighbourlist[k];
      if(body[i].phase == body[neigh].phase){
        dr0 = min_image_x( body[i].r[0] - body[neigh].r[0] );
        dr1 = body[i].r[1] - body[neigh].r[1];
        dr2 = body[i].r[2] - body[neigh].r[2];
        normdr = sqrt( dr0*dr0 + dr1*dr1 + dr2*dr2 );
        gi += body[neigh].mass/body[neigh].density*kernel_cubic( body[i].h, normdr );
        #pragma omp atomic
        gamma[neigh] += body[i].mass/body[i].density*kernel_cubic( body[neigh].h, normdr );
      }
    }
    #pragma omp atomic
    gamma[i] += gi;
  }
#pragma omp parallel for
  for(int i=0; i<N; i++){
    body[i].density = body[i].density/gamma[i];
  }

#pragma omp parallel for
  for(int i=0; i<N; i++) {
    body[i].pressure = body[i].rhozero*body[i].cvel*body[i].cvel*( pow((body[i].density/body[i].rhozero),body[i].Gamma) - 1.0 )/body[i].Gamma;//Monagan equation 10.1
    if(body[i].pressure < 0.0) body[i].pressure = 0.0;//forbid tensile pressure: avoids the SPH tensile instability at under-sampled (interface) particles
  }
}


void color_calc(Particle *body){
  //color calculations
#pragma omp parallel for
  for(int i=0; i<N; i++) {
    double dr0, dr1, dr2, normdr;
    //accumulate locally, flush once: see the note in forces()
    double ni0 = 0.0, ni1 = 0.0, ni2 = 0.0;
    for(size_t j=0; j<body[i].neighbourlist.size(); j++){
      int neigh = body[i].neighbourlist[j];

      if(body[i].phase == body[neigh].phase){
        dr0 = min_image_x( body[i].r[0] - body[neigh].r[0] );
        dr1 = body[i].r[1] - body[neigh].r[1];
        dr2 = body[i].r[2] - body[neigh].r[2];
        normdr = sqrt( dr0*dr0 + dr1*dr1 + dr2*dr2 );
        if(normdr < 1e-10) normdr = 1e-10;//avoid 0/0 for coincident particles

        dr0 =dr0/normdr;
        dr1 =dr1/normdr;
        dr2 =dr2/normdr;
        
        ni0 += body[neigh].mass*Dkernel_cubic(body[i].h,normdr)*dr0/body[neigh].density;
        ni1 += body[neigh].mass*Dkernel_cubic(body[i].h,normdr)*dr1/body[neigh].density;
        ni2 += body[neigh].mass*Dkernel_cubic(body[i].h,normdr)*dr2/body[neigh].density;
        
        //newtons third law

        #pragma omp atomic
        body[neigh].normal[0]-= body[i].mass*Dkernel_cubic(body[neigh].h,normdr)*dr0/body[i].density;
        #pragma omp atomic
        body[neigh].normal[1]-= body[i].mass*Dkernel_cubic(body[neigh].h,normdr)*dr1/body[i].density;
        #pragma omp atomic
        body[neigh].normal[2]-= body[i].mass*Dkernel_cubic(body[neigh].h,normdr)*dr2/body[i].density;
      }
    }
    #pragma omp atomic
    body[i].normal[0] += ni0;
    #pragma omp atomic
    body[i].normal[1] += ni1;
    #pragma omp atomic
    body[i].normal[2] += ni2;
  }
}


void forces(Particle *body, boundaries& limit){
  press_calc(body);
  color_calc(body);
#pragma omp parallel for
  for(int i=0; i<N; i++){
    double dr0, dr1, dr2, normdr, prhoj;
    double dv0, dv1, dv2;
    double normnormali, normnormalj;
    //accumulate particle i's own share locally and flush it once at the end of
    //the loop body. Writing body[i].a directly is a data race: another thread
    //handling a lower-indexed particle reaches i through its neighbourlist and
    //updates the same address atomically, and mixing atomic with non-atomic
    //access loses updates -- which breaks Newton's third law at random and
    //injects momentum and energy.
    double ai0 = 0.0, ai1 = 0.0, ai2 = 0.0;
    double prhoi = body[i].pressure/( body[i].density*body[i].density );
    for(size_t j=0; j<body[i].neighbourlist.size(); j++){
      int neigh = body[i].neighbourlist[j];

      prhoj = body[neigh].pressure/( body[neigh].density*body[neigh].density );
      
      dr0 = min_image_x( body[i].r[0] - body[neigh].r[0] );
      dr1 = body[i].r[1] - body[neigh].r[1];
      dr2 = body[i].r[2] - body[neigh].r[2];
      normdr = sqrt( dr0*dr0 + dr1*dr1 + dr2*dr2 );
      if(normdr < 1e-10) normdr = 1e-10;//avoid 0/0 for coincident particles

      dr0 = dr0/normdr;
      dr1 = dr1/normdr;
      dr2 = dr2/normdr;

      dv0 = body[neigh].v[0]-body[i].v[0];
      dv1 = body[neigh].v[1]-body[i].v[1];
      dv2 = body[neigh].v[2]-body[i].v[2];
      
      //Monaghan artificial viscosity.
      //dr0..2 are the UNIT vector rhat = (r_i - r_j)/|r_i - r_j| and dv0..2 are
      //(v_j - v_i), so with vdotrhat = (v_j - v_i).rhat the standard dot product
      //is v_ij . r_ij = (v_i - v_j).(r_i - r_j) = -normdr*vdotrhat. The pair is
      //approaching exactly when that is negative, i.e. when vdotrhat > 0.
      double hbar_ij = 0.5*(body[i].h + body[neigh].h);//also used by the viscosity below
      double pi_ij = 0.0;
      double vdotrhat = dv0*dr0 + dv1*dr1 + dv2*dr2;
      if(av_alpha > 0.0 && vdotrhat > 0.0){//approaching pairs only
        double cbar   = 0.5*(body[i].cvel + body[neigh].cvel);
        double rhobar = 0.5*(body[i].density + body[neigh].density);
        double mu = -hbar_ij*normdr*vdotrhat/( normdr*normdr + av_eps*hbar_ij*hbar_ij );//<0 here
        pi_ij = ( -av_alpha*cbar*mu + av_beta*mu*mu )/rhobar;//>=0
      }

      //pressure (+ artificial viscosity, which enters identically)
      double press_term = prhoi + prhoj + pi_ij;
      double dwi = Dkernel_cubic(body[i].h,normdr);
      double dwj = Dkernel_cubic(body[neigh].h,normdr);

      ai0 += -body[i].density*press_term*body[neigh].mass*dwi*dr0;
      ai1 += -body[i].density*press_term*body[neigh].mass*dwi*dr1;
      ai2 += -body[i].density*press_term*body[neigh].mass*dwi*dr2;

      #pragma omp atomic
      body[neigh].a[0] += body[neigh].density*press_term*body[i].mass*dwj*dr0;
      #pragma omp atomic
      body[neigh].a[1] += body[neigh].density*press_term*body[i].mass*dwj*dr1;
      #pragma omp atomic
      body[neigh].a[2] += body[neigh].density*press_term*body[i].mass*dwj*dr2;
      
      //Physical viscosity, Morris et al. (1997).
      //This replaces the Mueller form, which needed a THIRD kernel family
      //(DDkernel_vis, since removed) purely to supply a Laplacian. Morris builds the same
      //Laplacian out of the first derivative of the kernel already in use --
      //(r . grad W)/(r^2 + eta^2) is the Brookshaw estimator -- so the whole
      //solver now differentiates one and only one W.
      //Coefficient note: Morris carries (mu_i + mu_j), NOT the average
      //0.5*(mu_i + mu_j) the Mueller form used. That is not a doubling: the
      //Brookshaw Laplacian carries a compensating factor of 2, and for uniform
      //mu both forms reduce to (mu/rho) grad^2 v. The effective viscosity of
      //the run is unchanged.
      double visc = (body[i].mu + body[neigh].mu)*normdr
                    /( normdr*normdr + av_eps*hbar_ij*hbar_ij )/body[neigh].density;
      //dwi <= 0, so -visc*dwi >= 0 and the term drives v_i toward v_j: dissipative.
      ai0 += -visc*dwi*dv0;
      ai1 += -visc*dwi*dv1;
      ai2 += -visc*dwi*dv2;

      double viscj = (body[i].mu + body[neigh].mu)*normdr
                     /( normdr*normdr + av_eps*hbar_ij*hbar_ij )/body[i].density;
      #pragma omp atomic
      body[neigh].a[0] += viscj*dwj*dv0;
      #pragma omp atomic
      body[neigh].a[1] += viscj*dwj*dv1;
      #pragma omp atomic
      body[neigh].a[2] += viscj*dwj*dv2;
      
      //Surface tension
      
      if(body[i].phase == body[neigh].phase){
        normnormali = sqrt(body[i].normal[0]*body[i].normal[0]+body[i].normal[1]*body[i].normal[1]+body[i].normal[2]*body[i].normal[2]);
        normnormalj = sqrt(body[neigh].normal[0]*body[neigh].normal[0]+body[neigh].normal[1]*body[neigh].normal[1]+body[neigh].normal[2]*body[neigh].normal[2]);
        
        if(normnormali > body[i].threshold ){
          double ddk = DDkernel_cubic(body[i].h,normdr);
          double fx = -body[i].sigma*body[i].normal[0]*body[neigh].mass*ddk/(normnormali*body[neigh].density);
          double fy = -body[i].sigma*body[i].normal[1]*body[neigh].mass*ddk/(normnormali*body[neigh].density);
          double fz = -body[i].sigma*body[i].normal[2]*body[neigh].mass*ddk/(normnormali*body[neigh].density);
          double fmag = sqrt(fx*fx+fy*fy+fz*fz);
          if(fmag > surface_tension_amax){
            double scale = surface_tension_amax/fmag;
            fx *= scale; fy *= scale; fz *= scale;
          }
          ai0 += fx;
          ai1 += fy;
          ai2 += fz;
        }
        if(normnormalj > body[neigh].threshold ){
          double ddk = DDkernel_cubic(body[neigh].h,normdr);
          double fx = body[neigh].sigma*body[neigh].normal[0]*body[i].mass*ddk/(normnormalj*body[i].density);
          double fy = body[neigh].sigma*body[neigh].normal[1]*body[i].mass*ddk/(normnormalj*body[i].density);
          double fz = body[neigh].sigma*body[neigh].normal[2]*body[i].mass*ddk/(normnormalj*body[i].density);
          double fmag = sqrt(fx*fx+fy*fy+fz*fz);
          if(fmag > surface_tension_amax){
            double scale = surface_tension_amax/fmag;
            fx *= scale; fy *= scale; fz *= scale;
          }
          #pragma omp atomic
          body[neigh].a[0] += fx;
          #pragma omp atomic
          body[neigh].a[1] += fy;
          #pragma omp atomic
          body[neigh].a[2] += fz;
        }
      }
      else{
        //Interfacial repulsion between unlike phases: -Dkernel_cubic is >=0 for
        //normdr in [0,h] and vanishes smoothly at h, giving a short-range push
        //apart that keeps the two phases from interpenetrating/mixing
        ai0 += sigma_interphase*body[neigh].mass*(-Dkernel_cubic(body[i].h,normdr))*dr0/body[neigh].density;
        ai1 += sigma_interphase*body[neigh].mass*(-Dkernel_cubic(body[i].h,normdr))*dr1/body[neigh].density;
        ai2 += sigma_interphase*body[neigh].mass*(-Dkernel_cubic(body[i].h,normdr))*dr2/body[neigh].density;

        #pragma omp atomic
        body[neigh].a[0] += -sigma_interphase*body[i].mass*(-Dkernel_cubic(body[neigh].h,normdr))*dr0/body[i].density;
        #pragma omp atomic
        body[neigh].a[1] += -sigma_interphase*body[i].mass*(-Dkernel_cubic(body[neigh].h,normdr))*dr1/body[i].density;
        #pragma omp atomic
        body[neigh].a[2] += -sigma_interphase*body[i].mass*(-Dkernel_cubic(body[neigh].h,normdr))*dr2/body[i].density;
      }

    }
    
    //friction or damping. Everything accumulated in this loop is divided by
    //density below, so the density premultiply is what makes `friction` an
    //actual damping rate a = -friction*v. Without it the effective rate is
    //friction/density ~ 1e-3 1/s for water (tau ~ 1000 s, longer than the whole
    //run) and the two phases damp at different rates.
    ai0 += -friction*body[i].v[0]*body[i].density;
    ai1 += -friction*body[i].v[1]*body[i].density;
    ai2 += -friction*body[i].v[2]*body[i].density;

    //gravity
    ai2 += -9.82*body[i].density;

    //buoyancy (disabled: buoyancy==0, see the constant's definition)
    ai2 += body[i].buoyancy*(body[i].density - body[i].rhozero)*(-9.82);

    //boundary repulsion: smooth push-back as a particle nears a wall, integrated
    //through the same leapfrog scheme as gravity/pressure (note the density
    //premultiply here matches gravity/buoyancy above, since everything in this
    //loop gets divided by density below)
    if(!periodic_x)//no x walls when x is periodic
      ai0 += ( wall_repulsion(body[i].r[0]-limit.boundaryX[0], boundary_r0, boundary_D)
              - wall_repulsion(limit.boundaryX[1]-body[i].r[0], boundary_r0, boundary_D) )*body[i].density;
    ai1 += ( wall_repulsion(body[i].r[1]-limit.boundaryY[0], boundary_r0, boundary_D)
            - wall_repulsion(limit.boundaryY[1]-body[i].r[1], boundary_r0, boundary_D) )*body[i].density;
    ai2 += ( wall_repulsion(body[i].r[2]-limit.boundaryZ[0], boundary_r0, boundary_D)
            - wall_repulsion(limit.boundaryZ[1]-body[i].r[2], boundary_r0, boundary_D) )*body[i].density;
    double frock[3];
    rock_force(body[i].r, frock);//rocks repel like the walls, see rocks.h
    ai0 += frock[0]*body[i].density;
    ai1 += frock[1]*body[i].density;
    ai2 += frock[2]*body[i].density;

    #pragma omp atomic
    body[i].a[0] += ai0;
    #pragma omp atomic
    body[i].a[1] += ai1;
    #pragma omp atomic
    body[i].a[2] += ai2;
  }
  for(int i=0; i<N; i++){
    for(int j=0; j<3; j++){
      body[i].a[j] /= body[i].density;
    }
  }
}

//Monaghan-style boundary repulsion, bounded (no singularity): ramps from 0 at
//d=r0 up to D at d=0, so it stays numerically safe under a fixed dt.
//Past the wall (d<=0) it holds at D rather than dropping back to 0 -- returning
//0 there left a particle that overshot with nothing pushing it back in.
double wall_repulsion(double d, double r0, double D){
  if(d >= r0) return 0.0;
  if(d <= 0.0) return D;
  double x = 1.0 - d/r0;
  return D*x*x;
}

void useboundaries(Particle& body, boundaries& limit){
  /*
  //sphere
  double vecS[3], vecI[3], vecN[3], normvecN, normvel;
  if(body.r[0]*body.r[0] + body.r[1]*body.r[1] + body.r[2]*body.r[2] >= 1.0*1.0){
    normvel = sqrt(body.v[0]*body.v[0] + body.v[1]*body.v[1] + body.v[2]*body.v[2]);
    vecI[0]  = body.v[0]/normvel;
    vecI[1]  = body.v[1]/normvel;
    vecI[2]  = body.v[2]/normvel;

    vecN[0] = body.r[0] - 0.0;
    vecN[1] = body.r[1] - 0.0;
    vecN[2] = body.r[2] - 0.0;
    normvecN = sqrt(vecN[0]*vecN[0]+vecN[1]*vecN[1]+vecN[2]*vecN[2]);
    vecN[0] /= normvecN;
    vecN[1] /= normvecN;
    vecN[2] /= normvecN;

    vecS[0] = -1*(2*(vecN[0]*vecI[0] + vecN[1]*vecI[1] + vecN[2]*vecI[2])*vecN[0] - vecI[0]);
    vecS[1] = -1*(2*(vecN[0]*vecI[0] + vecN[1]*vecI[1] + vecN[2]*vecI[2])*vecN[1] - vecI[1]);
    vecS[2] = -1*(2*(vecN[0]*vecI[0] + vecN[1]*vecI[1] + vecN[2]*vecI[2])*vecN[2] - vecI[2]);

    body.v[0] = 0.97*normvel*vecS[0];
    body.v[1] = 0.97*normvel*vecS[1];
    body.v[2] = 0.97*normvel*vecS[2];
  }
  */

  //plane z=mx+b
  //m=atan(angle)
  /*
  double angle = 0.0;
  double vecS[3], vecI[3], vecN[3], normvecN, normvel, m=atan(angle*M_PI/180.0), b=0.0;
  vecN[0] = sin(angle*M_PI/180.0);
  vecN[1] = 0.0;
  vecN[2] = cos(angle*M_PI/180.0);
  if(body.r[2] <= -m*body.r[0]+b){
    normvel = sqrt(body.v[0]*body.v[0] + body.v[1]*body.v[1] + body.v[2]*body.v[2]);
    vecI[0]  = body.v[0]/normvel;
    vecI[1]  = body.v[1]/normvel;
    vecI[2]  = body.v[2]/normvel;
    
    vecS[0] = -1*(2*(vecN[0]*vecI[0] + vecN[1]*vecI[1] + vecN[2]*vecI[2])*vecN[0] - vecI[0]);
    vecS[1] = -1*(2*(vecN[0]*vecI[0] + vecN[1]*vecI[1] + vecN[2]*vecI[2])*vecN[1] - vecI[1]);
    vecS[2] = -1*(2*(vecN[0]*vecI[0] + vecN[1]*vecI[1] + vecN[2]*vecI[2])*vecN[2] - vecI[2]);

    body.v[0] = 0.97*normvel*vecS[0];
    body.v[1] = 0.97*normvel*vecS[1];
    body.v[2] = 0.97*normvel*vecS[2];
  }
  if(body.r[0] <= limit.boundaryX[0])    
    body.v[0] *= -0.97;//-1.0;
  if(body.r[0] >= limit.boundaryX[1])
    body.v[0] *= -0.97;//-1.0;
  if(body.r[1] <= limit.boundaryY[0])
    body.v[1] *= -0.97;//-1.0;
  if(body.r[1] >= limit.boundaryY[1])
    body.v[1] *= -0.97;//-1.0;
  */

  //cube
  //Safety net for anything that still tunnels through a wall. This clamps the
  //particle onto the wall and removes the inward normal velocity, so it is
  //strictly dissipative. The previous version mirrored the particle to +d and
  //flipped v, which teleported it from outside (where wall_repulsion was 0)
  //into the repulsion field, creating D*r0/3 ~ 25 J/kg out of nothing -- about
  //7 m/s of free velocity per event. The smooth wall_repulsion above does the
  //real work; this should almost never fire.
  if(periodic_x){
    body.r[0] = wrap_x(body.r[0]);//periodic: wrap instead of clamping
  }
  else if(body.r[0] < limit.boundaryX[0]){
    body.r[0] = limit.boundaryX[0];
    if(body.v[0] < 0.0) body.v[0] = 0.0;
  }
  else if(body.r[0] > limit.boundaryX[1]){
    body.r[0] = limit.boundaryX[1];
    if(body.v[0] > 0.0) body.v[0] = 0.0;
  }
  if(body.r[1] < limit.boundaryY[0]){
    body.r[1] = limit.boundaryY[0];
    if(body.v[1] < 0.0) body.v[1] = 0.0;
  }
  if(body.r[1] > limit.boundaryY[1]){
    body.r[1] = limit.boundaryY[1];
    if(body.v[1] > 0.0) body.v[1] = 0.0;
  }
  if(body.r[2] < limit.boundaryZ[0]){
    body.r[2] = limit.boundaryZ[0];
    if(body.v[2] < 0.0) body.v[2] = 0.0;
  }
  if(body.r[2] > limit.boundaryZ[1]){
    body.r[2] = limit.boundaryZ[1];
    if(body.v[2] > 0.0) body.v[2] = 0.0;
  }
  rock_clamp(body.r, body.v);//same safety net for the rocks, see rocks.h
}

void initialize(Particle *body, boundaries& limit, int phases, string shape){
  const double dim = box_dim;
  //Ceiling at 2*dim: the fluid occupies the lower half, so it has a real free
  //surface to relieve the inward push of the wall repulsion. At zheight==dim it
  //exactly fills the box and the boundary force can only compress it; at
  //dropZ+0.5 it had 5.5x its own volume to spread into.
  double zheight = 2.0*dim;
  if(shape == "file" && Ndrop > 0){
    //raise the ceiling clear of the drop: its own top, plus a few smoothing
    //lengths so the ceiling repulsion is not touching it at release
    zheight = dropZ + drop_radius_cells*s_water + 4.0*h_water_ref;
  }
  limit.set_boundaries("cube", dim, zheight);
  setup_rocks(dim);//after set_boundaries: needs box_Lx for the periodic distances
  initial_cond(shape, body, phases);
}

void verlet(Particle *body, boundaries& limit, double dt, unsigned int cycle){
  if(cycle == 1){
#pragma omp parallel for
    for(int i=0;i<N;i++){
      for(int j=0;j<3;j++){
	body[i].v[j] += 0.5*dt*body[i].a[j];
	body[i].r[j] += dt*body[i].v[j];
      }
      body[i].r[0] = wrap_x(body[i].r[0]);//periodic x: keep positions inside one period for the grid
    }
  }
  if(cycle == 2) {
#pragma omp parallel for
    for(int i=0;i<N;i++){
      for(int j=0;j<3;j++){
	      body[i].v[j] += 0.5*dt*body[i].a[j];
      }
      useboundaries(body[i], limit);
      //density is not reset here: press_calc overwrites it at the start of the
      //next step anyway, and zeroing it made the value unusable for diagnostics
    }
  }
}

void kinetic_E(Particle *body, int count){
  string name;
  stringstream convert; // stringstream used for the conversion
  ofstream file;
  double kinetic=0.0;
  double zsum1=0.0, zsum2=0.0;
  double rhosum=0.0;//mean density: should sit near rhozero once settled
  int n1=0, n2=0;

  convert<<setfill('0')<<setw(3)<<count<<".dat";//add the value of Number to the characters in the stream
  name = convert.str();

  file.open(name.c_str());
  for(int i=0;i<N;i++){
    kinetic += body[i].kin_energy();
    rhosum += body[i].density/body[i].rhozero;
    if(body[i].phase == 1){
      zsum1 += body[i].r[2];
      n1++;
    }
    else{
      zsum2 += body[i].r[2];
      n2++;
    }
    if(i<N-1)
      file<<body[i].r[0]<<", "<<body[i].r[1]<<", "<<body[i].r[2]<<", "<<body[i].v[0]<<", "<<body[i].v[1]<<", "<<body[i].v[2]<<", "<<body[i].phase<<endl;
    else
      file<<body[i].r[0]<<", "<<body[i].r[1]<<", "<<body[i].r[2]<<", "<<body[i].v[0]<<", "<<body[i].v[1]<<", "<<body[i].v[2]<<", "<<body[i].phase<<endl;
  }
  file<<endl;
  file.close();
  cout<<count<<"   KE="<<kinetic<<"   rho/rho0="<<rhosum/N
      <<"   water_z="<<(n1?zsum1/n1:0.0)<<"   petroleum_z="<<(n2?zsum2/n2:0.0)<<endl;
}


int main(){
  srand(0);  //srand((unsigned)time(0)) 
  vector<Particle> particles(N);
  Particle *liquid = particles.data();
  boundaries limit;
  unsigned int phases = 1;
  int count=0;
  //"cube" places all N particles uniformly in the box, whose volume is exactly
  //VolT, so the run starts at rest density everywhere. "file" reads 000.data,
  //a settled N1+N2 state (see initial_cond).
  initialize(liquid, limit, phases, "file"); //cube, sphere, file or layered
  for(int ts=0;ts<TS; ts++){
    double cvel_frac = (ts < cvel_ramp_steps) ? (double)ts/cvel_ramp_steps : 1.0;
#pragma omp parallel for
    for(int i=0; i<N; i++) liquid[i].ramp_cvel(cvel_frac);
    rock_grow(ts);//rocks grow in over the first steps, see rocks.h
    verlet(liquid, limit, dt, 1);
    forces(liquid, limit);
    verlet(liquid, limit, dt, 2);
    if(ts%frame_every ==0){
      count++;
      kinetic_E(liquid,count);
    }
  }
}

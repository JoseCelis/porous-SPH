//rocks.h -- static polyhedral obstacles (rocks) for sph.cpp.
//
//Everything about the rocks lives here: their shape, how many there are and
//where (rock_mode, rocks.csv or n_rocks), the geometry, and the hooks sph.cpp calls:
//  setup_rocks(dim)   in initialize():    builds the rocks, writes rocks.txt
//  rock_grow(ts)      in main():          the growth ramp at the start of a run
//  rock_force(r, f)   in forces():        wall-style repulsion from every rock
//  rock_clamp(r, v)   in useboundaries(): dissipative safety net, like the walls
//
//Not self-contained, unlike kernel.h: it uses the box/wall definitions in
//sph.cpp (wall_repulsion, boundary_r0, boundary_D, periodic_x, min_image_x,
//wrap_x, frame_every), so sph.cpp includes it right after those.
#ifndef ROCKS_H
#define ROCKS_H

#include <random>
#include <array>

/******Static polyhedral obstacles (rocks)******/
//Convex polyhedra fixed in the fluid, treated exactly like the box walls:
//wall_repulsion() pushes particles off their faces in forces(), and
//useboundaries() clamps anything that still gets inside (dissipative).
//
//Every rock is a copy of the template shape in rock_shape_file
//(rock_shape.csv), moved to its own centre, scaled by its own size factor and
//rotated. The template is the set of points p with  nhat_k.p <= d_k  for every
//face k; each row of the file is nx, ny, nz, d: an OUTWARD normal (any length,
//it is normalised) and the face's distance from the rock centre in metres at
//size 1. Any convex shape can be written this way.
//
//How many rocks, and where:
//  rock_mode  : "file"   -> one rock per row of rock_file (rocks.csv): centre,
//                           size and rotation. The number of rows is the
//                           number of rocks. Read at startup, so editing the
//                           file needs no rebuild.
//               "random" -> n_rocks rocks with random centre, size and
//                           orientation (seeded by rock_seed, so a run is
//                           reproducible), none overlapping another rock or
//                           a wall. See the rock_random_* settings.
//               "none"   -> no rocks.
//The final geometry is written to rocks.txt at startup; render_pyvista.py draws
//the rocks from that file.
const char rock_mode[] = "file"; //"file", "random" or "none"
const char rock_file[] = "rocks.csv"; //"file" mode, read from the current directory like 000.data
const int n_rocks = 8;                //"random" mode only

const char rock_shape_file[] = "rock_shape.csv"; //template shape, read from the current directory
vector< array<double,4> > rock_faces; //{nx, ny, nz, d} per face, filled by load_rock_shape()

//"random" mode settings
const unsigned int rock_seed = 1;
const double rock_random_size[2] = {0.4, 0.8}; //size factor range
const double rock_random_z[2] = {0.0, 0.63};   //range for the rock centres' z (the pool is 0.63 deep)
//minimum clearance between two rocks, and between a rock and a wall. Keep it
//to a few particle spacings (s ~ 0.045 m), or the fluid cannot get through.
const double rock_min_gap = 0.10;

//000.data has fluid where the rocks sit. The rocks therefore grow from a point
//to full size over these steps, so they push the fluid aside gently instead of
//crushing every overlapping particle onto their surfaces at step 0.
const int rock_ramp_steps = 5000;
double rock_scale = 1.0; //0..1, set each step in main()

//one rock in world orientation, at full size: p is inside when
//  n[k].(p - c) <= d[k]  for every face k  (n unit and already rotated)
struct Rock {
  double c[3];
  vector< array<double,3> > n;
  vector<double> d;
  double radius; //bounding radius about c, used to keep rocks apart
};
vector<Rock> rocks; //filled by setup_rocks()

//Distance from r to rock R's surface (negative inside) and, in n_out, the
//outward unit normal of the face that sets it. For a convex polyhedron the
//largest plane distance is exact wherever the nearest feature is a face; near
//edges and corners it slightly underestimates, which only makes the push there
//a little stronger.
inline double rock_distance(const Rock& R, const double r[3], double n_out[3]){
  double p[3] = { min_image_x(r[0]-R.c[0]), r[1]-R.c[1], r[2]-R.c[2] };
  double best = -1e300;
  n_out[0] = 0.0; n_out[1] = 0.0; n_out[2] = 1.0;//overwritten by the first face
  for(size_t k=0; k<R.d.size(); k++){
    double d = R.n[k][0]*p[0] + R.n[k][1]*p[1] + R.n[k][2]*p[2] - rock_scale*R.d[k];
    if(d > best){
      best = d;
      n_out[0] = R.n[k][0]; n_out[1] = R.n[k][1]; n_out[2] = R.n[k][2];
    }
  }
  return best;
}

//Reads the template shape from rock_shape_file into rock_faces: CSV rows
//nx, ny, nz, d. Lines starting with # and blank lines are skipped, and so is a
//header row (any row whose first field is not a number).
void load_rock_shape(){
  rock_faces.clear();
  ifstream file(rock_shape_file);
  if(!file.is_open()){
    cerr<<"load_rock_shape: could not open "<<rock_shape_file<<endl;
    exit(1);
  }
  string line;
  int lineno = 0;
  while( getline(file, line) ){
    lineno++;
    size_t first = line.find_first_not_of(" \t\r");
    if(first == string::npos || line[first] == '#') continue;
    for(size_t k=0; k<line.size(); k++) if(line[k] == ',') line[k] = ' ';
    stringstream ss(line);
    array<double,4> f;
    int nread = 0;
    while( nread < 4 && ss>>f[nread] ) nread++;
    if(nread == 0) continue;//header row
    string extra;
    if(nread < 4 || ss>>extra){
      cerr<<"load_rock_shape: "<<rock_shape_file<<" line "<<lineno<<": expected 4 numbers "
          <<"(nx, ny, nz, d), got: "<<line<<endl;
      exit(1);
    }
    if( f[0]*f[0] + f[1]*f[1] + f[2]*f[2] < 1e-20 ){
      cerr<<"load_rock_shape: "<<rock_shape_file<<" line "<<lineno<<": the normal is zero"<<endl;
      exit(1);
    }
    if(f[3] <= 0.0){
      cerr<<"load_rock_shape: "<<rock_shape_file<<" line "<<lineno
          <<": d must be > 0 (the rock centre must be inside every face)"<<endl;
      exit(1);
    }
    rock_faces.push_back(f);
  }
  if(rock_faces.size() < 4){
    cerr<<"load_rock_shape: "<<rock_shape_file<<" has "<<rock_faces.size()
        <<" faces; a closed polyhedron needs at least 4"<<endl;
    exit(1);
  }
  //The faces must close the shape: the rock is infinite along a direction u
  //exactly when no face points even partly along u (nhat_k.u <= 0 for all k).
  //If such a u exists, one also exists along an edge direction, i.e. along
  //n_a x n_b for some pair of faces, so testing every pair (both signs) is a
  //complete check. Sampling directions is not: it misses a shape that is open
  //along a single line, e.g. a column with no bottom face.
  const int nf = rock_faces.size();
  bool any_pair = false;
  for(int a=0; a<nf; a++)
    for(int b=a+1; b<nf; b++){
      const array<double,4>& A = rock_faces[a];
      const array<double,4>& B = rock_faces[b];
      double c[3] = { A[1]*B[2]-A[2]*B[1], A[2]*B[0]-A[0]*B[2], A[0]*B[1]-A[1]*B[0] };
      double cn = sqrt( c[0]*c[0] + c[1]*c[1] + c[2]*c[2] );
      if(cn < 1e-12) continue;//parallel faces: no edge direction
      any_pair = true;
      for(int sgn=-1; sgn<=1; sgn+=2){
        double u[3] = { sgn*c[0]/cn, sgn*c[1]/cn, sgn*c[2]/cn };
        double best = -1.0;
        for(int k=0; k<nf; k++){
          const array<double,4>& f = rock_faces[k];
          double nn = sqrt( f[0]*f[0] + f[1]*f[1] + f[2]*f[2] );
          best = max(best, ( f[0]*u[0] + f[1]*u[1] + f[2]*u[2] )/nn);
        }
        if(best <= 1e-9){
          cerr<<"load_rock_shape: the faces in "<<rock_shape_file<<" do not close the shape: "
              <<"it is open towards ("<<u[0]<<", "<<u[1]<<", "<<u[2]<<")"<<endl;
          exit(1);
        }
      }
    }
  if(!any_pair){
    cerr<<"load_rock_shape: all faces in "<<rock_shape_file<<" are parallel; the shape is not closed"<<endl;
    exit(1);
  }
}

//Builds rock i from the template: rotation matrix Q (row-major), centre, size.
Rock make_rock(const double c[3], double size, const double Q[3][3]){
  Rock R;
  const int nf = rock_faces.size();
  for(int j=0; j<3; j++) R.c[j] = c[j];
  R.n.resize(nf);
  R.d.resize(nf);
  for(int k=0; k<nf; k++){
    const double *f = rock_faces[k].data();
    double nn = sqrt( f[0]*f[0] + f[1]*f[1] + f[2]*f[2] );
    for(int j=0; j<3; j++)
      R.n[k][j] = ( Q[j][0]*f[0] + Q[j][1]*f[1] + Q[j][2]*f[2] )/nn;
    R.d[k] = size*f[3];
  }
  //bounding radius = farthest template corner. Corners are where three face
  //planes meet and that lie inside all the other planes.
  double rmax = 0.0;
  for(int a=0; a<nf; a++)
    for(int b=a+1; b<nf; b++)
      for(int e=b+1; e<nf; e++){
        const double *A = R.n[a].data(), *B = R.n[b].data(), *E = R.n[e].data();
        //solve [A;B;E] x = [d_a;d_b;d_e] with Cramer's rule
        double det = A[0]*(B[1]*E[2]-B[2]*E[1]) - A[1]*(B[0]*E[2]-B[2]*E[0]) + A[2]*(B[0]*E[1]-B[1]*E[0]);
        if(fabs(det) < 1e-9) continue;
        double da = R.d[a], db = R.d[b], de = R.d[e];
        double x[3];
        x[0] = ( da*(B[1]*E[2]-B[2]*E[1]) - A[1]*(db*E[2]-B[2]*de) + A[2]*(db*E[1]-B[1]*de) )/det;
        x[1] = ( A[0]*(db*E[2]-B[2]*de) - da*(B[0]*E[2]-B[2]*E[0]) + A[2]*(B[0]*de-db*E[0]) )/det;
        x[2] = ( A[0]*(B[1]*de-db*E[1]) - A[1]*(B[0]*de-db*E[0]) + da*(B[0]*E[1]-B[1]*E[0]) )/det;
        bool inside = true;
        for(int k=0; k<nf && inside; k++)
          if( R.n[k][0]*x[0] + R.n[k][1]*x[1] + R.n[k][2]*x[2] > R.d[k] + 1e-9 ) inside = false;
        if(inside) rmax = max(rmax, sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] ));
      }
  R.radius = rmax;
  return R;
}

//rotation: yaw about z, then tilt about x
void yaw_tilt(double yaw_deg, double tilt_deg, double Q[3][3]){
  double cy = cos(yaw_deg*M_PI/180.0), sy = sin(yaw_deg*M_PI/180.0);
  double ct = cos(tilt_deg*M_PI/180.0), st = sin(tilt_deg*M_PI/180.0);
  double Rz[3][3] = {{cy,-sy,0},{sy,cy,0},{0,0,1}};
  double Rx[3][3] = {{1,0,0},{0,ct,-st},{0,st,ct}};
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++){
      Q[i][j] = 0.0;
      for(int k=0; k<3; k++) Q[i][j] += Rx[i][k]*Rz[k][j];
    }
}

//uniformly random rotation, from a random unit quaternion
void random_rotation(mt19937& gen, double Q[3][3]){
  normal_distribution<double> g(0.0, 1.0);
  double q[4], nq = 0.0;
  for(int i=0; i<4; i++){ q[i] = g(gen); nq += q[i]*q[i]; }
  nq = sqrt(nq);
  double w = q[0]/nq, x = q[1]/nq, y = q[2]/nq, z = q[3]/nq;
  double M[3][3] = {
    {1-2*(y*y+z*z), 2*(x*y-w*z),   2*(x*z+w*y)},
    {2*(x*y+w*z),   1-2*(x*x+z*z), 2*(y*z-w*x)},
    {2*(x*z-w*y),   2*(y*z+w*x),   1-2*(x*x+y*y)}};
  memcpy(Q, M, sizeof(M));
}

//Fills `rocks` and writes rocks.txt. dim = box half-width in x and y.
void setup_rocks(double dim){
  rocks.clear();
  if( strcmp(rock_mode, "none") != 0 ) load_rock_shape();
  if( strcmp(rock_mode, "none") == 0 ){
    //nothing to build
  }
  else if( strcmp(rock_mode, "file") == 0 ){
    //CSV, one rock per row: x, y, z, size, yaw_deg, tilt_deg. Lines starting
    //with # and blank lines are skipped, and so is a header row (any row
    //whose first field is not a number).
    ifstream file(rock_file);
    if(!file.is_open()){
      cerr<<"setup_rocks: could not open "<<rock_file<<" (rock_mode = \"file\")"<<endl;
      exit(1);
    }
    string line;
    int lineno = 0;
    while( getline(file, line) ){
      lineno++;
      size_t first = line.find_first_not_of(" \t\r");
      if(first == string::npos || line[first] == '#') continue;
      for(size_t k=0; k<line.size(); k++) if(line[k] == ',') line[k] = ' ';
      stringstream ss(line);
      double row[6];
      int nread = 0;
      while( nread < 6 && ss>>row[nread] ) nread++;
      if(nread == 0) continue;//header row
      string extra;
      if(nread < 6 || ss>>extra){
        cerr<<"setup_rocks: "<<rock_file<<" line "<<lineno<<": expected 6 numbers "
            <<"(x, y, z, size, yaw_deg, tilt_deg), got: "<<line<<endl;
        exit(1);
      }
      if(row[3] <= 0.0){
        cerr<<"setup_rocks: "<<rock_file<<" line "<<lineno<<": size must be > 0"<<endl;
        exit(1);
      }
      double Q[3][3];
      yaw_tilt(row[4], row[5], Q);
      rocks.push_back( make_rock(row, row[3], Q) );
    }
  }
  else if( strcmp(rock_mode, "random") == 0 ){
    mt19937 gen(rock_seed);//own generator: does not shift the rand() stream used elsewhere
    uniform_real_distribution<double> u(0.0, 1.0);
    const int max_tries = 100000;
    int tries = 0;
    while( (int)rocks.size() < n_rocks ){
      if(++tries > max_tries){
        cerr<<"setup_rocks: could only place "<<rocks.size()<<" of "<<n_rocks
            <<" rocks without overlap. Lower n_rocks, rock_random_size or rock_min_gap."<<endl;
        exit(1);
      }
      double Q[3][3];
      random_rotation(gen, Q);
      double size = rock_random_size[0] + u(gen)*(rock_random_size[1]-rock_random_size[0]);
      double c[3] = { 0.0, 0.0, 0.0 };
      Rock R = make_rock(c, size, Q);
      double margin = R.radius + rock_min_gap;
      //x is periodic when periodic_x, so only y and the floor need a wall margin
      double xlo = periodic_x ? -dim : -dim + margin, xhi = periodic_x ? dim : dim - margin;
      double ylo = -dim + margin, yhi = dim - margin;
      double zlo = max(rock_random_z[0], margin), zhi = rock_random_z[1];
      if(xlo >= xhi || ylo >= yhi || zlo >= zhi) continue;//this rock is too big for the box
      R.c[0] = xlo + u(gen)*(xhi-xlo);
      R.c[1] = ylo + u(gen)*(yhi-ylo);
      R.c[2] = zlo + u(gen)*(zhi-zlo);
      bool clear = true;
      for(size_t j=0; j<rocks.size() && clear; j++){
        double dx = min_image_x(R.c[0]-rocks[j].c[0]), dy = R.c[1]-rocks[j].c[1], dz = R.c[2]-rocks[j].c[2];
        if( sqrt(dx*dx+dy*dy+dz*dz) < R.radius + rocks[j].radius + rock_min_gap ) clear = false;
      }
      if(clear) rocks.push_back(R);
    }
  }
  else{
    cerr<<"setup_rocks: unknown rock_mode \""<<rock_mode<<"\" (use \"file\", \"random\" or \"none\")"<<endl;
    exit(1);
  }

  //rocks.txt: what render_pyvista.py draws. Not *.dat, so the frame globs
  //(and `make fresh`) leave it alone.
  ofstream f("rocks.txt");
  f<<"# rocks used by sph.cpp. For each rock: its centre, then one line per face,"<<endl
   <<"# \"nx ny nz d\" = outward unit normal and distance from the centre at full size."<<endl
   <<"# The rocks grow from a point to full size over ramp_steps steps."<<endl;
  f<<setprecision(10);
  f<<"ramp_steps "<<rock_ramp_steps<<endl;
  f<<"steps_per_frame "<<frame_every<<endl;
  for(size_t i=0; i<rocks.size(); i++){
    f<<"rock "<<rocks[i].c[0]<<" "<<rocks[i].c[1]<<" "<<rocks[i].c[2]<<endl;
    for(size_t k=0; k<rocks[i].d.size(); k++)
      f<<rocks[i].n[k][0]<<" "<<rocks[i].n[k][1]<<" "<<rocks[i].n[k][2]<<" "<<rocks[i].d[k]<<endl;
  }
  cout<<"rocks: "<<rocks.size()<<" ("<<rock_mode<<"), written to rocks.txt"<<endl;
}

//Acceleration on a particle at r from all rocks: the same wall_repulsion() as
//the box walls, along the nearest face's outward normal. Returned in f (the
//caller multiplies by density, like the wall terms).
inline void rock_force(const double r[3], double f[3]){
  f[0] = f[1] = f[2] = 0.0;
  for(size_t ir=0; ir<rocks.size(); ir++){
    double nr[3];
    double frock = wall_repulsion( rock_distance(rocks[ir], r, nr), boundary_r0, boundary_D );
    f[0] += frock*nr[0];
    f[1] += frock*nr[1];
    f[2] += frock*nr[2];
  }
}

//Same safety net as the walls: put a particle that got inside a rock back onto
//the nearest face and remove the velocity component pointing into the rock.
inline void rock_clamp(double r[3], double v[3]){
  for(size_t ir=0; ir<rocks.size(); ir++){
    double nr[3];
    double d = rock_distance(rocks[ir], r, nr);
    if(d < 0.0){
      for(int j=0; j<3; j++) r[j] -= d*nr[j];
      r[0] = wrap_x(r[0]);
      double vn = v[0]*nr[0] + v[1]*nr[1] + v[2]*nr[2];
      if(vn < 0.0)
        for(int j=0; j<3; j++) v[j] -= vn*nr[j];
    }
  }
}

//growth ramp: rock size goes 0 -> 1 over the first rock_ramp_steps steps
inline void rock_grow(int ts){
  rock_scale = (ts < rock_ramp_steps) ? (double)ts/rock_ramp_steps : 1.0;
}

#endif

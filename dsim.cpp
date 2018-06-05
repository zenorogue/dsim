#include <stdio.h>
#define USDL
#ifdef EMSCRIPTEN
#include <emscripten.h>
#else
#define UGFX
#define UGD
#define UTTF
#endif
#include <array>
#include "graph2.h"
#include "mvec.h"

#define RLEIN

#ifdef EMSCRIPTEN
SDL_Surface *emscreen;
#endif

const ld elasticity = .1; // .55;

static const int SX = 800, SY = 800;
static const int ADD = 300;

int toss_results[80];

bool checkpixel = false;

bitmap rmap;

vector<color> sidecolors;

struct shape {
  int sides;
  string filename;
  ld scale;
  };

ld flat;

vector<shape> shapes = {
  {6, "basin6.png", 160.},
  {8, "basin8.png", 120.},
  {12, "basin12.png", 500.},
  {20, "basin20.png", 200.},
  {4, "basin4.png", 160.}
  };

int weighted = 0;

int current_shape;

auto& csh() { return shapes[current_shape]; }

vector<vec3 > corners;

vector<vec3 > weights;

vector<vec3 > asides, sides;

struct edge {
  style s;
  vec3 c1;
  vec3 c2;
  vector<vec3 > frontcheck;
  };

vector<edge> edges;

mat3 rotmatrix;

ld cy = -4;
ld vy = 0;

vec3 mom; 

int mx, my;

style linestyle(0xFFFFFFFF, 0xFFFFFFFF, 1);

font ff;
style tstyle(0xFFFFFFFF, 0xFFFFFFFF, 1);

vec where(vec3 m) {
  // for(int i=0; i<3; i++) printf("%lf ", double(m[i])); printf("=> ");
  m = rotmatrix * m;
  // for(int i=0; i<3; i++) printf("%lf ", double(m[i])); printf("\n");
  vec scr(m[checkpixel ? 2 : 0], m[1]);
  scr.y += cy;
  scr *= 80;
  scr.y += SY - 2;
  scr.x += ADD ? (SX+ADD/2) : SX/2;
  scr.x += 1e-3;
  return scr;
  }

bool infront(vec3 m) {
  int zdir = checkpixel ? 0 : 2;
  return (rotmatrix * m)[zdir] > -1e-3;
  }

void setpix(int x, int y, color col) {
#ifdef EMSCRIPTEN
  *((Uint32*)emscreen->pixels + y * (SX + ADD) + x) = SDL_MapRGBA(screen.s->format, part(col,2), part(col,1), part(col,0), 255);
#else
  screen[y][x] = col;
#endif
  }

void wizualizacja() {

#ifdef EMSCRIPTEN
  if (SDL_MUSTLOCK(emscreen)) SDL_LockSurface(emscreen);    
  #endif
  for(int y=0; y<SY; y++)
  for(int x=0; x<SX; x++)
    setpix(x, y, rmap[y][x]);
  for(int y=0; y<SY; y++)
  for(int x=0; x<ADD; x++)
    setpix(SX+x, y, 0);
    
  for(int f=0; f<2; f++) {

    for(auto e: edges) {
  #ifdef EMSCRIPTEN
      vec a1 = where(e.c1);
      vec a2 = where(e.c2);
      color col = e.s.stroke;
      bool front = false;
      for(auto& f: e.frontcheck) front = front || infront(f);
      if(!front) col = 0x404040;
      if(f == 0 && front) continue;
      if(f == 1 && !front) continue;
      for(int i=0; i<100; i++) {
        int x = int((a1.x * i + a2.x * (100-i))/100);
        int y = int((a1.y * i + a2.y * (100-i))/100);
        setpix(x, y, col);
        }
  #else
      line(e.s, where(e.c1), where(e.c2)).drawSDL(screen);
  #endif
      }
  
    int sid = 0;
    for(auto si: sides) {
      tstyle.fill = sidecolors[sid];
  #ifndef EMSCRIPTEN  
      text(tstyle, ff, where(si), center, 10, its(++sid)).drawSDL(screen);
  #else
      vec v = where(si);
      int x = int(v.x), y = int(v.y);
      if(x > SX+ADD-3) x = SX+ADD-3;
      if(y > SY-3) y = SY-3;
      bool front = infront(si);
      int b = front ? 2 : 1;
      if(f == 0 && front) { sid++; continue; }
      if(f == 1 && !front) { sid++; continue; }
      for(int dy=y-b; dy<=y+b; dy++) for(int dx=x-b; dx<=x+b; dx++)
        setpix(dx, dy, sidecolors[sid]);
      sid++;
  #endif
      }
    }
   
#if EMSCRIPTEN
  if (SDL_MUSTLOCK(emscreen)) SDL_UnlockSurface(emscreen);
  SDL_Flip(emscreen); 
#else
  screen.draw();
#endif
  }

void ttosses(int q); 

int lt;

vector<vec3> lastcorners;

ld inertia;
  
void apply(ld delta) {
  cy += vy * delta;
  vy += delta * 2;
  
  // vy = vy * exp(-delta / 10);
  // for(int z=0; z<3; z++) mom[z] *= exp(-delta / 10);
  
  mat3 aspin;
  
  vec3 mom2 = mom / inertia;
  
  if(inertia == 0) {
    for(int z=0; z<3; z++) {
      auto Z = addrot((z+2)%3, (z+1)%3, M_PI/2, 0);
      
      vec3 M = zero;
      
      for(auto pt: weights) {
        auto rpt = rotmatrix * pt;
        // for(int i=0; i<3; i++) printf("%lf ", (double) rpt[i]); printf(" X ");
        auto rpt2 = Z * rpt;
        // for(int i=0; i<3; i++) printf("%lf ", (double) rpt2[i]); printf("\n");    
        M = M + cross(rpt, rpt2);
        }
  
      for(int i=0; i<3; i++) M[i] /= Size(weights);
      
      // for(int i=0; i<3; i++) printf("%lf ", (double) M[i] + 1e-7); printf(" [%d] # ", z);
      
      for(int i=0; i<3; i++) aspin[z][i] = M[i];
      } // printf(" %d\n", Size(weights));
    
    bool symmetric = true;
    for(int i=0; i<3; i++) for(int j=0; j<3; j++) {
      if(i!=j && abs(aspin[i][j]) > 1e-6) symmetric = false;
      if(i==j && abs(aspin[i][j]-aspin[0][0]) > 1e-6) symmetric = false;
      }
    
    if(symmetric) {
      inertia = aspin[0][0];
      printf("the die is symmetric, inertia = %lf\n", double(inertia));
      }

    mom2 = inverse(aspin) * mom;
    }
  
  // for(int z=0; z<3; z++) printf("%lf ", double(mom[z])); printf("=> ");
  // for(int z=0; z<3; z++) printf("%lf ", double(mom2[z])); printf("\n");

  auto oldrot = rotmatrix;
  
  for(int z=0; z<3; z++) 
    rotmatrix = addrot((z+2)%3, (z+1)%3, mom2[z] * delta) * rotmatrix;
  
  if(false) { // test
    auto Z = (rotmatrix - oldrot) / delta;
    
    // for(int y=0; y<3; y++) { for(int x=0; x<3; x++) printf("%lf ", double(Z[y][x])); } printf("\n");
    
    vec3 M = zero;

    for(auto pt: weights) {
      auto rpt = rotmatrix * pt;
      // for(int i=0; i<3; i++) printf("%lf ", (double) rpt[i]); printf(" X ");
      auto rpt2 = Z * pt;
      // for(int i=0; i<3; i++) printf("%lf ", (double) rpt2[i]); printf("\n");    
      M = M + cross(rpt, rpt2) / Size(weights);
      }

    for(int i=0; i<3; i++) printf("%lf ", (double) M[i] + 1e-7); printf(" vs");
    for(int i=0; i<3; i++) printf("%lf ", (double) mom[i] + 1e-7); printf("\n");
    exit(1);
    }

  vector<vec3> newcorners;
  
  int i = 0;

  for(auto pt: corners) {
    auto pt2 = rotmatrix * pt;
    pt2[1] += cy;

    newcorners.push_back(pt2);
    auto& last = lastcorners[i++];

    if(pt2[1] > 0) {
      
      ld fval = (pt2[1] < last[1]) ? 1000*elasticity : 1000;
      // printf("%lf %lf fval = %lf\n", double(pt2[1]), double(last[1]), double(fval));

      ld force = -pt2[1] * delta * fval;
      vy += force;
      auto fvec = v3(0, force, 0);
      // for(int z=0; z<3; z++) printf("%lf ", double(pt2[z])); printf(" x ");
      // for(int z=0; z<3; z++) printf("%lf ", double(fvec[z])); printf("\n");
      
      mom = mom + cross(pt2, fvec);
      }
    }
  
  lastcorners = newcorners;
  }

bool stopped() {
  if(cy < flat) return false;
  if(abs(vy) > 1e-3) return false;
  for(int i=0; i<3; i++) if(i != 1 && abs(mom[i]) > 1e-3) return false;
  return true;
  }

int result() {
  int i = 0;
  ld least = 100;
  int at = 0;
  for(auto pt: sides) {
    i++;
    auto pt2 = rotmatrix * pt;
    if(pt2[1] < least) least = pt2[1], at = i;
    }
  return at;
  }

int tosses;

void icosa_colors() {
  ld cmin = 10, cmax = -10;
  for(auto& s: sides) cmin = min(cmin, s[0]), cmax = max(cmax, s[0]);
  sidecolors.clear();
  for(auto& s: sides) {
    color res;
    for(int i=0; i<3; i++)
      part(res, i) = 1 + (s[i] - cmin) / (cmax - cmin) * 254;
    sidecolors.push_back(res);
    }
  }

ld edgelength_squared;

ld inerval(vec3 m) { 
  return m[0] * m[0] + m[1] * m[1];
  }

ld uniform_inertia;

void build_edges_sides() {
  edges.clear(); sides.clear();
  edgelength_squared = 100;
  for(auto& c1: corners) for(auto& c2: corners) if(&c1 < &c2) 
    edgelength_squared = min(edgelength_squared, vdist(c1, c2));

  printf("edgelength_squared = %lf\n", double(edgelength_squared));
  
  ld tst = edgelength_squared + 1e-6;

  for(auto& c1: corners) for(auto& c2: corners) if(&c1 < &c2) if(vdist(c1, c2) < tst) {
    edges.push_back(edge{style(0xFFFFFFFF, 0xFFFFFFFF, 1), c1, c2 }); 
    }
  
  for(auto& c1: corners) for(auto& c2: corners) for(auto& c3: corners) if(&c1!=&c2 && &c1!=&c3 && &c2!=&c3) {
    vec3 sum = zero;
    bool bad = false;
    
    int cpf = 0;
    
    for(auto& co: corners) {
      ld v = inface(c1, c2, c3, co);
      if(v < -1e-5) bad = true;
      else if(v < 1e-5) sum = sum + co, cpf++;
      }
    
    if(bad) continue;
    
    sum = sum * ld(1./cpf);
    
    vec3 *curside = NULL;
    
    for(auto& si: sides) if(vdist(si, sum) < 1e-6)
      curside = &si;
    
    if(!curside) { sides.push_back(sum); curside = &sides.back(); }

    for(auto& e: edges)
      if(abs(inface(c1,c2,c3,e.c1)) < 1e-5 && abs(inface(c1,c2,c3,e.c2)) < 1e-5)
        e.frontcheck.push_back(sum);    
    }
  
  sort(sides.begin(), sides.end(), [] (const vec3& c1, const vec3& c2) {
    for(int i=0; i<3; i++)
      if(abs(c1[i] - c2[i]) > 1e-6) return c1[i] < c2[i];
    return false;
    });
  
  asides = sides;
  printf("sides = %d, edges = %d\n", Size(sides), Size(edges));
  
  ld sidecorner = 100;
  for(auto& si: sides) for(auto& c: corners) sidecorner = min(sidecorner, vdist(si, c));
  printf("sidecorner = %lf\n", double(sidecorner));
  
  ld total_inertia = 0;
  ld total_mass = 0;

  for(auto& si: sides) {
    vector<vec3> mycorners;
    for(auto& c: corners) if(vdist(si, c) < sidecorner + 1e-3) mycorners.push_back(c);
    auto M = rot_all(si);
    sort(mycorners.begin(), mycorners.end(), [=] (const vec3& v1, const vec3& v2) {
      auto x1 = M * v1;
      auto x2 = M * v2;
      return atan2(x1[2], x1[0]) < atan2(x2[2], x2[0]);
      });
    /* printf("corners = %d\n", Size(mycorners));
    for(auto& mc: mycorners) {
      for(int i=0; i<3; i++) printf("%10.7lf", double(mc[i])); printf("\n");
      } */
    for(int i=2; i<Size(mycorners); i++) {
      auto c1 = mycorners[0];
      auto c2 = mycorners[i-1];
      auto c3 = mycorners[i];
      ld mass = detv(c1, c2, c3) / 6;
      // printf("det = %lf\n", double(mass));
      total_mass += mass;
      
      int allabc = 0;
      
      int M = 40;
      ld Mi = M;
      
      /* for(int a=0; a<=M; a++) for(int b=0; b<=M; b++) for(int c=0; c<=M; c++) if(a+b+c<=M) allabc++,
        total_inertia += inerval(c1*ld(a/Mi)+cx2*ld(b/Mi)+c3*ld(c/Mi)) * mass / 12341;
      
      printf("allabc = %d\n", allabc); */
        
      total_inertia += (0 + inerval(c1) + inerval(c2) + inerval(c3) + inerval((c1+c2+c3)*ld(1./4)) * 16) * mass / 20;
      }
    }
  
  uniform_inertia = total_inertia / total_mass;
  printf("shape computed\n");
  }

void rotate_to_top() {
  auto M = rot_all(sides[0]);
  for(auto& c: corners) c = M * c;
  for(auto& c: sides) c = M * c;
  for(auto& c: asides) c = M * c;
  for(auto& c: edges) { c.c1 = M * c.c1, c.c2 = M * c.c2;
    for(auto& e: c.frontcheck) e = M * e;
    }
  }

void create_die() {

  corners.clear(); edges.clear(); sides.clear();
  
  printf("creating a %d-sided die\n", csh().sides);
  
  ld wiki_inertia;
  
  switch(csh().sides) {
  
    case 4:
      corners.push_back(v3(1,1,1));
      corners.push_back(v3(1,-1,-1));
      corners.push_back(v3(-1,1,-1));
      corners.push_back(v3(-1,-1,1));
      build_edges_sides();
      sides = corners;
      sidecolors = {0xFFFFFF, 0xFF0000, 0x00FF00, 0x0000FF};
      wiki_inertia = edgelength_squared / 20;
      break;
    
    case 6:

      for(int a=-1; a<=1; a+=2)
      for(int b=-1; b<=1; b+=2)
      for(int c=-1; c<=1; c+=2) 
        corners.push_back(v3(1.0 * a, b, c));
      
      build_edges_sides();
    
      // standard Western Rubik colors
      sidecolors = {
        0xC0C0C0, 0x00FF00, 0xFF0000, 0xFF8000, 0xFFFF00, 0x0000FF
        };
      
      wiki_inertia = edgelength_squared / 6;
      break;
  
    case 8:
      corners.push_back(v3(+2,0,0));
      corners.push_back(v3(-2,0,0));
      corners.push_back(v3(0,+2,0));
      corners.push_back(v3(0,-2,0));
      corners.push_back(v3(0,0,+2));
      corners.push_back(v3(0,0,-2));
      build_edges_sides();
      icosa_colors();
      wiki_inertia = edgelength_squared / 10;
      break;
    
    case 20:
    case 12:
      ld phi = (1 + sqrt(5)) / 2;
      corners.push_back(v3(0,+1,+phi));
      corners.push_back(v3(0,+1,-phi));
      corners.push_back(v3(0,-1,+phi));
      corners.push_back(v3(0,-1,-phi));
      corners.push_back(v3(+1,+phi,0));
      corners.push_back(v3(+1,-phi,0));
      corners.push_back(v3(-1,+phi,0));
      corners.push_back(v3(-1,-phi,0));
      corners.push_back(v3(+phi,0,+1));
      corners.push_back(v3(-phi,0,+1));
      corners.push_back(v3(+phi,0,-1));
      corners.push_back(v3(-phi,0,-1));
      
      build_edges_sides();
      wiki_inertia = phi * phi / 10 * edgelength_squared;
      if(csh().sides == 12) { 
        corners = sides; build_edges_sides(); 
        wiki_inertia = (39 * phi + 28) / 150 * edgelength_squared;
        }
      icosa_colors();
      break;
    
      
    }
  
  rotate_to_top();
  printf("sides[0] is: "); for(int i=0; i<3; i++) printf("%10.7lf", (double)sides[0][i]); printf("\n");
  flat = -sides[0][1]-1e-3;

  printf("inertia: computed = %lf, wiki = %lf\n", double(uniform_inertia), double(wiki_inertia));

  inertia = 0; weights.clear();
  switch(weighted&1) {
    case 0:
      inertia = uniform_inertia;
      printf("weight distribution: uniform\n");
      break;
    
    case 1:
      printf("weight distribution: all in corners\n");
      weights = corners;
      break;
    
    case 2:
      printf("weight distribution: face centers\n");
      for(auto& si: sides) weights.push_back(si);
      break;
    
    case 3:
      printf("weight distribution: edge centers\n");
      for(auto& ed: edges) weights.push_back((ed.c1 + ed.c2) * ld(.5));
      break;
    
    case 4:
      printf("weight distribution: unbalanced\n");
      for(int a=-1; a<=1; a+=2)
      for(int b=-1; b<=1; b+=2)
      for(int c=-1; c<=1; c+=2) 
        weights.push_back(v3(a * .3, b * .3, c));
    
      for(int i=0; i<8; i++) for(int j=0; j<3; j++) {
        if(i & (1<<j)) continue;
        edge ed = { style(0xFF0000FF, 0, 1), weights[i], weights[i ^ (1<<j)] };
        edges.push_back(ed);
        }
      break;
    }    

  for(int i=0; i<Size(sides); i++) toss_results[i] = 0;
  }

void new_toss() {
  checkpixel = false;
  tosses++;
  rotmatrix = unit();
  for(int i=0; i<1000; i++) {
    rotmatrix = addrot(i%3, (i+1)%3, randf() * 2 * M_PI) * rotmatrix;
    }
  
  for(int z=0; z<3; z++) mom[z] = 0;
  mom[0] = randf();
  // mom[1][0] = randf();
  mom[2] = randf();
  
  cy = -4; vy = 0;
  }
  
void may_reroll() {
  if(!checkpixel)
  if(stopped()) {
    int r = result();
    printf(" result = %d", result());
    toss_results[result()]++;
    for(int z=1; z<=Size(sides); z++) printf("%5d", toss_results[z]); printf("\n");
    new_toss();
    }
  }

void phys() {
  int ct = SDL_GetTicks();
  
  while(ct > lt) lt++, apply(.001), apply(.001), apply(.001);
  
  may_reroll();
  }

void ttosses(int q) { 
  int ctosses = tosses;
  while(tosses < ctosses + q) {
    for(int i=0; i<100; i++) apply(1e-3);
    may_reroll();
    }  
  }

void fractal_set(ld x, ld y) {
  checkpixel = true;
  cy = -4; vy = 0;
  rotmatrix = unit();
  mom = zero;
  mom[0] = x / csh().scale;
  mom[2] = y / csh().scale;
  }

#ifndef EMSCRIPTEN
void recolor() {
  vector<color> oldsidecolors = {
    0x0080FF, 0x00FF80, 0x80FF00, 0xFF8000, 0xFF0080, 0x8000FF,
    0x00C0C0, 0xC000C0, 0xC0C000, 0xFF0000, 0x00FF00, 0x0000FF,
    0x202020, 0x404040, 0x808080, 0xC0C0C0, 0xFFFFFF, 0x010101,
    0xFFD500, 0x305070
    };
  icosa_colors();
  rmap = readPng("d20-fractal-oldcolor.png");
  for(int y=0; y<SY; y++)
  for(int x=0; x<SX; x++) {
    rmap[y][x] &= 0xFFFFFF;
    color z = rmap[y][x];
    rmap[y][x] = 0;
    for(int i=0; i<20; i++) if(z == oldsidecolors[i]) rmap[y][x] = sidecolors[i];
    }
  }    
#endif

#ifdef EMSCRIPTEN
template<class A, class B> void writePng(A&,B&) {}
#endif

ld cdist(pair<int, int> p) {
  return hypot(p.first - (SX-1)/2., p.second - (SY-1)/2.);
  }

#ifndef EMSCRIPTEN
void genbasin(int id, bool vis) {
  current_shape = id;
  create_die();
  vector<pair<int, int> > pix;
  for(int y=0; y<SY; y++)
  for(int x=0; x<SX; x++)
    pix.push_back(make_pair(x,y));
    
  sort(pix.begin(), pix.end(), [] (auto p1, auto p2) { return cdist(p1) < cdist(p2); });
  
  int steps = 0;
  
  ll t = getVa();
  
  int allpix = Size(pix);
  
  for(auto [x, y]: pix) {
    steps ++;
    if(rmap[y][x]) continue;
    
    if(vis) if(x != SX/2 && y != SY/2 && x != y && x+y != SX) continue;
    // if((x > SX/2 || x < SX/2-1) && (y > SY/2 || y < SY/2-1)) continue;
    
    fractal_set(x - (SX-1)/2., y - (SY-1)/2.);
    int zsteps = 10000;
    while(!stopped()) { 
      for(int i=0; i<10; i++) apply(1e-3);
      zsteps--;
      if(zsteps < 0) break;
      }
    if(zsteps < 0)
      rmap[y][x] = 0xFFFF9090;
    else
      rmap[y][x] = 0xFF000000 | sidecolors[result()-1];
    
    if(vis && getVa() > t + 1000000) wizualizacja(), t = getVa();
    
    if(steps % (allpix / 100) == 0) {
      writePng(csh().filename, rmap);
      for(int i=0; i<100; i++) printf("%c", i * allpix / 100 < steps ? '#' : '.');
      printf(" %d\n", csh().sides);
      }
    }

  writePng(csh().filename, rmap);
  
  string rle = "rle" + its(csh().sides) + ".txt";
  FILE *f = fopen(rle.c_str(), "wt");
  fprintf(f, "255,%d,", csh().sides);
  int last = -1, qty = 0;
  for(int y=0; y<SY; y++) for(int x=0; x<SX; x++) {
    int c = rmap[y][x] & 0xFFFFFF;
    if(c == last && qty < 255) qty++;
    else {
      for(int a=0; a<Size(sides); a++) if(sidecolors[a] == last) fprintf(f, "%d,%d,", a, qty);
      last = c; qty = 1;
      }
    }
  for(int a=0; a<Size(sides); a++) if(sidecolors[a] == last) fprintf(f, "%d,%d,", a, qty);
  fprintf(f, "\n");
  fclose(f);
  }
#endif

void load_rle() {
#ifdef RLEIN
  unsigned char rlemap[] = {
    #include "rle-all.txt"
    };
  int cur = 0, left = 0, id = 0;
  while(!(rlemap[id] == 255 && rlemap[id+1] == csh().sides)) id += 2;
  id += 2;
  for(int y=0; y<SY; y++) for(int x=0; x<SX; x++) {
    if(left == 0) { cur = sidecolors[rlemap[id++]]; left = rlemap[id++]; }
    left--;
    rmap[y][x] = cur;
    }
#endif
  }

void klawisze() {
  SDL_Event event;
#ifndef EMSCRIPTEN
  SDL_Delay(1);
#endif
  bool ev;
  while((ev = SDL_PollEvent(&event))) switch (event.type) {
    case SDL_QUIT:
      exit(1);
      return;

    case SDL_MOUSEBUTTONDOWN: {
      if(mx < SX)
        fractal_set(mx - (SX-1)/2., my - (SY-1)/2.);
      else {
        current_shape = (1 + current_shape) % Size(shapes); create_die(); new_toss(); load_rle();
        }
      break;
      }
    
    case SDL_MOUSEMOTION: {
      mx = event.motion.x;
      my = event.motion.y;
      break;
      }
    
    case SDL_KEYDOWN: {
      int key = event.key.keysym.sym;
      int uni = event.key.keysym.unicode;
      
      printf("uni = %c key = %d\n", uni, key);
      if(uni == 'q') exit(1);
      
      if(uni == 'n') new_toss();
      
      if(uni == '1') ttosses(100);
      if(uni == '2') ttosses(600);
      if(uni == '3') ttosses(1000);
      if(uni == '4') ttosses(6000);
      
      if(uni == 'c') { current_shape = (1 + current_shape) % Size(shapes); create_die(); new_toss(); load_rle(); }
      if(uni == 'w') { weighted = (1 + weighted) % 5; create_die(); new_toss(); }
            
#ifndef EMSCRIPTEN      
      if(uni == 'f') genbasin(current_shape, true);
      if(uni == 'l') {
        rmap = readPng(csh().filename);
        for(int y=0; y<SY; y++)
        for(int x=0; x<SX; x++)
          rmap[y][x] |= 0xFF000000;
        }
#endif
      
      break;
      }
    
    }
  }

void mainloop() {
  wizualizacja();
  klawisze();
  phys();
  }

int main(int argc, char **argv) {
   
  current_shape = 3;
  create_die();
  rmap = emptyBitmap(SX, SY);
  lastcorners = corners;

#ifndef EMSCRIPTEN  
  if(argc >= 3 && argv[1] == string("-gb")) {
    genbasin(atoi(argv[2]), false);
    exit(0);
    }
#endif

  load_rle();

  srand(time(NULL));

#ifdef EMSCRIPTEN
  SDL_Init(SDL_INIT_VIDEO);
  emscreen = SDL_SetVideoMode(SX+ADD, SY, 32, SDL_SWSURFACE);
#else
  initGraph(SX+ADD, SY, "dsim", SDL_SWSURFACE);
#endif

  ff = makefont("DejaVuSans-Bold.ttf");
  
  new_toss();
  
  #ifdef EMSCRIPTEN
  emscripten_set_main_loop(mainloop, 0, false);
  #else

  while(true) {
    wizualizacja();
    klawisze();
    phys();
    }  
  #endif
  
  return 0;
  }


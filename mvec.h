typedef long long ll;

typedef array<ld, 3> vec3;
typedef array<array<ld, 3>, 3> mat3;

vec3 v3(ld x, ld y, ld z) {
  vec3 ret;
  ret[0] = x;
  ret[1] = y;
  ret[2] = z;
  return ret;
  }


vec3 zero = v3(0, 0, 0);

mat3 unit() {
  mat3 res;
  for(int i=0; i<3; i++) for(int k=0; k<3; k++) res[i][k] = (i == k);
  return res;
  }

vec3 operator + (const vec3& a, const vec3& b) { 
  return v3(a[0]+b[0], a[1]+b[1], a[2]+b[2]);
  }

vec3 operator - (const vec3& a, const vec3& b) { 
  return v3(a[0]-b[0], a[1]-b[1], a[2]-b[2]);
  }

vec3 operator * (const vec3& a, ld z) { 
  return v3(a[0]*z, a[1]*z, a[2]*z);
  }

vec3 operator / (const vec3& a, ld z) { 
  return v3(a[0]/z, a[1]/z, a[2]/z);
  }

mat3 operator - (const mat3& a, const mat3& b) {
  mat3 res;
  for(int i=0; i<3; i++)
  for(int j=0; j<3; j++)
    res[i][j] = a[i][j] - b[i][j];
  return res;
  }

mat3 operator * (const mat3& a, ld x) {
  mat3 res;
  for(int i=0; i<3; i++)
  for(int j=0; j<3; j++)
    res[i][j] = a[i][j] * x;
  return res;
  }

mat3 operator / (const mat3& a, ld x) {
  mat3 res;
  for(int i=0; i<3; i++)
  for(int j=0; j<3; j++)
    res[i][j] = a[i][j] / x;
  return res;
  }

vec3 operator * (const mat3& m, const vec3& a) {
  vec3 res;
  for(int i=0; i<3; i++) {
    ld r = 0;
    for(int j=0; j<3; j++) r += m[i][j] * a[j];
    res[i] = r;
    }
  return res;
  }

mat3 operator * (const mat3& a, const mat3& b) {
  mat3 res;
  for(int i=0; i<3; i++) for(int k=0; k<3; k++) {
    ld r = 0;
    for(int j=0; j<3; j++) r += a[i][j] * b[j][k];
    res[i][k] = r;
    }
  return res;
  }

vec3 cross(vec3 v1, vec3 v2) {
  vec3 ret;
  ret[0] = v1[1] * v2[2] - v1[2] * v2[1];
  ret[1] = v1[2] * v2[0] - v1[0] * v2[2];
  ret[2] = v1[0] * v2[1] - v1[1] * v2[0];
  return ret;
  }

ld sca(vec3 v1, vec3 v2) {
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
  }

double randf() {
  return (rand() % 1000000000) / 1000000000.;
  }

mat3 addrot(int i, int j, ld r, ld rem=1) {
  mat3 res;
  for(int i=0; i<3; i++) for(int j=0; j<3; j++) res[i][j] = 0;
  res[i][i] = cos(r);
  res[j][j] = cos(r);
  res[i][j] = sin(r);
  res[j][i] = -sin(r);
  res[3-i-j][3-i-j] = rem;
  return res;
  }

mat3 rot_from_to(int i, int j, vec3 what) {
  mat3 res;
  for(int i=0; i<3; i++) for(int j=0; j<3; j++) res[i][j] = i==j;
  ld x = what[i];
  ld y = what[j];
  if(abs(x) < 1e-6 && abs(y) < 1e-6) return res;
  ld z = sqrt(x*x+y*y);
  res[i][i] = y / z;
  res[i][j] = -x / z;
  res[j][i] = x / z;
  res[j][j] = y / z;
  return res;
  }

ld det(const mat3& T) {
  ld det = 0;
  for(int i=0; i<3; i++) 
    det += T[0][i] * T[1][(i+1)%3] * T[2][(i+2)%3];
  for(int i=0; i<3; i++) 
    det -= T[0][i] * T[1][(i+2)%3] * T[2][(i+1)%3];
  return det;
  }
  
mat3 inverse(const mat3& T) {
  
  ld d = det(T);
  mat3 T2;
  
  for(int i=0; i<3; i++) 
  for(int j=0; j<3; j++) 
    T2[j][i] = (T[(i+1)%3][(j+1)%3] * T[(i+2)%3][(j+2)%3] - T[(i+1)%3][(j+2)%3] * T[(i+2)%3][(j+1)%3]) / d;

  return T2;
  }

ld vdist(const vec3& c1, const vec3& c2) {
  auto s = c1 - c2;
  return sca(s, s);
  }

ld detv(const vec3& a, const vec3& b, const vec3& c) {
  mat3 z;
  for(int i=0; i<3; i++) {
    z[i][0] = a[i];
    z[i][1] = b[i];
    z[i][2] = c[i];
    }
  return det(z);
  }

ld inface(vec3& c1, vec3& c2, vec3& c3, vec3& co) {
  return detv(c2-c1, c3-c1, co-c1);
  }

mat3 rot_all(vec3 top) {
  auto M = rot_from_to(2, 0, top);
  top = M * top;
  auto M1 = rot_from_to(0, 1, top);
  return M1 * M;
  }

template<class T> int size(const T& t) { return t.size(); }

string its(int i) { return to_string(i); }

#include <sys/time.h>
ll getVa() {
  struct timeval tval;
  gettimeofday(&tval, NULL);
  return tval.tv_sec * 1000000 + tval.tv_usec;
  }

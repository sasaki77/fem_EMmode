#ifndef CALCEF
#define CALCEF

#include "common.h"

//============================================================
// 関数名：getDet()
// 引数　：調査する点(x0,y0),線分((始点)(x1,y1),(終点)(x2,y2))
// 戻り値：行列式の計算値
//---概要------------------------------------------------------
// 調査する点が線分の左右にあるか判断するための行列式を返す．
// 調査する点は戻り値が0ならば線分上，0より大きいなら左側，
// ０より小さいなら右側に位置する
//============================================================
double getDet(double x0,double y0,double x1,double y1,double x2,double y2)
{
  return ((x1-x0)*(y2-y0) - (y1-y0)*(x2-x0));
}

bool isInclude(Node p0,Node p1,Node p2, Node p3){
  Node nt[3] = {p1,p2,p3};
  bool isleft = true;
  for(int i=0;i<3;i++){
    double det = getDet( p0.z, p0.r, nt[i].z, nt[i].r, nt[(i+1)%3].z, nt[(i+1)%3].r );
    if(det<0){
      isleft = false;
      break;
    }
  }
  return isleft;
}

double Element::N(int i,double xi,double eta)
{
  double L1 = xi;
  double L2 = eta;
  double L3 = 1-xi-eta;
  switch(i){
  case 1: return L1*(2*L1 - 1);
  case 2: return L2*(2*L2 - 1);
  case 3: return L3*(2*L3 - 1);
  case 4: return 4*L1*L2;
  case 5: return 4*L2*L3;
  case 6: return 4*L3*L1;
  }
  return -1;
}

double Element::N_xi(int i,double xi,double eta)
{
  double L1 = xi;
  double L2 = eta;
  double L3 = 1-xi-eta;
  switch(i){
  case 1: return 4*L1 - 1;
  case 2: return 0;
  case 3: return -4*L3 + 1;
  case 4: return 4*L2;
  case 5: return -4*L2;
  case 6: return 4*L3 - 4*L1;
  }
  return -1;
}

double Element::N_eta(int i,double xi,double eta)
{
  double L1 = xi;
  double L2 = eta;
  double L3 = 1-xi-eta;
  switch(i){
  case 1: return 0;
  case 2: return 4*L2-1;
  case 3: return -4*L3 + 1;
  case 4: return 4*L1;
  case 5: return 4*L3 - 4*L2;
  case 6: return -4*L1;
  }
  return -1;
}

double Element::z_xi(double xi,double eta)
{
  double sum = 0;
  for(int i=0;i<6;i++)
    sum += p[i]->z*N_xi(i+1,xi,eta);
  return sum;
}
double Element::z_eta(double xi,double eta)
{
  double sum = 0;
  for(int i=0;i<6;i++)
    sum += p[i]->z*N_eta(i+1,xi,eta);
  return sum;
}
double Element::r_xi(double xi,double eta)
{
  double sum = 0;
  for(int i=0;i<6;i++)
    sum += p[i]->r*N_xi(i+1,xi,eta);
  return sum;
}
double Element::r_eta(double xi,double eta)
{
  double sum = 0;
  for(int i=0;i<6;i++)
    sum += p[i]->r*N_eta(i+1,xi,eta);
return sum;
}
double Element::jcb(double xi,double eta)
{
  return ( z_xi(xi,eta)*r_eta(xi,eta) - z_eta(xi,eta)*r_xi(xi,eta) );
}
double Element::N_z(int i,double xi,double eta)
{
  return (( r_eta(xi,eta)*N_xi(i,xi,eta) - r_xi(xi,eta)*N_eta(i,xi,eta))/jcb(xi,eta));
}
double Element::N_r(int i,double xi,double eta)
{
  return ((-z_eta(xi,eta)*N_xi(i,xi,eta) + z_xi(xi,eta)*N_eta(i,xi,eta))/jcb(xi,eta));
}
void Params::getElectroVector()
{
  double maxz,maxr;
  double minz,minr;
  double rz,rr;
  double square;
  maxz = maxr = DBL_MIN;
  minz = minr = DBL_MAX;

  for(int i=0;i<node.size();i++){
    if(minz > node[i].z) minz = node[i].z;
    if(minr > node[i].r) minr = node[i].r;
    if(maxz < node[i].z) maxz = node[i].z;
    if(maxr < node[i].r) maxr = node[i].r;
  }
  rz = maxz-minz;
  rr = maxr-minr;
  square = (rz>rr)?rz:rr;

  double d = square / ev_dev;
  double margin = square/100;
  for(double z = minz + margin; z<maxz; z+=d){
    for(double r = minr + margin; r<maxr; r+=d){
      Node p;
      p.set(z,r);
      for(int i=0;i<elem.size();i++){
        if( isInclude(p,*(elem[i].p[0]),*(elem[i].p[1]),*(elem[i].p[2])) ){
          Element e = elem[i];
          for(int j=0;j<6;j++){
            double xi,eta;
            switch(j+1){
            case 1: xi = 1.0; eta = 0.0; break;
            case 2: xi = 0.0; eta = 1.0; break;
            case 3: xi = 0.0; eta = 0.0; break;
            case 4: xi = 0.5; eta = 0.5; break;
            case 5: xi = 0.0; eta = 0.5; break;
            case 6: xi = 0.5; eta = 0.0; break;
            }
            for(int k=0;k<6;k++){
              p.Er -= e.p[k]->u * e.N_z(k+1,xi,eta)/(EPS_0*omega);
              p.Ez += e.p[k]->u * e.N(k+1,xi,eta)/p.r/(EPS_0*omega)
                     +e.p[k]->u * e.N_r(k+1,xi,eta)/(EPS_0*omega);
            }
          }
          ev.push_back(p);
        }
      }
    }
  }    
}

#endif

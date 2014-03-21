#ifndef CGM_CPP
#define CGM_CPP
// 配列でcgmを実装

#include <cstdlib>
#include <ctime>
#include <queue>

#include "common.h"

// ==================== 外部変数 ====================
// 以下の外部変数の初期化は FEM::cgm() で行われます．
double **a2;
double **b2;
double *u2;
int n_size;
// =================================================

void cV(double c, double *V)
{
  for(int i=0;i<n_size;i++)
    V[i] *= c;
}
void cM(double c, double **M)
{
  for(int i=0;i<n_size;i++)
    for(int j=0;j<n_size;j++)
      M[i][j] *= c;
}
void add(double *v1, double *v2, double *v3)
{
  for(int i=0;i<n_size;i++)
    v3[i] = v1[i] + v2[i];
}
void sub(double *v1, double *v2, double *v3)
{
  for(int i=0;i<n_size;i++)
    v3[i] = v1[i] - v2[i];
}
void MxV(double **M,double *V,double *MV)
{  
  for(int i=0;i<n_size;i++){
    MV[i] = 0;
    for(int j=0;j<n_size;j++)
      MV[i] += M[i][j] * V[j];
  }
}
double dot(double *v,double *w)
{
  double sum = 0;
  for(int i=0;i<n_size;i++) sum += v[i] * w[i];
  return sum;
}
double cgmdot(double *v,double **m, double *w)
{
  double *z = (double *)calloc(n_size,sizeof(double));
  MxV(m,w,z);
  double tmp = dot(v,z);
  free(z);
  return tmp;
}
void copy(double *v1,double *v2)
{
  for(int i=0;i<n_size;i++)
    v2[i] = v1[i];
}
double solve(double a,double b,double c)
{
  return ( -b + sqrt(b*b-4*a*c))/(2*a);
}
double ray(double *u)
{
  return( cgmdot(u,a2,u)/cgmdot(u,b2,u) );
}
void calc_r(double *u,double *got_r)
{
  double *au,*bu;
  double lambda = ray(u);
  au = (double *)calloc(n_size,sizeof(double));
  bu = (double *)calloc(n_size,sizeof(double));

  MxV(a2,u,au);
  MxV(b2,u,bu);
  cV(lambda,bu);
  sub(au,bu,got_r);
  cV(2/dot(u,bu),got_r); // <-- calculation of r is finished

  free(au);
  free(bu);
}

double cgm_main()
{
  double lambda;
  double *r_old,*r_new;
  double *u2_new,*p;
  double alpha,beta;
  double eps = 1e-16; // 収束判定子
  double ei; // 相対残差
  double *lx,*ax,*la,*lbx;
  int recalc_cnt = 0;
  bool recalc_flg = true;

  // ---------- メモリ確保 ----------
  r_old = (double *)calloc(n_size,sizeof(double));
  r_new = (double *)calloc(n_size,sizeof(double));
  u2_new = (double *)calloc(n_size,sizeof(double));
  p = (double *)calloc(n_size,sizeof(double));
  lx = (double *)calloc(n_size,sizeof(double));
  ax = (double *)calloc(n_size,sizeof(double));
  la = (double *)calloc(n_size,sizeof(double));
  lbx = (double *)calloc(n_size,sizeof(double));
  // ------------------------------

  srand((unsigned)time(NULL));

  while(recalc_flg){
    recalc_flg = false;
    recalc_cnt = 0;
    
    for(int i=0;i<n_size;i++)
      u2[i] = rand();
    
    calc_r(u2,r_old);
    copy(r_old,p);
    
    //cgm
    do{
      recalc_cnt++;
      double pap = cgmdot(p,a2,p);
      double pau = cgmdot(p,a2,u2);
      double uau = cgmdot(u2,a2,u2);
      double pbp = cgmdot(p,b2,p);
      double pbu = cgmdot(p,b2,u2);
      double ubu = cgmdot(u2,b2,u2);

      double c = (pap*pbu) - (pbp*pau);
      double d = (pap*ubu) - (pbp*uau);
      double e = (pau*ubu) - (pbu*uau);

      alpha = solve(c,d,e);
      
      if(isnan(alpha)){
	cout << "alpha IS NOT A NUMBER" << endl;
	exit(0);
      }

      cV(alpha,p);
      add(u2,p,u2_new);
      lambda = ray(u2_new);
      calc_r(u2_new,r_new);
      beta = dot(r_new,r_new)/dot(r_old,r_old);
      cV(beta/alpha,p);
      add(r_new,p,p);
      
      // 誤差評価
      copy(u2_new,lx);
      cV(lambda,lx);
      MxV(b2,lx,lbx);
      MxV(a2,u2_new,ax);
      sub(lbx,ax,la);
      ei = dot(la,la)/dot(lbx,lbx);
      
      cout << "prec = " << ei << endl;

      copy(r_new,r_old);
      copy(u2_new,u2);

    }while( ei > eps );
  }

  // ---------- メモリ解放 ----------
  free(r_old);
  free(r_new);
  free(u2_new);
  free(p);
  free(lx);free(ax);free(la);free(lbx);
  // ------------------------------

  return lambda;
}

void FEM::setbc()
{
  int line = 0,colum = 0;
  
  for(int i=0;i<n_size;i++){
    while( line < node.size() && node[line].known ) line++;
    for(int j=0;j<n_size;j++){
      while( colum < node.size() && node[colum].known ) colum++;
      a2[i][j] = A[line][colum];
      b2[i][j] = B[line][colum];
      colum++;
    }
    colum = 0;
    line++;
  }
}

void FEM::cgm()
{
  int n_known = 0;
  for(int i=0;i<node.size();i++)
    if(node[i].known)
      n_known++;

  n_size = node.size() - n_known;
  
  // メモリ確保
  a2 = (double **)calloc(n_size,sizeof(double*));
  for(int i=0;i<n_size;i++)
    a2[i] = (double *)calloc(n_size,sizeof(double));
  b2 = (double **)calloc(n_size,sizeof(double*));
  for(int i=0;i<n_size;i++)
    b2[i] = (double *)calloc(n_size,sizeof(double));
  u2 = (double *)calloc(n_size,sizeof(double));

  setbc();
  lambda = cgm_main();
  omega = LIGHT_SPEED* sqrt(lambda);

  // 固有ベクトルの更新
  int k=0;
  for(int i=0;i<node.size();i++){
    if(!node[i].known){
      u[i] = u2[k];
      k++;
    }
  }  

  // メモリ解放
  for(int i=0;i<n_size;i++)
    free(a2[i]);
  free(a2);
  for(int i=0;i<n_size;i++)
    free(b2[i]);
  free(b2);
  free(u2);
}

#endif // CGM_CPP

#ifndef IO_CPP
#define IO_CPP

#include "common.h"

// ==================== Maxwell's eq. より電場の計算 ====================

void FEM::calcEField()
{
  for(int i=0;i<elem.size();i++)
    elem[i]->calcEfield(&u,omega);
}

// ==================== 入射電力と磁場分布より実際の磁場の値を計算 ====================
void FEM::calcMField()
{
  cout << "input power [J] >> " ;
  cin >> P_in;
 
  double sum = 0;
  for(int i=0;i<elem.size();i++){
    sum += elem[i]->calcPsum(&u);
  }
  cout << "H^2 r dr dz = " << sum << endl;
  double normalizer = sqrt( P_in / (MU0*M_PI*sum) );
  cout << "normalizer = " << normalizer << endl;
  for(int i=0;i<u.size();i++) u[i] = normalizer * u[i];
}

// ==================== 線積分により壁面損失を計算 ====================
void FEM::calcPowerLoss()
{
  double _rs;
  P_wall = 0;

  for(int i=0;i<mt_edge.size();i++)
    P_wall += mt_edge[i]->calcPowerLoss(&u,omega);
  
  // 線積分の確かめ
  cout << "liner integration = " << P_wall << endl;
  cout << "surface resist = " << mt_edge[0]->getSurfResist(omega) << endl;

  // Q値の計算
  Q = omega * P_in / P_wall;

  // シャントインピーダンスの計算
  // Z軸上の平均電場を得る
  int n=0;
  E0=0;
  for(int i=0;i<node.size();i++){
    if( node[i].r < CALC_EPS ){
      n++;
      E0 += sqrt(node[i].Er*node[i].Er + 
		  node[i].Ez*node[i].Ez);
    } 
  }
  E0 /= n;
  R = LenZ*E0*E0/P_wall;
}

// ================================================================

#endif

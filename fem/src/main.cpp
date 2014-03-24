#ifndef MAIN_CPP
#define MAIN_CPP

#include "common.h"

FEM fem;

int main( int argv, char* argc[] )
{
  cout << "==================================================" << endl;
  cout << " fem " << endl;
  cout << "--------------------------------------------------" << endl;

  bool useCalculatedFile = false;
  for(int i=0;i<argv;i++)
    if((string)argc[i] == "-set")
      useCalculatedFile = true;

  string ifname = "";
  if( argv >= 2 && argc[1][0] != '-' ) ifname = (string)argc[1];

  cout << "fem.input()..." << endl;
  fem.input(ifname);
  cout << "end." << endl << endl;

  clock_t start_c=0,end_c=0;
  time_t start_t, end_t;

  if(useCalculatedFile){
    // 計算済みのファイルを使う
    cout << "fem.inpure_re()..." << endl;
    fem.input_re(ifname);
    cout << "end." << endl << endl ;
  }else{
    
    // 時間計測と記録
    // 1秒以下ならclock_tで，それ以上の長時間ならtime_tで時間を計測する
    start_c = clock();
    time(&start_t);

    cout << "fem.setLargeMatrix()..." << endl;
    fem.setLargeMatrix();
    cout << "end." << endl << endl ;

    cout << "fem.cgm()..." << endl;
    fem.cgm();
    cout << "end." << endl << endl;

    end_c = clock();
    time(&end_t);
  }

  cout << "fem.output()..." << endl;
  fem.calcMField();
  fem.calcEField();
  fem.calcPowerLoss();

  double timer;
  if( (timer = (double)(end_c-start_c)/CLOCKS_PER_SEC) > 1.0 )
    timer = difftime(end_t,start_t);

  fem.output(ifname,timer);
  fem.outputEonRAxis(ifname);
  cout << "end." << endl << endl;

  cout << "----------------------------------------" << endl;
  cout << " nodenum    = " << fem.node.size() << endl;
  cout << " edge size  = " << fem.getAveragedLength() << endl;
  cout << " area size  = " << fem.getAveragedArea() << endl;
  cout << " lambda     = " << fem.lambda << endl;
  cout << " frequency  = " << fem.omega/(2*M_PI) << endl;
  cout << " power loss = " << fem.P_wall << endl;
  cout << " Q          = " << fem.Q << endl;
  cout << " E0         = " << fem.E0 << endl;
  cout << " R          = " << fem.R << endl;
  cout << "----------------------------------------" << endl;

//   double timer;
//   if( (timer = (double)(end_c-start_c)/CLOCKS_PER_SEC) > 1.0 )
//     timer = difftime(end_t,start_t);
  //fem.memo(timer);

  cout << "==================================================" << endl;

  return 0;
}

#endif

//----------------------------------------------------------------------------------||
//-------------------                commute.h                   -------------------||
//----------------------------------------------------------------------------------||
//                                                                                  ||
//               __        ___    _   _ ____        ____ ___ ____                   ||
//               \ \      / / \  | \ | |  _ \      |  _ \_ _/ ___|                  ||
//                \ \ /\ / / _ \ |  \| | | | |_____| |_) | | |                      ||
//                 \ V  V / ___ \| |\  | |_| |_____|  __/| | |___                   ||
//                  \_/\_/_/   \_\_| \_|____/      |_|  |___\____|                  ||
//                                                                                  ||
//----------------------------------------------------------------------------------||
//--  (W)akefield (A)cceleration a(n)d (D)LA - (P)article (i)n (C)ell Simulation  --||
//----------------------------------------------------------------------------------||
//---Author-----------           : Tianhong Wang                --------------------||
//---Starting---------           : Feb-05-2019                  --------------------||
//---Email------------           : tw474@cornell.edu            --------------------||
//---Group------------           : Dr. Gennady Shvets' Group    --------------------||
//---Copyright--------           : (C) 2019 by Tianhong Wang    --------------------||
//----------------------------------------------------------------------------------||
//----------------------------------------------------------------------------------||



#ifndef H_COMM
#define H_COMM

#include<vector>

//---------------------------- Mesh class -----------------------
class Commute
   {
   friend class Domain;
   friend class Mesh;

private:
   Domain *p_domain() {return Domain::p_Domain;}; // pointer to domain class.


   int bufsize;

   double *SendSourceXm;  //Array to send for the sources;
   double *SendSourceYm;  //Array to send for the sources;
   double *SendSourceXp;  //Array to send for the sources;
   double *SendSourceYp;  //Array to send for the sources;

   double *ReceSourceXm;  //Array for Rece 
   double *ReceSourceYm;  //Array for Rece 
   double *ReceSourceXp;  //Array for Rece 
   double *ReceSourceYp;  //Array for Rece 



   // diagonal
   double *SendSourcemm;  //Array to send for the sources;
   double *SendSourcemp;  //Array to send for the sources;
   double *SendSourcepm;  //Array to send for the sources;
   double *SendSourcepp;  //Array to send for the sources;

   double *ReceSourcemm;  //Array for Rece 
   double *ReceSourcemp;  //Array for Rece 
   double *ReceSourcepm;  //Array for Rece 
   double *ReceSourcepp;  //Array for Rece 

   //double *SendFields;  //Array to send for the fields;

   int SendSouSizeX;
   int SendSouSizeY;
   int SendFieSizeX;

   int Rank;

   int RankIdx_X;
   int RankIdx_Y;

   int Xpa;
   int Ypa;

   int XmPE; 
   int XpPE;
   int YmPE;
   int YpPE; 

   int mmPE; 
   int mpPE;
   int pmPE;
   int ppPE; 

   int GridX;
   int GridY;

   int kold;

   std::vector<double> CellAccX; 
   std::vector<double> CellAccY; 

public:

   
   void DoCommuteT(exchange what, std::vector<int> SendN);
   void    UnPackT(exchange what, std::vector<int> ReceN);

   void DoCommute(exchange what, int k);
   void    DoPack(exchange what, int k);
   void    UnPack(exchange what, int k);

   int Get_bufsize() {return bufsize;};



   Commute(int XGridN, int YGridN);
   ~Commute();


};

#endif

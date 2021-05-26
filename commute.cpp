//----------------------------------------------------------------------------------||
//-------------------                commute.cpp                 -------------------||
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




#include "wand_PIC.h"



Commute::Commute(int XGridN, int YGridN)
{

	bufsize = ceil((p_domain()->Get_Buffersize())
		*p_domain()->p_Partition()->GetXpart());

	SendSourceXm = NULL;
	SendSourceXp = NULL;
	SendSourceYm = NULL;
	SendSourceYp = NULL;

	ReceSourceXm = NULL;
	ReceSourceXp = NULL;
	ReceSourceYm = NULL;
	ReceSourceYp = NULL;

	Rank = p_domain()->p_Partition()->Get_Rank();
	//Big Coordinates of the Rank
	RankIdx_X = p_domain()->p_Partition()->RankIdx_X(); 
	RankIdx_Y = p_domain()->p_Partition()->RankIdx_Y(); 

	//Defind the size of SendArray
	Xpa  = p_domain()->p_Partition()->GetXpart();
	Ypa  = p_domain()->p_Partition()->GetYpart();

	XmPE = p_domain()->p_Partition()->GetXmPE(); 
	XpPE = p_domain()->p_Partition()->GetXpPE(); 
	YmPE = p_domain()->p_Partition()->GetYmPE(); 
	YpPE = p_domain()->p_Partition()->GetYpPE(); 

	mmPE = p_domain()->p_Partition()->GetmmPE(); 
	mpPE = p_domain()->p_Partition()->GetmpPE(); 
	pmPE = p_domain()->p_Partition()->GetpmPE(); 
	ppPE = p_domain()->p_Partition()->GetppPE(); 

	GridX = XGridN;
	GridY = YGridN;

	kold=-10;

	if( (p_domain()->Get_Nbeam()) > 0 )
	{ 
		SendSouSizeX = YGridN * (SOU_DIM +BEA_DIM ) * 2;   //send in X directions: left and right
		SendSouSizeY = XGridN * (SOU_DIM +BEA_DIM ) * 2;   //send in Y directions: up   and down
	}
	else
	{ 
		//SOU_DIM = denn jx jy jxx jyy jxy
		SendSouSizeX = YGridN * SOU_DIM * 2;   //send in X directions: left and right
		SendSouSizeY = XGridN * SOU_DIM * 2;   //send in Y directions: up   and down	
	}

	SendSourceXm = new double[SendSouSizeX*bufsize];
	SendSourceXp = new double[SendSouSizeX*bufsize];

	ReceSourceXm = new double[SendSouSizeX*bufsize];
	ReceSourceXp = new double[SendSouSizeX*bufsize];

	SendSourceYm = new double[SendSouSizeY*bufsize];
	SendSourceYp = new double[SendSouSizeY*bufsize];

	ReceSourceYm = new double[SendSouSizeY*bufsize];
	ReceSourceYp = new double[SendSouSizeY*bufsize];

	// diagonal direction
	SendSourcemm = new double[SendSouSizeX/YGridN*bufsize];
	SendSourcemp = new double[SendSouSizeX/YGridN*bufsize];

	ReceSourcemm = new double[SendSouSizeX/YGridN*bufsize];
	ReceSourcemp = new double[SendSouSizeX/YGridN*bufsize];

	SendSourcepm = new double[SendSouSizeY/XGridN*bufsize];
	SendSourcepp = new double[SendSouSizeY/XGridN*bufsize];

	ReceSourcepm = new double[SendSouSizeY/XGridN*bufsize];
	ReceSourcepp = new double[SendSouSizeY/XGridN*bufsize];


	//Cell Position Accumulative
	CellAccX = std::vector<double> (XGridN+3,0.0);
	CellAccY = std::vector<double> (YGridN+3,0.0);

	for(int i=0;i<XGridN+2;i++)
	{	
		Cell &ccc = p_domain()->p_Mesh()->GetCell(i,0,0);
		CellAccX[i]= ccc.Xcord-ccc.dx*0.5;
	}
	Cell &c1 = p_domain()->p_Mesh()->GetCell(XGridN+1,0,0);
	CellAccX[XGridN+2]=c1.Xcord+c1.dx*0.5;

	for(int i=0;i<YGridN+2;i++)
	{	
		Cell &ccc = p_domain()->p_Mesh()->GetCell(0,i,0);
		CellAccY[i]= ccc.Ycord-ccc.dy*0.5;
	}
	Cell &c2 = p_domain()->p_Mesh()->GetCell(0,YGridN+1,0);
	CellAccY[YGridN+2]=c2.Ycord+c2.dy*0.5;


};

void Commute::DoCommute(int what, int k)
{


	//===============================================================
	//=====================           Pack         ==================
	//===============================================================
	//pack the fields or source all together in order to send.
	DoPack(what, k);
	int ssx;
	int ssy;

	int ssxd;
	int ssyd;

	MultiGrid *p_Multi = NULL;

	switch(what)
	{

		case COMMU_S:
		case COMMU_SO:
			ssx = SendSouSizeX;
			ssy = SendSouSizeY;
			ssxd= ssx/GridY;
			ssyd= ssy/GridX;
		break;

		case COMMU_F:
			ssx = GridY * WAK_DIM2;
			ssy = GridX * WAK_DIM2;
			ssxd= WAK_DIM2;
			ssyd= WAK_DIM2;
		break;

		case COMMU_A:
			ssx = GridY * 4;
			ssy = GridX * 4;
		break;

		case COMMU_MG_P:
		case COMMU_MG_R:
			p_Multi = p_domain()->p_MG();
			ssx = p_Multi->GetLayerGridY(k);
			ssy = p_Multi->GetLayerGridX(k);
			ssxd= 1;
			ssyd= 1;
		break;

		case COMMU_MG_P_C:
		case COMMU_MG_R_C:
		p_Multi = p_domain()->p_MG();
			ssx = p_Multi->GetLayerGridY(k)*2;
			ssy = p_Multi->GetLayerGridX(k)*2;
			ssxd= 1;
			ssyd= 1;
		break;
	}
		//===============================================================
		//=====================   ISend and IRecev     ==================
		//===============================================================
    	MPI_Request Request[16];
    	MPI_Status 	 Status[16];

 		MPI_Irecv(ReceSourceXm, ssx,  MPI_DOUBLE, XmPE, 0, MPI_COMM_WORLD, &Request[0]);
		MPI_Irecv(ReceSourceXp, ssx,  MPI_DOUBLE, XpPE, 1, MPI_COMM_WORLD, &Request[1]);
		MPI_Irecv(ReceSourceYm, ssy,  MPI_DOUBLE, YmPE, 2, MPI_COMM_WORLD, &Request[2]);
		MPI_Irecv(ReceSourceYp, ssy,  MPI_DOUBLE, YpPE, 3, MPI_COMM_WORLD, &Request[3]);

		MPI_Irecv(ReceSourcemm, ssxd, MPI_DOUBLE, mmPE, 4, MPI_COMM_WORLD, &Request[4]);
		MPI_Irecv(ReceSourcemp, ssxd, MPI_DOUBLE, mpPE, 5, MPI_COMM_WORLD, &Request[5]);
		MPI_Irecv(ReceSourcepm, ssyd, MPI_DOUBLE, pmPE, 6, MPI_COMM_WORLD, &Request[6]);
		MPI_Irecv(ReceSourcepp, ssyd, MPI_DOUBLE, ppPE, 7, MPI_COMM_WORLD, &Request[7]);

		MPI_Isend(SendSourcepp, ssxd, MPI_DOUBLE, ppPE, 4, MPI_COMM_WORLD, &Request[8]);
		MPI_Isend(SendSourcepm, ssxd, MPI_DOUBLE, pmPE, 5, MPI_COMM_WORLD, &Request[9]);
 		MPI_Isend(SendSourcemp, ssyd, MPI_DOUBLE, mpPE, 6, MPI_COMM_WORLD, &Request[10]);
		MPI_Isend(SendSourcemm, ssyd, MPI_DOUBLE, mmPE, 7, MPI_COMM_WORLD, &Request[11]);

		MPI_Isend(SendSourceXp, ssx,  MPI_DOUBLE, XpPE, 0, MPI_COMM_WORLD, &Request[12]);
		MPI_Isend(SendSourceXm, ssx,  MPI_DOUBLE, XmPE, 1, MPI_COMM_WORLD, &Request[13]);
 		MPI_Isend(SendSourceYp, ssy,  MPI_DOUBLE, YpPE, 2, MPI_COMM_WORLD, &Request[14]);
		MPI_Isend(SendSourceYm, ssy,  MPI_DOUBLE, YmPE, 3, MPI_COMM_WORLD, &Request[15]);

 		int ierr;
		ierr = MPI_Waitall(16, Request,Status);

		MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_RETURN);
		if (ierr!=0) 
		{
   			char errtxt[200];
			for (int i=0; i<16; i++) 
			{
				int err = Status[i].MPI_ERROR; 
				int len=200;
				MPI_Error_string(err,errtxt,&len);
				printf("%s; \n",errtxt);
			}
		MPI_Abort(MPI_COMM_WORLD,0);
		}

	//===============================================================
	//=============== For Conduction Boundary Condition =============
	//=============== This Part May Change for Periodic BC ==========
	//===============================================================
	switch(what)
	{
		case COMMU_F:

		int i = 0;
		if (RankIdx_X == 1) 
		{ 
			for (i = 0; i < ssx; i++)  {ReceSourceXm[i] = 0.0;}
			for (i = 0; i < ssxd; i++) {ReceSourcemm[i] = 0.0;ReceSourcemp[i] = 0.0;}
		};

		if (RankIdx_X == Xpa) 
		{ 
			for (i = 0; i < ssx; i++)  {ReceSourceXp[i] = 0.0;}
			for (i = 0; i < ssxd; i++) {ReceSourcepm[i] = 0.0;ReceSourcepm[i] = 0.0;}
		};

		if (RankIdx_Y == 1) 
		{ 
			for (i = 0; i < ssy; i++)  {ReceSourceYm[i] = 0.0;}
			for (i = 0; i < ssyd; i++) {ReceSourcemm[i] = 0.0;ReceSourcepm[i] = 0.0;}
		};

		if (RankIdx_Y == Ypa) 
		{
			for (i = 0; i < ssy; i++)  {ReceSourceYp[i] = 0.0;}
			for (i = 0; i < ssyd; i++) {ReceSourcemp[i] = 0.0;ReceSourcepp[i] = 0.0;}
		};
		break;

	}

	//===============================================================
	//=====================           UnPack       ==================
	//===============================================================
	//unpack the fields or source.
	UnPack(what, k);



	return;

}


void Commute::DoCommuteT(int what, std::vector<int> SendN)
{

	//===============================================================
	//=====================   ISend and IRecev     ==================
	//===============================================================
	int SendDIM;

	std::vector<int> ReceN(8,0);

    MPI_Request Request[16];
    MPI_Status 	 Status[16];

 	MPI_Irecv(&(*(ReceN.begin()+0)), 1, MPI_INT, mmPE, 0, MPI_COMM_WORLD, &Request[0]);
	MPI_Irecv(&(*(ReceN.begin()+1)), 1, MPI_INT, mpPE, 1, MPI_COMM_WORLD, &Request[1]);
	MPI_Irecv(&(*(ReceN.begin()+2)), 1, MPI_INT, pmPE, 2, MPI_COMM_WORLD, &Request[2]);
	MPI_Irecv(&(*(ReceN.begin()+3)), 1, MPI_INT, ppPE, 3, MPI_COMM_WORLD, &Request[3]);

	MPI_Irecv(&(*(ReceN.begin()+4)), 1, MPI_INT, XmPE, 4, MPI_COMM_WORLD, &Request[4]);
	MPI_Irecv(&(*(ReceN.begin()+5)), 1, MPI_INT, XpPE, 5, MPI_COMM_WORLD, &Request[5]);
	MPI_Irecv(&(*(ReceN.begin()+6)), 1, MPI_INT, YmPE, 6, MPI_COMM_WORLD, &Request[6]);
	MPI_Irecv(&(*(ReceN.begin()+7)), 1, MPI_INT, YpPE, 7, MPI_COMM_WORLD, &Request[7]);
	//-------
	MPI_Isend(&(*(SendN.begin()+3)), 1, MPI_INT, ppPE, 0, MPI_COMM_WORLD, &Request[8]);
	MPI_Isend(&(*(SendN.begin()+2)), 1, MPI_INT, pmPE, 1, MPI_COMM_WORLD, &Request[9]);
 	MPI_Isend(&(*(SendN.begin()+1)), 1, MPI_INT, mpPE, 2, MPI_COMM_WORLD, &Request[10]);
	MPI_Isend(&(*(SendN.begin()+0)), 1, MPI_INT, mmPE, 3, MPI_COMM_WORLD, &Request[11]);

	MPI_Isend(&(*(SendN.begin()+5)), 1, MPI_INT, XpPE, 4, MPI_COMM_WORLD, &Request[12]);
	MPI_Isend(&(*(SendN.begin()+4)), 1, MPI_INT, XmPE, 5, MPI_COMM_WORLD, &Request[13]);
 	MPI_Isend(&(*(SendN.begin()+7)), 1, MPI_INT, YpPE, 6, MPI_COMM_WORLD, &Request[14]);
	MPI_Isend(&(*(SendN.begin()+6)), 1, MPI_INT, YmPE, 7, MPI_COMM_WORLD, &Request[15]);


 	int ierr;
	ierr = MPI_Waitall(16, Request,Status);

	//==========test=========================
	//if(Rank == 9) printf("%d=%d\n",Rank,Sendym);
	//if(Rank == 5) printf("%d=%d\n",Rank,Receyp);
	//==========test=========================

    MPI_Request Request2[16];
    MPI_Status 	 Status2[16];

    switch(what)
    {
    	case COMMU_T:
    	SendDIM = SDT_DIM;
    	break;
    	case COMMU_P:
    	SendDIM = SDP_DIM;
    	break;
    }

    MPI_Irecv(ReceSourcemm, ReceN[0]*SendDIM, MPI_DOUBLE, mmPE, 0, MPI_COMM_WORLD, &Request2[0]);
	MPI_Irecv(ReceSourcemp, ReceN[1]*SendDIM, MPI_DOUBLE, mpPE, 1, MPI_COMM_WORLD, &Request2[1]);
	MPI_Irecv(ReceSourcepm, ReceN[2]*SendDIM, MPI_DOUBLE, pmPE, 2, MPI_COMM_WORLD, &Request2[2]);
	MPI_Irecv(ReceSourcepp, ReceN[3]*SendDIM, MPI_DOUBLE, ppPE, 3, MPI_COMM_WORLD, &Request2[3]);
 	
 	MPI_Irecv(ReceSourceXm, ReceN[4]*SendDIM, MPI_DOUBLE, XmPE, 4, MPI_COMM_WORLD, &Request2[4]);
	MPI_Irecv(ReceSourceXp, ReceN[5]*SendDIM, MPI_DOUBLE, XpPE, 5, MPI_COMM_WORLD, &Request2[5]);
	MPI_Irecv(ReceSourceYm, ReceN[6]*SendDIM, MPI_DOUBLE, YmPE, 6, MPI_COMM_WORLD, &Request2[6]);
	MPI_Irecv(ReceSourceYp, ReceN[7]*SendDIM, MPI_DOUBLE, YpPE, 7, MPI_COMM_WORLD, &Request2[7]);

	MPI_Isend(SendSourcepp, SendN[3]*SendDIM, MPI_DOUBLE, ppPE, 0, MPI_COMM_WORLD, &Request2[8]);
	MPI_Isend(SendSourcepm, SendN[2]*SendDIM, MPI_DOUBLE, pmPE, 1, MPI_COMM_WORLD, &Request2[9]);
 	MPI_Isend(SendSourcemp, SendN[1]*SendDIM, MPI_DOUBLE, mpPE, 2, MPI_COMM_WORLD, &Request2[10]);
	MPI_Isend(SendSourcemm, SendN[0]*SendDIM, MPI_DOUBLE, mmPE, 3, MPI_COMM_WORLD, &Request2[11]);
	
	MPI_Isend(SendSourceXp, SendN[5]*SendDIM, MPI_DOUBLE, XpPE, 4, MPI_COMM_WORLD, &Request2[12]);
	MPI_Isend(SendSourceXm, SendN[4]*SendDIM, MPI_DOUBLE, XmPE, 5, MPI_COMM_WORLD, &Request2[13]);
 	MPI_Isend(SendSourceYp, SendN[7]*SendDIM, MPI_DOUBLE, YpPE, 6, MPI_COMM_WORLD, &Request2[14]);
	MPI_Isend(SendSourceYm, SendN[6]*SendDIM, MPI_DOUBLE, YmPE, 7, MPI_COMM_WORLD, &Request2[15]);

	ierr = MPI_Waitall(16, Request2,Status2);

	if (ierr!=0) 
	{
		MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_RETURN);
   		char errtxt[200];
		for (int i=0; i<16; i++) 
		{
			int err = Status2[i].MPI_ERROR; 
			int len=200;
			MPI_Error_string(err,errtxt,&len);
			printf("%s; \n",errtxt);
		}
		MPI_Abort(MPI_COMM_WORLD,0);
	}


	UnPackT(what, ReceN);
	return;
}

		//      ___________
		//     |mp | yp| pp|
		//     |___|___|___|
		//     |xm |   | xp|
		//     |___|___|___|
		//     |mm | ym| pm|
		//     |___|___|___|
		//

void Commute::DoPack(int what, int k)

{

	Mesh *p_Meshs  = p_domain()->p_Mesh();
	int i,j,n,m;
	int LayerGridX;
	int LayerGridY;
	MultiGrid *p_Multi = NULL;

	switch (what)
	{

		//===============================================================
		//===================== Pack Plasma Source ======================
		//===============================================================
		case COMMU_S:
		case COMMU_SO:

		int NSource=SOU_DIM;
		if(p_domain()->Get_Nbeam()>0&&kold!=k) NSource=SOU_DIM+BEA_DIM;
		
		// Put the Sources at the Overlapping Cells into Send Array;
		// Send Direction: Y: up and down
		for (m = 0 ; m<=1; m++)
		{
			Cell &mp = p_Meshs->GetCell(m,GridY+1-m,k); 
			Cell &pp = p_Meshs->GetCell(GridX+1-m,GridY+1-m,k);

			for (i=0; i < GridX; i++)
			{

				Cell &cm = p_Meshs->GetCell(i+1,m,k);
				Cell &cp = p_Meshs->GetCell(i+1,GridY+1-m,k);

				for (n = 0; n < NSource; n++)
				{
					SendSourceYm[ (GridX*m+i)*NSource + n ] = cm.W_Source[n];
					SendSourceYp[ (GridX*m+i)*NSource + n ] = cp.W_Source[n];

					//diagonal
					if(i==0)
					{
						SendSourcemp[ m*NSource + n ] = mp.W_Source[n];
						SendSourcepp[ m*NSource + n ] = pp.W_Source[n];
					}
				}

			}
		}
		// Put the Sources at the Overlapping Cells into Send Array;
		// Send Direction: X: left and right
		for (m = 0 ; m<=1; m++)
		{

			Cell &mm = p_Meshs->GetCell(m,m,k);
			Cell &pm = p_Meshs->GetCell(GridX+1-m,m,k);

			for (j=0; j < GridY; j++)
			{

				Cell &cm = p_Meshs->GetCell(m,j+1,k);
				Cell &cp = p_Meshs->GetCell(GridX+1-m,j+1,k);
				for (n = 0; n < NSource; n++)
				{
					SendSourceXm[ (GridY*m+j)*NSource + n] = cm.W_Source[n];
					SendSourceXp[ (GridY*m+j)*NSource + n] = cp.W_Source[n];

					//diagonal
					if(j==0)
					{
						SendSourcemm[ m*NSource + n ] = mm.W_Source[n];
						SendSourcepm[ m*NSource + n ] = pm.W_Source[n];
					}
				}

			}
		}

		break;



		//===============================================================
		//===================== Pack Wakefields    ======================
		//===============================================================
		case COMMU_F:

			Cell &mp = p_Meshs->GetCell(1,GridY,k);
			Cell &pp = p_Meshs->GetCell(GridX,GridY,k);

			for (i=1; i <= GridX; i++)
			{
				Cell &cm = p_Meshs->GetCell(i,	  1,k);
				Cell &cp = p_Meshs->GetCell(i,GridY,k);

				for (n = 0; n < WAK_DIM2; n++)
				{
					SendSourceYm[ (i-1)*WAK_DIM2 + n ] = cm.W_Fields[n+5];
					SendSourceYp[ (i-1)*WAK_DIM2 + n ] = cp.W_Fields[n+5];

					if(i==1)
					{
						SendSourcemp[n] = mp.W_Source[n+5];
						SendSourcepp[n] = pp.W_Source[n+5];
					}

				}

			}

			Cell &mm = p_Meshs->GetCell(1,1,k);
			Cell &pm = p_Meshs->GetCell(GridX,1,k);

			for (j=1; j <= GridY; j++)
			{
				Cell &cm = p_Meshs->GetCell(1,	  j,k);
				Cell &cp = p_Meshs->GetCell(GridX,j,k);
				for (n = 0; n < WAK_DIM2; n++)
				{
					SendSourceXm[ (j-1)*WAK_DIM2 + n ] = cm.W_Fields[n+5];
					SendSourceXp[ (j-1)*WAK_DIM2 + n ] = cp.W_Fields[n+5];

					if(j==1)
					{
						SendSourcemm[n] = mm.W_Source[n+5];
						SendSourcepm[n] = pm.W_Source[n+5];
					}
				}

			}

		break;


		//===============================================================
		//===================== Pack Vector Potential  ==================
		//===============================================================




		//===============================================================
		//===================== Pack Multigrid Field  ===================
		//===============================================================
		//Exchange multigrid potential
		case COMMU_MG_P:

		p_Multi = p_domain()->p_MG();

		LayerGridX = p_Multi->GetLayerGridX(k);
		LayerGridY = p_Multi->GetLayerGridY(k);

		for (j = 1; j<= LayerGridY; j++)
		{
			MG_Cell &cxm = p_Multi->GetMGCell(1,		  j, k);
			MG_Cell &cxp = p_Multi->GetMGCell(LayerGridX, j, k);
			SendSourceXm[j-1] = cxm.M_value[0];
			SendSourceXp[j-1] = cxp.M_value[0];
		}

		for (i = 1; i<= LayerGridX; i++)
		{
			MG_Cell &cym = p_Multi->GetMGCell(i,		  1, k);
			MG_Cell &cyp = p_Multi->GetMGCell(i, LayerGridY, k);


			SendSourceYm[i-1] = cym.M_value[0];
			SendSourceYp[i-1] = cyp.M_value[0];
		}
		break;


		case COMMU_MG_P_C:

		p_Multi = p_domain()->p_MG();

		LayerGridX = p_Multi->GetLayerGridX(k);
		LayerGridY = p_Multi->GetLayerGridY(k);

		for (j = 1; j<= LayerGridY; j++)
		{
			MG_Cell &cxm = p_Multi->GetMGCell(1,		  j, k);
			MG_Cell &cxp = p_Multi->GetMGCell(LayerGridX, j, k);

			SendSourceXm[(j-1)*2+0] = (cxm.C_value[0]).real();
			SendSourceXm[(j-1)*2+1] = (cxm.C_value[0]).imag();

			SendSourceXp[(j-1)*2+0] = (cxp.C_value[0]).real();
			SendSourceXp[(j-1)*2+1] = (cxp.C_value[0]).imag();
		}

		for (i = 1; i<= LayerGridX; i++)
		{
			MG_Cell &cym = p_Multi->GetMGCell(i,		  1, k);
			MG_Cell &cyp = p_Multi->GetMGCell(i, LayerGridY, k);

			SendSourceYm[(i-1)*2+0] = (cym.C_value[0]).real();
			SendSourceYm[(i-1)*2+1] = (cym.C_value[0]).imag();
			SendSourceYp[(i-1)*2+0] = (cyp.C_value[0]).real();
			SendSourceYp[(i-1)*2+1] = (cyp.C_value[0]).imag();
		}
		break;


		//===============================================================
		//===================== Pack Multigrid Source  ==================
		//===============================================================
		//Exchange multigrid residual
		case COMMU_MG_R:

		p_Multi = p_domain()->p_MG();

		LayerGridX = p_Multi->GetLayerGridX(k);
		LayerGridY = p_Multi->GetLayerGridY(k);

		for (j = 1; j<= LayerGridY; j++)
		{
			MG_Cell &cxm = p_Multi->GetMGCell(1,		  j, k);
			MG_Cell &cxp = p_Multi->GetMGCell(LayerGridX, j, k);
			SendSourceXm[j-1] = cxm.M_value[2];
			SendSourceXp[j-1] = cxp.M_value[2];
		}

		for (i = 1; i<= LayerGridX; i++)
		{
			MG_Cell &cym = p_Multi->GetMGCell(i,		  1, k);
			MG_Cell &cyp = p_Multi->GetMGCell(i, LayerGridY, k);


			SendSourceYm[i-1] = cym.M_value[2];
			SendSourceYp[i-1] = cyp.M_value[2];
		}
		break;



		case COMMU_MG_R_C:

		p_Multi = p_domain()->p_MG();

		LayerGridX = p_Multi->GetLayerGridX(k);
		LayerGridY = p_Multi->GetLayerGridY(k);

		for (j = 1; j<= LayerGridY; j++)
		{
			MG_Cell &cxm = p_Multi->GetMGCell(1,		  j, k);
			MG_Cell &cxp = p_Multi->GetMGCell(LayerGridX, j, k);

			SendSourceXm[(j-1)*2+0] = (cxm.C_value[2]).real();
			SendSourceXm[(j-1)*2+1] = (cxm.C_value[2]).imag();
			SendSourceXp[(j-1)*2+0] = (cxp.C_value[2]).real();
			SendSourceXp[(j-1)*2+1] = (cxp.C_value[2]).imag();
		}

		for (i = 1; i<= LayerGridX; i++)
		{
			MG_Cell &cym = p_Multi->GetMGCell(i,		  1, k);
			MG_Cell &cyp = p_Multi->GetMGCell(i, LayerGridY, k);

			SendSourceYm[(i-1)*2+0] = (cym.C_value[2]).real();
			SendSourceYm[(i-1)*2+1] = (cym.C_value[2]).imag();
			SendSourceYp[(i-1)*2+0] = (cyp.C_value[2]).real();
			SendSourceYp[(i-1)*2+1] = (cyp.C_value[2]).imag();
		}
		break;


	}



	return;
}


void Commute::UnPack(int what, int k)
{
	Mesh *p_Meshs  = p_domain()->p_Mesh();

	int i,j,n,m;
	MultiGrid *p_Multi = NULL;
	int LayerGridX;
	int LayerGridY;

	switch (what)
	{
		//===============================================================
		//===================== Unpack Plasma Source  ===================
		//===============================================================
		case COMMU_S:
		case COMMU_SO:
		int NSource=SOU_DIM;
		if(p_domain()->Get_Nbeam()>0&&kold!=k) 
		{
			NSource=SOU_DIM+BEA_DIM;
			kold=k;
		}
		// Pull sources from the Rece Array, and add on the edging cells ;
		// Receive Direction: Y
		for (m = 0; m<=1; m++)
		{

			Cell &mm = p_Meshs->GetCell(1-m,1-m,k);
			Cell &pm = p_Meshs->GetCell(GridX+m,1-m,k);

			for (i = 0; i < GridX; i++)
			{

				Cell &cm = p_Meshs->GetCell(i+1,1-m,k);
				Cell &cp = p_Meshs->GetCell(i+1,GridY+m,k);
				for (n = 0; n < NSource; n++)
				{
					cm.W_Source[n] += ReceSourceYm[ (GridX*m+i)*NSource+ n ];
					cp.W_Source[n] += ReceSourceYp[ (GridX*m+i)*NSource+ n ];
				
					if(i==0)
					{
						mm.W_Source[n] += ReceSourcemm[ m*NSource+ n ];
						pm.W_Source[n] += ReceSourcepm[ m*NSource+ n ];
					}
				}

			}
		}

		for (m = 0; m<=1; m++)
		{

			Cell &mp = p_Meshs->GetCell(1-m,GridY+m,k);
			Cell &pp = p_Meshs->GetCell(GridX+m,GridY+m,k);


			for (j = 0; j < GridY; j++)
			{

				Cell &cm = p_Meshs->GetCell(1-m,j+1,k);
				Cell &cp = p_Meshs->GetCell(GridX+m,j+1,k);
				for (n = 0; n < NSource; n++)
				{
					cm.W_Source[n] += ReceSourceXm[ (GridY*m+j)*NSource + n];
					cp.W_Source[n] += ReceSourceXp[ (GridY*m+j)*NSource + n];

					if(j==0)
					{
						mp.W_Source[n] += ReceSourcemp[ m*NSource+ n ];
						pp.W_Source[n] += ReceSourcepp[ m*NSource+ n ];
					}
				}



			}
		}

		//======== if Conduction Boundary Condition ======
		if(p_domain()->Get_BC()==1)
		{	
			if (RankIdx_X == 1) 
			{	
				for (i=0; i<=GridY+1; i++)
				{	Cell &c = p_Meshs->GetCell(0, i,k); c.W_Denn  = c.W_Deni;
					for (n = 1; n < NSource; n++) {c.W_Source[n]==0;};

				}
			}

			if (RankIdx_X == Xpa) 
			{	
				for (i=0; i<=GridY+1; i++)
				{	Cell &c = p_Meshs->GetCell(GridX+1,i,k); c.W_Denn  = c.W_Deni;
					for (n = 1; n < NSource; n++) {c.W_Source[n]==0;};
				}
			
			}

			if (RankIdx_Y == 1) 
			{	
				for (i=0; i<=GridX+1; i++)
				{	Cell &c = p_Meshs->GetCell(i,0,k); c.W_Denn  = c.W_Deni;
					for (n = 1; n < NSource; n++) {c.W_Source[n]==0;};
				}
			}

			if (RankIdx_Y == Ypa) 
			{	
				for (i=0; i<=GridX+1; i++)
				{	Cell &c = p_Meshs->GetCell(i,GridY+1,k); c.W_Denn  = c.W_Deni;
					for (n = 1; n < NSource; n++) {c.W_Source[n]==0;};
				}
			}

		}

		break;


		//===============================================================
		//===================== Unpack Wakefields     ===================
		//===============================================================
		case COMMU_F:

			Cell &mm = p_Meshs->GetCell(0,0,k);
			Cell &pm = p_Meshs->GetCell(GridX+1,0,k);

			for (i=1; i <= GridX; i++)
			{
				Cell &cm = p_Meshs->GetCell(i,	    0,k);
				Cell &cp = p_Meshs->GetCell(i,GridY+1,k);
				for (n = 0; n < WAK_DIM2; n++)
				{
					cm.W_Fields[n+5] = ReceSourceYm[ (i-1)*WAK_DIM2 + n ];
					cp.W_Fields[n+5] = ReceSourceYp[ (i-1)*WAK_DIM2 + n ];

					if(i==1)
					{
						mm.W_Fields[n+5] = ReceSourcemm[ n ];
						pm.W_Fields[n+5] = ReceSourcepm[ n ];
					}

				}

			}

			Cell &mp = p_Meshs->GetCell(0,GridY+1,k);
			Cell &pp = p_Meshs->GetCell(GridX+1,GridY+1,k);

			for (j=1; j <= GridY; j++)
			{
				Cell &cm = p_Meshs->GetCell(0,	    j,k);
				Cell &cp = p_Meshs->GetCell(GridX+1,j,k);
				for (n = 0; n < WAK_DIM2; n++)
				{
					cm.W_Fields[n+5] = ReceSourceXm[ (j-1)*WAK_DIM2 + n ];
					cp.W_Fields[n+5] = ReceSourceXp[ (j-1)*WAK_DIM2 + n ];

					if(j==1)
					{
						mp.W_Fields[n+5] = ReceSourcemp[ n ];
						pp.W_Fields[n+5] = ReceSourcepp[ n ];
					}
				}

			}
		break;

		//===============================================================
		//===================== Unpack Vector Potential   ===============
		//===============================================================



		//===============================================================
		//===================== Unpack Multigird Field ==================
		//===============================================================
		case COMMU_MG_P:

		p_Multi = p_domain()->p_MG();

		LayerGridX = p_Multi->GetLayerGridX(k);
		LayerGridY = p_Multi->GetLayerGridY(k);

		for (j = 1; j<= LayerGridY; j++)
		{
			MG_Cell &cxm = p_Multi->GetMGCell(0,			 j, k);
			MG_Cell &cxp = p_Multi->GetMGCell(LayerGridX+1,  j, k);
			cxm.M_value[0] = ReceSourceXm[j-1];
			cxp.M_value[0] = ReceSourceXp[j-1];
		}

		for (i = 1; i<= LayerGridX; i++)
		{
			MG_Cell &cym = p_Multi->GetMGCell(i,			 0, k);
			MG_Cell &cyp = p_Multi->GetMGCell(i,  LayerGridY+1, k);
			cym.M_value[0] = ReceSourceYm[i-1];
			cyp.M_value[0] = ReceSourceYp[i-1];
		}


		if (RankIdx_X == 1) 
		{ for (i = 0; i <= LayerGridY+1; i++) { (p_Multi->GetMGCell(0, 			  i, k)).M_value[0] = 0.0; } };

		if (RankIdx_X == Xpa) 
		{ for (i = 0; i <= LayerGridY+1; i++) { (p_Multi->GetMGCell(LayerGridX+1, i, k)).M_value[0] = 0.0; } };

		if (RankIdx_Y == 1) 
		{ for (i = 0; i <= LayerGridX+1; i++) { (p_Multi->GetMGCell(i, 			  0, k)).M_value[0] = 0.0; } };

		if (RankIdx_Y == Ypa) 
		{ for (i = 0; i <= LayerGridX+1; i++) { (p_Multi->GetMGCell(i, LayerGridY+1, k)).M_value[0] = 0.0; } };

		break;


		case COMMU_MG_P_C:

		p_Multi = p_domain()->p_MG();

		LayerGridX = p_Multi->GetLayerGridX(k);
		LayerGridY = p_Multi->GetLayerGridY(k);

		for (j = 1; j<= LayerGridY; j++)
		{
			MG_Cell &cxm = p_Multi->GetMGCell(0,			 j, k);
			MG_Cell &cxp = p_Multi->GetMGCell(LayerGridX+1,  j, k);
			cxm.C_value[0] = ReceSourceXm[(j-1)*2+0]+ci*ReceSourceXm[(j-1)*2+1];
			cxp.C_value[0] = ReceSourceXp[(j-1)*2+0]+ci*ReceSourceXp[(j-1)*2+1];
		}

		for (i = 1; i<= LayerGridX; i++)
		{
			MG_Cell &cym = p_Multi->GetMGCell(i,			 0, k);
			MG_Cell &cyp = p_Multi->GetMGCell(i,  LayerGridY+1, k);
			cym.C_value[0] = ReceSourceYm[(i-1)*2+0]+ci*ReceSourceYm[(i-1)*2+1];
			cyp.C_value[0] = ReceSourceYp[(i-1)*2+0]+ci*ReceSourceYp[(i-1)*2+1];
		}

		if (RankIdx_X == 1) 
		{ for (i = 0; i <= LayerGridY+1; i++) { (p_Multi->GetMGCell(0, 			  i, k)).C_value[0] = 0.0; } };

		if (RankIdx_X == Xpa) 
		{ for (i = 0; i <= LayerGridY+1; i++) { (p_Multi->GetMGCell(LayerGridX+1, i, k)).C_value[0] = 0.0; } };

		if (RankIdx_Y == 1) 
		{ for (i = 0; i <= LayerGridX+1; i++) { (p_Multi->GetMGCell(i, 			  0, k)).C_value[0] = 0.0; } };

		if (RankIdx_Y == Ypa) 
		{ for (i = 0; i <= LayerGridX+1; i++) { (p_Multi->GetMGCell(i, LayerGridY+1, k)).C_value[0] = 0.0; } };

		break;

		//===============================================================
		//===================== Unpack Multigird Source =================
		//===============================================================
		case COMMU_MG_R:

		p_Multi = p_domain()->p_MG();

		LayerGridX = p_Multi->GetLayerGridX(k);
		LayerGridY = p_Multi->GetLayerGridY(k);

		for (j = 1; j<= LayerGridY; j++)
		{
			MG_Cell &cxm = p_Multi->GetMGCell(0,			 j, k);
			MG_Cell &cxp = p_Multi->GetMGCell(LayerGridX+1,  j, k);
			cxm.M_value[2] = ReceSourceXm[j-1];
			cxp.M_value[2] = ReceSourceXp[j-1];
		}

		for (i = 1; i<= LayerGridX; i++)
		{
			MG_Cell &cym = p_Multi->GetMGCell(i,			 0, k);
			MG_Cell &cyp = p_Multi->GetMGCell(i,  LayerGridY+1, k);
			cym.M_value[2] = ReceSourceYm[i-1];
			cyp.M_value[2] = ReceSourceYp[i-1];
		}



		if (RankIdx_X == 1) 
		{ for (i = 0; i <= LayerGridY+1; i++) 
		{ (p_Multi->GetMGCell(0, i, k)).M_value[2]			  = (p_Multi->GetMGCell(1, i, k)).M_value[2]; } };

		if (RankIdx_X == Xpa) 
		{ for (i = 0; i <= LayerGridY+1; i++) 
		{ (p_Multi->GetMGCell(LayerGridX+1, i, k)).M_value[2] = (p_Multi->GetMGCell(LayerGridX, i, k)).M_value[2]; } };

		if (RankIdx_Y == 1) 
		{ for (i = 0; i <= LayerGridX+1; i++) 
		{ (p_Multi->GetMGCell(i, 0, k)).M_value[2] 			  = (p_Multi->GetMGCell(i, 1, k)).M_value[2]; } };

		if (RankIdx_Y == Ypa) 
		{ for (i = 0; i <= LayerGridX+1; i++)
		{ (p_Multi->GetMGCell(i, LayerGridY+1, k)).M_value[2] = (p_Multi->GetMGCell(i, LayerGridY, k)).M_value[2]; } };

		break;


		case COMMU_MG_R_C:

		p_Multi = p_domain()->p_MG();

		LayerGridX = p_Multi->GetLayerGridX(k);
		LayerGridY = p_Multi->GetLayerGridY(k);

		for (j = 1; j<= LayerGridY; j++)
		{
			MG_Cell &cxm = p_Multi->GetMGCell(0,			 j, k);
			MG_Cell &cxp = p_Multi->GetMGCell(LayerGridX+1,  j, k);
			cxm.C_value[2] = ReceSourceXm[(j-1)*2+0]+ci*ReceSourceXm[(j-1)*2+1];
			cxp.C_value[2] = ReceSourceXp[(j-1)*2+0]+ci*ReceSourceXp[(j-1)*2+1];
		}

		for (i = 1; i<= LayerGridX; i++)
		{
			MG_Cell &cym = p_Multi->GetMGCell(i,			 0, k);
			MG_Cell &cyp = p_Multi->GetMGCell(i,  LayerGridY+1, k);
			cym.C_value[2] = ReceSourceYm[(i-1)*2+0]+ci*ReceSourceYm[(i-1)*2+1];
			cyp.C_value[2] = ReceSourceYp[(i-1)*2+0]+ci*ReceSourceYp[(i-1)*2+1];
		}



		if (RankIdx_X == 1) 
		{ for (i = 0; i <= LayerGridY+1; i++) 
		{ (p_Multi->GetMGCell(0, i, k)).C_value[2]			  = (p_Multi->GetMGCell(1, i, k)).C_value[2]; } };

		if (RankIdx_X == Xpa) 
		{ for (i = 0; i <= LayerGridY+1; i++) 
		{ (p_Multi->GetMGCell(LayerGridX+1, i, k)).C_value[2] = (p_Multi->GetMGCell(LayerGridX, i, k)).C_value[2]; } };

		if (RankIdx_Y == 1) 
		{ for (i = 0; i <= LayerGridX+1; i++) 
		{ (p_Multi->GetMGCell(i, 0, k)).C_value[2] 			  = (p_Multi->GetMGCell(i, 1, k)).C_value[2]; } };

		if (RankIdx_Y == Ypa) 
		{ for (i = 0; i <= LayerGridX+1; i++)
		{ (p_Multi->GetMGCell(i, LayerGridY+1, k)).C_value[2] = (p_Multi->GetMGCell(i, LayerGridY, k)).C_value[2]; } };

		break;


	}




	return;
}



void Commute::UnPackT(int what, std::vector<int> ReceN)
{

	int n;
	double x0,y0,z0,px,py,pz;
	double xt,yt, sx, sy;
	double vx,vy,old_x,old_y,old_vx,old_vy; 
	double q2m, weight;
	double Ex0,Ey0,Ez0;
	int type;

	Trajectory *p =NULL;
	Particle  *pp =NULL;


	int TpCellx =  p_domain()->p_Mesh()->Get_TpCellx();
	int TpCelly =  p_domain()->p_Mesh()->Get_TpCelly();


	switch(what)
	{
		
		case COMMU_T: // 
		

		for(int dir=0; dir<8; dir++)
		{
			double* Re=NULL;

			if(dir==0) Re=ReceSourcemm;
			if(dir==1) Re=ReceSourcemp;
			if(dir==2) Re=ReceSourcepm;
			if(dir==3) Re=ReceSourcepp;
			if(dir==4) Re=ReceSourceXm;
			if(dir==5) Re=ReceSourceXp;
			if(dir==6) Re=ReceSourceYm;
			if(dir==7) Re=ReceSourceYp;

			for(n=0; n<ReceN[dir]; n++)
			{
				xt = *Re; Re++; 
				yt = *Re; Re++;

				x0 = *Re; Re++; 
				y0 = *Re; Re++; 
				z0 = *Re; Re++; 

				vx = *Re; Re++; 
				vy = *Re; Re++; 

				old_x  = *Re; Re++; 
				old_y  = *Re; Re++; 
				old_vx = *Re; Re++; 
				old_vy = *Re; Re++; 

				sx = *Re; Re++; 
				sy = *Re; Re++; 

				p = new Trajectory(x0, y0, z0, TpCellx, TpCelly,sx,sy);

				p->x = xt;
				p->y = yt;
				p->Vx= vx;
				p->Vy= vy;

				p->old_x  = old_x;
				p->old_y  = old_y;
				p->old_vx = old_vx;
				p->old_vy = old_vy;

				p->Vxx = (p->Vx)*(p->Vx);
				p->Vyy = (p->Vy)*(p->Vy);
				p->Vxy = (p->Vx)*(p->Vy);

				auto upper=std::upper_bound(CellAccX.begin(),CellAccX.end(),xt);
				p->idx_i= (upper-CellAccX.begin()-1);
				upper=std::upper_bound(CellAccY.begin(),CellAccY.end(),yt);
				p->idx_j=(upper-CellAccY.begin()-1);
			}

		}
		break; //COMMU_T: // 



		case COMMU_P:
		
		// for(n=0; n<Recexm; n++)
		// {
		// 	x0 = ReceSourceXm[n*SDP_DIM + 3];
		// 	y0 = ReceSourceXm[n*SDP_DIM + 4];
		// 	z0 = ReceSourceXm[n*SDP_DIM + 5];
		// 	px = ReceSourceXm[n*SDP_DIM + 6];
		// 	py = ReceSourceXm[n*SDP_DIM + 7];
		// 	pz = ReceSourceXm[n*SDP_DIM + 8];
		// 	Ex0= ReceSourceXm[n*SDP_DIM + 9];
		// 	Ey0= ReceSourceXm[n*SDP_DIM +10];
		// 	Ez0= ReceSourceXm[n*SDP_DIM +11];

		// 	type=(int)ReceSourceXm[n*SDP_DIM +12];
		// 	q2m    =  ReceSourceXm[n*SDP_DIM +13];
		// 	weight =  ReceSourceXm[n*SDP_DIM +14];

		// 	switch(type)
		// 	{
		// 	case ELECTRON:
		// 	pp = new Electron(x0, y0, z0, px, py, pz,
		// 					Ex0, Ey0, Ez0, q2m, weight);
		// 	break;

		// 	case ION:
		// 	pp = new Ion(x0, y0, z0, px, py, pz,
		// 					Ex0, Ey0, Ez0, q2m, weight);
		// 	break;
		// 	}

		// 	pp->x = ReceSourceXm[n*SDP_DIM + 0];
		// 	pp->y = ReceSourceXm[n*SDP_DIM + 1];
		// 	pp->z = ReceSourceXm[n*SDP_DIM + 2];

		// 	pp->Wxw = ReceSourceXm[n*SDP_DIM + 15];
		// 	pp->Wyw = ReceSourceXm[n*SDP_DIM + 16];
		// 	pp->Wzw = ReceSourceXm[n*SDP_DIM + 17];
		// 	pp->Wxl = ReceSourceXm[n*SDP_DIM + 18];
		// 	pp->Wyl = ReceSourceXm[n*SDP_DIM + 19];
		// 	pp->Wzl = ReceSourceXm[n*SDP_DIM + 20];


		// }


		// for(n=0; n<Recexp; n++)
		// {
		// 	x0 = ReceSourceXp[n*SDP_DIM + 3];
		// 	y0 = ReceSourceXp[n*SDP_DIM + 4];
		// 	z0 = ReceSourceXp[n*SDP_DIM + 5];
		// 	px = ReceSourceXp[n*SDP_DIM + 6];
		// 	py = ReceSourceXp[n*SDP_DIM + 7];
		// 	pz = ReceSourceXp[n*SDP_DIM + 8];
		// 	Ex0= ReceSourceXp[n*SDP_DIM + 9];
		// 	Ey0= ReceSourceXp[n*SDP_DIM +10];
		// 	Ez0= ReceSourceXp[n*SDP_DIM +11];

		// 	type=(int)ReceSourceXp[n*SDP_DIM +12];
		// 	q2m=      ReceSourceXp[n*SDP_DIM +13];
		// 	weight =  ReceSourceXp[n*SDP_DIM +14];
		// 	switch(type)
		// 	{
		// 	case ELECTRON:
		// 	pp = new Electron(x0, y0, z0, px, py, pz,
		// 					Ex0, Ey0, Ez0, q2m, weight);
		// 	break;

		// 	case ION:
		// 	pp = new Ion(x0, y0, z0, px, py, pz,
		// 					Ex0, Ey0, Ez0, q2m, weight);
		// 	break;
		// 	}

		// 	pp->x = ReceSourceXp[n*SDP_DIM + 0];
		// 	pp->y = ReceSourceXp[n*SDP_DIM + 1];
		// 	pp->z = ReceSourceXp[n*SDP_DIM + 2];

		// 	pp->Wxw = ReceSourceXp[n*SDP_DIM + 15];
		// 	pp->Wyw = ReceSourceXp[n*SDP_DIM + 16];
		// 	pp->Wzw = ReceSourceXp[n*SDP_DIM + 17];
		// 	pp->Wxl = ReceSourceXp[n*SDP_DIM + 18];
		// 	pp->Wyl = ReceSourceXp[n*SDP_DIM + 19];
		// 	pp->Wzl = ReceSourceXp[n*SDP_DIM + 20];
		// }

		// for(n=0; n<Receym; n++)
		// {
		// 	x0 = ReceSourceYm[n*SDP_DIM + 3];
		// 	y0 = ReceSourceYm[n*SDP_DIM + 4];
		// 	z0 = ReceSourceYm[n*SDP_DIM + 5];
		// 	px = ReceSourceYm[n*SDP_DIM + 6];
		// 	py = ReceSourceYm[n*SDP_DIM + 7];
		// 	pz = ReceSourceYm[n*SDP_DIM + 8];
		// 	Ex0= ReceSourceYm[n*SDP_DIM + 9];
		// 	Ey0= ReceSourceYm[n*SDP_DIM +10];
		// 	Ez0= ReceSourceYm[n*SDP_DIM +11];

		// 	type=(int)ReceSourceYm[n*SDP_DIM +12];
		// 	q2m =     ReceSourceYm[n*SDP_DIM +13];
		// 	weight =  ReceSourceYm[n*SDP_DIM +14];
		// 	switch(type)
		// 	{
		// 	case ELECTRON:
		// 	pp = new Electron(x0, y0, z0, px, py, pz,
		// 					Ex0, Ey0, Ez0, q2m, weight);
		// 	break;

		// 	case ION:
		// 	pp = new Ion(x0, y0, z0, px, py, pz,
		// 					Ex0, Ey0, Ez0, q2m, weight);
		// 	break;
		// 	}

		// 	pp->x = ReceSourceYm[n*SDP_DIM + 0];
		// 	pp->y = ReceSourceYm[n*SDP_DIM + 1];
		// 	pp->z = ReceSourceYm[n*SDP_DIM + 2];

		// 	pp->Wxw = ReceSourceYm[n*SDP_DIM + 15];
		// 	pp->Wyw = ReceSourceYm[n*SDP_DIM + 16];
		// 	pp->Wzw = ReceSourceYm[n*SDP_DIM + 17];
		// 	pp->Wxl = ReceSourceYm[n*SDP_DIM + 18];
		// 	pp->Wyl = ReceSourceYm[n*SDP_DIM + 19];
		// 	pp->Wzl = ReceSourceYm[n*SDP_DIM + 20];
		// }

		// for(n=0; n<Receyp; n++)
		// {
		// 	x0 = ReceSourceYp[n*SDP_DIM + 3];
		// 	y0 = ReceSourceYp[n*SDP_DIM + 4];
		// 	z0 = ReceSourceYp[n*SDP_DIM + 5];
		// 	px = ReceSourceYp[n*SDP_DIM + 6];
		// 	py = ReceSourceYp[n*SDP_DIM + 7];
		// 	pz = ReceSourceYp[n*SDP_DIM + 8];
		// 	Ex0= ReceSourceYp[n*SDP_DIM + 9];
		// 	Ey0= ReceSourceYp[n*SDP_DIM +10];
		// 	Ez0= ReceSourceYp[n*SDP_DIM +11];

		// 	type=(int)ReceSourceYp[n*SDP_DIM +12];
		// 	q2m=      ReceSourceYp[n*SDP_DIM +13];
		// 	weight =  ReceSourceYp[n*SDP_DIM +14];
			
		// 	switch(type)
		// 	{
		// 	case ELECTRON:
		// 	pp = new Electron(x0, y0, z0, px, py, pz,
		// 					Ex0, Ey0, Ez0, q2m, weight);
		// 	break;

		// 	case ION:
		// 	pp = new Ion(x0, y0, z0, px, py, pz,
		// 					Ex0, Ey0, Ez0, q2m, weight);
		// 	break;
		// 	}

		// 	pp->x = ReceSourceYp[n*SDP_DIM + 0];
		// 	pp->y = ReceSourceYp[n*SDP_DIM + 1];
		// 	pp->z = ReceSourceYp[n*SDP_DIM + 2];

		// 	pp->Wxw = ReceSourceYp[n*SDP_DIM + 15];
		// 	pp->Wyw = ReceSourceYp[n*SDP_DIM + 16];
		// 	pp->Wzw = ReceSourceYp[n*SDP_DIM + 17];
		// 	pp->Wxl = ReceSourceYp[n*SDP_DIM + 18];
		// 	pp->Wyl = ReceSourceYp[n*SDP_DIM + 19];
		// 	pp->Wzl = ReceSourceYp[n*SDP_DIM + 20];
		// }

		break;

	}

	return;
}


Commute::~Commute()
{
	delete[] SendSourceXm;  //Array to send for the sources;
	delete[] SendSourceYm;  //Array to send for the sources;
   	delete[] SendSourceXp;  //Array to send for the sources;
  	delete[] SendSourceYp;  //Array to send for the sources;

   	delete[] ReceSourceXm;  //Array for Rece 
  	delete[] ReceSourceYm;  //Array for Rece 
	delete[] ReceSourceXp;  //Array for Rece 
	delete[] ReceSourceYp;  //Array for Rece 


}





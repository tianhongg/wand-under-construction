//----------------------------------------------------------------------------------||
//-------------------                pushtrajs.cpp               -------------------||
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
//---Starting---------           : Feb-22-2019                  --------------------||
//---Email------------           : tw474@cornell.edu            --------------------||
//---Group------------           : Dr. Gennady Shvets' Group    --------------------||
//---Copyright--------           : (C) 2019 by Tianhong Wang    --------------------||
//----------------------------------------------------------------------------------||
//----------------------------------------------------------------------------------||

#include "wand_PIC.h"


void Mesh::PushTrajectory(double k0, int k, int step)
{

	// temporary remove
	return;

}

void Mesh::PushTrajectory_Half()
{

	Trajectory *p = NULL;

	int Xpa  = p_domain()->p_Partition()->GetXpart();
	int Ypa  = p_domain()->p_Partition()->GetYpart();
	
	double xt, yt, Vx, Vy, xtp, ytp;
	double dztmp;

	p = p_Trajectory;

	int i;
	int j;

	while (p)
	{

		xt = p-> x;
		yt = p-> y;

		Vx = p-> Vx;
		Vy = p-> Vy;

		i=p->idx_i;
		j=p->idx_j;

		//=================================================
		//============Trajectory Outside Boundary =========
		//=================================================
		if(RankIdx_X ==1	&& i==0)
		{ p->Vx = p->Vy = p->Vxx = p->Vyy = p->Vxy=0; p = p->p_PrevTraj; continue;}
		if(RankIdx_X == Xpa && i==GridX)
		{ p->Vx = p->Vy = p->Vxx = p->Vyy = p->Vxy=0; p = p->p_PrevTraj; continue;}
		if(RankIdx_Y == 1	&& j==0)
		{ p->Vx = p->Vy = p->Vxx = p->Vyy = p->Vxy=0; p = p->p_PrevTraj; continue;}
		if(RankIdx_Y == Ypa && j==GridY)
		{ p->Vx = p->Vy = p->Vxx = p->Vyy = p->Vxy=0; p = p->p_PrevTraj; continue;}
		//==================================================

		// Push  position only......
		dztmp=dzz;
		xtp = p-> old_x  + Vx*dztmp;
		ytp = p-> old_y  + Vy*dztmp;

		p-> x = xtp;
		p-> y = ytp;

		p-> old_x = xtp;
		p-> old_y = ytp;

		// update the cell index;
		Cell &c=GetCell(i,j,0);
		while(xtp>c.Xcord+c.dx*0.5&&i<GridX+1)  { p->idx_i++; i++; c=GetCell(i,j,0); }
		while(xtp<c.Xcord-c.dx*0.5&&i>0) 		{ p->idx_i--; i--; c=GetCell(i,j,0); }
		while(ytp>c.Ycord+c.dy*0.5&&j<GridY+1) 	{ p->idx_j++; j++; c=GetCell(i,j,0); }
		while(ytp<c.Ycord-c.dy*0.5&&j>0) 		{ p->idx_j--; j--; c=GetCell(i,j,0); }
		p = p->p_PrevTraj;

	}
	//============================================
	//=========== Exchange Particles =============
	ExchangeT();
	//============================================
	return;

}

void Mesh::PushTrajectory_HalfE(int k) 
{

	Trajectory *p = NULL;
	double ddx;
	double ddy;
	double Ex, Ey, Ez, Psi, Pondx, Pondy, Asq;
	double Fx, Fy, gamma;

	double xt, yt, Vx, Vy, Vxp, Vyp, xtp, ytp;
	double wmm, wmp, wpm, wpp;
	double dxdy = dx*dy;
	double dztmp;
	int i,j, im, jm;

	int Xpa  = p_domain()->p_Partition()->GetXpart();
	int Ypa  = p_domain()->p_Partition()->GetYpart();
	double Xmax = Offset_X+GridX*dx;
	double Ymax = Offset_Y+GridY*dy;

	p = p_Trajectory;
	Vmax = 0.0;
	while (p)
	{

		xt = p-> x;
		yt = p-> y;
		Vx = p-> Vx;
		Vy = p-> Vy;

		i=p->idx_i;
		j=p->idx_j;

		im=i;
		jm=j;
		//=================================================
		//============Trajectory Outside Boundary =========
		//=================================================
		if(RankIdx_X ==1	&& i==0)
		{ p->Vx = p->Vy = p->Vxx = p->Vyy = p->Vxy=0; p = p->p_PrevTraj; continue;}
		if(RankIdx_X == Xpa && i==GridX)
		{ p->Vx = p->Vy = p->Vxx = p->Vyy = p->Vxy=0; p = p->p_PrevTraj; continue;}
		if(RankIdx_Y == 1	&& j==0)
		{ p->Vx = p->Vy = p->Vxx = p->Vyy = p->Vxy=0; p = p->p_PrevTraj; continue;}
		if(RankIdx_Y == Ypa && j==GridY)
		{ p->Vx = p->Vy = p->Vxx = p->Vyy = p->Vxy=0; p = p->p_PrevTraj; continue;}
		//==================================================

		Cell &c=GetCell(i,j,k);
		if(xt<c.Xcord) im--;
		if(yt<c.Ycord) jm--;

		Cell &cmm = GetCell(im,  jm,  k);
		Cell &cmp = GetCell(im,  jm+1,k);
		Cell &cpm = GetCell(im+1,jm,  k);
		Cell &cpp = GetCell(im+1,jm+1,k);

		ddx= (xt-cmm.Xcord)/(cmm.dx*0.5+cpm.dx*0.5);
		ddy= (yt-cmm.Ycord)/(cmm.dy*0.5+cmp.dy*0.5);

		wmm = (1-ddx)*(1-ddy);
		wmp = (1-ddx)*(ddy);
		wpm = (ddx)*(1-ddy);
		wpp = (ddx)*(ddy);

		Ex  = wmm*cmm.W_Ex  + wmp*cmp.W_Ex  + wpm*cpm.W_Ex  + wpp*cpp.W_Ex;
		Ey  = wmm*cmm.W_Ey  + wmp*cmp.W_Ey  + wpm*cpm.W_Ey  + wpp*cpp.W_Ey;
		Ez  = wmm*cmm.W_Ez  + wmp*cmp.W_Ez  + wpm*cpm.W_Ez  + wpp*cpp.W_Ez;

		Psi   = wmm*cmm.W_Psi  + wmp*cmp.W_Psi  + wpm*cpm.W_Psi  + wpp*cpp.W_Psi;
		Pondx = wmm*cmm.W_Ponx + wmp*cmp.W_Ponx + wpm*cpm.W_Ponx + wpp*cpp.W_Ponx;
		Pondy = wmm*cmm.W_Pony + wmp*cmp.W_Pony + wpm*cpm.W_Pony + wpp*cpp.W_Pony;
		
		Asq   = wmm*cmm.W_Asq  + wmp*cmp.W_Asq  + wpm*cpm.W_Asq  + wpp*cpp.W_Asq;

		gamma = 0.5*(1+Psi)*(Vx*Vx+Vy*Vy+1)+0.5*(1.0+0.5*Asq)/(1+Psi);
		Fx = ((gamma*Ex-Pondx*0.25)/(1+Psi) - Vx*(Vx*Ex+Vy*Ey+Ez))/(1+Psi);
		Fy = ((gamma*Ey-Pondy*0.25)/(1+Psi) - Vy*(Vx*Ex+Vy*Ey+Ez))/(1+Psi);

		dztmp=dzz;
		Vxp = p-> old_vx + Fx*dztmp*0.5;
		Vyp = p-> old_vy + Fy*dztmp*0.5;

		//==========================================
		//=========== Adapteive Z Step =============
		double Vr = sqrt(Vxp*Vxp+Vyp*Vyp);
		if(AdaptiveStep>0 & Vr >=Vlim*AdaptiveStep)
		{
			Vxp = Vlim*AdaptiveStep*Vxp/Vr;
			Vyp = Vlim*AdaptiveStep*Vyp/Vr;
		}

		p-> Vx = Vxp;
		p-> Vy = Vyp;

		p-> old_vx = Vxp;
		p-> old_vy = Vyp;
		
		p-> Vxx = Vxp*Vxp;
		p-> Vxy = Vxp*Vyp;
		p-> Vyy = Vyp*Vyp;

		p = p->p_PrevTraj;

		double Vrr=sqrt(Vxp*Vxp+Vyp*Vyp);
		Vmax = std::max(Vmax, Vrr*dz/sqrt(cmm.dx*cmm.dx+cmm.dy*cmm.dy)); // how many grids it can cross

	}
	return;

}


void Mesh::PushTrajectory_HalfB(int k)
{

	Trajectory *p = NULL;
	double ddx;
	double ddy;
	double Psi, Bx, By, Bz;
	double Fx, Fy, gamma;

	double xt, yt, Vx, Vy, Vxp, Vyp, xtp, ytp;
	double wmm, wmp, wpm, wpp;
	double dxdy = dx*dy;
	double dztmp;
	int i,j, im, jm;

	int Xpa  = p_domain()->p_Partition()->GetXpart();
	int Ypa  = p_domain()->p_Partition()->GetYpart();
	double Xmax = Offset_X+GridX*dx;
	double Ymax = Offset_Y+GridY*dy;

	p = p_Trajectory;
	//====== for adpative z step========
	Vmax = 0.0;
	while (p)
	{

		xt = p-> x;
		yt = p-> y;
		Vx = p-> Vx;
		Vy = p-> Vy;

		i=p->idx_i;
		j=p->idx_j;

		im=i;
		jm=j;
		//=================================================
		//============Trajectory Outside Boundary =========
		//=================================================
		if(RankIdx_X ==1	&& i==0)
		{ p->Vx = p->Vy = p->Vxx = p->Vyy = p->Vxy=0; p = p->p_PrevTraj; continue;}
		if(RankIdx_X == Xpa && i==GridX)
		{ p->Vx = p->Vy = p->Vxx = p->Vyy = p->Vxy=0; p = p->p_PrevTraj; continue;}
		if(RankIdx_Y == 1	&& j==0)
		{ p->Vx = p->Vy = p->Vxx = p->Vyy = p->Vxy=0; p = p->p_PrevTraj; continue;}
		if(RankIdx_Y == Ypa && j==GridY)
		{ p->Vx = p->Vy = p->Vxx = p->Vyy = p->Vxy=0; p = p->p_PrevTraj; continue;}
		//==================================================

		Cell &c=GetCell(i,j,k);
		if(xt<c.Xcord) im--;
		if(yt<c.Ycord) jm--;

		Cell &cmm = GetCell(im,  jm,  k);
		Cell &cmp = GetCell(im,  jm+1,k);
		Cell &cpm = GetCell(im+1,jm,  k);
		Cell &cpp = GetCell(im+1,jm+1,k);

		ddx= (xt-cmm.Xcord)/(cmm.dx*0.5+cpm.dx*0.5);
		ddy= (yt-cmm.Ycord)/(cmm.dy*0.5+cmp.dy*0.5);

		wmm = (1-ddx)*(1-ddy);
		wmp = (1-ddx)*(ddy);
		wpm = (ddx)*(1-ddy);
		wpp = (ddx)*(ddy);

		Bx  = wmm*cmm.W_Bx  + wmp*cmp.W_Bx  + wpm*cpm.W_Bx  + wpp*cpp.W_Bx;
		By  = wmm*cmm.W_By  + wmp*cmp.W_By  + wpm*cpm.W_By  + wpp*cpp.W_By;
		Bz  = wmm*cmm.W_Bz  + wmp*cmp.W_Bz  + wpm*cpm.W_Bz  + wpp*cpp.W_Bz;
		Psi = wmm*cmm.W_Psi + wmp*cmp.W_Psi + wpm*cpm.W_Psi + wpp*cpp.W_Psi;
		
		Fx = (- Vy*Bz - By )/(1+Psi);
		Fy = (  Vx*Bz + Bx )/(1+Psi);

		dztmp=dzz;
		Vxp = p-> old_vx + Fx*dztmp;
		Vyp = p-> old_vy + Fy*dztmp;

		p-> old_vx = Vxp;
		p-> old_vy = Vyp;
		
		p-> Vxx = Vxp*Vxp;
		p-> Vxy = Vxp*Vyp;
		p-> Vyy = Vyp*Vyp;

		p = p->p_PrevTraj;


	}

	return;

}

Trajectory* Mesh::Reconnect(Trajectory* p_Traj)
{

	Trajectory* p_temp;

	if(p_Traj->p_PrevTraj)
	{
		p_Traj->p_PrevTraj->p_NextTraj = p_Traj->p_NextTraj;
		p_temp = p_Traj->p_PrevTraj;
	}
	else
	{
		p_temp = NULL;
	}

	if(p_Traj->p_NextTraj)
	{
		p_Traj->p_NextTraj->p_PrevTraj = p_Traj->p_PrevTraj;
	}
	else
	{
		p_Trajectory = p_Traj->p_PrevTraj;

	}

	delete p_Traj;
	return p_temp;

}




void Mesh::PackT(Trajectory* p_Traj, int Sendn, int where)
{
	Commute *p_COMM = p_domain()->p_Com();

	switch(where)
	{

		case 0:
		p_COMM->SendSourceXm[(Sendn-1)*SDT_DIM + 0] = p_Traj-> x;
		p_COMM->SendSourceXm[(Sendn-1)*SDT_DIM + 1] = p_Traj-> y;
		p_COMM->SendSourceXm[(Sendn-1)*SDT_DIM + 2] = p_Traj-> x0;
		p_COMM->SendSourceXm[(Sendn-1)*SDT_DIM + 3] = p_Traj-> y0;
		p_COMM->SendSourceXm[(Sendn-1)*SDT_DIM + 4] = p_Traj-> z0;
		p_COMM->SendSourceXm[(Sendn-1)*SDT_DIM + 5] = p_Traj-> Vx;
		p_COMM->SendSourceXm[(Sendn-1)*SDT_DIM + 6] = p_Traj-> Vy;
		p_COMM->SendSourceXm[(Sendn-1)*SDT_DIM + 7] = p_Traj-> old_x;
		p_COMM->SendSourceXm[(Sendn-1)*SDT_DIM + 8] = p_Traj-> old_y;
		p_COMM->SendSourceXm[(Sendn-1)*SDT_DIM + 9] = p_Traj-> old_vx;
		p_COMM->SendSourceXm[(Sendn-1)*SDT_DIM +10] = p_Traj-> old_vy;
		p_COMM->SendSourceXm[(Sendn-1)*SDT_DIM +11] = p_Traj-> sx;
		p_COMM->SendSourceXm[(Sendn-1)*SDT_DIM +12] = p_Traj-> sy;

		break;
	
		case 1:
		p_COMM->SendSourceXp[(Sendn-1)*SDT_DIM + 0] = p_Traj-> x;
		p_COMM->SendSourceXp[(Sendn-1)*SDT_DIM + 1] = p_Traj-> y;
		p_COMM->SendSourceXp[(Sendn-1)*SDT_DIM + 2] = p_Traj-> x0;
		p_COMM->SendSourceXp[(Sendn-1)*SDT_DIM + 3] = p_Traj-> y0;
		p_COMM->SendSourceXp[(Sendn-1)*SDT_DIM + 4] = p_Traj-> z0;
		p_COMM->SendSourceXp[(Sendn-1)*SDT_DIM + 5] = p_Traj-> Vx;
		p_COMM->SendSourceXp[(Sendn-1)*SDT_DIM + 6] = p_Traj-> Vy;
		p_COMM->SendSourceXp[(Sendn-1)*SDT_DIM + 7] = p_Traj-> old_x;
		p_COMM->SendSourceXp[(Sendn-1)*SDT_DIM + 8] = p_Traj-> old_y;
		p_COMM->SendSourceXp[(Sendn-1)*SDT_DIM + 9] = p_Traj-> old_vx;
		p_COMM->SendSourceXp[(Sendn-1)*SDT_DIM +10] = p_Traj-> old_vy;
		p_COMM->SendSourceXp[(Sendn-1)*SDT_DIM +11] = p_Traj-> sx;
		p_COMM->SendSourceXp[(Sendn-1)*SDT_DIM +12] = p_Traj-> sy;
		break;

		case 2:
		p_COMM->SendSourceYm[(Sendn-1)*SDT_DIM + 0] = p_Traj-> x;
		p_COMM->SendSourceYm[(Sendn-1)*SDT_DIM + 1] = p_Traj-> y;
		p_COMM->SendSourceYm[(Sendn-1)*SDT_DIM + 2] = p_Traj-> x0;
		p_COMM->SendSourceYm[(Sendn-1)*SDT_DIM + 3] = p_Traj-> y0;
		p_COMM->SendSourceYm[(Sendn-1)*SDT_DIM + 4] = p_Traj-> z0;
		p_COMM->SendSourceYm[(Sendn-1)*SDT_DIM + 5] = p_Traj-> Vx;
		p_COMM->SendSourceYm[(Sendn-1)*SDT_DIM + 6] = p_Traj-> Vy;
		p_COMM->SendSourceYm[(Sendn-1)*SDT_DIM + 7] = p_Traj-> old_x;
		p_COMM->SendSourceYm[(Sendn-1)*SDT_DIM + 8] = p_Traj-> old_y;
		p_COMM->SendSourceYm[(Sendn-1)*SDT_DIM + 9] = p_Traj-> old_vx;
		p_COMM->SendSourceYm[(Sendn-1)*SDT_DIM +10] = p_Traj-> old_vy;
		p_COMM->SendSourceYm[(Sendn-1)*SDT_DIM +11] = p_Traj-> sx;
		p_COMM->SendSourceYm[(Sendn-1)*SDT_DIM +12] = p_Traj-> sy;
		break;

		case 3:
		p_COMM->SendSourceYp[(Sendn-1)*SDT_DIM + 0] = p_Traj-> x;
		p_COMM->SendSourceYp[(Sendn-1)*SDT_DIM + 1] = p_Traj-> y;
		p_COMM->SendSourceYp[(Sendn-1)*SDT_DIM + 2] = p_Traj-> x0;
		p_COMM->SendSourceYp[(Sendn-1)*SDT_DIM + 3] = p_Traj-> y0;
		p_COMM->SendSourceYp[(Sendn-1)*SDT_DIM + 4] = p_Traj-> z0;
		p_COMM->SendSourceYp[(Sendn-1)*SDT_DIM + 5] = p_Traj-> Vx;
		p_COMM->SendSourceYp[(Sendn-1)*SDT_DIM + 6] = p_Traj-> Vy;
		p_COMM->SendSourceYp[(Sendn-1)*SDT_DIM + 7] = p_Traj-> old_x;
		p_COMM->SendSourceYp[(Sendn-1)*SDT_DIM + 8] = p_Traj-> old_y;
		p_COMM->SendSourceYp[(Sendn-1)*SDT_DIM + 9] = p_Traj-> old_vx;
		p_COMM->SendSourceYp[(Sendn-1)*SDT_DIM +10] = p_Traj-> old_vy;	
		p_COMM->SendSourceYp[(Sendn-1)*SDT_DIM +11] = p_Traj-> sx;
		p_COMM->SendSourceYp[(Sendn-1)*SDT_DIM +12] = p_Traj-> sy;
		break;
	}

	return;
}



void Mesh::ExchangeT()
{
	Trajectory *p = NULL;
	
	//=========Send and Receive Buf Size===============
	double bufsize = p_domain()->p_Com()->Get_bufsize();
	bufsize *= (GridX*SOU_DIM*2.0/SDT_DIM);
	//=================================================

	double xtp, ytp;

	int Xpa  = p_domain()->p_Partition()->GetXpart();
	int Ypa  = p_domain()->p_Partition()->GetYpart();

	int Sendxm, Sendxp, Sendym, Sendyp;
	int S_SUM, A_SUM;
	int i;
	int j;

	while(1)
	{

		Sendxm = Sendxp = Sendym = Sendyp = 0; 
		p = p_Trajectory;
		while (p)
		{
			xtp = p-> x;
			ytp = p-> y;

			i=p->idx_i;
			j=p->idx_j;

		//=================================================
		//============Trajectory Outside Boundary =========
		//=================================================
		if(p_domain()->Get_BC()==1)
		{
		if(RankIdx_X ==1	&& i==0)
		{ p = p->p_PrevTraj; continue;}
		if(RankIdx_X == Xpa && i==GridX+1)
		{ p = p->p_PrevTraj; continue;}
		if(RankIdx_Y == 1	&& j==0)
		{ p = p->p_PrevTraj; continue;}
		if(RankIdx_Y == Ypa && j==GridY+1)
		{ p = p->p_PrevTraj; continue;}
		}
		//==================================================


		//====================================
		//====== Send to left Neighbor =======
		//====================================
		if(i==0)
		{
			Sendxm +=1;
			PackT(p, Sendxm, 0);
			p = Reconnect(p);
		}
		//====================================
		//====== Send to right Neighbor ======
		//====================================
		else if(i==GridX+1)
		{
			Sendxp +=1;
			PackT(p, Sendxp, 1);
			p = Reconnect(p);
		}
		//====================================
		//====== Send to Neighbor Below ======
		//====================================
		else if(j==0)
		{
			Sendym +=1;
			PackT(p, Sendym, 2);
			p = Reconnect(p);
		}
		//====================================
		//====== Send to Neighbor Above ======
		//====================================
		else if(j==GridY+1)
		{
			Sendyp +=1;
			PackT(p, Sendyp, 3);
			p = Reconnect(p);
		}
		else
		{
			p = p->p_PrevTraj;
		}

		}

		//=====================================================================
		//======== Exchange Trajectoryies with Neighboring Processors =========
		//=====================================================================
		if(bufsize<Sendxm || bufsize<Sendxp || bufsize<Sendym || bufsize<Sendyp)
		{	
			printf("==== Mesh: At Rank: %5d. ==================\n",Rank);
			std::cout << "==== Mesh: Send Too Many Trajectoryies.  ====\n";
			std::cout << "==== Mesh: May Cause Memory Problems.    ====\n";
			std::cout << "==== Mesh: Try to Increase the Buf Size. ====\n";
		}

		S_SUM = Sendxm+Sendxp+Sendym+Sendyp;
		MPI_Allreduce(&S_SUM, &A_SUM, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		if( A_SUM == 0 ) {break;};
		p_domain()->p_Com()->DoCommuteT(COMMU_T, Sendxm, Sendxp, Sendym, Sendyp);
	}
	return;
}

void Mesh::AdjustZstep(double k0, int k, double &dz2dz)
{
	//==========================================
	//=========== Adapteive Z Step =============
	if(AdaptiveStep>0)
	{	
		double A_Vmax;
		MPI_Allreduce(&Vmax, &A_Vmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		if(A_Vmax>Vlim) 
 		{
			dzz = dz/(A_Vmax/Vlim);
		}
		else
		{
			dzz = dz;
		}
		if(k0 + dzz/dz > k &dzz<dz &(k!=k0)) dzz=(k-k0)*dz;
	}

	dz2dz= dzz/dz;
	//==========================================
	return;
}

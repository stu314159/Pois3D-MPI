//C++ includes
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

//MPI includes
#include <mpi.h>

// My includes
#include "lattice_vars.h"
#include "vtk_lib.h"
#include "lbm_lib.h"

#define HALO 1

using namespace std;

int main(int argc, char * argv[]){

  int rank, size, rc;
  MPI_Status stat;
  MPI_Request rq_in1, rq_in2;
  MPI_Request rq_out1,rq_out2;
  int tag_d = 666;
  int tag_u = 999;

  double time_start, time_end, ex_time, LPU_sec, gNumLP;


  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  int nd_m, nd_p; //node minus and node plus
  nd_m = rank-1;
  if(nd_m<0)
    nd_m=(size-1);

  nd_p = rank+1;
  if(nd_p==size)
    nd_p=0;

 // declare input parameters
  int LatticeType;
  int Dynamics;
  int Entropic;
  int Initialization;
  int Num_ts;
  int ts_rep_freq;
  int plot_freq;
  int obst_type;
  float obst_param1;
  float obst_param2;
  float obst_param3;
  float obst_param4;
  float rho_lbm;
  float umax_lbm;
  float omega;
  int Nx;
  int Ny;
  int Nz;
  float t_conv_fact;
  float l_conv_fact;

  // read the input file
  ifstream input_params("params.lbm",ios::in);
  input_params >> LatticeType;
  input_params >> Dynamics;
  input_params >> Entropic;
  input_params >> Initialization;
  input_params >> Num_ts;
  input_params >> ts_rep_freq;
  input_params >> plot_freq;
  input_params >> obst_type;
  input_params >> obst_param1;
  input_params >> obst_param2;
  input_params >> obst_param3;
  input_params >> obst_param4;
  input_params >> rho_lbm;
  input_params >> umax_lbm;
  input_params >> omega;
  input_params >> Nx;
  input_params >> Ny;
  input_params >> Nz;
  input_params >> t_conv_fact;
  input_params >> l_conv_fact;

  int numSpd;
  switch (LatticeType){
  case(1):
    numSpd = 15;
    break;

  case(2):
    numSpd = 19;
    break;

  case(3):
    numSpd=27;
    break;

  }
  float * ex; float * ey; float * ez; float * w; int * bb_spd;

  int numPspeeds;
  int * Pspeeds;
  int numMspeeds;
  int * Mspeeds;
 switch (numSpd){
    //lattice parameters are defined in lattice_vars.h
  case (15):
    ex = ex15;
    ey = ey15;
    ez = ez15;
    w = w15;
    bb_spd = bb15;
    // speeds for sharing streaming data.
    //Pspeeds stream from lower ranks to higher (p_out to m_in)
    //Mspeeds stream from higher ranks to lower (m_out to p_in)
    numPspeeds=numPspeedsD3Q15;
    Pspeeds = PspeedsD3Q15;//<--zero-based
    numMspeeds=numMspeedsD3Q15;
    Mspeeds = MspeedsD3Q15;//<-- zeros-based
    
    //M = M15;
    break;

  case(19):
    ex = ex19;
    ey = ey19;
    ez = ez19;
    w = w19;
    bb_spd=bb19;
    //M = M19;
    break;

  case(27):
    ex = ex27;
    ey = ey27;
    ez = ez27;
    w = w27;
    bb_spd = bb27;
    //M = NULL; //none implemented for 27-speed model yet.
    break;

  }


  // prepare local partition variables fIn,fOut,inl,onl,snl
 // sort out geometric partition
  int numMySlices = Nz/size;
  if(rank<(Nz%size))//<-- floor Nz/size
    numMySlices+=1;//<-- add one starting at rank 0 for remainder


  int numMyNodes = Nx*Ny*numMySlices;

  int firstSlice, lastSlice;
  firstSlice=(Nz/size)*rank; //<-- starting point

  firstSlice=(Nz/size)*rank;//<--starting point
  if((Nz%size)<rank){ //<-- add the minimum of the rank or Nz%size
    firstSlice+=(Nz%size);
  }else{
    firstSlice+=rank;
  }

  lastSlice=firstSlice+numMySlices-1;
  int totalSlices=numMySlices+2*HALO;//<-- this can be used like Nz for a local partition
  int nnodes = totalSlices*Nx*Ny; //local value for nnodes.

 
  // declare and allocate fIn and fOut
  float * fEven = new float[nnodes*numSpd];
  float * fOdd = new float[nnodes*numSpd];
  // declare and allocate snl, inl and onl
  int * snl = new int[nnodes];
  int * inl = new int[nnodes];
  int * onl = new int[nnodes];

  initialize_lattice_partition(fOdd,rho_lbm,w,Nx,Ny,totalSlices,numSpd);
  initialize_lattice_partition(fEven,rho_lbm,w,Nx,Ny,totalSlices,numSpd);
  //initialize_snl_partition(snl,Nx,Ny,totalSlices,numSpd);



 
  //-- for inl and onl, this works.  Stay with it. ----

  //initialize inl and onl to zero
  int tid_l,z,y;
  for(int nd=0;nd<nnodes;nd++){
    inl[nd]=0;
    onl[nd]=0;
    snl[nd]=0;
  }

  //y==0, snl=1
  y=0;
  for(int z=0;z<totalSlices;z++){
    for(int x=0;x<Nx;x++){
      tid_l=x+y*Nx+z*Nx*Ny;
      snl[tid_l]=1;
    }
  }

  y=(Ny-1);
  for(int z=0;z<totalSlices;z++){
    for(int x=0;x<Nx;x++){
      tid_l=x+y*Nx+z*Nx*Ny;
      snl[tid_l]=1;
    }
  }

  //determine if I have a slice on the inlet.
  //stop being cute... second slice of rank 0 is on the inlet.  last slice of rank (size-1) is on the inlet.
  if(rank==0){
    z=1;
    for(int y=1;y<(Ny-1);y++){//<-- skip top and bottom
      for(int x=0;x<Nx;x++){
  	tid_l = x+y*Nx+z*Nx*Ny;
  	inl[tid_l]=1;
      }
    }
  }
  if(rank==(size-1)){
    z=totalSlices-1;
    for(int y=1;y<(Ny-1);y++){
      for(int x=0;x<Nx;x++){
  	tid_l=x+y*Nx+z*Nx*Ny;
  	inl[tid_l]=1;
      }
    }
  }


  //first slice of rank 0 is on the outlet.  second to last slice of rank (size-1) is on the outlet
  if(rank==0){
    z=0;
    for(int y=1;y<(Ny-1);y++){
      for(int x=0;x<Nx;x++){
  	tid_l=x+y*Nx+z*Nx*Ny;
  	onl[tid_l]=1;
      }
    }
  }
  if(rank==(size-1)){
    z=totalSlices-2;
    for(int y=1;y<(Ny-1);y++){
      for(int x=0;x<Nx;x++){
  	tid_l=x+y*Nx+z*Nx*Ny;
  	onl[tid_l]=1;
      }
    }

  }

  float * u_bc = new float[nnodes];
  initialize_ubc_pois3D(u_bc,umax_lbm,inl,onl,Nx,Ny,totalSlices);

  int firstNdm = Nx*Ny*HALO;
  int lastNdm = Nx*Ny*(HALO+1);
  
  int firstNdp = nnodes-2*(Nx*Ny*HALO);
  int lastNdp = nnodes-(Nx*Ny*HALO);

//get pointers to the start of the ghost array data
// not in this case, the "out" pointers point to the HALO area
// where streamed information to be communicated is stored.
// the "in" pointers point to boundary slice data where
// the streamed data is to be communicated to.
  float * ghost_out_odd_p = fOdd+(Nx*Ny*(totalSlices-HALO)*numSpd);
  float * ghost_out_odd_m = fOdd;
  float * ghost_in_odd_p = fOdd+(Nx*Ny*numSpd*numMySlices);
  float * ghost_in_odd_m = fOdd+(Nx*Ny*numSpd*HALO);

  float * ghost_out_even_p = fEven+(Nx*Ny*(totalSlices-HALO)*numSpd);
  float * ghost_out_even_m = fEven;
  float * ghost_in_even_p = fEven+(Nx*Ny*numSpd*numMySlices);
  float * ghost_in_even_m = fEven+(Nx*Ny*numSpd*HALO);

  int numHALO = (Nx*Ny*numPspeeds*HALO);//<--number of values in each HALO buffer.

  //ghost buffers for streamed data
  float * ghost_in_m; ghost_in_m = new float[Nx*Ny*numMspeeds*HALO];
  float * ghost_out_m;ghost_out_m=new float[Nx*Ny*numMspeeds*HALO];
  float * ghost_in_p;ghost_in_p=new float[Nx*Ny*numPspeeds*HALO];
  float * ghost_out_p;ghost_out_p=new float[Nx*Ny*numPspeeds*HALO];





  // visualization stuff...
  string densityFileStub("density");
  string fileSuffix(".vtk");
  stringstream ts_ind;
  string ts_ind_str;
  int vtk_ts = 0;
  string fileName;
  string dataName("densityMPI");
  int dims[3];
  dims[0]=Nx; dims[1]=Ny; dims[2]=Nz;
  float origin[3];
  origin[0]=0.; origin[1]=0.; origin[2]=0.;
  float spacing[3];
  spacing[0]=l_conv_fact; spacing[1]=l_conv_fact; spacing[2]=l_conv_fact;

  float * rho_l = new float[Nx*Ny*numMySlices];
  float * rho_g;

  if(rank==0){
    rho_g = new float[Nx*Ny*Nz];
  }

  //temporary variables for visualization
  float tmp_rho,tmp_ux;
  int tid_g;

  // save initial data
  //save density
  for(int z=HALO;z<(totalSlices-HALO);z++){
    for(int y=0;y<Ny;y++){
      for(int x=0;x<Nx;x++){
  	tid_l=x+y*Nx+(z-HALO)*Nx*Ny;
  	tmp_rho=0; tid_g=x+y*Nx+z*Nx*Ny;
  	tmp_ux=0.;
  	for(int spd=0;spd<numSpd;spd++){
  	  tmp_rho+=fEven[spd+tid_g*numSpd];
  	}
  	rho_l[tid_l]=tmp_rho;
      }
    }
  }
  MPI_Gather(rho_l,numMySlices*Nx*Ny,MPI_FLOAT,rho_g,
  	     numMySlices*Nx*Ny,MPI_FLOAT,0,MPI_COMM_WORLD);

  if(rank==0){
    ts_ind << vtk_ts;
	
    fileName=densityFileStub+ts_ind.str()+fileSuffix;
    ts_ind.str("");
    SaveVTKImageData_ascii(rho_g,fileName,dataName,origin,spacing,dims);
  }

  //write files for other data (like velocity) if desired...

  //after other files are written, increment the vtk_ts.
  if(rank==0){
    vtk_ts+=1;
  }



  time_start=MPI_Wtime();




  for(int ts=0;ts<Num_ts;ts++){
    //say something comforting about my progress
    if((ts+1)%ts_rep_freq==0){
      if(rank==0){
	cout << "Executing time step number " << ts+1 << endl;
      }
    }


    if(ts%2==0){
      //even time step 
    //collide/stream lower boundary slices
    ts_pois3D_D3Q15_LBGK_r2(fEven,fOdd,snl,inl,onl,u_bc,omega,
		    Nx,Ny,HALO,HALO+1);
    // collect data in ghost_m_out
    stream_out_collect(ghost_out_odd_m,ghost_out_m,numSpd,numMspeeds,
		       Mspeeds,Nx,Ny,totalSlices,HALO);
    //begin communication to ghost_p_in
    MPI_Isend(ghost_out_m,numHALO,MPI_FLOAT,nd_m,tag_d,MPI_COMM_WORLD,&rq_out1);
    MPI_Irecv(ghost_in_p,numHALO,MPI_FLOAT,nd_p,tag_d,MPI_COMM_WORLD,&rq_in1);


    //collide/stream upper boundary slices
    ts_pois3D_D3Q15_LBGK_r2(fEven,fOdd,snl,inl,onl,u_bc,omega,
		    Nx,Ny,totalSlices-2*HALO,totalSlices-HALO);

    //collect data in ghost_p_out
    stream_out_collect(ghost_out_odd_p,ghost_out_p,numSpd,numPspeeds,
		       Pspeeds,Nx,Ny,totalSlices,HALO);
    //begin communication to ghost_m_in
    MPI_Isend(ghost_out_p,numHALO,MPI_FLOAT,nd_p,tag_u,MPI_COMM_WORLD,&rq_out2);
    MPI_Irecv(ghost_in_m,numHALO,MPI_FLOAT,nd_m,tag_u,MPI_COMM_WORLD,&rq_in2);

    //collide/stream interior lattice points
    ts_pois3D_D3Q15_LBGK_r2(fEven,fOdd,snl,inl,onl,u_bc,omega,
		    Nx,Ny,HALO+1,totalSlices-2*HALO);
    //ensure communication of boundary lattice points is complete
    MPI_Wait(&rq_in1,&stat);
    MPI_Wait(&rq_in2,&stat);
    //copy data from ghost_in_p to Mspeeds on P boundary points
    stream_in_distribute(ghost_in_odd_p,ghost_in_p,numSpd,
    			 numMspeeds,Mspeeds,Nx,Ny,totalSlices,HALO);
    //copy data from ghost_m_in to Pspeeds on M boundary points
    stream_in_distribute(ghost_in_odd_m,ghost_in_m,numSpd,
    			 numPspeeds,Pspeeds,Nx,Ny,totalSlices,HALO);

    }else{
      //odd time step
      //collide/stream lower boundary slices
      ts_pois3D_D3Q15_LBGK_r2(fOdd,fEven,snl,inl,onl,u_bc,omega,
			   Nx,Ny,HALO,HALO+1);
      // collect data in ghost_m_out
      stream_out_collect(ghost_out_even_m,ghost_out_m,numSpd,numMspeeds,
			 Mspeeds,Nx,Ny,totalSlices,HALO);
      //begin communication to ghost_p_in
      MPI_Isend(ghost_out_m,numHALO,MPI_FLOAT,nd_m,tag_d,MPI_COMM_WORLD,&rq_out1);
      MPI_Irecv(ghost_in_p,numHALO,MPI_FLOAT,nd_p,tag_d,MPI_COMM_WORLD,&rq_in1);


      //collide/stream upper boundary slices
      ts_pois3D_D3Q15_LBGK_r2(fOdd,fEven,snl,inl,onl,u_bc,omega,
			   Nx,Ny,totalSlices-2*HALO,totalSlices-HALO);

      //collect data in ghost_p_out
      stream_out_collect(ghost_out_even_p,ghost_out_p,numSpd,numPspeeds,
			 Pspeeds,Nx,Ny,totalSlices,HALO);
      //begin communication to ghost_m_in
      MPI_Isend(ghost_out_p,numHALO,MPI_FLOAT,nd_p,tag_u,MPI_COMM_WORLD,&rq_out2);
      MPI_Irecv(ghost_in_m,numHALO,MPI_FLOAT,nd_m,tag_u,MPI_COMM_WORLD,&rq_in2);

      //collide/stream interior lattice points
      ts_pois3D_D3Q15_LBGK_r2(fOdd,fEven,snl,inl,onl,u_bc,omega,
			   Nx,Ny,HALO+1,totalSlices-2*HALO);
      //ensure communication of boundary lattice points is complete
      MPI_Wait(&rq_in1,&stat);
      MPI_Wait(&rq_in2,&stat);
      //copy data from ghost_p_in to Mspeeds on P boundary points
      stream_in_distribute(ghost_in_even_p,ghost_in_p,numSpd,
      			   numMspeeds,Mspeeds,Nx,Ny,totalSlices,HALO);
      //copy data from ghost_m_in to Pspeeds on M boundary points
      stream_in_distribute(ghost_in_even_m,ghost_in_m,numSpd,
      			   numPspeeds,Pspeeds,Nx,Ny,totalSlices,HALO);


    }


    //at selected time step intervals, gather and write data to vtk file for
    //visualization
    if((ts+1)%plot_freq==0){

      //save density
      if(ts%2==1){
	for(int z=HALO;z<(totalSlices-HALO);z++){
	  for(int y=0;y<Ny;y++){
	    for(int x=0;x<Nx;x++){
	      tid_l=x+y*Nx+(z-HALO)*Nx*Ny;
	      tmp_rho=0; tid_g=x+y*Nx+z*Nx*Ny;
	      tmp_ux=0.;
	      for(int spd=0;spd<numSpd;spd++){
		tmp_rho+=fEven[spd+tid_g*numSpd];
	      }
	      rho_l[tid_l]=tmp_rho;
	    }
	  }
	}
      }else{
	for(int z=HALO;z<(totalSlices-HALO);z++){
	  for(int y=0;y<Ny;y++){
	    for(int x=0;x<Nx;x++){
	      tid_l=x+y*Nx+(z-HALO)*Nx*Ny;
	      tmp_rho=0; tid_g=x+y*Nx+z*Nx*Ny;
	      tmp_ux=0.;
	      for(int spd=0;spd<numSpd;spd++){
		tmp_rho+=fOdd[spd+tid_g*numSpd];
	      }
	      rho_l[tid_l]=tmp_rho;
	    }
	  }
	}

      }
      MPI_Gather(rho_l,numMySlices*Nx*Ny,MPI_FLOAT,rho_g,
    		 numMySlices*Nx*Ny,MPI_FLOAT,0,MPI_COMM_WORLD);

      if(rank==0){
    	ts_ind << vtk_ts;
	
    	fileName=densityFileStub+ts_ind.str()+fileSuffix;
    	ts_ind.str("");
    	SaveVTKImageData_ascii(rho_g,fileName,dataName,origin,spacing,dims);
      }

      //write other files...

      //after other files are written, increment the vtk_ts.
      if(rank==0){
    	vtk_ts+=1;
      }

    }


  }//end time step

  time_end = MPI_Wtime();

  if(rank==0){

    ex_time = time_end-time_start;
    gNumLP = Nx*Ny*Nz;
    LPU_sec = ((double)gNumLP* (double)Num_ts)/ex_time;
    cout << "Estimated LPU/sec = " << LPU_sec << endl;


  }



  //clean up dynamically allocated memory
  delete [] fEven;
  delete [] fOdd;
  delete [] snl;
  delete [] inl;
  delete [] onl;
  delete [] u_bc;
  delete [] rho_l;

  delete [] ghost_in_m;
  delete [] ghost_out_m;
  delete [] ghost_in_p;
  delete [] ghost_out_p;

  if(rank==0){
    delete [] rho_g;
  }

  MPI_Finalize();
  return 0;

}

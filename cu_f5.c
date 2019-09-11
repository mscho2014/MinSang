/* 
 Cu case calculation  
	Copper hole population calculation
	- modification after Hamburg's visiting (190701) 

	[190704]
	almost final version of copper calculation! 
	put (1) collisional process  (2) electron-photon coupling  (3) n4 electron excitation at 000 fs 


	[190812]
	professor's idea! 


*/ 

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<gsl-2.5/gsl/gsl_complex.h>
#include<gsl-2.5/gsl/gsl_complex_math.h>
#include<gsl-2.5/gsl/gsl_linalg.h>
#include<gsl-2.5/gsl/gsl_matrix.h>
#include<gsl-2.5/gsl/gsl_blas.h>
#include<string.h> 


int main(void){
FILE *fp0;
FILE *fp1;
FILE *fp2; 
FILE *fp3;
 
char NAME[50];

/* Set-up condition variables */
        int nam ; // set the filename

int i, t;// i is energy grid and tm is time grid
int j;
int mxsize = 6; // matrix size
int vxsize; 
	vxsize = mxsize-1; 
int s; 

double dt = 1e-2;// time grid [fs]
double dt_r ;
	dt_r = (1e-15)*dt;
double tm; 
int tma; 
double tmax; 
	
double n[vxsize];
double ne[vxsize]; // electron population 
double nm[vxsize]; // maximun hole population 

double n_temp[vxsize]; // update hole population 
double n_temp2[vxsize];


/* New algorithm */ 
double pex; 
double pex_temp;
double jj;
	


/* pumping laser condition */
double convert = 1.6e-19; // leV = 1.6e-19 J
double JR = 0.1;// (0.1)  //graph laser fluence : 0.015 J/cm2 / real laser fluence : 0.23 J/cm2 (absorbed) <- real incident fluence is 0.46 J/cm2 ///  Absorbed Laser fluence [J/cm2]
double pulse = 80; // laser pulse duration [fs]
double p_center = 50; //Gaussian laser pulse center [fs]
double Ep ; // laser wavelength unit change [nm] -> [eV] 
	Ep = 1239.84/400;

//	JR = jj; 
double Jv;
	Jv = JR/(pulse*1e-15)/Ep/convert;

/* copper material coefficient */ 
double Nav = 6.02e23 ; // Avogadro's number 
double Md = 1/63.546; // Cu [mole/g]
double cumd = 8.96; // Cu [g/cm3]
double ceff_ab = 1.43e5; //1.43 for 70nm, 7.09 for 14nm cu absorption coefficient at 400 nm 
double sgm_ab; // absorption cross-section = (absorption coefficient)/(atomic number density) 
	sgm_ab = ceff_ab/(Nav*Md*cumd);

	
double ntot = 8; // total electron number : ntot = n1 + n2 + n3

double tau1;
double tau2;
double tau3; 
double tauc; 
double tauep; 

int tau_1; // input parameter for tau_1
int tau_2; // input parameter for tau_2
int tau_3; // input parameter for tau_3
double tau_c; 
double tau_ep;
 
pex_temp = Jv*sgm_ab;


/* Define initial distribution of Hole and make an initial condition file f_0fs */


/* Maximum value of Hole population condition */
double n0 ; // max hole in excited states // max value : it was an infinity in principle, but n0 = n1 + n2 + n3 will be fine. 
double n1 = 3.2 ; // 3.2 / 2 (value : 0.5)  -1.5 ~ 0 eV max value  this value is related with the shape!!!! 
	// n1 = 1.5 ; which will be discussed 
double n2 = 3.4; // -2eV max value
double n3 = 3.2; // -3eV max value 
double n4 = 1.97;

	n0 = n1 + n2 + n3 + n4; 

  nm[0] = n0;
  nm[1] = n1;
  nm[2] = n2;
  nm[3] = n3;
  nm[4] = n4; 

/* Define initial distribution of Hole and make an initial condition file f_0fs */
  n[0] = n0; // holes in excited states  
  n[1] = 0; // -1.5eV ~ 0 eV  hole #
  n[2] = 0; // -2.5 ~ -1.5 eV hole #
  n[3] = 0; // -3.5 ~ -2.5 eV hole #
  n[4] = 0; 

  ne[0] = nm[0]-n[0];
  ne[1] = nm[1]-n[1];
  ne[2] = nm[2]-n[2];
  ne[3] = nm[3]-n[3];
  ne[4] = nm[4]-n[4];

 
double netot;
double rv_netot; 
int col_p; // collision process control parameter  (0 : no collision, 0.5 : 1/2 n1 population participates the collision, 1 : whole n1 participates the collision process) 
double vep;
double v_ep;
        printf("*****************************************************************\n");	
	printf("(0/3)Let's setup basic simulation parameter \n0-1)Filename please \n");
	scanf("%d",&nam);
	printf("Result file name will be 'rf_h%d' for hole and 'rf_e%d' for e\n",nam,nam); 
	printf("calculation time?[fs] \n");
        
	//scaf("%d",&tma);
	tma = 4000;
	tmax = tma/dt; 
	printf("calculation time is ");     
	printf("%d fs \n",tma);
	

        printf("what is tau2?[fs] : this is hole-lifetime from -2eV to conduction (col off : 500 fs / col on : 300 fs) \n");
	scanf("%d",&tau_2);
	//tau_2 = 800;
	tau2 = tau_2*1e-15; 
        printf("tau2 = ");
        printf("%d fs \n",tau_2);

        printf("what is tau3?[fs] : this is hole-lifetime from -3eV to -2eV (col on/off : 300 fs) \n");
        scanf("%d",&tau_3);
	//tau_3 = 300;
	tau3 = tau_3*1e-15; 
	printf("tau3 = ");
	printf("%d fs \n",tau_3);
        


        printf("*****************************************************************\n");
	printf("(1/3)Do you want to turn on collision process? (Yes : 1, No :0) \n");  
	scanf("%d",&col_p);
	//col_p = 1; 
	
	if (col_p ==1){
         printf("what is tau_c (collisional lifetime)? [fs] (10,000 fs or 30,000 fs?) \n" );
         scanf("%le", &tau_c);
	 //tau_c = 10000;
         tauc = tau_c*1e-15;
         }
         else{
         tau_c = 0.1;
         tauc = 0.1;
         }

        printf("*****************************************************************\n");
	printf("(2/3)Do you want to turn on electron-photon coupling? (Yes : 1, No :0) \n");
	//scanf("%le",&vep);
	vep = 0;
	//v_ep = vep; 
	
	if (vep == 0 ){
	 tau_ep = 0.1;
	 tauep = tau_ep*1e-15;
	 printf("electron-photon coupling turn off!!\n"); 
	}
	else{
	printf("what is the tau_ep (electron-photon lifetime)? [fs] (1000 fs?)\n");
        scanf("%le", &tau_ep);
	//tau_ep = 1000; 
        tauep = tau_ep*1e-15;
	} 

double ebeam; 
double t4st, t4ed;
double delt; 
double tau4, tau_4; 
int n4ch; 
double n4m; 
double n4p; 

        printf("*****************************************************************\n");
	printf("(3/3)n4 exciting on? (Yes:1, no :0)\n");
	scanf("%d",&n4ch);
	//n4ch = 0; 

	if(n4ch==1){
	printf("when n4 excitation on? [fs] (700 fs?)\n");
	//scanf("%le",&t4st);
	t4st = 700;
	printf("how long you want to turn on? [fs] (200 fs?) -> delt = 200 fs \n");
	//scanf("%le",&delt);
	delt = 200; 
	t4ed = t4st + delt;
	printf("what is tau4?[fs] (50 fs?) yes\n");
	scanf("%le",&tau_4);
	//tau_4 = 50; 
	tau4 = tau_4*1e-15;
	printf("what is n4 hole max-value? (col off : 0.15 / col on : 0.12) \n");
	scanf("%le",&n4m);
 	//n4m = 0.15;  
	}
	

	else{
		printf("n4 exciting turn off!!\n");
		ebeam = 0;
		t4st = 10000;
		t4ed = 20000;
		tau_4 = 50;
		tau4 = tau_4*1e-15; 
		n4m = 0;  
	}
	n4p=n4m*1e13; 

double n_1;
//	printf("n1\n");
//	scanf("%le",&n_1);
//	n1 = n_1; 

		

	sprintf(NAME,"Init_cond");
	fp0 = fopen(NAME,"w");
	
/* print out E-distribution at t=0 */
                       
        fprintf(fp0,"hole in excited : %le \n",n[0]);
        fprintf(fp0,"-1.5 ~ 0 eV (or sp band) : %le \n",n[1]);
	fprintf(fp0,"-2 : %le \n",n[2]);
	fprintf(fp0,"-3 : %le \n",n[3]);
	fprintf(fp0,"-4 : %le \n",n[4]); 
	fprintf(fp0,"tau2 : %le, tau3 = %le, calculation time = %d [fs], tauc = %le \n",tau2, tau3, tma, tauc);
	fclose(fp0);

/************** calculation **************/
	/* H : total / P : photo-excitation / C: collisional - deexcitation */ 
	gsl_matrix *P = gsl_matrix_alloc (mxsize,mxsize);

	gsl_vector *F = gsl_vector_alloc (mxsize);
	gsl_vector *F_t = gsl_vector_alloc (mxsize);
	gsl_vector *F_tt = gsl_vector_alloc (mxsize);

	gsl_vector *Fe = gsl_vector_alloc (mxsize); 
	
	/* Vectors in non-linear term */  
	gsl_vector *Bvec = gsl_vector_alloc (mxsize);
	gsl_vector *poad = gsl_vector_alloc (mxsize); // photo-excitation add
	gsl_vector *coad = gsl_vector_alloc (mxsize); // collisional excitation add 
	gsl_vector *epad = gsl_vector_alloc (mxsize); // electron-phonon coupling add
		
	/* matrix for linear algebra */ 
	gsl_permutation *p = gsl_permutation_alloc (mxsize);
	gsl_vector *FNDM = gsl_vector_alloc (mxsize);
	
	sprintf(NAME,"rf_h%d",nam);
	fp1 = fopen(NAME,"w");
	fprintf(fp1,"time(fs) hex h1 h2 h3 h4\n");
	fprintf(fp1,"0.000 %le %le %le %le %le \n",n[0],n[1],n[2],n[3],n[4]);

	sprintf(NAME,"rf_e%d",nam);
	fp2 = fopen(NAME,"w");
	fprintf(fp2,"time(fs) eex e1 e2 e3 e4\n");
	fprintf(fp2,"0.000 %le %le %le %le %le \n",ne[0],ne[1],ne[2],ne[3],ne[4]); 

	sprintf(NAME,"rf_tot%d",nam);
	fp3 = fopen(NAME,"w");
	fprintf(fp3,"time(fs) electron(ex) hole(con)\n");
	fprintf(fp3,"0.000 %le %le \n",ne[0],-(n[1]+n[2]+n[3]+n[4]));

	/* number test parameter : boundary condition */ 
	double e0;
	double e1;
	double e2;
	double e3;
	double e4; 

	double ee1;
	double ee2;
	double ee3;	
	double ee4; 


for(t = 1; t <= tmax ; ++t){
	tm = (double)(t)*dt; //

	pex = pex_temp*(pulse/40.5013)*exp(-pow((tm-p_center),2)*log(2)/pow((pulse/2),2));

	/* Initialization of the population in each time step t */ 
	for (i = 0; i <mxsize; ++i){
                        gsl_vector_set(F,i,n[i]);
			ne[i] = nm[i] - n[i];
			gsl_vector_set(Fe,i,ne[i]); 
        }
		//netot = n[1]+n[2]+n[3]+n[4]; // # of electrons in excited band and also # of holes below Fermi energy  

	/* number limitation */
     if (n[0] >= n0){e0 = 0;}
     else{e0 = 1;}

     if (n[1] >= n1){e1 = 0;}
     else{e1 = 1;}
     if (n[2] >= n2){e2 = 0;}	
     else{e2 = 1;}
     if (n[3] >= n3){e3 = 0;}
     else{e3 = 1;}
     if (n[4] >= n4){e4 = 0;}
     else{e4 = 1;} 

     if (n4ch == 1 && tm > t4st && tm < t4ed){ebeam = 1;}
     else {ebeam = 0;}

     if (n[1]<0){ee1 =0;}
     else{ee1=1;}
     if (n[2]<0){ee2=0;}
     else{ee2=1;}
     if (n[3]<0){ee3=0;}
     else{ee3=1;}
     if (n[4]<0){ee4=0;}
     else{ee4=1;}




	
	/* [Linear matrix calculation) : declare the calculation matrix - Decay term */ 

	gsl_matrix_set(P,0,0,1);
	gsl_matrix_set(P,0,1,0);
	gsl_matrix_set(P,0,2,(dt_r)*(-1/tau2)*e1);
	gsl_matrix_set(P,0,3,(dt_r)*(-1/tau3)*e3*e2);
	gsl_matrix_set(P,0,4,(dt_r)*(-1/tau4)*e3*e1);

	gsl_matrix_set(P,1,0,0);
	gsl_matrix_set(P,1,1,1);
        gsl_matrix_set(P,1,2,(dt_r)*(2/tau2)*e1);
	gsl_matrix_set(P,1,3,(dt_r)*(1/tau3)*e2*e1);
	gsl_matrix_set(P,1,4,(dt_r)*(1/tau4)*e3*e1);

        gsl_matrix_set(P,2,0,0);
        gsl_matrix_set(P,2,1,0);
        gsl_matrix_set(P,2,2,1+(dt_r)*(-1/tau2)*e1);	     
	gsl_matrix_set(P,2,3,(dt_r)*(1/tau3)*e2*e1);
	gsl_matrix_set(P,2,4,0);

	gsl_matrix_set(P,3,0,0);
        gsl_matrix_set(P,3,1,0);
        gsl_matrix_set(P,3,2,0);
        gsl_matrix_set(P,3,3,1+(dt_r)*(-1/tau3)*e2*e1);
	gsl_matrix_set(P,3,4,(dt_r)*(1/tau4)*e3*e1);

        gsl_matrix_set(P,4,0,0);
        gsl_matrix_set(P,4,1,0);
        gsl_matrix_set(P,4,2,0);
        gsl_matrix_set(P,4,3,0);
	gsl_matrix_set(P,4,4,1+(dt_r)*(-1/tau4)*e3*e1);

        gsl_blas_dgemv(CblasNoTrans, 1.0, P,F,0.0,F_t); // photo-excitation calculation [matrix multiplication]
 	
 
	netot = n[1]+n[2]+n[3]+n[4];
	
	if (netot<1e-50)
	{ rv_netot = 0; }
	else{rv_netot = 1/netot;}
	
	//printf("%le \n", rv_netot);
	/* Non-linear calculation : (1) Photo-excitation + n4 population pop-up + (2) collision + (3) el-ph coupling  */ 
	
	/* (1) Photo-excitation and n4 pop-up input*/ 
	gsl_vector_set(poad,0,-(dt_r)*pex*e3-(dt_r)*(ebeam*n4p)*e4);
	gsl_vector_set(poad,1,0);
	gsl_vector_set(poad,2,0);
	gsl_vector_set(poad,3,(dt_r)*(pex)*e3);
	gsl_vector_set(poad,4,(dt_r)*(ebeam*n4p)*e4); 

	gsl_vector_add(F_t,poad);


	/* (2) collisional excitation add */ 
	gsl_vector_set(coad,0,(dt_r)*(col_p/tauc)*netot*(ne[2]*e2*e0*ee1+ne[3]*e3*e0*ee2*ee1+ne[4]*e4*e0*ee3*ee1));
        gsl_vector_set(coad,1,(dt_r)*(-col_p/tauc)*netot*(2*ne[2]*e2*e0*ee1+ne[3]*e3*e0*ee2*ee1+ne[4]*e4*e0*ee3*ee1));
        gsl_vector_set(coad,2,(dt_r)*(col_p/tauc)*netot*ne[2]*e2*e0*ee1+(dt_r)*(-col_p/tauc)*netot*ne[3]*e3*e0*ee2*ee1);
        gsl_vector_set(coad,3,0+(dt_r)*(col_p/tauc)*netot*(ne[3]*e3*e0*ee2*ee1-ne[4]*e4*e0*ee3*ee1));
        gsl_vector_set(coad,4,0+0+(dt_r)*(col_p/tauc)*netot*ne[4]*e4*e0*ee3*ee1);

        gsl_vector_add(F_t,coad);

	
	/* (3) el-ph coupling term */   
	gsl_vector_set(epad,0,(dt_r)*(v_ep/tauep)*netot*e0);
        gsl_vector_set(epad,1,(dt_r)*(-v_ep/tauep)*n[1]*ee1);
        gsl_vector_set(epad,2,(dt_r)*(-v_ep/tauep)*n[2]*ee2);
        gsl_vector_set(epad,3,(dt_r)*(-v_ep/tauep)*n[3]*ee3);
        gsl_vector_set(epad,4,(dt_r)*(-v_ep/tauep)*n[4]*ee4);

        gsl_vector_add(F_t,epad);
	

         /* update n_temp[i] : distribution after photo-excitation */
         for (i=0 ; i <mxsize; ++i){
                                  double a;
                                  n_temp[i] = gsl_vector_get(F_t,i);
				  ne[i] = nm[i]-n_temp[i]; 
         }

	fprintf(fp1,"%le %le %le %le %le %le \n",tm,n_temp[0],n_temp[1],n_temp[2],n_temp[3],n_temp[4]);
	fprintf(fp2,"%le %le %le %le %le %le \n",tm,ne[0],ne[1],ne[2],ne[3],ne[4]);
	fprintf(fp3,"%le %le %le \n",tm,-ne[0],netot);  
  
	for (i = 0; i <mxsize; ++i){
		n[i] = n_temp[i];
	}	
	
	
	} // The end of the "t"
	
	fclose(fp1);
	fclose(fp2);
	fclose(fp3); 
	printf("calculation is finished\n");

return 0;
} // The end of the "main" 



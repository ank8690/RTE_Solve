#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<math.h>
#define sigma 5.67e-8
main()
{
	FILE *fp1,*fp2,*fp3;
	int i,j,l,m,ni,nj,nt,np;
	double Int_old[50][50][10][20],Int_new[50][50][10][20],*D_x,*D_y,*D_omega,dx,dy,dtheta,dphi,theta,phi,rms,Sum_n;
	double eps_b,eps_t,eps_l,eps_r,kappa,Tb,Tt,Tl,Tr,Tg,I_w,Nume,Deno,Source,Sum_o,flux_bot[50],flux_mid[50],q_w,q_u,q_l, PI;
	
	PI = 4.0*atan(1.0);
	printf("Enter the number of cells in X and Y directions\n");
	scanf("%d%d",&ni,&nj);
	printf("\nEnter the number of cells in Theta and Phi directions\n");
	scanf("%d%d",&nt,&np);
	printf("\nEnter the emissivity of Bottom, Top, Left and Right surfaces respectivily\n");
	scanf("%lf%lf%lf%lf",&eps_b,&eps_t,&eps_l,&eps_r);
	printf("\nEnter the Temperature of Bottom, Top, Left and Right surfaces in Kelvin respectivily\n");
	scanf("%lf%lf%lf%lf",&Tb,&Tt,&Tl,&Tr);
	printf("\nEnter the Temperature of Gas in kelvin\n");
	scanf("%lf",&Tg);
	printf("\nEnter the absorption coefficient of Gas\n");
	scanf("%lf",&kappa);
	fp1=fopen("intensity_Abs_Emit.dat","w");
	fp2=fopen("flux_Abs_Emit_mid.dat","w");
	fp3=fopen("flux_Abs_Emit_bot.dat","w");
	if(fp1==NULL)
	printf("Unable to open intensity file\n");
	if(fp2==NULL)
	printf("Unable to open flux Mid file\n");
	if(fp3==NULL)
	printf("Unable to open flux Bot file\n");
	dx=1.0/ni;
	dy=1.0/nj;
	dtheta=PI/(nt);
	dphi=2.0*PI/np;
		
		D_x=(double *)malloc((nt*np)*sizeof(double));
		D_y=(double *)malloc((nt*np)*sizeof(double));
		D_omega=(double *)malloc((nt*np)*sizeof(double));
		for(l=0;l<nt;l++)
			{
			 theta=dtheta*(l+0.5);
			 for(m=0;m<np;m++)
			 	{
				 phi=dphi*(m+0.5);
				 *(D_x+np*l+m)=cos(phi)*sin(dphi/2.0)*(dtheta-cos(2.0*theta)*sin(dtheta));
				 *(D_y+np*l+m)=sin(phi)*sin(dphi/2.0)*(dtheta-cos(2.0*theta)*sin(dtheta));
				 *(D_omega+np*l+m)=2.0*dphi*sin(theta)*sin(dtheta/2.0);
				 printf("%e\t%e\t%e\n",*(D_x+np*l+m)*dy,*(D_y+np*l+m)*dx,*(D_omega+np*l+m));
				  }
			}
			/*double x=0.0,y=0.0,z=0.0;
			for(l=nt/2;l<nt;l++)
				for(m=0;m<np;m++)
					{
					x=x+*(D_y+np*l+m);
					y=y+*(D_z+np*l+m);
					z=z+*(D_omega+np*l+m);
					}
			printf ("%e\t%e\t%e",x,y,z);
			int abc;
			scanf("%d",&abc);*/
					
			for(l=0;l<nt;l++)
				for(m=0;m<np;m++)
					for(j=0;j<=(nj+1);j++)
						for(i=0;i<=(ni+1);i++)
							Int_new[i][j][l][m]=0.0;
			
			do
			 {
			  for(l=0;l<nt;l++)
			  	for(m=0;m<np;m++)
			  		for(j=0;j<=(nj+1);j++)
						for(i=0;i<=(ni+1);i++)
							Int_old[i][j][l][m]=Int_new[i][j][l][m];
							
			
			for(l=0;l<nt;l++)
				for(m=0;m<np;m++)
					{
					 if(*(D_x+np*l+m)>0&&*(D_y+np*l+m)>0)
					 	{
						 for(j=1;j<=nj;j++)
						 	for(i=1;i<=ni;i++)
								{
								 Nume=dy*(fabs(*(D_x+np*l+m)))*(Int_new[i-1][j][l][m])+dx*(fabs(*(D_y+np*l+m)))*(Int_new[i][j-1][l][m]);
								 Source=kappa*sigma*pow(Tg,4)*dx*dy*(*(D_omega+np*l+m))/PI;
								 Deno=dy*(fabs(*(D_x+np*l+m)))+dx*(fabs(*(D_y+np*l+m)))+kappa*dx*dy*(*(D_omega+np*l+m));
								 (Int_new[i][j][l][m])=(Nume+Source)/Deno;
								 	if(Int_new[i][j][l][m]<0)
								 		Int_new[i][j][l][m]=0.0;
								 }
						 }
					if(*(D_x+np*l+m)>0&&*(D_y+np*l+m)<0)
						{
						 for(j=nj;j>=1;j--)
						 	for(i=1;i<=ni;i++)
								{
								 Nume=dy*(fabs(*(D_x+np*l+m)))*(Int_new[i-1][j][l][m])+dx*(fabs(*(D_y+np*l+m)))*(Int_new[i][j+1][l][m]);
								 Source=kappa*sigma*pow(Tg,4)*dx*dy*(*(D_omega+np*l+m))/PI;
								 Deno=dy*(fabs(*(D_x+np*l+m)))+dx*(fabs(*(D_y+np*l+m)))+kappa*dx*dy*(*(D_omega+np*l+m));
								 (Int_new[i][j][l][m])=(Nume+Source)/Deno;
								 	if(Int_new[i][j][l][m]<0)
								 		Int_new[i][j][l][m]=0.0;
								 }
						 }
					
					if(*(D_x+np*l+m)<0&&*(D_y+np*l+m)>0)
						{
						 for(j=1;j<=nj;j++)
						 	for(i=ni;i>=1;i--)
								{
								 Nume=dy*(fabs(*(D_x+np*l+m)))*(Int_new[i+1][j][l][m])+dx*(fabs(*(D_y+np*l+m)))*(Int_new[i][j-1][l][m]);
								 Source=kappa*sigma*pow(Tg,4)*dx*dy*(*(D_omega+np*l+m))/PI;
								 Deno=dy*(fabs(*(D_x+np*l+m)))+dx*(fabs(*(D_y+np*l+m)))+kappa*dx*dy*(*(D_omega+np*l+m));
								 (Int_new[i][j][l][m])=(Nume+Source)/Deno;
								 	if(Int_new[i][j][l][m]<0)
								 		Int_new[i][j][l][m]=0.0;
								 }
						 }
					
					if(*(D_x+np*l+m)<0&&*(D_y+np*l+m)<0)
					 	{
						 for(j=nj;j>=1;j--)
						 	for(i=ni;i>=1;i--)
								{
								 Nume=dy*(fabs(*(D_x+np*l+m)))*(Int_new[i+1][j][l][m])+dx*(fabs(*(D_y+np*l+m)))*(Int_new[i][j+1][l][m]);
								 Source=kappa*sigma*pow(Tg,4)*dx*dy*(*(D_omega+np*l+m))/PI;
								 Deno=dy*(fabs(*(D_x+np*l+m)))+dx*(fabs(*(D_y+np*l+m)))+kappa*dx*dy*(*(D_omega+np*l+m));
								 (Int_new[i][j][l][m])=(Nume+Source)/Deno;
								 	if(Int_new[i][j][l][m]<0)
								 		Int_new[i][j][l][m]=0.0;
								 }
						 }
					 }
			
			for(i=1;i<=ni;i++)
				{
				 I_w=0.0;
				for(l=0;l<nt;l++)
					for(m=np/2;m<np;m++)
						I_w=I_w-(1-eps_b)*(Int_new[i][1][l][m])*(*(D_y+np*l+m))/PI;
						 
				for(l=0;l<nt;l++)
					for(m=0;m<np/2;m++)
						(Int_new[i][0][l][m])=2.0*(I_w+eps_b*sigma*pow(Tb,4)/PI)-(Int_new[i][1][l][m]);
							if(Int_new[i][0][l][m]<0)
								 Int_new[i][0][l][m]=0.0;
				}
			
			for(i=1;i<=ni;i++)
				{
				 I_w=0.0;
				for(l=0;l<nt;l++)
					for(m=0;m<np/2;m++)
						I_w=I_w-(1-eps_t)*(Int_new[i][ni][l][m])*(*(D_y+np*l+m))/PI;
						 
				for(l=0;l<nt;l++)
					for(m=np/2;m<np;m++)
						(Int_new[i][nj+1][l][m])=2.0*(I_w+eps_t*sigma*pow(Tt,4)/PI)-(Int_new[i][nj][l][m]);
							if(Int_new[i][nj+1][l][m]<0)
								 Int_new[i][nj+1][l][m]=0.0;
				}
			
			for(j=1;j<=nj;j++)
				{
				 I_w=0.0;
				for(l=0;l<nt;l++)
					for(m=3*np/4;m>np/4;m--)
						I_w=I_w-(1-eps_l)*(Int_new[1][j][l][m])*(*(D_x+np*l+m))/PI;
						
				for(l=0;l<nt;l++)
					for(m=np/4;m<3*np/4;m++)
						(Int_new[0][j][l][m])=2.0*(I_w+eps_l*sigma*pow(Tl,4)/PI)-(Int_new[1][j][l][m]);
							if(Int_new[0][j][l][m]<0)
								 Int_new[0][j][l][m]=0.0;
				}
			
			for(j=1;j<=nj;j++)
				{
				 I_w=0.0;
				for(l=0;l<nt;l++)
					for(m=np/4;m<3*np/4;m++)
						 I_w=I_w-(1-eps_r)*(Int_new[ni][j][l][m])*(*(D_x+np*l+m))/PI;
						 
				for(l=0;l<nt;l++)
					for(m=3*np/4;m>np/4;m--)
						(Int_new[ni+1][j][l][m])=2.0*(I_w+eps_r*sigma*pow(Tr,4)/PI)-(Int_new[ni][j][l][m]);
							if(Int_new[ni+1][j][l][m]<0)
								 Int_new[ni+1][j][l][m]=0.0;
				}
			         Sum_n=0.0;
				 Sum_o=0.0;
			for(l=0;l<nt;l++)
				for(m=0;m<np;m++)
					for(j=0;j<=(nj+1);j++)
						for(i=0;i<=(ni+1);i++)
							{
							Sum_n=Sum_n+(Int_new[i][i][l][m]);
							Sum_o=Sum_o+(Int_old[i][i][l][m]);
							 }
							rms=fabs((Sum_n-Sum_o))/fabs(Sum_n);
							
			printf("%e\n",rms);
			
			}while(rms>1.0e-7);
			
		for(l=0;l<nt;l++)
				{
			  	for(m=0;m<np;m++)
					{
			  		for(j=1;j<=(nj);j++)
						{
						for(i=1;i<=(ni);i++)
							fprintf(fp1,"%e\t%e\t%e\t%e\t%e\n",dx*(2.0*i-1.0)/2.0,dy*(2.0*j-1.0)/2.0,dtheta*(l+0.5),dphi*(m+0.5),Int_new[i][j][l][m]);
						}
					}
				}
		
		for(i=1;i<(ni+1);i++)
			{
			 q_w=0.0;
			for(l=0;l<nt;l++)
				for(m=np/2;m<np;m++)
					q_w=q_w+Int_new[i][1][l][m]*(*(D_y+np*l+m));
				
			flux_bot[i]=eps_l*(sigma*pow(Tl,4)-q_w)/(sigma*pow(Tg,4));
			 }
		for(i=1;i<(ni+1);i++)
			{
			 q_u=q_l=0.0;
			for(l=0;l<nt;l++)
				{
				for(m=0;m<np;m++)
					{
					 q_u=q_u+Int_new[i][10][l][m]*(*(D_y+np*l+m));
					 q_l=q_l+Int_new[i][11][l][m]*(*(D_y+np*l+m));
					 }
				}
			flux_mid[i]=eps_l*(sigma*pow(Tl,4)-(q_u+q_l))/(sigma*pow(Tg,4));
			 }
		for(i=1;i<=ni;i++)
			{
			fprintf(fp2,"%e\t%e\n",dx*(2.0*i-1.0)/2.0,flux_mid[i]);
			fprintf(fp3,"%e\t%e\n",dx*(2.0*i-1.0)/2.0,flux_bot[i]);
			}
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	return 0;
}

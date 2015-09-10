#include "blocknufft3d.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <QTime>
#include "omp.h"

double eval_kernel(double x,double y,double z) {
    //return 1+0*x+0*y+0*z;
    //return exp(-x*x)*exp(-y*y)*exp(-z*z);
    return x*y*z;
}

double *s_exp_lookup1=0;
double *s_exp_lookup2=0;
double *s_exp_lookup3=0;

void exp_lookup_init(double tau1,double tau2,double tau3) {
    if (s_exp_lookup1) free(s_exp_lookup1);
    if (s_exp_lookup2) free(s_exp_lookup2);
    if (s_exp_lookup3) free(s_exp_lookup3);
    int NN=100;
    s_exp_lookup1=(double *)malloc(sizeof(double)*NN);
    s_exp_lookup2=(double *)malloc(sizeof(double)*NN);
    s_exp_lookup3=(double *)malloc(sizeof(double)*NN);
    for (int k0=0; k0<NN; k0++) {
        s_exp_lookup1[k0]=exp(-k0*k0/(4*tau1));
        s_exp_lookup2[k0]=exp(-k0*k0/(4*tau2));
        s_exp_lookup3[k0]=exp(-k0*k0/(4*tau3));
    }
}

bool blockspread3d(const BlockNufft3DOptions &opts,double *out,double *x,double *y,double *z,double *d) {
	printf("blockspread3d.\n");

    exp_lookup_init(opts.tau1,opts.tau2,opts.tau3);

	int M=opts.M;
	int K1=opts.K1,K2=opts.K2,K3=opts.K3;
	int N1=opts.N1,N2=opts.N2,N3=opts.N3; int N1N2N3=N1*N2*N3; int N1N2=N1*N2;
	int R1=opts.R1,R2=opts.R2,R3=opts.R3;
	if ((R1>K1)||(R2>K2)||(R3>K3)) {
		printf("R cannot be bigger than K\n");
		return false;
	}
	for (int m=0; m<M; m++) {
		if ((x[m]<0)||(x[m]>=opts.N1)) {
			printf("x is out of bounds: x[%d]=%g\n",m,x[m]);
			return false;
		}
		if ((y[m]<0)||(y[m]>=opts.N2)) {
			printf("y is out of bounds: y[%d]=%g\n",m,y[m]);
			return false;
		}
		if ((z[m]<0)||(z[m]>=opts.N3)) {
			printf("z is out of bounds: z[%d]=%g\n",m,z[m]);
			return false;
		}
	}

	omp_set_num_threads(opts.num_threads);

	for (int n=0; n<N1N2N3; n++) out[n]=0;

	printf("setting up input blocks... ");
	QTime timer0; timer0.start();
	int num_blocks_x=ceil(opts.N1*1.0/K1);
	int num_blocks_y=ceil(opts.N2*1.0/K2);
	int num_blocks_z=ceil(opts.N3*1.0/K3);
	int num_blocks=num_blocks_x*num_blocks_y*num_blocks_z;
	int *input_block_counts=(int *)malloc(sizeof(int)*num_blocks); for (int i=0; i<num_blocks; i++) input_block_counts[i]=0;
    #pragma omp parallel
	{
        int local_input_block_counts[num_blocks]; for (int i=0; i<num_blocks; i++) local_input_block_counts[i]=0;

        #pragma omp for //causes a crash
		for (int m=0; m<M; m++) {
            int b1=(int)(x[m]/K1);
            int c1=(int)(x[m]-b1*K1+0.5);
            int b2=(int)(y[m]/K2);
            int c2=(int)(y[m]-b2*K2+0.5);
            int b3=(int)(z[m]/K3);
            int c3=(int)(z[m]-b3*K3+0.5);
			bool left_side=false,right_side=false;
			bool top_side=false,bottom_side=false;
			bool back_side=false,front_side=false;
			if ((c1>=K1-R1/2)&&(b1+1<num_blocks_x)) right_side=true;
			if ((c1<=-R1/2+R1-1)&&(b1-1>=0)) left_side=true;
			if ((c2>=K2-R2/2)&&(b2+1<num_blocks_y)) top_side=true;
			if ((c2<=-R2/2+R2-1)&&(b2-1>=0)) bottom_side=true;
			if ((c3>=K3-R3/2)&&(b3+1<num_blocks_z)) front_side=true;
			if ((c3<=-R3/2+R3-1)&&(b3-1>=0)) back_side=true;

			int inds1[27]; int inds2[27]; int inds3[27];
			for (int ii=0; ii<27; ii++) {inds1[ii]=inds2[ii]=inds3[ii]=-99;}
			int ee=0;
			for (int i1=-1; i1<=1; i1++) {
				if ( ((i1==-1)&&(left_side)) || (i1==0) || ((i1==1)&&(right_side)) ) {
					for (int i2=-1; i2<=1; i2++) {
						if ( ((i2==-1)&&(bottom_side)) || (i2==0) || ((i2==1)&&(top_side)) ) {
							for (int i3=-1; i3<=1; i3++) {
								if ( ((i3==-1)&&(back_side)) || (i3==0) || ((i3==1)&&(front_side)) ) {
									inds1[ee]=i1; inds2[ee]=i2; inds3[ee]=i3; ee++;
								}
							}
						}
					}
				}
			}
			for (int ii=0; ii<ee; ii++) {
                local_input_block_counts[(b1+inds1[ii])+num_blocks_x*(b2+inds2[ii])+num_blocks_x*num_blocks_y*(b3+inds3[ii])]++;
			}
		}

        #pragma omp critical
        {
            for (int i=0; i<num_blocks; i++) input_block_counts[i]+=local_input_block_counts[i];
        }
	}
	int input_size=0;
	int *input_block_indices=(int *)malloc(sizeof(int)*num_blocks);
	for (int i=0; i<num_blocks; i++) {
		input_block_indices[i]=input_size;
		input_size+=input_block_counts[i];
	}
	printf("Elapsed: %d ms\n",timer0.elapsed());

	QTime timerA; timerA.start();
	printf("setting up output blocks... ");
	int *output_block_counts=(int *)malloc(sizeof(int)*num_blocks); for (int i=0; i<num_blocks; i++) output_block_counts[i]=0;
	for (int n3=0; n3<N3; n3++) {
		for (int n2=0; n2<N2; n2++) {
			for (int n1=0; n1<N1; n1++) {
				int b1=n1/K1;
				int b2=n2/K2;
				int b3=n3/K3;
				output_block_counts[b1+num_blocks_x*b2+num_blocks_x*num_blocks_y*b3]++;
			}
		}
	}
	int output_size=0;
	int *output_block_indices=(int *)malloc(sizeof(int)*num_blocks);
	for (int i=0; i<num_blocks; i++) {
		output_block_indices[i]=output_size;
		output_size+=output_block_counts[i];
	}
	double *output_tmp=(double *)malloc(sizeof(double)*output_size);
	for (int ii=0; ii<output_size; ii++) output_tmp[ii]=0;
	printf("Elapsed: %d ms\n",timerA.elapsed());

	printf("setting input... ");
	QTime timerB; timerB.start();
	double *input_x=(double *)malloc(sizeof(double)*input_size);
	double *input_y=(double *)malloc(sizeof(double)*input_size);
	double *input_z=(double *)malloc(sizeof(double)*input_size);
	double *input_d=(double *)malloc(sizeof(double)*input_size);
	int *input_ii=(int *)malloc(sizeof(int)*num_blocks);
	for (int ii=0; ii<num_blocks; ii++) input_ii[ii]=0;
    omp_lock_t block_locks[num_blocks];
    for (int ii=0; ii<num_blocks; ii++) omp_init_lock(&block_locks[ii]);
    //#pragma omp parallel -- in the future we should parallelize this step. For now there is a problem with needing to lock/unlock around input_ii[bb]++
	{
        //#pragma omp for
		for (int m=0; m<M; m++) {
			int b1=(int)(x[m]/K1);
            int c1=(int)(x[m]-b1*K1+0.5);
			int b2=(int)(y[m]/K2);
            int c2=(int)(y[m]-b2*K2+0.5);
			int b3=(int)(z[m]/K3);
            int c3=(int)(z[m]-b3*K3+0.5);
			bool left_side=false,right_side=false;
			bool top_side=false,bottom_side=false;
            bool back_side=false,front_side=false;
			if ((c1>=K1-R1/2)&&(b1+1<num_blocks_x)) right_side=true;
			if ((c1<=-R1/2+R1-1)&&(b1-1>=0)) left_side=true;
			if ((c2>=K2-R2/2)&&(b2+1<num_blocks_y)) top_side=true;
			if ((c2<=-R2/2+R2-1)&&(b2-1>=0)) bottom_side=true;
			if ((c3>=K3-R3/2)&&(b3+1<num_blocks_z)) front_side=true;
			if ((c3<=-R3/2+R3-1)&&(b3-1>=0)) back_side=true;

			int inds1[27]; int inds2[27]; int inds3[27];
			for (int ii=0; ii<27; ii++) {inds1[ii]=inds2[ii]=inds3[ii]=-99;}
			int ee=0;
			for (int i1=-1; i1<=1; i1++) {
				if ( ((i1==-1)&&(left_side)) || (i1==0) || ((i1==1)&&(right_side)) ) {
					for (int i2=-1; i2<=1; i2++) {
						if ( ((i2==-1)&&(bottom_side)) || (i2==0) || ((i2==1)&&(top_side)) ) {
							for (int i3=-1; i3<=1; i3++) {
								if ( ((i3==-1)&&(back_side)) || (i3==0) || ((i3==1)&&(front_side)) ) {
									inds1[ee]=i1; inds2[ee]=i2; inds3[ee]=i3; ee++;
								}
							}
						}
					}
				}
			}
            for (int ii=0; ii<ee; ii++) {
                int bb=(b1+inds1[ii])+num_blocks_x*(b2+inds2[ii])+num_blocks_x*num_blocks_y*(b3+inds3[ii]);
                //omp_set_lock(&block_locks[bb]); //very important to lock here because we are incrementing input_ii[bb] as part of this whole operation!
                int jj=input_block_indices[bb]+input_ii[bb];
                input_ii[bb]++;
                //omp_unset_lock(&block_locks[bb]);
                input_x[jj]=x[m];
                input_y[jj]=y[m];
                input_z[jj]=z[m];
                input_d[jj]=d[m];
            }
		}
	}
	printf("Elapsed: %d ms\n",timerA.elapsed());

	printf("spreading... ");
	QTime timerC; timerC.start();
	#pragma omp parallel
	{
		#pragma omp for
		for (int cc=0; cc<num_blocks; cc++) {
			int cc1=cc%num_blocks_x;
			int cc2=(cc%(num_blocks_x*num_blocks_y))/num_blocks_x;
			int cc3=cc/(num_blocks_x*num_blocks_y);
			int factor1=K1; if ((cc1+1)*K1>=N1) factor1=N1-cc1*K1;
			int factor2=K2; if ((cc2+1)*K2>=N2) factor2=N2-cc2*K2;
			//int factor3=K3; if ((cc3+1)*K3>=N3) factor3=N3-cc3*K3;
			int factor12=factor1*factor2;
            int block_xmin=cc1*K1,block_xmax=(cc1+1)*K1-1; if (block_xmax>=N1) block_xmax=N1-1;
            int block_ymin=cc2*K2,block_ymax=(cc2+1)*K2-1; if (block_ymax>=N2) block_ymax=N2-1;
            int block_zmin=cc3*K3,block_zmax=(cc3+1)*K3-1; if (block_zmax>=N3) block_zmax=N3-1;
			int jj=input_block_indices[cc];
			int tmp=jj+input_block_counts[cc];
			while (jj<tmp) {
				double x0=input_x[jj];
				double y0=input_y[jj];
				double z0=input_z[jj];
				double d0=input_d[jj];

                int x_integer=(int)(x0+0.5);
                double x_diff=x0-x_integer;
                double x_term1=exp(-x_diff*x_diff/(4*opts.tau1));
                double x_term2_factor=exp(-2*x_diff/(4*opts.tau1));
                int xmin=x_integer-R1/2;
                int xmax=x_integer-R1/2+R1-1;

                int y_integer=(int)(y0+0.5);
                double y_diff=y0-y_integer;
                double y_term1=exp(-y_diff*y_diff/(4*opts.tau2));
                double y_term2_factor=exp(-2*y_diff/(4*opts.tau2));
                int ymin=y_integer-R1/2;
                int ymax=y_integer-R1/2+R1-1;

                int z_integer=(int)(z0+0.5);
                double z_diff=z0-z_integer;
                double z_term1=exp(-z_diff*z_diff/(4*opts.tau3));
                double z_term2_factor=exp(-2*z_diff/(4*opts.tau3));
                int zmin=z_integer-R1/2;
                int zmax=z_integer-R1/2+R1-1;

                double x_term2[(R1+1)/2+1]; x_term2[0]=1;
                for (int aa=1; aa<(R1+1)/2+1; aa++) x_term2[aa]=x_term2[aa-1]*x_term2_factor;
                double y_term2[(R2+1)/2+1]; y_term2[0]=1;
                for (int aa=1; aa<(R2+1)/2+1; aa++) y_term2[aa]=y_term2[aa-1]*y_term2_factor;
                double z_term2[(R3+1)/2+1]; z_term2[0]=1;
                for (int aa=1; aa<(R3+1)/2+1; aa++) z_term2[aa]=z_term2[aa-1]*z_term2_factor;

                for (int aa=0; aa<(R1+1)/2+1; aa++) x_term2[aa]*=s_exp_lookup1[aa];
                for (int aa=0; aa<(R2+1)/2+1; aa++) y_term2[aa]*=s_exp_lookup2[aa];
                for (int aa=0; aa<(R3+1)/2+1; aa++) z_term2[aa]*=s_exp_lookup3[aa];

                double kernval0=d0*x_term1*y_term1*z_term1;
                for (int iz=zmin; iz<=zmax; iz++) {
                    if ((block_zmin<=iz)&&(iz<=block_zmax)) {
                        int kkk1=output_block_indices[cc]+factor12*(iz-block_zmin);
                        int iiz=abs(iz-z_integer);
                        double kernval1=kernval0*z_term2[iiz];
                        for (int iy=ymin; iy<=ymax; iy++) {
                            if ((block_ymin<=iy)&&(iy<=block_ymax)) {
                                int kkk2=kkk1+factor1*(iy-block_ymin);
                                int iiy=abs(iy-y_integer);
                                double kernval2=kernval1*y_term2[iiy];
                                int kkk3=kkk2+(xmin-block_xmin);
                                for (int ix=xmin; ix<=xmax; ix++) {
                                    if ((block_xmin<=ix)&&(ix<=block_xmax)) {
                                        //The following lines are time-consuming!!
                                        int iix=abs(ix-x_integer);
                                        output_tmp[kkk3]+=kernval2*x_term2[iix];
                                    }
                                    kkk3++;
                                }
                            }
                        }
                    }
                }
                //printf("%d,%d,%d,%d  %d,%d,%d,%d  %d,%d,%d,%d  %d\n",xmin,xmax,block_xmin,block_xmax,ymin,ymax,block_ymin,block_ymax,zmin,zmax,block_zmin,block_zmax,debug_ct);
                jj++;
            }
		}
	}
	printf("Elapsed: %d ms\n",timerC.elapsed());

	printf("setting output... ");
	QTime timerD; timerD.start();
    #pragma omp parallel
	{
        #pragma omp for
		for (int cc=0; cc<num_blocks; cc++) {
			int cc1=cc%num_blocks_x;
			int cc2=(cc%(num_blocks_x*num_blocks_y))/num_blocks_x;
			int cc3=cc/(num_blocks_x*num_blocks_y);
			int factor1=K1; if ((cc1+1)*K1>=N1) factor1=N1-cc1*K1;
			int factor2=K2; if ((cc2+1)*K2>=N2) factor2=N2-cc2*K2;
			int factor3=K3; if ((cc3+1)*K3>=N3) factor3=N3-cc3*K3;
			int dd1=cc1*K1;
			int dd2=cc2*K2;
			int dd3=cc3*K3;
			int kkk=output_block_indices[cc];
			for (int i3=0; i3<factor3; i3++) {
				for (int i2=0; i2<factor2; i2++) {
					for (int i1=0; i1<factor1; i1++) {
						out[(dd1+i1)+(dd2+i2)*N1+(dd3+i3)*N1N2]=output_tmp[kkk];
						kkk++;
					}
				}
			}
		}
	}
	printf("Elapsed: %d ms\n",timerD.elapsed());

	printf("freeing... ");
    for (int ii=0; ii<num_blocks; ii++) omp_destroy_lock(&block_locks[ii]);
	QTime timerE; timerE.start();
	free(input_block_counts);
	free(input_block_indices);
	free(output_block_counts);
	free(output_block_indices);
	free(output_tmp);
	free(input_x);
	free(input_y);
	free(input_z);
	free(input_d);
	free(input_ii);
	printf("Elapsed: %d ms\n",timerE.elapsed());

	return true;
}

void test_blockspread3d(BlockNufft3DOptions &opts)
{
	printf("test_blockspread2d: preparing.\n");
	double *x=(double *)malloc(sizeof(double)*opts.M);
	double *y=(double *)malloc(sizeof(double)*opts.M);
	double *z=(double *)malloc(sizeof(double)*opts.M);
	double *d=(double *)malloc(sizeof(double)*opts.M);
	double *out1=(double *)malloc(sizeof(double)*opts.N1*opts.N2*opts.N3);
	double *out2=(double *)malloc(sizeof(double)*opts.N1*opts.N2*opts.N3);

	printf("setting input data.\n");
	for (int m=0; m<opts.M; m++) {
        x[m]=1+(qrand()%(opts.N1-1))+sin(m);
        y[m]=1+(qrand()%(opts.N2-1))+cos(m);
        z[m]=1+(qrand()%(opts.N3-1))+sin(2*m);
		d[m]=1;
	}

    if (1) {
		QTime timer; timer.start();
		blockspread3d(opts,out1,x,y,z,d);
		for (int aa=0; aa<10; aa++) printf("%d:%g, ",aa,out1[aa]); printf("\n");
		printf("Time for blockspread3d: %d ms\n\n",timer.elapsed());
	}
	if (1) {
		QTime timer; timer.start();
        opts.K1=opts.N1; opts.K2=opts.N2; opts.K3=opts.N3;
		blockspread3d(opts,out2,x,y,z,d);
		for (int aa=0; aa<10; aa++) printf("%d:%g, ",aa,out2[aa]); printf("\n");
		printf("Time for blockspread3d: %d ms\n\n",timer.elapsed());

		bool ok=true;
        double maxdiff=0;
		for (int i=0; i<opts.N1*opts.N2*opts.N3; i++) {
			if (out1[i]!=out2[i]) ok=false;
            maxdiff=fmax(fabs(out1[i]-out2[i]),maxdiff);
		}
		if (ok) printf("Outputs match!\n");
        else printf("Outputs do not match. maxdiff=%g\n",maxdiff);
	}

	free(x);
	free(y);
	free(z);
	free(d);
}

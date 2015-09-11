#include "blocknufft3d.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <QTime>
#include "omp.h"
#include "fftw3.h"

#define LEFT_SIDE 1
#define RIGHT_SIDE 2
#define BOTTOM_SIDE 4
#define TOP_SIDE 8
#define BACK_SIDE 16
#define FRONT_SIDE 32

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
void exp_lookup_destroy() {
	free(s_exp_lookup1); s_exp_lookup1=0;
	free(s_exp_lookup2); s_exp_lookup2=0;
	free(s_exp_lookup3); s_exp_lookup3=0;
}

bool blockspread3d(const BlockSpread3DOptions &opts,double *out,double *x,double *y,double *z,double *d) {
	printf("blockspread3d.\n");

	QTime timer0;

    exp_lookup_init(opts.tau1,opts.tau2,opts.tau3);
	omp_set_num_threads(opts.num_threads);

	int M=opts.M;
	int K1=opts.K1,K2=opts.K2,K3=opts.K3;
	int N1=opts.N1o,N2=opts.N2o,N3=opts.N3o; int N1N2=N1*N2;
	int N1_times_2=N1*2; int N1N2_times_2=N1N2*2;
	int R1=opts.R1,R2=opts.R2,R3=opts.R3;
	if ((R1>K1)||(R2>K2)||(R3>K3)) {
		printf("R cannot be bigger than K\n");
		return false;
	}
	printf("Checking inputs... "); timer0.start();
	bool ok=true;
	#pragma omp parallel
	{
		bool local_ok=true;
		#pragma omp for
		for (int m=0; m<M; m++) {
			if ((x[m]<0)||(x[m]>=opts.N1o)) {
				local_ok=false;
			}
			if ((y[m]<0)||(y[m]>=opts.N2o)) {
				local_ok=false;
			}
			if ((z[m]<0)||(z[m]>=opts.N3o)) {
				local_ok=false;
			}
		}
		#pragma omp critical
		{
			if (!local_ok) ok=false;
		}
	}
	printf("Elapsed: %d ms\n",timer0.elapsed());

	if (!ok) {
		printf("Something is out of range.\n");
		for (int i=0; i<N1*N2*N3*2; i++) out[i]=0;
		return false;
	}

	int num_blocks_x=ceil(opts.N1o*1.0/K1);
	int num_blocks_y=ceil(opts.N2o*1.0/K2);
	int num_blocks_z=ceil(opts.N3o*1.0/K3);
	int num_blocks=num_blocks_x*num_blocks_y*num_blocks_z;

	printf("Defining block ids and location codes... "); timer0.start();
	int *location_codes=(int *)malloc(sizeof(int)*M);
	int *block_ids=(int *)malloc(sizeof(int)*M);
	#pragma omp parallel
	{
		#pragma omp for
		for (int m=0; m<M; m++) {
			int b1=(int)(x[m]/K1);
			int c1=(int)(x[m]-b1*K1+0.5);
			int b2=(int)(y[m]/K2);
			int c2=(int)(y[m]-b2*K2+0.5);
			int b3=(int)(z[m]/K3);
			int c3=(int)(z[m]-b3*K3+0.5);
			int code=0;
			if ((c1>=K1-R1/2)&&(b1+1<num_blocks_x)) code=code|RIGHT_SIDE;
			if ((c1<=-R1/2+R1-1)&&(b1-1>=0)) code=code|LEFT_SIDE;
			if ((c2>=K2-R2/2)&&(b2+1<num_blocks_y)) code=code|TOP_SIDE;
			if ((c2<=-R2/2+R2-1)&&(b2-1>=0)) code=code|BOTTOM_SIDE;
			if ((c3>=K3-R3/2)&&(b3+1<num_blocks_z)) code=code|FRONT_SIDE;
			if ((c3<=-R3/2+R3-1)&&(b3-1>=0)) code=code|BACK_SIDE;
			location_codes[m]=code;
			block_ids[m]=b1+num_blocks_x*b2+num_blocks_x*num_blocks_y*b3;
		}
	}
	printf("Elapsed: %d ms\n",timer0.elapsed());

	printf("Computing block counts... "); timer0.start();
	int input_block_counts[num_blocks]; for (int i=0; i<num_blocks; i++) input_block_counts[i]=0;
	#pragma omp parallel
	{
		int local_input_block_counts[num_blocks]; for (int i=0; i<num_blocks; i++) local_input_block_counts[i]=0;
		#pragma omp for
		for (int m=0; m<M; m++) {
			int code=location_codes[m];
			int bb0=block_ids[m];
			for (int i1=-1; i1<=1; i1++) {
				if ( ((i1==-1)&&(code&LEFT_SIDE)) || (i1==0) || ((i1==1)&&(code&RIGHT_SIDE)) ) {
					int bb1=bb0+i1;
					for (int i2=-1; i2<=1; i2++) {
						if ( ((i2==-1)&&(code&BOTTOM_SIDE)) || (i2==0) || ((i2==1)&&(code&TOP_SIDE)) ) {
							int bb2=bb1+i2*num_blocks_x;
							for (int i3=-1; i3<=1; i3++) {
								if ( ((i3==-1)&&(code&BACK_SIDE)) || (i3==0) || ((i3==1)&&(code&FRONT_SIDE)) ) {
									int bb3=bb2+i3*num_blocks_x*num_blocks_y;
									local_input_block_counts[bb3]++;
								}
							}
						}
					}
				}
			}
		}
		#pragma omp critical
		{
			for (int i=0; i<num_blocks; i++) input_block_counts[i]+=local_input_block_counts[i];
		}
	}
	printf("Elapsed: %d ms\n",timer0.elapsed());

	printf("Computing sizes... "); timer0.start();
	int input_size=0;
	int input_block_indices[num_blocks];
	int output_size=0;
	int output_block_indices[num_blocks];
	for (int i=0; i<num_blocks; i++) {
		int b1=i%num_blocks_x;
		int b2=(i/num_blocks_x)%num_blocks_y;
		int b3=(i/num_blocks_x/num_blocks_y);
		input_block_indices[i]=input_size;
		input_size+=input_block_counts[i];
		output_block_indices[i]=output_size;
		int F1=fmin(K1,N1-b1*K1);
		int F2=fmin(K2,N2-b2*K2);
		int F3=fmin(K3,N3-b3*K3);
		output_size+=F1*F2*F3;
	}
	int input_size_times_2=input_size*2;
	int output_size_times_2=output_size*2; //for complex
	double *output_tmp=(double *)malloc(sizeof(double)*output_size_times_2); //complex
	for (int ii=0; ii<output_size_times_2; ii++) output_tmp[ii]=0;
	printf("Elapsed: %d ms\n",timer0.elapsed());

	printf("setting input... "); timer0.start();
	int block_ii[num_blocks]; for (int i=0; i<num_blocks; i++) block_ii[i]=0;
	double *input_x=(double *)malloc(sizeof(double)*input_size);
	double *input_y=(double *)malloc(sizeof(double)*input_size);
	double *input_z=(double *)malloc(sizeof(double)*input_size);
	double *input_d=(double *)malloc(sizeof(double)*input_size_times_2); //times 2 because complex
	for (int m=0; m<M; m++) { //can this be parallelized? Not sure!
		int code=location_codes[m];
		int bb0=block_ids[m];
		for (int i1=-1; i1<=1; i1++) {
			if ( ((i1==-1)&&(code&LEFT_SIDE)) || (i1==0) || ((i1==1)&&(code&RIGHT_SIDE)) ) {
				int bb1=bb0+i1;
				for (int i2=-1; i2<=1; i2++) {
					if ( ((i2==-1)&&(code&BOTTOM_SIDE)) || (i2==0) || ((i2==1)&&(code&TOP_SIDE)) ) {
						int bb2=bb1+i2*num_blocks_x;
						for (int i3=-1; i3<=1; i3++) {
							if ( ((i3==-1)&&(code&BACK_SIDE)) || (i3==0) || ((i3==1)&&(code&FRONT_SIDE)) ) {
								int bb3=bb2+i3*num_blocks_x*num_blocks_y;
								int iii=input_block_indices[bb3]+block_ii[bb3];
								input_x[iii]=x[m];
								input_y[iii]=y[m];
								input_z[iii]=z[m];
								input_d[iii*2]=d[m*2];
								input_d[iii*2+1]=d[m*2+1];
								block_ii[bb3]++;
							}
						}
					}
				}
			}
		}
	}
	printf("Elapsed: %d ms\n",timer0.elapsed());

	printf("spreading... "); timer0.start();
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
			int factor1_times_2=factor1*2; int factor12_times_2=factor12*2;
            int block_xmin=cc1*K1,block_xmax=(cc1+1)*K1-1; if (block_xmax>=N1) block_xmax=N1-1;
            int block_ymin=cc2*K2,block_ymax=(cc2+1)*K2-1; if (block_ymax>=N2) block_ymax=N2-1;
            int block_zmin=cc3*K3,block_zmax=(cc3+1)*K3-1; if (block_zmax>=N3) block_zmax=N3-1;
			int jj=input_block_indices[cc];
			int tmp=jj+input_block_counts[cc];
			double x_term2[(R1+1)/2+1];
			double y_term2[(R2+1)/2+1];
			double z_term2[(R3+1)/2+1];
			double precomp_x_term2[R1+1]; //conservatively large
			while (jj<tmp) {
				double x0=input_x[jj],y0=input_y[jj],z0=input_z[jj],d0_re=input_d[jj*2],d0_im=input_d[jj*2+1];

                int x_integer=(int)(x0+0.5);
                double x_diff=x0-x_integer;
                double x_term1=exp(-x_diff*x_diff/(4*opts.tau1));
                double x_term2_factor=exp(-2*x_diff/(4*opts.tau1));
				int xmin=fmax(x_integer-R1/2,block_xmin);
				int xmax=fmin(x_integer-R1/2+R1-1,block_xmax);

                int y_integer=(int)(y0+0.5);
                double y_diff=y0-y_integer;
                double y_term1=exp(-y_diff*y_diff/(4*opts.tau2));
                double y_term2_factor=exp(-2*y_diff/(4*opts.tau2));
				int ymin=fmax(y_integer-R1/2,block_ymin);
				int ymax=fmin(y_integer-R1/2+R1-1,block_ymax);

                int z_integer=(int)(z0+0.5);
                double z_diff=z0-z_integer;
                double z_term1=exp(-z_diff*z_diff/(4*opts.tau3));
                double z_term2_factor=exp(-2*z_diff/(4*opts.tau3));
				int zmin=fmax(z_integer-R1/2,block_zmin);
				int zmax=fmin(z_integer-R1/2+R1-1,block_zmax);

				x_term2[0]=1;
				int aamax=(R1+1)/2;
				for (int aa=1; aa<=aamax; aa++) x_term2[aa]=x_term2[aa-1]*x_term2_factor;
				y_term2[0]=1;
				int bbmax=(R2+1)/2;
				for (int aa=1; aa<=bbmax; aa++) y_term2[aa]=y_term2[aa-1]*y_term2_factor;
				z_term2[0]=1;
				int ccmax=(R3+1)/2;
				for (int aa=1; aa<=ccmax; aa++) z_term2[aa]=z_term2[aa-1]*z_term2_factor;

				for (int aa=0; aa<=aamax; aa++) x_term2[aa]*=s_exp_lookup1[aa];
				for (int aa=0; aa<=bbmax; aa++) y_term2[aa]*=s_exp_lookup2[aa];
				for (int aa=0; aa<=ccmax; aa++) z_term2[aa]*=s_exp_lookup3[aa];

				int precomp_x_term2_sz=xmax-xmin+1;
				for (int ix=xmin; ix<=xmax; ix++) {
					int iix=abs(ix-x_integer);
					precomp_x_term2[ix-xmin]=x_term2[iix];
				}

				double kernval0=x_term1*y_term1*z_term1;
                for (int iz=zmin; iz<=zmax; iz++) {
					int kkk1=output_block_indices[cc]*2+factor12_times_2*(iz-block_zmin); //complex index
					int iiz=abs(iz-z_integer);
					double kernval1=kernval0*z_term2[iiz];
					for (int iy=ymin; iy<=ymax; iy++) {
						int kkk2=kkk1+factor1_times_2*(iy-block_ymin);
						int iiy=abs(iy-y_integer);
						double kernval2=kernval1*y_term2[iiy];
						int kkk3=kkk2+(xmin-block_xmin)*2; //times 2 because complex
						/* //since this inner-most loop matters a lot, we optimize the following by precomputing...
						for (int ix=xmin; ix<=xmax; ix++) {
							int iix=abs(ix-x_integer);
							output_tmp[kkk3]+=d0_re*kernval2*x_term2[iix];
							output_tmp[kkk3+1]+=d0_im*kernval2*x_term2[iix];
							kkk3+=2; //plus two because complex
						}
						*/
						//replaced by the following for efficiency -- note, we don't need to compute abs(ix-x_integer)
						for (int iii=0; iii<precomp_x_term2_sz; iii++) {
							//most of the time is spent within this code block!!!
							double tmp0=kernval2*precomp_x_term2[iii];
							output_tmp[kkk3]+=d0_re*tmp0;
							output_tmp[kkk3+1]+=d0_im*tmp0; //most of the time is spent on this line!!!
							kkk3+=2; //plus two because complex
						}
					}
                }

                //printf("%d,%d,%d,%d  %d,%d,%d,%d  %d,%d,%d,%d  %d\n",xmin,xmax,block_xmin,block_xmax,ymin,ymax,block_ymin,block_ymax,zmin,zmax,block_zmin,block_zmax,debug_ct);
                jj++;
            }
		}
	}
	printf("Elapsed: %d ms\n",timer0.elapsed());

	printf("setting output... "); timer0.start();
	//for (int n=0; n<N1N2N3_times_2; n++) out[n]=0; //not necessary?
    #pragma omp parallel
	{
        #pragma omp for
		for (int cc=0; cc<num_blocks; cc++) {
			int cc1=cc%num_blocks_x;
			int cc2=(cc%(num_blocks_x*num_blocks_y))/num_blocks_x;
			int cc3=cc/(num_blocks_x*num_blocks_y);
			int factor1=K1; if ((cc1+1)*K1>=N1) factor1=(N1-cc1*K1);
			int factor2=K2; if ((cc2+1)*K2>=N2) factor2=N2-cc2*K2;
			int factor3=K3; if ((cc3+1)*K3>=N3) factor3=N3-cc3*K3;
			int dd1=cc1*K1;
			int dd2=cc2*K2;
			int dd3=cc3*K3;
			int kkk=output_block_indices[cc]*2; //times 2 because complex
			for (int i3=0; i3<factor3; i3++) {
				for (int i2=0; i2<factor2; i2++) {
					for (int i1=0; i1<factor1; i1++) { //make this inner loop more efficient by not doing the multiplication here?
						int jjj=(dd1+i1)*2+(dd2+i2)*N1_times_2+(dd3+i3)*N1N2_times_2; //complex index
						out[jjj]=output_tmp[kkk];
						out[jjj+1]=output_tmp[kkk+1];
						kkk+=2; //add 2 because complex
					}
				}
			}
		}
	}
	printf("Elapsed: %d ms\n",timer0.elapsed());

	printf("freeing... "); timer0.start();
	exp_lookup_destroy();
	free(output_tmp);
	free(input_x);
	free(input_y);
	free(input_z);
	free(input_d);
	free(location_codes);
	free(block_ids);
	printf("Elapsed: %d ms\n",timer0.elapsed());

	return true;
}

void test_blockspread3d(BlockSpread3DOptions &opts)
{
	printf("test_blockspread3d: preparing.\n");
	double *x=(double *)malloc(sizeof(double)*opts.M);
	double *y=(double *)malloc(sizeof(double)*opts.M);
	double *z=(double *)malloc(sizeof(double)*opts.M);
	double *d=(double *)malloc(sizeof(double)*opts.M*2);
	double *out1=(double *)malloc(sizeof(double)*opts.N1o*opts.N2o*opts.N3o*2);
	double *out2=(double *)malloc(sizeof(double)*opts.N1o*opts.N2o*opts.N3o*2);

	printf("setting input data.\n");
	for (int m=0; m<opts.M; m++) {
		x[m]=1+(qrand()%(opts.N1o-1))+sin(m);
		y[m]=1+(qrand()%(opts.N2o-1))+cos(m);
		z[m]=1+(qrand()%(opts.N3o-1))+sin(2*m);
		d[m*2]=1;
		d[m*2+1]=0;
	}

    if (1) {
		QTime timer; timer.start();
		blockspread3d(opts,out1,x,y,z,d);
		for (int aa=0; aa<10; aa++) printf("%d:%g, ",aa,out1[aa]); printf("\n");
		printf("Time for blockspread3d: %d ms\n\n",timer.elapsed());
	}
	if (1) {
		QTime timer; timer.start();
		opts.K1=opts.N1o; opts.K2=opts.N2o; opts.K3=opts.N3o;
		blockspread3d(opts,out2,x,y,z,d);
		for (int aa=0; aa<10; aa++) printf("%d:%g, ",aa,out2[aa]); printf("\n");
		printf("Time for blockspread3d: %d ms\n\n",timer.elapsed());

		bool ok=true;
        double maxdiff=0;
		for (int i=0; i<opts.N1o*opts.N2o*opts.N3o; i++) {
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

bool do_fft_3d(int N1,int N2,int N3,double *out,double *in) {
	int N1N2N3=N1*N2*N3;
	fftw_complex *in2=(fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N1N2N3);
	fftw_complex *out2=(fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N1N2N3);
	printf("Setting in...\n");
	for (int ii=0; ii<N1N2N3; ii++) {
		in2[ii][0]=in[ii*2];
		in2[ii][1]=in[ii*2+1];
	}
	printf("fftw_plan_dft_3d...\n");
	fftw_plan p=fftw_plan_dft_3d(N1,N2,N3,in2,out2,FFTW_FORWARD,FFTW_ESTIMATE);
	printf("fftw_execute...\n");
	fftw_execute(p);
	printf("Setting out...\n");
	for (int ii=0; ii<N1N2N3; ii++) {
		out[ii*2]=out2[ii][0];
		out[ii*2+1]=out2[ii][1];
	}
	printf("fftw_free...\n");
	fftw_free(in2);
	fftw_free(out2);

	return true;
}
void do_fix_3d(const BlockNufft3DOptions &opts,double *out,double *out_oversamp_hat) {
	int N1N2N3_times_2=opts.N1*opts.N2*opts.N3*2;
	for (int ii=0; ii<N1N2N3_times_2; ii++) {
		out[ii]=out_oversamp_hat[ii]; //fix, haha
	}
}

//x,y,z should be in range [0,2pi)
bool blocknufft3d(const BlockNufft3DOptions &opts,double *out,double *x,double *y,double *z,double *d) {
	BlockSpread3DOptions sopts;
	sopts.N1o=opts.N1*opts.oversamp; sopts.N2o=opts.N2*opts.oversamp; sopts.N3o=opts.N3*opts.oversamp;
	sopts.K1=opts.K1; sopts.K2=opts.K2; sopts.K3=opts.K3;
	sopts.M=opts.M;
	sopts.num_threads=opts.num_threads;
	sopts.R1=opts.R1; sopts.R2=opts.R2; sopts.R3=opts.R3;
	sopts.tau1=opts.tau1; sopts.tau2=opts.tau2; sopts.tau3=opts.tau3;

	printf("Allocating...\n");
	double *out_oversamp=(double *)malloc(sizeof(double)*sopts.N1o*sopts.N2o*sopts.N3o*2);
	double *out_oversamp_hat=(double *)malloc(sizeof(double)*sopts.N1o*sopts.N2o*sopts.N3o*2);
	double *x2=(double *)malloc(sizeof(double)*opts.M);
	double *y2=(double *)malloc(sizeof(double)*opts.M);
	double *z2=(double *)malloc(sizeof(double)*opts.M);

	printf("Scaling x,y,z...\n");
	double factor_x=sopts.N1o/(2*M_PI);
	double factor_y=sopts.N2o/(2*M_PI);
	double factor_z=sopts.N3o/(2*M_PI);
	for (int ii=0; ii<opts.M; ii++) {
		x2[ii]=x[ii]*factor_x;
		y2[ii]=y[ii]*factor_y;
		z2[ii]=z[ii]*factor_z;
	}
	printf("Spreading...\n");
	blockspread3d(sopts,out_oversamp,x2,y2,z2,d);

	printf("fft...\n");
	if (!do_fft_3d(sopts.N1o,sopts.N2o,sopts.N3o,out_oversamp_hat,out_oversamp)) {
		printf("problem in do_fft_3d\n");
		free(out_oversamp);
		free(out_oversamp_hat);
		free(x2); free(y2); free(z2);
		return false;
	}

	printf("fix...\n");
	do_fix_3d(opts,out,out_oversamp_hat);

	printf("free...\n");
	free(out_oversamp_hat);
	free(out_oversamp);
	free(x2); free(y2); free(z2);
	printf("done.\n");
	return true;
}

void test_blocknufft3d(BlockNufft3DOptions &opts)
{
	printf("test_blocknufft3d: preparing.\n");
	double *x=(double *)malloc(sizeof(double)*opts.M);
	double *y=(double *)malloc(sizeof(double)*opts.M);
	double *z=(double *)malloc(sizeof(double)*opts.M);
	double *d=(double *)malloc(sizeof(double)*opts.M*2);
	double *out1=(double *)malloc(sizeof(double)*opts.N1*opts.N2*opts.N3*2);
	double *out2=(double *)malloc(sizeof(double)*opts.N1*opts.N2*opts.N3*2);

	printf("setting input data.\n");
	for (int m=0; m<opts.M; m++) {
		x[m]=(1+sin(m+1))*M_PI; //+1 is important so we don't get a 2pi
		y[m]=(1+sin(m+1))*M_PI;
		z[m]=(1+sin(m*0.67+1))*M_PI;
		d[m*2]=1;
		d[m*2+1]=0;
	}

	if (1) {
		QTime timer; timer.start();
		printf("blocknufft3d\n");
		blocknufft3d(opts,out1,x,y,z,d);
		for (int aa=0; aa<10; aa++) printf("%d:%g, ",aa,out1[aa]); printf("\n");
		printf("Time for blocknufft3d: %d ms\n\n",timer.elapsed());
	}
	if (1) {
		QTime timer; timer.start();
		opts.K1=opts.N1*opts.oversamp; opts.K2=opts.N2*opts.oversamp; opts.K3=opts.K3*opts.oversamp;
		printf("blocknufft3d\n");
		blocknufft3d(opts,out2,x,y,z,d);
		for (int aa=0; aa<10; aa++) printf("%d:%g, ",aa,out2[aa]); printf("\n");
		printf("Time for blockspread3d: %d ms\n\n",timer.elapsed());

		bool ok=true;
		double maxdiff=0;
		for (int i=0; i<opts.N1*opts.N2*opts.N3*2; i++) {
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

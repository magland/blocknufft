#include "blocknufft3d.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "qute.h"
#include "omp.h"
#include "fftw3.h"

#define LEFT_SIDE 1
#define RIGHT_SIDE 2
#define BOTTOM_SIDE 4
#define TOP_SIDE 8
#define BACK_SIDE 16
#define FRONT_SIDE 32

#define ROUND_2_INT(f) ((int)(f >= 0.0 ? (f + 0.5) : (f - 0.5)))

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
            int c1=ROUND_2_INT(x[m]-b1*K1);
			int b2=(int)(y[m]/K2);
            int c2=ROUND_2_INT(y[m]-b2*K2);
			int b3=(int)(z[m]/K3);
            int c3=ROUND_2_INT(z[m]-b3*K3);
			int code=0;
            int k1=K1,k2=K2,k3=K3;
            if (b1==num_blocks_x-1) k1=fmin(N1-(num_blocks_x-1)*K1,K1); //the min is needed in case K1>N1*2
            if (b2==num_blocks_y-1) k2=fmin(N2-(num_blocks_y-1)*K2,K2);
            if (b3==num_blocks_z-1) k3=fmin(N3-(num_blocks_z-1)*K3,K3);
            if ((k1-c1<=R1/2)) code=code|RIGHT_SIDE;
            if ((c1<=-R1/2+R1-1)) code=code|LEFT_SIDE;
            if ((k2-c2<=R2/2)) code=code|TOP_SIDE;
            if ((c2<=-R2/2+R2-1)) code=code|BOTTOM_SIDE;
            if ((k3-c3<=R3/2)) code=code|FRONT_SIDE;
            if ((c3<=-R3/2+R3-1)) code=code|BACK_SIDE;
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
            int bx=bb0%num_blocks_x;
            int by=(bb0/num_blocks_x)%num_blocks_y;
            int bz=(bb0/num_blocks_x/num_blocks_y);
            for (int i1=-1; i1<=1; i1++) {
                if ( ((i1==-1)&&(code&LEFT_SIDE)) || (i1==0) || ((i1==1)&&(code&RIGHT_SIDE)) ) {
                    int bbx=bx+i1;
                    if (bbx<0) {bbx+=num_blocks_x;}
                    if (bbx>=num_blocks_x) {bbx-=num_blocks_x;}
                    for (int i2=-1; i2<=1; i2++) {
                        if ( ((i2==-1)&&(code&BOTTOM_SIDE)) || (i2==0) || ((i2==1)&&(code&TOP_SIDE)) ) {
                            int bby=by+i2;
                            if (bby<0) {bby+=num_blocks_y;}
                            if (bby>=num_blocks_y) {bby-=num_blocks_y;}
                            for (int i3=-1; i3<=1; i3++) {
                                if ( ((i3==-1)&&(code&BACK_SIDE)) || (i3==0) || ((i3==1)&&(code&FRONT_SIDE)) ) {
                                    int bbz=bz+i3;
                                    if (bbz<0) {bbz+=num_blocks_z;}
                                    if (bbz>=num_blocks_z) {bbz-=num_blocks_z;}
                                    int bbb=bbx+num_blocks_x*bby+num_blocks_x*num_blocks_y*bbz;
                                    local_input_block_counts[bbb]++;
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

    printf("Computing sizes and initializing output array... "); timer0.start();
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
    #pragma omp parallel
    {
        #pragma omp for
        for (int ii=0; ii<output_size_times_2; ii++)
            output_tmp[ii]=0;
    }
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
        int bx=bb0%num_blocks_x;
        int by=(bb0/num_blocks_x)%num_blocks_y;
        int bz=(bb0/num_blocks_x/num_blocks_y);
        for (int i1=-1; i1<=1; i1++) {
            if ( ((i1==-1)&&(code&LEFT_SIDE)) || (i1==0) || ((i1==1)&&(code&RIGHT_SIDE)) ) {
                int bbx=bx+i1;
                int wrapx=0;
                if (bbx<0) {bbx+=num_blocks_x; wrapx=N1;}
                if (bbx>=num_blocks_x) {bbx-=num_blocks_x; wrapx=-N1;}
                for (int i2=-1; i2<=1; i2++) {
                    if ( ((i2==-1)&&(code&BOTTOM_SIDE)) || (i2==0) || ((i2==1)&&(code&TOP_SIDE)) ) {
                        int bby=by+i2;
                        int wrapy=0;
                        if (bby<0) {bby+=num_blocks_y; wrapy=N2;}
                        if (bby>=num_blocks_y) {bby-=num_blocks_y; wrapy=-N2;}
                        for (int i3=-1; i3<=1; i3++) {
                            if ( ((i3==-1)&&(code&BACK_SIDE)) || (i3==0) || ((i3==1)&&(code&FRONT_SIDE)) ) {
                                int bbz=bz+i3;
                                int wrapz=0;
                                if (bbz<0) {bbz+=num_blocks_z; wrapz=N3;}
                                if (bbz>=num_blocks_z) {bbz-=num_blocks_z; wrapz=-N3;}
                                int bbb=bbx+num_blocks_x*bby+num_blocks_x*num_blocks_y*bbz;
                                int iii=input_block_indices[bbb]+block_ii[bbb];
                                input_x[iii]=x[m]+wrapx;
                                input_y[iii]=y[m]+wrapy;
                                input_z[iii]=z[m]+wrapz;
                                input_d[iii*2]=d[m*2];
                                input_d[iii*2+1]=d[m*2+1];
                                block_ii[bbb]++;
                            }
                        }
                    }
                }
            }
        }
	}
	printf("Elapsed: %d ms\n",timer0.elapsed());

    printf("spreading []... "); timer0.start();
	#pragma omp parallel
	{
        if (omp_get_thread_num()==0) printf("#################### Using %d threads (%d prescribed)\n",omp_get_num_threads(),opts.num_threads);
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

                int x_integer=ROUND_2_INT(x0);
                double x_diff=x0-x_integer;
                double x_term1=exp(-x_diff*x_diff/(4*opts.tau1));
                double x_term2_factor=exp(-2*x_diff/(4*opts.tau1));
				int xmin=fmax(x_integer-R1/2,block_xmin);
				int xmax=fmin(x_integer-R1/2+R1-1,block_xmax);

                int y_integer=ROUND_2_INT(y0);
                double y_diff=y0-y_integer;
                double y_term1=exp(-y_diff*y_diff/(4*opts.tau2));
                double y_term2_factor=exp(-2*y_diff/(4*opts.tau2));
				int ymin=fmax(y_integer-R1/2,block_ymin);
				int ymax=fmin(y_integer-R1/2+R1-1,block_ymax);

                int z_integer=ROUND_2_INT(z0);
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
            int factor1=K1; if ((cc1+1)*K1>=N1) factor1=N1-cc1*K1;
			int factor2=K2; if ((cc2+1)*K2>=N2) factor2=N2-cc2*K2;
			int factor3=K3; if ((cc3+1)*K3>=N3) factor3=N3-cc3*K3;
			int dd1=cc1*K1;
			int dd2=cc2*K2;
			int dd3=cc3*K3;
			int kkk=output_block_indices[cc]*2; //times 2 because complex
			for (int i3=0; i3<factor3; i3++) {
				for (int i2=0; i2<factor2; i2++) {
					for (int i1=0; i1<factor1; i1++) { //make this inner loop more efficient by not doing the multiplication here?
                        int jjj=(dd1+i1)+(dd2+i2)*N1+(dd3+i3)*N1N2;
                        out[jjj*2]=output_tmp[kkk];
                        out[jjj*2+1]=output_tmp[kkk+1];
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

bool do_fft_3d(int N1,int N2,int N3,double *out,double *in,int num_threads=1) {
    if (num_threads>1) {
        fftw_init_threads();
        fftw_plan_with_nthreads(num_threads);
    }

	int N1N2N3=N1*N2*N3;
	fftw_complex *in2=(fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N1N2N3);
	fftw_complex *out2=(fftw_complex *)fftw_malloc(sizeof(fftw_complex)*N1N2N3);
	for (int ii=0; ii<N1N2N3; ii++) {
		in2[ii][0]=in[ii*2];
		in2[ii][1]=in[ii*2+1];
	}
	fftw_plan p=fftw_plan_dft_3d(N1,N2,N3,in2,out2,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(p);
	for (int ii=0; ii<N1N2N3; ii++) {
		out[ii*2]=out2[ii][0];
		out[ii*2+1]=out2[ii][1];
	}
	fftw_free(in2);
	fftw_free(out2);

    if (num_threads>1) {
        fftw_cleanup_threads();
    }

	return true;
}
void do_fix_3d(int N1,int N2,int N3,int oversamp,double tau,double *out,double *out_oversamp_hat) {
    double *correction_vals1=(double *)malloc(sizeof(double)*(N1/2+2));
    double *correction_vals2=(double *)malloc(sizeof(double)*(N2/2+2));
    double *correction_vals3=(double *)malloc(sizeof(double)*(N3/2+2));
    double tmp1=2*M_PI/(N1*oversamp);
    double tmp2=2*M_PI/(N2*oversamp);
    double tmp3=2*M_PI/(N3*oversamp);
    for (int i=0; i<N1/2+2; i++) {
        correction_vals1[i]=exp(-i*i*tau*tmp1*tmp1);
    }
    for (int i=0; i<N2/2+2; i++) {
        correction_vals2[i]=exp(-i*i*tau*tmp2*tmp2);
    }
    for (int i=0; i<N3/2+2; i++) {
        correction_vals3[i]=exp(-i*i*tau*tmp3*tmp3);
    }
    for (int i3=0; i3<N3; i3++) {
        int aa3=i3*N1*N2;
        int bb3=0;
        double correction3=1;
        if (i3<N3/2) { //should this be <= ?  It doesn't seem to make a significant difference. Not sure why
            bb3=i3*N1*oversamp*N2*oversamp;
            correction3=correction_vals3[i3];
        }
        else {
            bb3=(N3*oversamp-(N3-i3))*N1*oversamp*N2*oversamp;
            correction3=correction_vals3[N3-i3];
        }
        for (int i2=0; i2<N2; i2++) {
            int aa2=i2*N1;
            int bb2=0;
            double correction2=1;
            if (i2<N2/2) {
                bb2=i2*N1*oversamp;
                correction2=correction_vals2[i2];
            }
            else {
                bb2=(N2*oversamp-(N2-i2))*N1*oversamp;
                correction2=correction_vals2[N2-i2];
            }
            aa2+=aa3;
            bb2+=bb3;
            correction2*=correction3;
            for (int i1=0; i1<N1; i1++) {
                int aa1=i1;
                int bb1=0;
                double correction1=1;
                if (i1<N1/2) {
                    bb1=i1;
                    correction1=correction_vals1[i1];
                }
                else {
                    bb1=N1*oversamp-(N1-i1);
                    correction1=correction_vals1[N1-i1];
                }
                aa1+=aa2;
                bb1+=bb2;
                correction1*=correction2;
                out[aa1*2]=out_oversamp_hat[bb1*2]/correction1;
                out[aa1*2+1]=out_oversamp_hat[bb1*2+1]/correction1;
            }
        }
    }
    free(correction_vals1);
    free(correction_vals2);
    free(correction_vals3);
}

/*
void write_mda_re_im(const QString &path,int N1,int N2,int N3,double *X) {
    Mda A_re; A_re.allocate(N1,N2,N3);
    Mda A_im; A_im.allocate(N1,N2,N3);
    double *A_re_ptr=A_re.dataPtr();
    double *A_im_ptr=A_im.dataPtr();
    int N1N2N3=N1*N2*N3;
    for (int i=0; i<N1N2N3; i++) {
        A_re_ptr[i]=X[2*i];
        A_im_ptr[i]=X[2*i+1];
    }
    A_re.write(path+"_re.mda");
    A_im.write(path+"_im.mda");
}
*/

//x,y,z should be in range [0,2pi)
bool blocknufft3d(const BlockNufft3DOptions &opts,double *out,double *x,double *y,double *z,double *d) {
    QTime timer0; QTime timer_total; timer_total.start();

    /*
c     The oversampled regular mesh is defined by
c
c     nf1 = rat*ms  points, where rat is the oversampling ratio.
c     nf2 = rat*mt  points, where rat is the oversampling ratio.
c
c     For simplicity, we set
c
c         rat = 2 for eps > 1.0d-11
c         rat = 3 for eps <= 1.0d-11.
c
c     The Gaussian used for convolution is:
c
c        g(x) = exp(-x^2 / 4tau)
c
c     It can be shown [DR] that the precision eps is achieved when
c
c     nspread = int(-log(eps)/(pi*(rat-1d0)/(rat-.5d0)) + .5d0)
c     and tau is chosen as
c
c     tau = pi*lambda/(ms**2)
c     lambda = nspread/(rat(rat-0.5)).
c     Note that the Fourier transform of g(x) is
c
c     G(s) = exp(-s^2 tau) = exp(-pi*lambda s^2/ms^2)
    */

    double eps=opts.eps;
    int oversamp=3; if (eps<= 1e-11) oversamp=2;
    int R=(int)(-log(eps)/(M_PI*(oversamp-1)/(oversamp-.5)) + .5);
    double lambda=R/(oversamp*(oversamp-0.5));
    double tau=lambda/2; //need to ask Leslie about this...
    printf("Using oversamp=%d, nspread=%d, tau=%g\n",oversamp,R,tau);

	BlockSpread3DOptions sopts;
    sopts.N1o=opts.N1*oversamp; sopts.N2o=opts.N2*oversamp; sopts.N3o=opts.N3*oversamp;
	sopts.K1=opts.K1; sopts.K2=opts.K2; sopts.K3=opts.K3;
	sopts.M=opts.M;
	sopts.num_threads=opts.num_threads;
    sopts.R1=R; sopts.R2=R; sopts.R3=R;
    sopts.tau1=tau; sopts.tau2=tau; sopts.tau3=tau;

    printf("Allocating... "); timer0.start();
	double *out_oversamp=(double *)malloc(sizeof(double)*sopts.N1o*sopts.N2o*sopts.N3o*2);
	double *out_oversamp_hat=(double *)malloc(sizeof(double)*sopts.N1o*sopts.N2o*sopts.N3o*2);
    printf("Elapsed: %d ms\n",timer0.elapsed());

    printf("Scaling coordinates... "); timer0.start();
	double factor_x=sopts.N1o/(2*M_PI);
	double factor_y=sopts.N2o/(2*M_PI);
	double factor_z=sopts.N3o/(2*M_PI);
	for (int ii=0; ii<opts.M; ii++) {
        x[ii]*=factor_x;
        y[ii]*=factor_y;
        z[ii]*=factor_z;
	}
    printf("Elapsed: %d ms\n",timer0.elapsed());

    printf("Spreading...\n"); timer0.start();
    blockspread3d(sopts,out_oversamp,x,y,z,d);
    double spreading_time=timer0.elapsed();
    printf("Total time for spreading: %d ms\n",timer0.elapsed());

    printf("fft... "); timer0.start();
    if (!do_fft_3d(sopts.N1o,sopts.N2o,sopts.N3o,out_oversamp_hat,out_oversamp,opts.num_threads)) {
		printf("problem in do_fft_3d\n");
        for (int ii=0; ii<opts.M; ii++) {
            x[ii]/=factor_x;
            y[ii]/=factor_y;
            z[ii]/=factor_z;
        }
		free(out_oversamp);
		free(out_oversamp_hat);
		return false;
	}
    double fft_time=timer0.elapsed();
    printf("Elapsed: %d ms\n",timer0.elapsed());

    printf("fix... "); timer0.start();
    do_fix_3d(opts.N1,opts.N2,opts.N3,oversamp,tau,out,out_oversamp_hat);
    printf("Elapsed: %d ms\n",timer0.elapsed());

    printf("Restoring coordinates... "); timer0.start();
    for (int ii=0; ii<opts.M; ii++) {
        x[ii]/=factor_x;
        y[ii]/=factor_y;
        z[ii]/=factor_z;
    }
    printf("Elapsed: %d ms\n",timer0.elapsed());

    printf("free...\n"); timer0.start();
	free(out_oversamp_hat);
	free(out_oversamp);
    printf("Elapsed: %d ms\n",timer0.elapsed());

    double total_time=timer_total.elapsed();
    double other_time=total_time-spreading_time-fft_time;

    printf("Elapsed time: %.3f seconds: (%.3f spreading, %.3f fft, %.3f other) (%.1f%% spreading, %.1f%% fft, %.1f%% other)\n",
            total_time/1000,spreading_time/1000,fft_time/1000,other_time/1000,spreading_time/total_time*100,fft_time/total_time*100,other_time/total_time*100);

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
        d[m*2]=0;//sin(m+0.1);
		d[m*2+1]=0;
	}
    x[0]=20*M_PI/180; y[0]=10*M_PI/180; z[0]=5*M_PI/180; d[0]=1;
    double sum_d[2]={0,0};
    for (int i=0; i<opts.M; i++) {
        sum_d[0]+=d[i*2];
        sum_d[1]+=d[i*2+1];
    }

	if (1) {
		QTime timer; timer.start();
		printf("blocknufft3d\n");
		blocknufft3d(opts,out1,x,y,z,d);
		for (int aa=0; aa<10; aa++) printf("%d:%g, ",aa,out1[aa]); printf("\n");
        printf("sum_d = %g + %gi\n",sum_d[0],sum_d[1]);
		printf("Time for blocknufft3d: %d ms\n\n",timer.elapsed());
	}
	if (1) {
		QTime timer; timer.start();
        opts.K1=opts.N1*10; opts.K2=opts.N2*10; opts.K3=opts.K3*10;
		printf("blocknufft3d\n");
		blocknufft3d(opts,out2,x,y,z,d);
		for (int aa=0; aa<10; aa++) printf("%d:%g, ",aa,out2[aa]); printf("\n");
        printf("Time for blocknufft3d: %d ms\n\n",timer.elapsed());

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

void blocknufft3d(int N1,int N2,int N3,int M,double *out,double *xyz,double *d,double eps,int K1,int K2,int K3,int num_threads) {
    BlockNufft3DOptions opts;
    opts.eps=eps;
    opts.K1=K1; opts.K2=K2; opts.K3=K3;
    opts.N1=N1; opts.N2=N2; opts.N3=N3;
    opts.M=M;
    opts.num_threads=num_threads;

    double *x=&xyz[0];
    double *y=&xyz[M];
    double *z=&xyz[2*M];
    blocknufft3d(opts,out,x,y,z,d);

}

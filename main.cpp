#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "omp.h"
#include "blocknufft3d.h"
#include <chrono>
#include "qute.h"

int main(int argc, char *argv[])
{
	if (0) {
		BlockSpread3DOptions opts3;
        opts3.N1o=1e3; opts3.N2o=1e3; opts3.N3o=1e2;
		opts3.M=1e7;
		opts3.R1=5; opts3.R2=5; opts3.R3=5;
        opts3.K1=100; opts3.K2=100; opts3.K3=100;
		//opts3.K1=opts3.N1; opts3.K2=opts3.N2; opts3.K3=opts3.N3;
		opts3.num_threads=10;
		opts3.tau1=1; opts3.tau2=1; opts3.tau3=1;
		test_blockspread3d(opts3);
	}
    if (0) {
		BlockNufft3DOptions opts3;
        opts3.N1=256; opts3.N2=256; opts3.N3=16;
        opts3.M=1e6;
        //opts3.N1=66; opts3.N2=66; opts3.N3=66;
        //opts3.M=1;
        opts3.eps=1e-10;
        opts3.K1=1000; opts3.K2=100; opts3.K3=100;
        //opts3.K1=opts3.N1*2; opts3.K2=opts3.N2*2; opts3.K3=opts3.N3*2;
        opts3.num_threads=10;
		test_blocknufft3d(opts3);
	}

    {
        double *out=(double *)malloc(sizeof(double)*3*3*3*2);
        double *xyz=(double *)malloc(sizeof(double)*1*3);
        double *d=(double *)malloc(sizeof(double)*1*2);
        xyz[0]=xyz[1]=xyz[2]=0;
        d[0]=1; d[1]=0;
        blocknufft3d(3,3,3,1,out,xyz,d,1e-5,10,10,10,1);
        printf("%g,%g,%g,%g\n",out[0],out[2],out[8],out[22]);
        free(out);
        free(xyz);
        free(d);
    }

	return 0;
}

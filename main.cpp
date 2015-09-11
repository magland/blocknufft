#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "omp.h"
#include <QTime>
#include <QDebug>
#include "blocknufft1d.h"
#include "blocknufft2d.h"
#include "blocknufft3d.h"

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
	if (1) {
		BlockNufft3DOptions opts3;
		opts3.N1=1e2; opts3.N2=1e2; opts3.N3=1e2;
		opts3.M=100;
		opts3.R1=5; opts3.R2=5; opts3.R3=5;
		//opts3.K1=100; opts3.K2=100; opts3.K3=100;
		opts3.K1=opts3.N1*2; opts3.K2=opts3.N2*2; opts3.K3=opts3.N3*2;
		opts3.num_threads=1;
		opts3.tau1=1; opts3.tau2=1; opts3.tau3=1;
		opts3.oversamp=2;
		test_blocknufft3d(opts3);
	}

	return 0;
}

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
	BlockNufft1DOptions opts1;
	opts1.N=1e7;
	opts1.M=1e8;
	opts1.R=50;
	opts1.K=1e4;
	opts1.num_threads=6;
	//test_blockspread1d(opts1);

	BlockNufft2DOptions opts2;
	opts2.N1=1e3; opts2.N2=1e4;
	opts2.M=1e8;
	opts2.R1=6; opts2.R2=6;
	opts2.K1=200; opts2.K2=200;
	//opts2.K1=opts2.N1; opts2.K2=opts2.N2;
	opts2.num_threads=6;
	//test_blockspread2d(opts2);

	BlockNufft3DOptions opts3;
    opts3.N1=1e2; opts3.N2=1e2; opts3.N3=1e2;
    opts3.M=1e7;
    opts3.R1=10; opts3.R2=10; opts3.R3=10;
    opts3.K1=20; opts3.K2=20; opts3.K3=10;
	//opts3.K1=opts3.N1; opts3.K2=opts3.N2; opts3.K3=opts3.N3;
    opts3.num_threads=20;
    opts3.tau1=1; opts3.tau2=1; opts3.tau3=1;
	test_blockspread3d(opts3);

	return 0;
}

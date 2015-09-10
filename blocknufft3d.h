#ifndef BLOCKNUFFT3D_H
#define BLOCKNUFFT3D_H

struct BlockNufft3DOptions {
	int N1,N2,N3; //number of grid points
	int K1,K2,K3; //block size
	int R1,R2,R3; //kernel size
	int M; //number of non-uniform points
	int num_threads;
};

bool blockspread3d(const BlockNufft3DOptions &opts,double *out,double *x,double *y,double *z,double *d);
void test_blockspread3d(BlockNufft3DOptions &opts);

#endif // BLOCKNUFFT3D_H


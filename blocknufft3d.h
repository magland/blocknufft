#ifndef BLOCKNUFFT3D_H
#define BLOCKNUFFT3D_H

struct BlockNufft3DOptions {
	int N1,N2,N3; //uniform grid size
	int K1,K2,K3; //block size
	int R1,R2,R3; //kernel size
	int M; //number of non-uniform points
	int num_threads; //
	double tau1,tau2,tau3; //spreading factor
	int oversamp; //oversampling factor
};

//x,y,z should be in range [0,2pi)
bool blocknufft3d(const BlockNufft3DOptions &opts,double *out,double *x,double *y,double *z,double *d);
void test_blocknufft3d(BlockNufft3DOptions &opts);

struct BlockSpread3DOptions {
	int N1o,N2o,N3o; //size of grid
	int K1,K2,K3; //block size
	int R1,R2,R3; //kernel size
	int M; //number of non-uniform points
    int num_threads; //
    double tau1,tau2,tau3; //spreading factor
};

bool blockspread3d(const BlockSpread3DOptions &opts,double *out,double *x,double *y,double *z,double *d);
void test_blockspread3d(BlockSpread3DOptions &opts);

#endif // BLOCKNUFFT3D_H

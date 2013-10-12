
#include <stdio.h>
#include <stdlib.h>

#include <xmmintrin.h>		// SSE

#include "COclManager.h"

//#define FILENAME_MS3D "data/Dophi.ms3d"
#define KernelFunctionNameString	"transformVectorByMatrix"
#define KernelFileNameString		"MatrixMultPoint2D.cl"


#define OPTIMIZE_OPENCL			1

#define CENTER_ROTATE			0
#define SCALE					0.995f
#define ANGLE_THETA				50
#define PI 3.1415927
#define ANGLE_THETA_RAD			((ANGLE_THETA) * PI /180)

#define OPTIMIZE_UNROLL_LOOP	0 // 展开/合并循环
#define OPTIMIZE_INNER_LOOP		1 // 内层循环 并行化 
#define OPTIMIZE_SSE			1
#define OPTIMIZE_OMP			1
#define SIZE_POINT_PER_TIME		1

#define OPTIMIZE_SERIAL			1 // 串行优化
#define OUTPUT_TEXT_OR_IMAGE	0 // 1文本 ；0图形

#define ONE_MATRIX_EACH_POINT	0

#define SIZE_WIDTH			(1<<10)
#define SIZE_HEIGHT		(SIZE_WIDTH)

#if !ONE_MATRIX_EACH_POINT
#define SIZE_SCALE_M	1
#define SIZE_MATRIX		(1<<0)
#else
#define SIZE_SCALE_M	SIZE_WIDTH
#define SIZE_MATRIX		SIZE_HEIGHT
#endif

#define ELEMENT_COUNT_POINT		2
#define ELEMENT_COUNT_LINE		4
#define ELEMENT_LENGTH_LINE		4

#define ELEMENT_COUNT_MATIRX	(ELEMENT_COUNT_LINE*ELEMENT_LENGTH_LINE)

class CMatrixMultPoint2DWeightOCL
{
public:

	int mmp( bool bMulti );
	void Implement( bool bMulti );
	void Init( );
	void UnInit( );
	bool	verify();

	float* getOutput( );

	void	accumulate();

	void setEnvOpenCL(ENV_OPENCL* pEnv) { m_pEnvOpenCL = pEnv;}
	
	bool ExecuteKernel();
	void SetupKernel();

	CMatrixMultPoint2DWeightOCL();
	~CMatrixMultPoint2DWeightOCL();
protected:
	void mmpSerial( );

	void mmpParallel( );

	void mmpRef( );

	void kernel(  float imgIn[][SIZE_WIDTH][ELEMENT_COUNT_POINT], float imgOut[][SIZE_WIDTH][ELEMENT_COUNT_POINT], float m[][ELEMENT_COUNT_LINE], int i, float rad);
	void kernelElement(float* pIn, float* pOut, float* pMat);

	void kernelSSE(   float* imgIn, float* imgOut, int i, float rad=0.0f );
	void matrixMultiply(float mLeft[][ELEMENT_LENGTH_LINE], float mRight[][ELEMENT_LENGTH_LINE], float mResult[][ELEMENT_LENGTH_LINE]);
private:
	float (*m_imgIn)[SIZE_WIDTH][ELEMENT_COUNT_POINT],(*m_pOutRef)[SIZE_WIDTH][ELEMENT_COUNT_POINT],(*m_imgOut)[SIZE_WIDTH][ELEMENT_COUNT_POINT];
	float (*m_pMat)[ELEMENT_COUNT_LINE][ELEMENT_LENGTH_LINE];
	int	(*m_pIndex)[SIZE_WIDTH];
	int m_nSizePoint, m_nSizeMatrix;
	bool m_bMulti;

	float *m_pIn, *m_pOut;
	float *m_pMatrix;

	ENV_OPENCL*	m_pEnvOpenCL;

	struct  OCLKernelArguments
	{
		cl_mem m_pfInputBuffer ;
		cl_mem m_pfOCLOutputBuffer ;
		
		size_t globalWorkSize[2];
		size_t localWorkSize[2];
	};
	OCLKernelArguments m_argOCL;
};



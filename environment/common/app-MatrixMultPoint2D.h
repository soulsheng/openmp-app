
#include <stdio.h>
#include <stdlib.h>

#include <xmmintrin.h>		// SSE

#define OPTIMIZE_UNROLL_LOOP	0 // 展开/合并循环
#define OPTIMIZE_INNER_LOOP		0 // 内层循环 并行化 
#define OPTIMIZE_SSE			0
#define OPTIMIZE_OMP			0
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

class CMatrixMultPoint2D
{
public:

	int mmp( bool bMulti );
	void Implement( bool bMulti );
	void Init( );
	void UnInit( );
	bool	verify();

	float* getOutput( );

	void	accumulate();

	CMatrixMultPoint2D();
	~CMatrixMultPoint2D();
protected:
	void mmpSerial( );

	void mmpParallel( );

	void mmpRef( );

	void kernel(  float imgIn[][SIZE_WIDTH][ELEMENT_COUNT_POINT], float imgOut[][SIZE_WIDTH][ELEMENT_COUNT_POINT], float m[][ELEMENT_COUNT_LINE], int i );
	void kernelElement(float* pIn, float* pOut, float* pMat);

	void kernelSSE(   float* imgIn, float* imgOut, __m128& m0, __m128& m1, __m128& m2, int i );
	void matrixMultiply(float mLeft[][ELEMENT_LENGTH_LINE], float mRight[][ELEMENT_LENGTH_LINE], float mResult[][ELEMENT_LENGTH_LINE]);

private:
	float (*m_imgIn)[SIZE_WIDTH][ELEMENT_COUNT_POINT],(*m_pOutRef)[SIZE_WIDTH][ELEMENT_COUNT_POINT],(*m_imgOut)[SIZE_WIDTH][ELEMENT_COUNT_POINT];
	float (*m_pMat)[ELEMENT_COUNT_LINE][ELEMENT_LENGTH_LINE];
	int	(*m_pIndex)[SIZE_WIDTH];
	int m_nSizePoint, m_nSizeMatrix;
	bool m_bMulti;

	float *m_pIn, *m_pOut;
	float *m_pMatrix;
};



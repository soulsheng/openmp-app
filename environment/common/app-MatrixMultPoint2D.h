
#include <stdio.h>
#include <stdlib.h>

#define OUTPUT_TEXT_OR_IMAGE	1

#define ONE_MATRIX_EACH_POINT	0

#define SIZE_WIDTH			(1<<10)
#define SIZE_HEIGHT		(1<<10)

#if !ONE_MATRIX_EACH_POINT
#define SIZE_SCALE_M	1
#define SIZE_MATRIX		(1<<0)
#else
#define SIZE_SCALE_M	SIZE_WIDTH
#define SIZE_MATRIX		SIZE_HEIGHT
#endif

#define ELEMENT_COUNT_POINT		2
#define ELEMENT_COUNT_LINE		3
#define ELEMENT_LENGTH_LINE		3

#define ELEMENT_COUNT_MATIRX	(ELEMENT_COUNT_LINE*ELEMENT_LENGTH_LINE)

class CMatrixMultPoint2D
{
public:

	int mmp( bool bMulti );
	void Implement( bool bMulti );
	void Init( );
	void UnInit( );
	bool	verify();

	int* getOutput( );

	void	accumulate();

	CMatrixMultPoint2D();
	~CMatrixMultPoint2D();
protected:
	void mmpSerial( );

	void mmpParallel( );

	void mmpRef( );

	void kernel(int (*imgIn)[SIZE_WIDTH][ELEMENT_COUNT_POINT], int (*imgOut)[SIZE_WIDTH][ELEMENT_COUNT_POINT], float(*pMat)[ELEMENT_LENGTH_LINE], int (*pIndex)[SIZE_WIDTH]);
	void kernelElement(int* pIn, int* pOut, float* pMat);

private:
	int (*m_imgIn)[SIZE_WIDTH][ELEMENT_COUNT_POINT],(*m_pOutRef)[SIZE_WIDTH][ELEMENT_COUNT_POINT],(*m_imgOut)[SIZE_WIDTH][ELEMENT_COUNT_POINT];
	float (*m_pMat)[ELEMENT_COUNT_LINE][ELEMENT_LENGTH_LINE];
	int	(*m_pIndex)[SIZE_WIDTH];
	int m_nSizePoint, m_nSizeMatrix;
	bool m_bMulti;
};



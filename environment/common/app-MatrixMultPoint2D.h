
#include <stdio.h>
#include <stdlib.h>

class CMatrixMultPoint2D
{
public:

	int mmp( bool bMulti );
	void Implement( bool bMulti );
	void Init( );
	void UnInit( );
	bool	verify();

	CMatrixMultPoint2D();
	~CMatrixMultPoint2D();
protected:
	void mmpSerial( );

	void mmpParallel( );

	void mmpRef( );

	void kernel(int* pIn, int* pOut, float* pMat, int* pIndex);
	void kernelElement(int* pIn, int* pOut, float* pMat);

private:
	int *m_pIn,*m_pOutRef,*m_pOut;
	float *m_pMat;
	int	*m_pIndex;
	int m_nSizePoint, m_nSizeMatrix;
	bool m_bMulti;
};




#include <stdio.h>
#include <stdlib.h>

class CMatrixMultPoint3D
{
public:

	int mmp( bool bMulti );
	void Implement( bool bMulti );
	void Init( );
	void UnInit( );
	bool	verify();

	CMatrixMultPoint3D();
	~CMatrixMultPoint3D();
protected:
	void mmpSerial( );

	void mmpParallel( );

	void mmpRef( );

	void kernel(float* pIn, float* pOut, float* pMat, int* pIndex);
	void kernelElement(float* pIn, float* pOut, float* pMat);

private:
	float *m_pIn,*m_pMat,*m_pOut;
	float *m_pOutRef;
	int	*m_pIndex;
	int m_nSizePoint, m_nSizeMatrix;
	bool m_bMulti;
};



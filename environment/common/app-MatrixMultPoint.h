
#include <stdio.h>
#include <stdlib.h>

class CMatrixMultPoint
{
public:

	int mmp( bool bMulti );
	void Implement( bool bMulti );
	void Init( );
	void UnInit( );
	bool	verify();

	CMatrixMultPoint();
	~CMatrixMultPoint();
protected:
	void mmpSerial( );

	void mmpParallel( );

	void mmpRef( );

private:
	float *a,*b,*c;
	float *aRef;
	int m, n;
	bool m_bMulti;
};



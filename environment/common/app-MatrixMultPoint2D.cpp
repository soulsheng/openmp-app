
#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>

#include "app-MatrixMultPoint2D.h"

#define ONE_MATRIX_EACH_POINT	0

#define SIZE_SCALE_P			(1<<10)
#define SIZE_POINT		(1<<10)

#if !ONE_MATRIX_EACH_POINT
	#define SIZE_SCALE_M	1
	#define SIZE_MATRIX		(1<<0)
#else
	#define SIZE_SCALE_M	SIZE_SCALE_P
	#define SIZE_MATRIX		SIZE_POINT
#endif

#define ELEMENT_COUNT_POINT		2
#define ELEMENT_COUNT_LINE		3
#define ELEMENT_LENGTH_LINE		3

#define ELEMENT_COUNT_MATIRX	(ELEMENT_COUNT_LINE*ELEMENT_LENGTH_LINE)


int CMatrixMultPoint2D::mmp(  bool bMulti )
{
	m_bMulti = bMulti;

	Init();

	Implement( bMulti );

	UnInit();

   return(0);
}

void CMatrixMultPoint2D::mmpSerial( )
{
	kernel(m_pIn, m_pOut, m_pMat, m_pIndex);

}


void CMatrixMultPoint2D::mmpParallel( )
{
#pragma omp parallel for
	for (int i=0;i<m_nSizePoint;i++)
	{
		int *pInOne = m_pIn + i*ELEMENT_COUNT_POINT ;
		int *pOutOne = m_pOut + i*ELEMENT_COUNT_POINT ;
#if ONE_MATRIX_EACH_POINT
		float *pMatOne = m_pMat + i*ELEMENT_COUNT_MATIRX ;
#else
		float *pMatOne = m_pMat + m_pIndex[i]*ELEMENT_COUNT_MATIRX ;
#endif
		//kernelElement( pInOne, pOutOne, pMatOne );
		pOutOne[0] =
			pMatOne[0*ELEMENT_LENGTH_LINE+0] * pInOne[0] +
			pMatOne[1*ELEMENT_LENGTH_LINE+0] * pInOne[1] +
			pMatOne[2*ELEMENT_LENGTH_LINE+0] ;

		pOutOne[1] =
			pMatOne[0*ELEMENT_LENGTH_LINE+1] * pInOne[0] +
			pMatOne[1*ELEMENT_LENGTH_LINE+1] * pInOne[1] +
			pMatOne[2*ELEMENT_LENGTH_LINE+1] ;
	} /*-- End of omp parallel for --*/
}

void CMatrixMultPoint2D::Init()
{
	m_nSizePoint = SIZE_POINT * SIZE_SCALE_P;
	m_nSizeMatrix = SIZE_MATRIX * SIZE_SCALE_M;

	m_pIn=(int *)malloc( m_nSizePoint*sizeof(int) *ELEMENT_COUNT_POINT );
	m_pMat=(float *)malloc( m_nSizeMatrix*sizeof(float) *ELEMENT_COUNT_MATIRX );
	m_pOut=(int *)malloc( m_nSizePoint*sizeof(int) *ELEMENT_COUNT_POINT );
	m_pOutRef=(int *)malloc( m_nSizePoint*sizeof(int) *ELEMENT_COUNT_POINT );
	m_pIndex = (int*)malloc( m_nSizePoint*sizeof(int) );

	int i,j;
	for (j=0; j<m_nSizePoint*ELEMENT_COUNT_POINT; j++)
		m_pIn[j] = j;
	
	for (i=0; i<m_nSizeMatrix; i++)
		for (j=0; j<ELEMENT_COUNT_MATIRX; j++)
			m_pMat[i*ELEMENT_COUNT_MATIRX+j] = j/(float)ELEMENT_COUNT_MATIRX;

	for (j=0; j<m_nSizePoint; j++)
		m_pIndex[j] = j%m_nSizeMatrix;
}

void CMatrixMultPoint2D::UnInit()
{
	if (m_pIn)	{ free(m_pIn); m_pIn=NULL; }
	if (m_pMat)	{ free(m_pMat); m_pMat=NULL; }
	if (m_pOut)	{ free(m_pOut); m_pOut=NULL; }
	if (m_pOutRef){ free(m_pOutRef); m_pOutRef=NULL; }
	if (m_pIndex)	{ free(m_pIndex); m_pIndex=NULL; }
}

void CMatrixMultPoint2D::Implement( bool bMulti )
{
	m_bMulti = bMulti;

	if( m_bMulti )
		mmpParallel( );
	else
		mmpSerial( );
}

CMatrixMultPoint2D::CMatrixMultPoint2D()
{
	m_pIn = m_pOut = m_pOutRef = NULL;
	m_pMat = NULL;
	m_pIndex = NULL;
	m_bMulti = false;
}

CMatrixMultPoint2D::~CMatrixMultPoint2D()
{
	UnInit();
}
#define ABS(a, b)  ((a)>(b)?(a-b):(b-a))  
bool CMatrixMultPoint2D::verify()
{
	mmpRef( );

	for (int i=0; i<m_nSizePoint*ELEMENT_COUNT_POINT; i++)
	{
		if ( ABS(m_pOut[i] , m_pOutRef[i]) > 1 )
		{
			return false;
		}
	}
	
	return true;
}

void CMatrixMultPoint2D::mmpRef()
{
	kernel(m_pIn, m_pOutRef, m_pMat, m_pIndex);

}

void CMatrixMultPoint2D::kernel( int* pIn, int* pOut, float* pMat, int* pIndex )
{
	for (int i=0;i<m_nSizePoint;i++)
	{
		int *pInOne = pIn + i*ELEMENT_COUNT_POINT ;
		int *pOutOne = pOut + i*ELEMENT_COUNT_POINT ;
#if ONE_MATRIX_EACH_POINT
		float *pMatOne = m_pMat + i*ELEMENT_COUNT_MATIRX ;
#else
		float *pMatOne = m_pMat + m_pIndex[i]*ELEMENT_COUNT_MATIRX ;
#endif
		//kernelElement( pInOne, pOutOne, pMatOne );
		pOutOne[0] =
			pMatOne[0*ELEMENT_LENGTH_LINE+0] * pInOne[0] +
			pMatOne[1*ELEMENT_LENGTH_LINE+0] * pInOne[1] +
			pMatOne[2*ELEMENT_LENGTH_LINE+0] ;

		pOutOne[1] =
			pMatOne[0*ELEMENT_LENGTH_LINE+1] * pInOne[0] +
			pMatOne[1*ELEMENT_LENGTH_LINE+1] * pInOne[1] +
			pMatOne[2*ELEMENT_LENGTH_LINE+1] ;

	}
	
}

void CMatrixMultPoint2D::kernelElement( int* pIn, int* pOut, float* pMat )
{
	pOut[0] =
		pMat[0*ELEMENT_LENGTH_LINE+0] * pIn[0] +
		pMat[1*ELEMENT_LENGTH_LINE+0] * pIn[1] +
		pMat[2*ELEMENT_LENGTH_LINE+0] ;

	pOut[1] =
		pMat[0*ELEMENT_LENGTH_LINE+1] * pIn[0] +
		pMat[1*ELEMENT_LENGTH_LINE+1] * pIn[1] +
		pMat[2*ELEMENT_LENGTH_LINE+1] ;
}

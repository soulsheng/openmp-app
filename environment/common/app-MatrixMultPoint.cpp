
#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>

#include "app-MatrixMultPoint.h"

#define SIZE_SCALE_P			(1<<13)
#define SIZE_SCALE_M			1
#define SIZE_POINT		(1<<7)
#define SIZE_MATRIX		(1<<7)

#define ELEMENT_COUNT_POINT		3
#define ELEMENT_COUNT_LINE		3

#define ELEMENT_COUNT_MATIRX	(ELEMENT_COUNT_LINE*4)


int CMatrixMultPoint::mmp(  bool bMulti )
{
	m_bMulti = bMulti;

	Init();

	Implement( bMulti );

	UnInit();

   return(0);
}

void CMatrixMultPoint::mmpSerial( )
{
	kernel(m_pIn, m_pOut, m_pMat, m_pIndex);

}


void CMatrixMultPoint::mmpParallel( )
{
#pragma omp parallel for
	for (int i=0;i<m_nSizePoint;i++)
	{
		float *pInOne = m_pIn + i*ELEMENT_COUNT_POINT ;
		float *pOutOne = m_pOut + i*ELEMENT_COUNT_POINT ;
		float *pMatOne = m_pMat + m_pIndex[i]*ELEMENT_COUNT_MATIRX ;

		//kernelElement( pInOne, pOutOne, pMatOne );
		pOutOne[0] =
			(pMatOne[0*4+0] * pInOne[0] +
			pMatOne[1*4+0] * pInOne[1] +
			pMatOne[2*4+0] * pInOne[2] +
			pMatOne[3*4+0]) ;

		pOutOne[1] =
			(pMatOne[0*4+1] * pInOne[0] +
			pMatOne[1*4+1] * pInOne[1] +
			pMatOne[2*4+1] * pInOne[2] +
			pMatOne[3*4+1]) ;

		pOutOne[2] =
			(pMatOne[0*4+2] * pInOne[0] +
			pMatOne[1*4+2] * pInOne[1] +
			pMatOne[2*4+2] * pInOne[2] +
			pMatOne[3*4+2]) ;
	} /*-- End of omp parallel for --*/
}

void CMatrixMultPoint::Init()
{
	m_nSizePoint = SIZE_POINT * SIZE_SCALE_P;
	m_nSizeMatrix = SIZE_MATRIX * SIZE_SCALE_M;

	m_pIn=(float *)malloc( m_nSizePoint*sizeof(float) *ELEMENT_COUNT_POINT );
	m_pMat=(float *)malloc( m_nSizeMatrix*sizeof(float) *ELEMENT_COUNT_MATIRX );
	m_pOut=(float *)malloc( m_nSizePoint*sizeof(float) *ELEMENT_COUNT_POINT );
	m_pOutRef=(float *)malloc( m_nSizePoint*sizeof(float) *ELEMENT_COUNT_POINT );
	m_pIndex = (int*)malloc( m_nSizePoint*sizeof(int) );

	int i,j;
	for (j=0; j<m_nSizePoint*ELEMENT_COUNT_POINT; j++)
		m_pIn[j] = 2.0;
	
	for (i=0; i<m_nSizeMatrix; i++)
		for (j=0; j<ELEMENT_COUNT_MATIRX; j++)
			m_pMat[i*ELEMENT_COUNT_MATIRX+j] = i;

	for (j=0; j<m_nSizePoint; j++)
		m_pIndex[j] = j%ELEMENT_COUNT_MATIRX;
}

void CMatrixMultPoint::UnInit()
{
	if (m_pIn)	{ free(m_pIn); m_pIn=NULL; }
	if (m_pMat)	{ free(m_pMat); m_pMat=NULL; }
	if (m_pOut)	{ free(m_pOut); m_pOut=NULL; }
	if (m_pOutRef){ free(m_pOutRef); m_pOutRef=NULL; }
	if (m_pIndex)	{ free(m_pIndex); m_pIndex=NULL; }
}

void CMatrixMultPoint::Implement( bool bMulti )
{
	m_bMulti = bMulti;

	if( m_bMulti )
		mmpParallel( );
	else
		mmpSerial( );
}

CMatrixMultPoint::CMatrixMultPoint()
{
	m_pIn = m_pOut = m_pOutRef = NULL;
	m_pMat = NULL;
	m_pIndex = NULL;
	m_bMulti = false;
}

CMatrixMultPoint::~CMatrixMultPoint()
{
	UnInit();
}

bool CMatrixMultPoint::verify()
{
	mmpRef( );

	for (int i=0; i<m_nSizePoint*ELEMENT_COUNT_POINT; i++)
	{
		if ( fabs(m_pOut[i] - m_pOutRef[i]) > 1.0e-3 )
		{
			return false;
		}
	}
	
	return true;
}

void CMatrixMultPoint::mmpRef()
{
	kernel(m_pIn, m_pOutRef, m_pMat, m_pIndex);

}

void CMatrixMultPoint::kernel( float* pIn, float* pOut, float* pMat, int* pIndex )
{
	for (int i=0;i<m_nSizePoint;i++)
	{
		float *pInOne = pIn + i*ELEMENT_COUNT_POINT ;
		float *pOutOne = pOut + i*ELEMENT_COUNT_POINT ;
		float *pMatOne = pMat + pIndex[i]*ELEMENT_COUNT_MATIRX ;

		//kernelElement( pInOne, pOutOne, pMatOne );
		pOutOne[0] =
			(pMatOne[0*4+0] * pInOne[0] +
			pMatOne[1*4+0] * pInOne[1] +
			pMatOne[2*4+0] * pInOne[2] +
			pMatOne[3*4+0]) ;

		pOutOne[1] =
			(pMatOne[0*4+1] * pInOne[0] +
			pMatOne[1*4+1] * pInOne[1] +
			pMatOne[2*4+1] * pInOne[2] +
			pMatOne[3*4+1]) ;

		pOutOne[2] =
			(pMatOne[0*4+2] * pInOne[0] +
			pMatOne[1*4+2] * pInOne[1] +
			pMatOne[2*4+2] * pInOne[2] +
			pMatOne[3*4+2]) ;
	}
	
}

void CMatrixMultPoint::kernelElement( float* pIn, float* pOut, float* pMat )
{
	pOut[0] =
		(pMat[0*4+0] * pIn[0] +
		pMat[1*4+0] * pIn[1] +
		pMat[2*4+0] * pIn[2] +
		pMat[3*4+0]) ;

	pOut[1] =
		(pMat[0*4+1] * pIn[0] +
		pMat[1*4+1] * pIn[1] +
		pMat[2*4+1] * pIn[2] +
		pMat[3*4+1]) ;

	pOut[2] =
		(pMat[0*4+2] * pIn[0] +
		pMat[1*4+2] * pIn[1] +
		pMat[2*4+2] * pIn[2] +
		pMat[3*4+2]) ;
}

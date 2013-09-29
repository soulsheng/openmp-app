
#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>

#include "app-MatrixMultPoint2D.h"



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
	kernel(imgIn, imgOut, m_pMat[0], m_pIndex);

}


void CMatrixMultPoint2D::mmpParallel( )
{
	float(*m)[ELEMENT_LENGTH_LINE] = m_pMat[0];

#pragma omp parallel for
	for (int i=0;i<SIZE_HEIGHT;i++)
	{
		for (int j=0;j<SIZE_WIDTH;j++)
		{
			int x =	imgIn[i][j][0] ;
			int y =	imgIn[i][j][1] ;

			imgOut[i][j][0] = x * m[0][0] + y * m[1][0] + m[2][0] ;
			imgOut[i][j][1] = x * m[0][1] + y * m[1][1] + m[2][1] ;
		}
	} /*-- End of omp parallel for --*/
}

void CMatrixMultPoint2D::Init()
{
	m_nSizePoint = SIZE_HEIGHT * SIZE_WIDTH;
	m_nSizeMatrix = SIZE_MATRIX * SIZE_SCALE_M;
#if 0
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
#else
	imgIn = new int[SIZE_HEIGHT][SIZE_WIDTH][ELEMENT_COUNT_POINT];
	imgOut = new int[SIZE_HEIGHT][SIZE_WIDTH][ELEMENT_COUNT_POINT];
	m_pOutRef = new int[SIZE_HEIGHT][SIZE_WIDTH][ELEMENT_COUNT_POINT];
	
	m_pIndex = new int[SIZE_HEIGHT][SIZE_WIDTH];
	
	m_pMat=new float[m_nSizeMatrix][ELEMENT_COUNT_LINE][ELEMENT_LENGTH_LINE];

	for (int i=0;i<SIZE_HEIGHT;i++)
		for (int j=0;j<SIZE_WIDTH;j++)
			for (int k=0;k<ELEMENT_COUNT_POINT;k++)
				imgIn[i][j][k] = ((i*SIZE_WIDTH + j)*ELEMENT_COUNT_POINT ) +k;

	for (int i=0; i<m_nSizeMatrix; i++)
		for (int j=0;j<ELEMENT_COUNT_LINE;j++)
			for (int k=0;k<ELEMENT_LENGTH_LINE;k++)
				m_pMat[i][j][k] = (j*ELEMENT_COUNT_LINE + k)/(float)(ELEMENT_LENGTH_LINE*ELEMENT_COUNT_LINE) ;

	for (int i=0;i<SIZE_HEIGHT;i++)
		for (int j=0;j<SIZE_WIDTH;j++)
			m_pIndex[i][j] = j%m_nSizeMatrix;
#endif
}

void CMatrixMultPoint2D::UnInit()
{
	if (imgIn)	{ delete[](imgIn); imgIn=NULL; }
	if (m_pMat)	{ delete[](m_pMat); m_pMat=NULL; }
	if (imgOut)	{ delete[](imgOut); imgOut=NULL; }
	if (m_pOutRef){ delete[](m_pOutRef); m_pOutRef=NULL; }
	if (m_pIndex)	{ delete[](m_pIndex); m_pIndex=NULL; }
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
	imgIn = imgOut = m_pOutRef = NULL;
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
		if ( ABS(imgOut[i] , m_pOutRef[i]) > 1 )
		{
			return false;
		}
	}
	
	return true;
}

void CMatrixMultPoint2D::mmpRef()
{
	kernel(imgIn, m_pOutRef, m_pMat[0], m_pIndex);

}

void CMatrixMultPoint2D::kernel( int (*pIn)[SIZE_WIDTH][ELEMENT_COUNT_POINT], int (*pOut)[SIZE_WIDTH][ELEMENT_COUNT_POINT], float(*m)[ELEMENT_LENGTH_LINE], int (*pIndex)[SIZE_WIDTH] )
{
	for (int i=0;i<SIZE_HEIGHT;i++)
	{
		for (int j=0;j<SIZE_WIDTH;j++)
		{
			int x =	imgIn[i][j][0] ;
			int y =	imgIn[i][j][1] ;

			imgOut[i][j][0] = x * m[0][0] + y * m[1][0] + m[2][0] ;
			imgOut[i][j][1] = x * m[0][1] + y * m[1][1] + m[2][1] ;
		}

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

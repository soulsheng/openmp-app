
#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>

#include "app-MatrixMultPoint.h"

#define SIZE_SCALE_M			16
#define SIZE_SCALE_N			16
#define SIZE_MATRIX_M		100
#define SIZE_MATRIX_N		100

/*
   a[m][1]
   =
   b[m][n] *  c[n][1]
*/
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
   int i, j;

   for (i=0; i<m; i++)
   {
      a[i] = 0.0;
      for (j=0; j<n; j++)
         a[i] += b[i*n+j]*c[j];
   }
}


void CMatrixMultPoint::mmpParallel( )
{
	int i, j;

#pragma omp parallel for   private(i,j)
	for (i=0; i<m; i++)
	{
		a[i] = 0.0;
		for (j=0; j<n; j++)
			a[i] += b[i*n+j]*c[j];
	} /*-- End of omp parallel for --*/
}

void CMatrixMultPoint::Init()
{
	m = SIZE_MATRIX_M * SIZE_SCALE_M;
	n = SIZE_MATRIX_N * SIZE_SCALE_N;

	a=(float *)malloc(m*sizeof(float));
	b=(float *)malloc(m*n*sizeof(float));
	c=(float *)malloc(n*sizeof(float));
	aRef=(float *)malloc(m*sizeof(float));

	int i,j;
	for (j=0; j<n; j++)
		c[j] = 2.0;
	for (i=0; i<m; i++)
		for (j=0; j<n; j++)
			b[i*n+j] = i;

}

void CMatrixMultPoint::UnInit()
{
	if (a)	{ free(a); a=NULL; }
	if (b)	{ free(b); b=NULL; }
	if (c)	{ free(c); c=NULL; }
	if (aRef)	{ free(aRef); aRef=NULL; }
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
	a = b = c = NULL;
	aRef = NULL;
	m_bMulti = false;
}

CMatrixMultPoint::~CMatrixMultPoint()
{
	UnInit();
}

bool CMatrixMultPoint::verify()
{
	mmpRef( );

	for (int i=0; i<m; i++)
	{
		if ( fabs(a[i] - aRef[i]) > 1.0e-3 )
		{
			return false;
		}
	}
	
	return true;
}

void CMatrixMultPoint::mmpRef()
{
	int i, j;

	for (i=0; i<m; i++)
	{
		aRef[i] = 0.0;
		for (j=0; j<n; j++)
			aRef[i] += b[i*n+j]*c[j];
	}
}

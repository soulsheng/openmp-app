
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
	float m[ELEMENT_COUNT_LINE][ELEMENT_LENGTH_LINE];
	memcpy( m, m_pMat[0], sizeof(float) * ELEMENT_COUNT_LINE * ELEMENT_LENGTH_LINE ) ;

#if 1
	__m128 mat[3];
	mat[0] = _mm_load_ps( m_pMatrix ) ;
	mat[1] = _mm_load_ps( m_pMatrix+4 ) ;
	mat[2] = _mm_load_ps( m_pMatrix+8 ) ;
	for (int i=0;i<SIZE_HEIGHT;i++)
	{
#if OPTIMIZE_SSE
		kernelSSE(m_pIn, m_pOut, mat[0], mat[1], mat[2], i);
#else
		kernel(m_imgIn, m_imgOut, m, i);
#endif
	}
#else
	for( int i = 0; i < SIZE_HEIGHT; i ++ )
		for( int j = 0; j < SIZE_WIDTH; j++ )	{
			m_imgOut[i][j][0] 
				= m_imgIn[i][j][0]*m[0][0] + m_imgIn[i][j][1]*m[1][0] + m[2][0];
			m_imgOut[i][j][1] 
				= m_imgIn[i][j][0]*m[0][1] + m_imgIn[i][j][1]*m[1][1] + m[2][1];
		}
#endif

}


void CMatrixMultPoint2D::mmpParallel( )
{
	float m[ELEMENT_COUNT_LINE][ELEMENT_LENGTH_LINE];
	memcpy( m, m_pMat[0], sizeof(float) * ELEMENT_COUNT_LINE * ELEMENT_LENGTH_LINE ) ;

	float (*imgIn)[SIZE_WIDTH][ELEMENT_COUNT_POINT] ;
	imgIn = m_imgIn;
	float (*imgOut)[SIZE_WIDTH][ELEMENT_COUNT_POINT] ;
	imgOut = m_imgOut;

	__m128 mat[3];
	mat[0] = _mm_load_ps( m_pMatrix ) ;
	mat[1] = _mm_load_ps( m_pMatrix+4 ) ;
	mat[2] = _mm_load_ps( m_pMatrix+8 ) ;
#pragma omp parallel for 
	for (int i=0;i<SIZE_HEIGHT;i++)
	{
#if OPTIMIZE_SSE
		kernelSSE( m_pIn, m_pOut, mat[0], mat[1], mat[2], i);
#else
		kernel(m_imgIn, m_imgOut, m, i);
#endif
	}

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
	m_imgIn = new float[SIZE_HEIGHT][SIZE_WIDTH][ELEMENT_COUNT_POINT];
	m_imgOut = new float[SIZE_HEIGHT][SIZE_WIDTH][ELEMENT_COUNT_POINT];
	//m_pOutRef = new float[SIZE_HEIGHT][SIZE_WIDTH][ELEMENT_COUNT_POINT];
	
	m_pIndex = new int[SIZE_HEIGHT][SIZE_WIDTH];
	
	m_pMat=new float[m_nSizeMatrix][ELEMENT_COUNT_LINE][ELEMENT_LENGTH_LINE];
	
	m_pIn = (float*)_aligned_malloc( SIZE_HEIGHT*SIZE_WIDTH*ELEMENT_COUNT_POINT*sizeof(float), 16 );//new float[];
	m_pOut = (float*)_aligned_malloc( SIZE_HEIGHT*SIZE_WIDTH*ELEMENT_COUNT_POINT*sizeof(float), 16 );//new float[];
	m_pMatrix = (float*)_aligned_malloc( ELEMENT_COUNT_LINE*ELEMENT_LENGTH_LINE*sizeof(float), 16 );//new float[];

	for (int i=0;i<SIZE_HEIGHT;i++)
		for (int j=0;j<SIZE_WIDTH;j++)
		{
			m_imgIn[i][j][0] = i;
			m_imgIn[i][j][1] = j;

			m_imgOut[i][j][0] = i;
			m_imgOut[i][j][1] = j;
		}
	memcpy( m_pIn, &m_imgIn[0][0][0], SIZE_HEIGHT*SIZE_WIDTH*ELEMENT_COUNT_POINT*sizeof(float));
	memcpy( m_pOut, &m_imgOut[0][0][0], SIZE_HEIGHT*SIZE_WIDTH*ELEMENT_COUNT_POINT*sizeof(float));
	
	for (int i=0; i<m_nSizeMatrix; i++)
		for (int j=0;j<ELEMENT_COUNT_LINE;j++)
			for (int k=0;k<ELEMENT_LENGTH_LINE;k++)
				m_pMat[i][j][k] = 0;//(j*ELEMENT_COUNT_LINE + k)/(float)(ELEMENT_LENGTH_LINE*ELEMENT_COUNT_LINE) ;

	for (int i=0;i<SIZE_HEIGHT;i++)
		for (int j=0;j<SIZE_WIDTH;j++)
			m_pIndex[i][j] = j%m_nSizeMatrix;

	float theta = 15;
	float scale = 1.05;
#define PI 3.1415927
	float mr[ELEMENT_LENGTH_LINE][ELEMENT_LENGTH_LINE] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} } ; // 旋转变换矩阵初始化
	float ms[ELEMENT_LENGTH_LINE][ELEMENT_LENGTH_LINE] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} } ; // 缩放变换矩阵初始化
	float (*m)[ELEMENT_LENGTH_LINE] = m_pMat[0] ;	 // 统一变换矩阵初始化

	// 旋转矩阵
	float rad = theta /180 * PI ;
	mr[0][0] = cos( rad ) ;
	mr[0][1] = sin( rad ) ;
	mr[1][0] = -sin( rad ) ;
	mr[1][1] = cos( rad ) ;

	// 缩放矩阵
	ms[0][0] = scale ;
	ms[1][1] = scale ;

	// 获得统一变换矩阵
	for( int i = 0 ; i < ELEMENT_LENGTH_LINE ; i++ )
		for( int j = 0 ; j < ELEMENT_LENGTH_LINE ; j++ )
		{
			for( int k = 0 ; k < ELEMENT_LENGTH_LINE ; k++ )
			{
				m[i][j] += ms[i][k] * mr[k][j];
			}
			*(m_pMatrix+j*ELEMENT_LENGTH_LINE+i) = m[i][j];
		}
#endif
}

void CMatrixMultPoint2D::UnInit()
{
	if (m_imgIn)	{ delete[](m_imgIn); m_imgIn=NULL; }
	if (m_pMat)	{ delete[](m_pMat); m_pMat=NULL; }
	if (m_imgOut)	{ delete[](m_imgOut); m_imgOut=NULL; }
	if (m_pOutRef){ delete[](m_pOutRef); m_pOutRef=NULL; }
	if (m_pIndex)	{ delete[](m_pIndex); m_pIndex=NULL; }
	if (m_pMatrix)	{ _aligned_free(m_pMatrix); m_pMatrix=NULL; }
	if (m_pIn)	{ _aligned_free(m_pIn); m_pIn=NULL; }
	if (m_pOut)	{ _aligned_free(m_pOut); m_pOut=NULL; }
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
	m_imgIn = m_imgOut = m_pOutRef = NULL;
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
		//if ( ABS(m_imgOut[i] , m_pOutRef[i]) > 1 )
		{
			return false;
		}
	}
	
	return true;
}

void CMatrixMultPoint2D::mmpRef()
{
	//kernel(m_imgIn, m_pOutRef, m_pMat[0], m_pIndex);

}

void CMatrixMultPoint2D::kernel( float imgIn[][SIZE_WIDTH][ELEMENT_COUNT_POINT], float imgOut[][SIZE_WIDTH][ELEMENT_COUNT_POINT], float m[][ELEMENT_COUNT_LINE], int i )
{
	for (int j=0;j<SIZE_WIDTH ;j+=SIZE_POINT_PER_TIME)
	{
#if OPTIMIZE_SERIAL
		for (int k=0; k< SIZE_POINT_PER_TIME ; k++ )
		{
			float x =	imgIn[i][j+k][0] ;
			float y =	imgIn[i][j+k][1] ;

			imgOut[i][j+k][0] = x * m[0][0] + y * m[1][0] + m[2][0] ;
			imgOut[i][j+k][1] = x * m[0][1] + y * m[1][1] + m[2][1] ;
		}
#else
		int k=0;
		for (int k=0; k< SIZE_POINT_PER_TIME ; k++ )
		{
			imgOut[i][j+k][0] = imgIn[i][j+k][0] * m[0][0] + imgIn[i][j+k][1] * m[1][0] + m[2][0] ;
			imgOut[i][j+k][1] = imgIn[i][j+k][0] * m[0][1] + imgIn[i][j+k][1] * m[1][1] + m[2][1] ;
		}
#endif
	}

	
}

void CMatrixMultPoint2D::kernelElement( float* pIn, float* pOut, float* pMat )
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

float* CMatrixMultPoint2D::getOutput()
{	
#if !OPTIMIZE_SSE
	return &m_imgOut[0][0][0]; 
#else
	return m_pOut;
#endif
}

void CMatrixMultPoint2D::accumulate()
{
	//kernel(m_imgIn, m_imgOut, m_pMat[0], m_pIndex);
	Implement(true);
	memcpy( m_imgIn, m_imgOut, sizeof(float) * SIZE_HEIGHT * SIZE_WIDTH * ELEMENT_COUNT_POINT );
	memcpy( m_pIn, m_pOut, sizeof(float) * SIZE_HEIGHT * SIZE_WIDTH * ELEMENT_COUNT_POINT );
}


//#define _MM_SHUFFLE(fp3,fp2,fp1,fp0) (((fp3) << 6) | ((fp2) << 4) | \
//	((fp1) << 2) | ((fp0)))


#define __MM_SELECT(v, fp)                                                          \
	_mm_shuffle_ps((v), (v), _MM_SHUFFLE((fp),(fp),(fp),(fp)))

void CMatrixMultPoint2D::kernelSSE( float* imgIn, float* imgOut, __m128& m0, __m128& m1, __m128& m2, int i )
{
	__m128 *pSrcPos = (__m128*)imgIn;
	__m64 *pDestPos = (__m64*)imgOut;

	for (int j=0;j<SIZE_WIDTH/2 ;j++)
	{
		__m128 vIn0x, vIn0y, vIn1x, vIn1y;

		int index = j+i*SIZE_WIDTH/2;
		vIn0x = __MM_SELECT( pSrcPos[index], 0); 
		vIn0y = __MM_SELECT( pSrcPos[index], 1); 

		__m128 vOut0;
		__m128 tmp00, tmp01;
		tmp00 = _mm_mul_ps( m0, vIn0x );
		tmp01 = _mm_mul_ps( m1, vIn0y );

		vOut0 = _mm_add_ps( tmp00, tmp01 );
		vOut0 = _mm_add_ps( vOut0, m2 );
		_mm_storel_pi( (__m64*)(&pDestPos[index*2]), vOut0 );
		
		vIn1x = __MM_SELECT( pSrcPos[index], 2); 
		vIn1y = __MM_SELECT( pSrcPos[index], 3); 
		__m128 vOut1;
		__m128 tmp10, tmp11;
		tmp10 = _mm_mul_ps( m0, vIn1x );
		tmp11 = _mm_mul_ps( m1, vIn1y );

		vOut1 = _mm_add_ps( tmp10, tmp11 );
		vOut1 = _mm_add_ps(    vOut1, m2 );
		_mm_storel_pi( (__m64*)(&pDestPos[index*2+1]), vOut1 );
		
	}


}

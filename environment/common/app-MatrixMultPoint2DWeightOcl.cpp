
#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>

#include "app-MatrixMultPoint2DWeightOCL.h"
#include "utils.h"


int CMatrixMultPoint2DWeightOCL::mmp(  bool bMulti )
{
	m_bMulti = bMulti;

	Init();

	Implement( bMulti );

	UnInit();

   return(0);
}

void CMatrixMultPoint2DWeightOCL::mmpSerial( )
{
	float rad = ANGLE_THETA_RAD ;

#if OPTIMIZE_SSE

	__m128 mat[3], matPad[3];
	mat[0] = _mm_load_ps( m_pMatrix ) ;
	mat[1] = _mm_load_ps( m_pMatrix+4 ) ;
	mat[2] = _mm_load_ps( m_pMatrix+8 ) ;

	matPad[0] = _mm_shuffle_ps( mat[0], mat[0], _MM_SHUFFLE(1,0,1,0) ); // m00, m01, m00, m01
	matPad[1] = _mm_shuffle_ps( mat[1], mat[1], _MM_SHUFFLE(1,0,1,0) ); // m10, m11, m10, m11
	matPad[2] = _mm_shuffle_ps( mat[2], mat[2], _MM_SHUFFLE(1,0,1,0) ); // m20, m21, m20, m21

	for (int i=0;i<SIZE_HEIGHT;i++)
		kernelSSE(m_pIn, m_pOut, i, rad);
#else

	float m[ELEMENT_COUNT_LINE][ELEMENT_LENGTH_LINE];
	memcpy( m, m_pMat[0], sizeof(float) * ELEMENT_COUNT_LINE * ELEMENT_LENGTH_LINE ) ;

#if OPTIMIZE_UNROLL_LOOP
	float *pIn = &m_imgIn[0][0][0];
	float *pOut = &m_imgOut[0][0][0];

	for( int i = 0; i < SIZE_HEIGHT*SIZE_WIDTH; i ++ )
		{
			pOut[2*i]	= pIn[2*i]*m[0][0] + pIn[2*i+1]*m[1][0] + m[2][0];
			pOut[2*i+1]	= pIn[2*i]*m[0][1] + pIn[2*i+1]*m[1][1] + m[2][1];
		}
#else

#if !OPTIMIZE_SERIAL
	for( int i = 0; i < SIZE_HEIGHT; i ++ )
		for( int j = 0; j < SIZE_WIDTH; j++ )	{
			m_imgOut[i][j][0] 
				= m_imgIn[i][j][0]*m[0][0] + m_imgIn[i][j][1]*m[1][0] + m[2][0];
			m_imgOut[i][j][1] 
				= m_imgIn[i][j][0]*m[0][1] + m_imgIn[i][j][1]*m[1][1] + m[2][1];
		}
#else
	for (int i=0;i<SIZE_HEIGHT;i++)
		kernel(m_imgIn, m_imgOut, m, i, rad);

#endif //OPTIMIZE_SERIAL

#endif // OPTIMIZE_UNROLL_LOOP

#endif // OPTIMIZE_SSE

}


void CMatrixMultPoint2DWeightOCL::mmpParallel( )
{
	float rad = ANGLE_THETA_RAD ;

#if OPTIMIZE_SSE


#if OPTIMIZE_OMP && !OPTIMIZE_INNER_LOOP
#pragma omp parallel for //shared(nOffset)
#endif	
	for (int i=0;i<SIZE_HEIGHT;i++)
		kernelSSE(m_pIn, m_pOut, i, rad);
#else

	float m[ELEMENT_COUNT_LINE][ELEMENT_LENGTH_LINE];
	memcpy( m, m_pMat[0], sizeof(float) * ELEMENT_COUNT_LINE * ELEMENT_LENGTH_LINE ) ;

#if OPTIMIZE_UNROLL_LOOP
	float *pIn = &m_imgIn[0][0][0];
	float *pOut = &m_imgOut[0][0][0];

#if OPTIMIZE_OMP
#pragma omp parallel for
#endif
	for( int i = 0; i < SIZE_HEIGHT*SIZE_WIDTH; i ++ )
	{
		pOut[2*i]	= pIn[2*i]*m[0][0] + pIn[2*i+1]*m[1][0] + m[2][0];
		pOut[2*i+1]	= pIn[2*i]*m[0][1] + pIn[2*i+1]*m[1][1] + m[2][1];
	}
#else

#if !OPTIMIZE_SERIAL
	for( int i = 0; i < SIZE_HEIGHT; i ++ )
		for( int j = 0; j < SIZE_WIDTH; j++ )	{
			float x = m_imgIn[i][j][0];
			float y = m_imgIn[i][j][1];

			float offset = 4.0f;
			float distance = sqrt( x*x + y*y + offset) ;
			float weight = 1.0f / sqrt(distance)  ;

			float radOne = rad * weight;
			float mOne[2][2];

			mOne[0][0] = cos( radOne ) ;
			mOne[0][1] = sin( radOne ) ;
			mOne[1][0] = -sin( radOne ) ;
			mOne[1][1] = cos( radOne ) ;

			// ���ž��� ƽ�ƾ���
			//m[0][0] *= SCALE ;
			//m[1][1] *= SCALE ;

			m_imgOut[i][j][0] = x * mOne[0][0] + y * mOne[1][0] ;
			m_imgOut[i][j][1] = x * mOne[0][1] + y * mOne[1][1] ;
		}
#else
	for (int i=0;i<SIZE_HEIGHT;i++)
		kernel(m_imgIn, m_imgOut, m, i, rad);

#endif //OPTIMIZE_SERIAL

#endif // OPTIMIZE_UNROLL_LOOP

#endif // OPTIMIZE_SSE

}

void CMatrixMultPoint2DWeightOCL::Init()
{
	m_nSizePoint = SIZE_HEIGHT * SIZE_WIDTH;
	m_nSizeMatrix = SIZE_MATRIX * SIZE_SCALE_M;

#if !OPTIMIZE_SSE
	m_imgIn = new float[SIZE_HEIGHT][SIZE_WIDTH][ELEMENT_COUNT_POINT];
	m_imgOut = new float[SIZE_HEIGHT][SIZE_WIDTH][ELEMENT_COUNT_POINT];
	//m_pOutRef = new float[SIZE_HEIGHT][SIZE_WIDTH][ELEMENT_COUNT_POINT];
	
	m_pIndex = new int[SIZE_HEIGHT][SIZE_WIDTH];
	
	m_pMat=new float[m_nSizeMatrix][ELEMENT_COUNT_LINE][ELEMENT_LENGTH_LINE];
	for (int i=0;i<SIZE_HEIGHT;i++)
		for (int j=0;j<SIZE_WIDTH;j++)
		{
			m_imgIn[i][j][0] = i;
			m_imgIn[i][j][1] = j;

			m_imgOut[i][j][0] = i;
			m_imgOut[i][j][1] = j;
		}

	for (int i=0; i<m_nSizeMatrix; i++)
		for (int j=0;j<ELEMENT_COUNT_LINE;j++)
			for (int k=0;k<ELEMENT_LENGTH_LINE;k++)
				m_pMat[i][j][k] = 0;//(j*ELEMENT_COUNT_LINE + k)/(float)(ELEMENT_LENGTH_LINE*ELEMENT_COUNT_LINE) ;

	for (int i=0;i<SIZE_HEIGHT;i++)
		for (int j=0;j<SIZE_WIDTH;j++)
			m_pIndex[i][j] = j%m_nSizeMatrix;
#else
	m_pIn = (float*)_aligned_malloc( SIZE_HEIGHT*SIZE_WIDTH*ELEMENT_COUNT_POINT*sizeof(float), 16 );//new float[];
	m_pOut = (float*)_aligned_malloc( SIZE_HEIGHT*SIZE_WIDTH*ELEMENT_COUNT_POINT*sizeof(float), 16 );//new float[];
	m_pMatrix = (float*)_aligned_malloc( ELEMENT_COUNT_LINE*ELEMENT_LENGTH_LINE*sizeof(float), 16 );//new float[];
	memset( m_pMatrix, 0, ELEMENT_COUNT_LINE*ELEMENT_LENGTH_LINE*sizeof(float));

	for (int i=0;i<SIZE_HEIGHT;i++)
		for (int j=0;j<SIZE_WIDTH;j++)
		{
			*(m_pIn + (i*SIZE_WIDTH+j)*2+0) = i;
			*(m_pIn + (i*SIZE_WIDTH+j)*2+1) = j;

			*(m_pOut + (i*SIZE_WIDTH+j)*2+0) = i;
			*(m_pOut + (i*SIZE_WIDTH+j)*2+1) = j;
		}

#endif

	
	

	float offset = 0.0f;
	float theta = ANGLE_THETA;
	float scale = SCALE;
	float mr[ELEMENT_LENGTH_LINE][ELEMENT_LENGTH_LINE] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} } ; // ��ת�任�����ʼ��
	float ms[ELEMENT_LENGTH_LINE][ELEMENT_LENGTH_LINE] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} } ; // ���ű任�����ʼ��
	float (*m)[ELEMENT_LENGTH_LINE] = m_pMat[0] ;	 // ͳһ�任�����ʼ��

	// ��ת����
	float rad = ANGLE_THETA_RAD ;
	mr[0][0] = cos( rad ) ;
	mr[0][1] = sin( rad ) ;
	mr[1][0] = -sin( rad ) ;
	mr[1][1] = cos( rad ) ;

	// ���ž��� ƽ�ƾ���
	ms[0][0] = scale ;
	ms[1][1] = scale ;

	ms[2][0] = offset;
	ms[2][1] = -offset;

	// ���ͳһ�任����
	for( int i = 0 ; i < ELEMENT_LENGTH_LINE ; i++ )
		for( int j = 0 ; j < ELEMENT_LENGTH_LINE ; j++ )
		{
			for( int k = 0 ; k < ELEMENT_LENGTH_LINE ; k++ )
			{
#if !OPTIMIZE_SSE
				m[i][j] += ms[i][k] * mr[k][j];
#else
				*(m_pMatrix + i*ELEMENT_LENGTH_LINE + j) += ms[i][k] * mr[k][j];
#endif
			}
		}

}

void CMatrixMultPoint2DWeightOCL::UnInit()
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

void CMatrixMultPoint2DWeightOCL::Implement( bool bMulti )
{
#if  OPTIMIZE_OPENCL
	ExecuteKernel();
#else
	m_bMulti = bMulti;

	if( m_bMulti )
		mmpParallel( );
	else
		mmpSerial( );
#endif
}

CMatrixMultPoint2DWeightOCL::CMatrixMultPoint2DWeightOCL()
{
	m_imgIn = m_imgOut = m_pOutRef = NULL;
	m_pMat = NULL;
	m_pIndex = NULL;
	m_bMulti = false;
}

CMatrixMultPoint2DWeightOCL::~CMatrixMultPoint2DWeightOCL()
{
	UnInit();
}
#define ABS(a, b)  ((a)>(b)?(a-b):(b-a))  
bool CMatrixMultPoint2DWeightOCL::verify()
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

void CMatrixMultPoint2DWeightOCL::mmpRef()
{
	//kernel(m_imgIn, m_pOutRef, m_pMat[0], m_pIndex);

}

void CMatrixMultPoint2DWeightOCL::kernel( float imgIn[][SIZE_WIDTH][ELEMENT_COUNT_POINT], float imgOut[][SIZE_WIDTH][ELEMENT_COUNT_POINT], float m[][ELEMENT_COUNT_LINE], int i, float rad )
{
#if OPTIMIZE_OMP && OPTIMIZE_INNER_LOOP
#pragma omp parallel for //shared(i)
#endif
		for (int j=0;j<SIZE_WIDTH ;j++)
		{
			float x = imgIn[i][j][0];
			float y = imgIn[i][j][1];
#if  CENTER_ROTATE
			float offset = 4.0f;
			float distance = sqrt( (x-SIZE_WIDTH/2.0f)*(x-SIZE_WIDTH/2.0f) + (y-SIZE_HEIGHT/2.0f)*(y-SIZE_HEIGHT/2.0f) + offset) ;
			float weight = 1.0f / sqrt(distance)  ;

			float radOne = rad * weight;

			float mOne[ELEMENT_LENGTH_LINE][ELEMENT_LENGTH_LINE] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} } ; // ��ϱ任�����ʼ��
			float mt[ELEMENT_LENGTH_LINE][ELEMENT_LENGTH_LINE] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} } ; // �ֲ��任�����ʼ�� ƽ��
			float mr[ELEMENT_LENGTH_LINE][ELEMENT_LENGTH_LINE] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} } ; // �ֲ��任�����ʼ�� ��ת
			float ms[ELEMENT_LENGTH_LINE][ELEMENT_LENGTH_LINE] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} } ; // �ֲ��任�����ʼ�� ����
			// ����ƽ�Ƶ�ԭ��
			mt[2][0] = SIZE_WIDTH/2.0f ;
			mt[2][1] = SIZE_HEIGHT/2.0f ;
			matrixMultiply( mt, mOne, mOne );

			mr[0][0] = cos( radOne ) ;
			mr[0][1] = sin( radOne ) ;
			mr[1][0] = -sin( radOne ) ;
			mr[1][1] = cos( radOne ) ;
			matrixMultiply( mr, mOne, mOne );
			
			// ���ž��� ƽ�ƾ���
			ms[0][0] = SCALE ;
			ms[1][1] = SCALE ;
			matrixMultiply( ms, mOne, mOne );
			
			// ƽ�ƻ�ԭ����
			mt[2][0] = -SIZE_WIDTH/2.0f ;
			mt[2][1] = -SIZE_HEIGHT/2.0f ;
			matrixMultiply( mt, mOne, mOne );

			imgOut[i][j][0] = x * mOne[0][0] + y * mOne[1][0] + mOne[2][0];
			imgOut[i][j][1] = x * mOne[0][1] + y * mOne[1][1] + mOne[2][1];
#else
			float offset = 4.0f;
			float distance = sqrt( x*x + y*y + offset) ;
			float weight = 1.0f / sqrt(distance)  ;

			float radOne = rad * weight;
			float mOne [3][3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} } ; // ��ϱ任�����ʼ��

			// ��ת ����
			mOne[0][0] = cos( radOne ) * SCALE;
			mOne[0][1] = sin( radOne ) * SCALE ;
			mOne[1][0] = -sin( radOne ) * SCALE ;
			mOne[1][1] = cos( radOne ) * SCALE ;

#endif
			imgOut[i][j][0] = x * mOne[0][0] + y * mOne[1][0] + mOne[2][0];
			imgOut[i][j][1] = x * mOne[0][1] + y * mOne[1][1] + mOne[2][1];
		}

}


void CMatrixMultPoint2DWeightOCL::kernelElement( float* pIn, float* pOut, float* pMat )
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

float* CMatrixMultPoint2DWeightOCL::getOutput()
{	
#if !OPTIMIZE_SSE
	return &m_imgOut[0][0][0]; 
#else
	return m_pOut;
#endif
}

void CMatrixMultPoint2DWeightOCL::accumulate()
{
	//kernel(m_imgIn, m_imgOut, m_pMat[0], m_pIndex);
	Implement(m_bMulti);
#if OPTIMIZE_SSE
	memcpy( m_pIn, m_pOut, sizeof(float) * SIZE_HEIGHT * SIZE_WIDTH * ELEMENT_COUNT_POINT );
#else
	memcpy( m_imgIn, m_imgOut, sizeof(float) * SIZE_HEIGHT * SIZE_WIDTH * ELEMENT_COUNT_POINT );
#endif
}


//#define _MM_SHUFFLE(fp3,fp2,fp1,fp0) (((fp3) << 6) | ((fp2) << 4) | \
//	((fp1) << 2) | ((fp0)))


#define __MM_SELECT(v, fp)                                                          \
	_mm_shuffle_ps((v), (v), _MM_SHUFFLE((fp),(fp),(fp),(fp)))

void CMatrixMultPoint2DWeightOCL::kernelSSE( float* imgIn, float* imgOut, int i, float rad )
{
#if 0
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
#else
	__m128 *pSrcPos = (__m128*)imgIn;
	__m128 *pDestPos = (__m128*)imgOut;
	int nOffset = i*SIZE_WIDTH/2;

#if OPTIMIZE_OMP && OPTIMIZE_INNER_LOOP
#pragma omp parallel for //shared(nOffset)
#endif
	for (int j=0;j<SIZE_WIDTH/2 ;j++)
	{
		// (x1,x1,x2,x2)  (y1,y1,y2,y2)
		__m128 vIn = pSrcPos[j+nOffset];

		__m128 vInx = _mm_shuffle_ps( vIn, vIn, _MM_SHUFFLE(2,2,0,0) );
		__m128 vIny = _mm_shuffle_ps( vIn, vIn, _MM_SHUFFLE(3,3,1,1) );

		__m128 xySqure = _mm_add_ps( _mm_mul_ps( vInx, vInx ) , _mm_mul_ps( vIny, vIny ) );
		xySqure = _mm_add_ps( xySqure, _mm_set_ps1(4.0f) );
		__m128 weight = _mm_sqrt_ps( _mm_sqrt_ps(xySqure) );
		weight = _mm_rcp_ps(weight);
		__m128 radWeight = _mm_load_ps1(&rad);
		radWeight = _mm_mul_ps( radWeight, weight);

		 __declspec(align(16)) float radLast[4];
		_mm_store_ps( radLast, radWeight );

		 __declspec(align(16)) float mOne[2][4];

		mOne[0][0] = cos( radLast[0] ) ;
		mOne[0][1] = sin( radLast[0] ) ;
		mOne[1][0] = -mOne[0][1] ;
		mOne[1][1] = mOne[0][0] ;

		mOne[0][2] = cos( radLast[2] ) ;
		mOne[0][3] = sin( radLast[2] ) ;
		mOne[1][2] = -mOne[0][3] ;
		mOne[1][3] = mOne[0][2] ;

		__m128 mat[2];
		mat[0] = _mm_load_ps( mOne[0] ) ;
		mat[1] = _mm_load_ps( mOne[1] ) ;

		__m128 tmp00 = _mm_mul_ps( mat[0], vInx );
		__m128 tmp01 = _mm_mul_ps( mat[1], vIny );

		__m128 vOut0 = _mm_add_ps( tmp00, tmp01 );
		vOut0 = _mm_mul_ps( vOut0, _mm_set_ps1( SCALE ) );

		pDestPos[j+nOffset] = vOut0 ;


	}
#endif

}

void CMatrixMultPoint2DWeightOCL::matrixMultiply( float mLeft[][ELEMENT_LENGTH_LINE], float mRight[][ELEMENT_LENGTH_LINE], float mResult[][ELEMENT_LENGTH_LINE] )
{
	double mMix[ELEMENT_LENGTH_LINE][ELEMENT_LENGTH_LINE] = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} } ; // ��ϱ任�����ʼ��
	for( int i = 0 ; i < ELEMENT_LENGTH_LINE ; i++ )
		for( int j = 0 ; j < ELEMENT_LENGTH_LINE ; j++ )
		{
			for( int k = 0 ; k < ELEMENT_LENGTH_LINE ; k++ )
			{
				mMix[i][j] += mLeft[i][k] * mRight[k][j];
			}
		}

		for( int i = 0 ; i < ELEMENT_LENGTH_LINE ; i++ )
			for( int j = 0 ; j < ELEMENT_LENGTH_LINE ; j++ )
				mResult[i][j] = mMix[i][j];
}

void CMatrixMultPoint2DWeightOCL::SetupKernel()
{
	const cl_mem_flags INFlags  = CL_MEM_COPY_HOST_PTR | CL_MEM_READ_ONLY; 
	const cl_mem_flags OUTFlags = CL_MEM_COPY_HOST_PTR | CL_MEM_READ_WRITE;

	m_argOCL.m_pfInputBuffer = clCreateBuffer( m_pEnvOpenCL->g_context, INFlags, 
		 SIZE_HEIGHT*SIZE_WIDTH*ELEMENT_COUNT_POINT*sizeof(float),	m_pIn,	NULL); 

	m_argOCL.m_pfOCLOutputBuffer = clCreateBuffer( m_pEnvOpenCL->g_context, INFlags, 
		SIZE_HEIGHT*SIZE_WIDTH*ELEMENT_COUNT_POINT*sizeof(float),	m_pOut,	NULL); 

	// thread size
	m_argOCL.localWorkSize[0] = LocalWorkX;
	m_argOCL.localWorkSize[1] = LocalWorkX;

	size_t  workGroupSizeMaximum;
	clGetKernelWorkGroupInfo( m_pEnvOpenCL->g_kernel, m_pEnvOpenCL->g_device_ID, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), (void *)&workGroupSizeMaximum, NULL);

	int nElementSizePadding = roundToPowerOf2( SIZE_HEIGHT*SIZE_WIDTH/2 );
	if ( nElementSizePadding > workGroupSizeMaximum )
	{
		m_argOCL.globalWorkSize[0] = workGroupSizeMaximum;
		m_argOCL.globalWorkSize[1] = nElementSizePadding / workGroupSizeMaximum;
	}
}

bool CMatrixMultPoint2DWeightOCL::ExecuteKernel()
{
	cl_int   err;

	err = clEnqueueWriteBuffer( m_pEnvOpenCL->g_cmd_queue, m_argOCL.m_pfInputBuffer, CL_TRUE, 0, 
		SIZE_HEIGHT*SIZE_WIDTH*ELEMENT_COUNT_POINT*sizeof(float) , m_pIn, 0, NULL, NULL);

	if (err != CL_SUCCESS)
	{
		printf("ERROR: Failed to clEnqueueReadBuffer...\n");
		return false;
	}

	//Set kernel arguments
	clSetKernelArg( m_pEnvOpenCL->g_kernel, 0, sizeof(cl_mem), (void *) &m_argOCL.m_pfInputBuffer );
	clSetKernelArg( m_pEnvOpenCL->g_kernel, 1, sizeof(cl_mem), (void *) &m_argOCL.m_pfOCLOutputBuffer );

	float rad = ANGLE_THETA_RAD ;
	clSetKernelArg( m_pEnvOpenCL->g_kernel, 2, sizeof(float), &rad);

	float scale = SCALE ;
	clSetKernelArg( m_pEnvOpenCL->g_kernel, 3, sizeof(float), &scale);
	
	int sizeMax = SIZE_HEIGHT*SIZE_WIDTH/2 ;
	clSetKernelArg( m_pEnvOpenCL->g_kernel, 4, sizeof(int), &sizeMax);



	//execute opencl
	cl_event g_perf_event = NULL;
#if LocalWorkSizeDef
	err= clEnqueueNDRangeKernel( m_pEnvOpenCL->g_cmd_queue, m_pEnvOpenCL->g_kernel, 2, NULL, m_argOCL.globalWorkSize, NULL, 0, NULL, &g_perf_event);
#else
	err= clEnqueueNDRangeKernel( m_pEnvOpenCL->g_cmd_queue, m_pEnvOpenCL->g_kernel, 2, NULL, m_argOCL.globalWorkSize, m_argOCL.localWorkSize, 0, NULL, &g_perf_event);
#endif
	if (err != CL_SUCCESS)
	{
		printf("ERROR: Failed to execute kernel...\n");
		return false;
	}
	err = clWaitForEvents(1, &g_perf_event);
	if (err != CL_SUCCESS)
	{
		printf("ERROR: Failed to clWaitForEvents...\n");
		return false;
	}

	void* tmp_ptr = NULL;
	err = clEnqueueReadBuffer( m_pEnvOpenCL->g_cmd_queue, m_argOCL.m_pfOCLOutputBuffer, CL_TRUE, 0, 
		SIZE_HEIGHT*SIZE_WIDTH*ELEMENT_COUNT_POINT*sizeof(float) , m_pOut, 0, NULL, NULL);

	if (err != CL_SUCCESS)
	{
		printf("ERROR: Failed to clEnqueueReadBuffer...\n");
		return false;
	}

	clFinish( m_pEnvOpenCL->g_cmd_queue );

	return true;
}

/*

   DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS HEADER.

   Copyright 2009 Sun Microsystems, Inc. All rights reserved.

   The contents of this file are subject to the terms of the BSD License("BSD")(the "License"). 
   You can obtain a copy of the License at: http://www.opensparc.net/pubs/t1/licenses/BSD+_License.txt

   The BSD License

   Redistribution and use in source and binary forms, with or without 
   modification, are permitted provided that the following conditions are met:

       * Redistribution of source code must retain the above copyright 
         notice, this list of conditions and the following disclaimer.
       * Redistribution in binary form must reproduce the above copyright 
         notice, this list of conditions and the following disclaimer in 
         the documentation and/or other materials provided with the 
         distribution.
       * Neither the name of Sun Microsystems, Inc. or the names of 
         contributors may be used to endorse or promote products derived 
         from this software without specific prior written permission.

   This software is provided "AS IS," without a warranty of any kind. ALL 
   EXPRESS OR IMPLIED CONDITIONS, REPRESENTATIONS AND WARRANTIES, INCLUDING ANY 
   IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE OR 
   NON-INFRINGEMENT, ARE HEREBY EXCLUDED. SUN MICROSYSTEMS, INC. ("SUN") AND 
   ITS LICENSORS SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY LICENSEE AS A 
   RESULT OF USING, MODIFYING OR DISTRIBUTING THIS SOFTWARE OR ITS DERIVATIVES. 
   IN NO EVENT WILL SUN OR ITS LICENSORS BE LIABLE FOR ANY LOST REVENUE, PROFIT 
   OR DATA, OR FOR DIRECT, INDIRECT, SPECIAL, CONSEQUENTIAL, INCIDENTAL OR 
   PUNITIVE DAMAGES, HOWEVER CAUSED AND REGARDLESS OF THE THEORY OF LIABILITY, 
   ARISING OUT OF THE USE OF OR INABILITY TO USE THIS SOFTWARE, EVEN IF SUN HAS 
   BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.

   You acknowledge that this software is not designed, licensed or intended for 
   use in the design, construction, operation or maintenance of any nuclear facility. 

*/
#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>

#include "app-MatrixMultVector.h"

#define SIZE_MATRIX_M		1000
#define SIZE_MATRIX_N		1000

/*
   a[m][1]
   =
   b[m][n] *  c[n][1]
*/
int CMatrixMultVector::mxv( )
{

	mxvInit();

	mxvImplement();

	mxvUnInit();

   return(0);
}

void CMatrixMultVector::mxvSerial(int m, int n, double *  a, double *  b,
         double *  c)
{
   int i, j;

   for (i=0; i<m; i++)
   {
      a[i] = 0.0;
      for (j=0; j<n; j++)
         a[i] += b[i*n+j]*c[j];
   }
}


void CMatrixMultVector::mxvParallel(int m, int n, double *  a, double *  b,
	double *  c)
{
	int i, j;

#pragma omp parallel for default(none) \
	shared(m,n,a,b,c) private(i,j)
	for (i=0; i<m; i++)
	{
		a[i] = 0.0;
		for (j=0; j<n; j++)
			a[i] += b[i*n+j]*c[j];
	} /*-- End of omp parallel for --*/
}

void CMatrixMultVector::mxvInit()
{
	m = SIZE_MATRIX_M;
	n = SIZE_MATRIX_N;

	if ( (a=(double *)malloc(m*sizeof(double))) == NULL )
		perror("memory allocation for a");
	if ( (b=(double *)malloc(m*n*sizeof(double))) == NULL )
		perror("memory allocation for b");
	if ( (c=(double *)malloc(n*sizeof(double))) == NULL )
		perror("memory allocation for c");
	int i,j;
	for (j=0; j<n; j++)
		c[j] = 2.0;
	for (i=0; i<m; i++)
		for (j=0; j<n; j++)
			b[i*n+j] = i;

}

void CMatrixMultVector::mxvUnInit()
{
	if (a)	{ free(a); a=NULL; }
	if (b)	{ free(b); b=NULL; }
	if (c)	{ free(c); c=NULL; }
}

void CMatrixMultVector::mxvImplement()
{
#if ENABLE_PARALLEL
	(void) mxvParallel(m, n, a, b, c);
#else
	(void) mxvSerial(m, n, a, b, c);
#endif
}

CMatrixMultVector::CMatrixMultVector()
{
	a = b = c = NULL;
}

CMatrixMultVector::~CMatrixMultVector()
{
	mxvUnInit();
}

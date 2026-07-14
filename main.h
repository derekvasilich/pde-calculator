/*
 *  main.h
 *  Assignment Five
 *
 *  Created by Derek Williams on 10-01-30.
 *  Copyright 2010 Derek Williams. All rights reserved.
 *
 */

#define SM_FORWARD_SPACE		0
#define SM_BACKWARD_SPACE		1
#define SM_CENTRAL_SPACE		2
#define	SM_LAX_FREDRICHS		3
#define SM_LEAPFROG				4
#define SM_EQUILIBRIUM			5
#define SM_LAX_WENDROFF			6
#define SM_CRANK_NICOLSON		7
#define SM_SOR_METHOD			8
#define SM_CG_METHOD			9
#define NUM_SCHM				10

#define SM_FORWARD_SPACE_NAME	"Forward-Time Forward-Space"
#define SM_BACKWARD_SPACE_NAME	"Forward-Time Backward-Space"
#define SM_CENTRAL_SPACE_NAME	"Forward-Time Central-Space"
#define	SM_LAX_FREDRICHS_NAME	"Lax Fredrichs"
#define SM_LEAPFROG_NAME		"Leapfrog"
#define SM_EQUILIBRIUM_NAME		"Equilibrium"
#define SM_LAX_WENDROFF_NAME	"Lax Wendroff"
#define SM_CRANK_NICOLSON_NAME	"Crank Nicolson (with ADI)"
#define SM_SOR_METHOD_NAME		"Successive Overrelaxation"
#define SM_CG_METHOD_NAME		"Conjugate Gradient"

#define A_CONST					1000
#define A_SPACE					1001

#define A_CONST_NAME			"constant"
#define A_SPACE_NAME			"1 + 1/4(3 - x)(1 + x)"

#define INIT_ONE				0
#define INIT_TWO				1
#define INIT_THREE				2
#define INIT_FOUR				3
#define INIT_FIVE				4
#define INIT_SIX				5
#define INIT_SEVEN				6
#define INIT_EIGHT				7
#define INIT_NINE				8
#define NUM_INIT				9

#define H_10	10
#define H_20	20
#define H_30	30
#define H_40	40
#define H_80	80

#define LAMBDA_8				100
#define LAMBDA_10				101
#define LAMBDA_16				102

#define INIT_ONE_NAME			"cos^2(π x), |x| <= 1/2"
#define INIT_TWO_NAME			"1 - |x|, |x| <= 1"
#define INIT_THREE_NAME			"sin(2π x)"
#define INIT_FOUR_NAME			"cos(2π x)"
#define INIT_FIVE_NAME			"cos(ξπ x)cos^2(π/2 x)"
#define INIT_SIX_NAME			"1 - |x| if |x| < 1/2, 1/4 if |x| = 1/2"
#define INIT_SEVEN_NAME			"sin(1.2(x - y))cosh(x + 2y)"
#define INIT_EIGHT_NAME			"exp(-(x^2 + y^2)"
#define INIT_NINE_NAME			"cos(x)sin(y)"

#define MAX_VM_N				3000
#define SAVE_ID					3001

#define F_ZERO					4000
#define F_ONE					4001

#define F_ZERO_NAME				"F(t, x) = 0"
#define F_ONE_NAME				"F(t, x) = sin(π (x-t))"

#define SETTINGS_ID				5000

#define BOUND_ONE				0
#define BOUND_TWO				1
#define BOUND_THREE				2
#define BOUND_FOUR				3
#define NUM_BOUND				4

#define BOUND_ONE_NAME			"v(n+1,M) = v(n,M-1)"
#define BOUND_TWO_NAME			"v(n+1,M) - v(n,M) + λ(v(n+1,M) - v(n+1,M-1)) = kv(n+1,M)"
#define BOUND_THREE_NAME		"v(n+1,M) = 2v(n+1,M-1) - v(n+1,M-2)"
#define BOUND_FOUR_NAME			"Dirichlet and Neumann"

#define TIME_START				0
	
#define Y_CONST 1
#define X_CONST 2

#define BLOWUP_VAL				5.0

#define ALIGN_CENTER			0
#define ALIGN_RIGHT				1
#define ALIGN_LEFT				-1

#define TINY					-1e100
	
#define BUFLEN					1024

#define CONVERG					10e-7
#define	SUM_MAX					1e6

#define CTRL_SCHEME				1
#define CTRL_H					2
#define CTRL_TIME				3
#define CTRL_UTERM				4
#define CTRL_LAMBDA				5
#define CTRL_XI					6
#define CTRL_IC					7
#define CTRL_TOP				8
#define CTRL_BOTTOM				9
#define CTRL_RIGHT				10
#define CTRL_LEFT				11
#define CTRL_REPEAT				12
#define CTRL_BC					13
#define CTRL_SCALE				14
#define CTRL_MU					15
#define CTRL_A					16
#define CTRL_B					17
#define CTRL_DIV				18
#define CTRL_SAVE				19
#define CTRL_LOAD				20

#define CF_SCHEME				0
#define CF_LAMBDA				1
#define CF_MU					2
#define CF_TOP					3
#define CF_BOTTOM				4
#define CF_LEFT					5
#define CF_RIGHT				6
#define CF_BOUNDARY				7
#define CF_TIME					8
#define CF_THREEDIM				9
#define CF_WIRE					10
#define CF_ACONST				11
#define CF_BCONST				12
#define CF_A_FUNC				13
#define CF_B_FUNC				14
#define CF_INIT_FUNC			15
#define CF_BOUNDS_FUNC			16
#define CF_H					17
#define CF_PBC					18
#define CF_RHO					19
#define CF_THETA				20
#define CF_PHI					21
#define CF_R					22
#define CF_SIGMA				23
#define CF_K					24

/* handy vector operations */

#define VSUB3(A, B, C)				\
	A[0] = B[0] - C[0];				\
	A[1] = B[1] - C[1];				\
	A[2] = B[2] - C[2];

#define VCROSS3(A, B, C)			\
	A[0] =   B[1]*C[2] - B[2]*C[1];	\
	A[1] = - B[0]*C[2] + B[2]*C[0];	\
	A[2] =   B[0]*C[1] - B[1]*C[0];

#define VEQ3(A, B)					\
	A[0] = B[0];					\
	A[1] = B[1];					\
	A[2] = B[2];

#define VSET2(A, X, Y)				\
	A[0] = X;						\
	A[1] = Y;

#define VSET3(A, X, Y, Z)			\
	A[0] = X;						\
	A[1] = Y;						\
	A[2] = Z;

/* default values */

#define DEFAULT_A				&aConst
#define DEFAULT_B				&bConst

#define DEFAULT_LAMBDA			0.8

#define BOUNDS_BOTTOM			0
#define BOUNDS_TOP				1

#define	DEFAULT_THREEDIM		1

#define TIME_END				1

#define BOUNDS_LEFT				0
#define BOUNDS_RIGHT			1

#define DEFAULT_SCHEME			SM_CRANK_NICOLSON
#define DEFAULT_BOUNDARY		&BZero
#define DEFAULT_INIT			&initOne
#define DEFAULT_MU				0.0

#define DEFAULT_INIT_NUM		INIT_NINE
#define DEFAULT_BOUND			BOUND_ONE

//#define	DEFAULT_CONFIG			"sor.cfg"
#define	DEFAULT_CONFIG			"crank3d.cfg"
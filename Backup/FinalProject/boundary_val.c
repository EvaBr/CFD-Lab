#include "init.h"
#include "boundary_val.h"
#include "helper.h"

void boundaryvalues_no_slip(
			int i,
			int j,
			int k,
			double ***U,
			double ***V,
			double ***W,
			int ***Flag
			){
        // No-slip boundary conditions for U, V and W.
	// 26 (6 + 12 + 8) cases in total.
	//printf( "cell type:  %d \n", getcelltype(Flag[i  ][j  ][k  ]));
	//int flags = (~(2730&Flag[i  ][j  ][k  ]))&getwallbit(0);
	//printf("flags: %d\n", flags);
	//int sth = ((flags&2) >> 1);
	//printf("sth: %d\n", sth);
	//sth += ((flags >> 2)&2) + ((flags >> 3)&4) + ((flags >> 4)&8) + ((flags >> 5)&16) + ((flags >> 6)&32);

	switch(getboundarytype(Flag[i  ][j  ][k  ])){
		case B_O:
			U[i  ][j  ][k  ]   = 0.0;
			V[i  ][j-1][k  ] = -V[i+1][j-1][k ];
			V[i  ][j  ][k  ] = -V[i+1][j  ][k ];
			W[i  ][j  ][k-1] = -W[i+1][j  ][k-1];
			W[i  ][j  ][k  ] = -W[i+1][j  ][k ];
			break;
		case B_W:
			U[i-1][j  ][k  ] = 0.0;
			V[i  ][j-1][k  ] = -V[i-1][j-1][k  ];
			V[i  ][j  ][k  ] = -V[i-1][j  ][k  ];
			W[i  ][j  ][k-1] = -W[i-1][j  ][k-1];
			W[i  ][j  ][k  ] = -W[i-1][j  ][k  ];
			break;
		case B_N:
			V[i  ][j  ][k  ] = 0.0;
			U[i-1][j  ][k  ] = -U[i-1][j+1][k  ];
			U[i  ][j  ][k  ] = -U[i  ][j+1][k  ];
			W[i  ][j  ][k-1] = -W[i  ][j+1][k-1];
			W[i  ][j  ][k  ] = -W[i  ][j+1][k  ];
			break;
		case B_S:
			V[i  ][j-1][k  ] = 0.0;
			U[i-1][j  ][k  ] = -U[i-1][j-1][k  ];
			U[i  ][j  ][k  ] = -U[i  ][j-1][k  ];
			W[i  ][j  ][k-1] = -W[i  ][j-1][k-1];
			W[i  ][j  ][k  ] = -W[i  ][j-1][k  ];
			break;
		case B_U:
			W[i  ][j  ][k  ] = 0.0;
			U[i-1][j  ][k  ] = -U[i-1][j  ][k+1];
			U[i  ][j  ][k  ] = -U[i  ][j  ][k+1];
			V[i  ][j-1][k  ] = -V[i  ][j-1][k+1];
			V[i  ][j  ][k  ] = -V[i  ][j  ][k+1];
			break;
		case B_D:
			W[i  ][j  ][k-1] = 0.0;
			U[i-1][j  ][k  ] = -U[i-1][j  ][k-1];
			U[i  ][j  ][k  ] = -U[i  ][j  ][k-1];
			V[i  ][j-1][k  ] = -V[i  ][j-1][k-1];
			V[i  ][j  ][k  ] = -V[i  ][j  ][k-1];
			break;

		case B_NO:
			U[i  ][j  ][k  ] = 0.0;
			U[i-1][j  ][k  ] = -U[i-1][j+1][k  ];
			V[i  ][j  ][k  ] = 0.0;
			V[i  ][j-1][k  ] =  -V[i+1][j-1][k  ];
			W[i  ][j  ][k  ] = -(W[i  ][j+1][k  ] + W[i+1][j  ][k  ]) * 0.5;
			W[i  ][j  ][k-1] = -(W[i  ][j+1][k-1] + W[i+1][j  ][k-1]) * 0.5;
			break;
		case B_NW:
			U[i-1][j  ][k  ] = 0.0;
			U[i  ][j  ][k  ] = -U[i  ][j+1][k  ];
			V[i  ][j  ][k  ] = 0.0;
			V[i  ][j-1][k  ] =  -V[i-1][j-1][k  ];
			W[i  ][j  ][k  ] = -(W[i  ][j+1][k  ] + W[i-1][j  ][k  ]) * 0.5;
			W[i  ][j  ][k-1] = -(W[i  ][j+1][k-1] + W[i-1][j  ][k-1]) * 0.5;
			break;
		case B_NU:
			V[i  ][j  ][k  ] = 0.0;
			V[i  ][j-1][k  ] = -V[i  ][j-1][k+1];
			W[i  ][j  ][k  ] = 0.0;
			W[i  ][j  ][k-1] =  -W[i  ][j+1][k-1];
			U[i  ][j  ][k  ] = -(U[i  ][j+1][k  ] + U[i  ][j  ][k+1]) * 0.5;
			U[i-1][j  ][k  ] = -(U[i-1][j+1][k  ] + U[i-1][j  ][k+1]) * 0.5;
			break;
		case B_ND:
			V[i  ][j  ][k  ] = 0.0;
			V[i  ][j-1][k  ] = -V[i  ][j-1][k-1];
			W[i  ][j  ][k-1] = 0.0;
			W[i  ][j  ][k  ] =  -W[i  ][j+1][k  ];
			U[i  ][j  ][k  ] = -(U[i  ][j+1][k  ] + U[i  ][j  ][k-1]) * 0.5;
			U[i-1][j  ][k  ] = -(U[i-1][j+1][k  ] + U[i-1][j  ][k-1]) * 0.5;
			break;
		case B_SO:
			U[i  ][j  ][k  ] = 0.0;
			U[i-1][j  ][k  ] = -U[i-1][j-1][k  ];
			V[i  ][j-1][k  ] = 0.0;
			V[i  ][j  ][k  ] =  -V[i+1][j  ][k  ];
			W[i  ][j  ][k  ] = -(W[i  ][j-1][k  ] + W[i+1][j  ][k  ]) * 0.5;
			W[i  ][j  ][k-1] = -(W[i  ][j-1][k-1] + W[i+1][j  ][k-1]) * 0.5;
			break;
		case B_SW:
			U[i-1][j  ][k  ] = 0.0;
			U[i  ][j  ][k  ] = -U[i  ][j-1][k  ];
			V[i  ][j-1][k  ] = 0.0;
			V[i  ][j  ][k  ] =  -V[i-1][j  ][k  ];
			W[i  ][j  ][k  ] = -(W[i  ][j-1][k  ] + W[i-1][j  ][k  ]) * 0.5;
			W[i  ][j  ][k-1] = -(W[i  ][j-1][k-1] + W[i-1][j  ][k-1]) * 0.5;
			break;
		case B_SU:
			V[i  ][j-1][k  ] = 0.0;
			V[i  ][j  ][k  ] = -V[i  ][j  ][k+1];
			W[i  ][j  ][k  ] = 0.0;
			W[i  ][j  ][k-1] =  -W[i  ][j-1][k-1];
			U[i  ][j  ][k  ] = -(U[i  ][j-1][k  ] + U[i  ][j  ][k+1]) * 0.5;
			U[i-1][j  ][k  ] = -(U[i-1][j-1][k  ] + U[i-1][j  ][k+1]) * 0.5;
			break;
		case B_SD:
			V[i  ][j-1][k  ] = 0.0;
			V[i  ][j  ][k  ] = -V[i  ][j  ][k-1];
			W[i  ][j  ][k-1] = 0.0;
			W[i  ][j  ][k  ] =  -W[i  ][j-1][k  ];
			U[i  ][j  ][k  ] = -(U[i  ][j-1][k  ] + U[i  ][j  ][k-1]) * 0.5;
			U[i-1][j  ][k  ] = -(U[i-1][j-1][k  ] + U[i-1][j  ][k-1]) * 0.5;
			break;
		case B_OU:
			U[i  ][j  ][k  ] = 0.0;
			U[i-1][j  ][k  ] = -U[i-1][j  ][k+1];
			W[i  ][j  ][k  ] = 0.0;
			W[i  ][j  ][k-1] =  -W[i+1][j  ][k-1];
			V[i  ][j  ][k  ] = -(V[i+1][j  ][k  ] + V[i  ][j  ][k+1]) * 0.5;
			V[i  ][j-1][k  ] = -(V[i+1][j-1][k  ] + V[i  ][j-1][k+1]) * 0.5;
			break;
		case B_WU:
			U[i-1][j  ][k  ] = 0.0;
			U[i  ][j  ][k  ] = -U[i  ][j  ][k+1];
			W[i  ][j  ][k  ] = 0.0;
			W[i  ][j  ][k-1] =  -W[i-1][j  ][k-1];
			V[i  ][j  ][k  ] = -(V[i-1][j  ][k  ] + V[i  ][j  ][k+1]) * 0.5;
			V[i  ][j-1][k  ] = -(V[i-1][j-1][k  ] + V[i  ][j-1][k+1]) * 0.5;
			break;
		case B_OD:
			U[i  ][j  ][k  ] = 0.0;
			U[i-1][j  ][k  ] = -U[i-1][j  ][k-1];
			W[i  ][j  ][k-1] = 0.0;
			W[i  ][j  ][k  ] = -W[i+1][j  ][k  ];
			V[i  ][j  ][k  ] = -(V[i+1][j  ][k  ] + V[i  ][j  ][k-1]) * 0.5;
			V[i  ][j-1][k  ] = -(V[i+1][j-1][k  ] + V[i  ][j-1][k-1]) * 0.5;
			break;
		case B_WD:
			U[i-1][j  ][k  ] = 0.0;
			U[i  ][j  ][k  ] = -U[i  ][j  ][k-1];
			W[i  ][j  ][k-1] = 0.0;
			W[i  ][j  ][k  ] = -W[i-1][j  ][k  ];
			V[i  ][j  ][k  ] = -(V[i-1][j  ][k  ] + V[i  ][j  ][k-1]) * 0.5;
			V[i  ][j-1][k  ] = -(V[i-1][j-1][k  ] + V[i  ][j-1][k-1]) * 0.5;
			break;
		case B_NOU:
			U[i  ][j  ][k  ] = 0.0;
			U[i-1][j  ][k  ] = -(U[i-1][j+1][k  ] + U[i-1][j  ][k+1]) * 0.5;
			V[i  ][j  ][k  ] = 0.0;
			V[i  ][j-1][k  ] = -(V[i+1][j-1][k  ] + V[i  ][j-1][k+1]) * 0.5;
			W[i  ][j  ][k  ] = 0.0;
			W[i  ][j  ][k-1] = -(W[i+1][j  ][k-1] + W[i  ][j+1][k-1]) * 0.5;
			break;
		case B_NWU:
			U[i-1][j  ][k  ] = 0.0;
			U[i  ][j  ][k  ] = -(U[i  ][j+1][k  ] + U[i  ][j  ][k+1]) * 0.5;
			V[i  ][j  ][k  ] = 0.0;
			V[i  ][j-1][k  ] = -(V[i-1][j-1][k  ] + V[i  ][j-1][k+1]) * 0.5;
			W[i  ][j  ][k  ] = 0.0;
			W[i  ][j  ][k-1] = -(W[i-1][j  ][k-1] + W[i  ][j+1][k-1]) * 0.5;
			break;
		case B_NOD:
			U[i  ][j  ][k  ] = 0.0;
			U[i-1][j  ][k  ] = -(U[i-1][j+1][k  ] + U[i-1][j  ][k-1]) * 0.5;
			V[i  ][j  ][k  ] = 0.0;
			V[i  ][j-1][k  ] = -(V[i+1][j-1][k  ] + V[i  ][j-1][k-1]) * 0.5;
			W[i  ][j  ][k-1] = 0.0;
			W[i  ][j  ][k  ] = -(W[i+1][j  ][k  ] + W[i  ][j+1][k  ]) * 0.5;
			break;
		case B_NWD:
			U[i-1][j  ][k  ] = 0.0;
			U[i  ][j  ][k  ] = -(U[i  ][j+1][k  ] + U[i  ][j  ][k-1]) * 0.5;
			V[i  ][j  ][k  ] = 0.0;
			V[i  ][j-1][k  ] = -(V[i-1][j-1][k  ] + V[i  ][j-1][k-1]) * 0.5;
			W[i  ][j  ][k-1] = 0.0;
			W[i  ][j  ][k  ] = -(W[i-1][j  ][k  ] + W[i  ][j+1][k  ]) * 0.5;
			break;
		case B_SOU:
			U[i  ][j  ][k  ] = 0.0;
			U[i-1][j  ][k  ] = -(U[i-1][j-1][k  ] + U[i-1][j  ][k+1]) * 0.5;
			V[i  ][j-1][k  ] = 0.0;
			V[i  ][j  ][k  ] = -(V[i+1][j  ][k  ] + V[i  ][j  ][k+1]) * 0.5;
			W[i  ][j  ][k  ] = 0.0;
			W[i  ][j  ][k-1] = -(W[i+1][j  ][k-1] + W[i  ][j-1][k-1]) * 0.5;
			break;
		case B_SWU:
			U[i-1][j  ][k  ] = 0.0;
			U[i  ][j  ][k  ] = -(U[i  ][j-1][k  ] + U[i  ][j  ][k+1]) * 0.5;
			V[i  ][j-1][k  ] = 0.0;
			V[i  ][j  ][k  ] = -(V[i-1][j  ][k  ] + V[i  ][j  ][k+1]) * 0.5;
			W[i  ][j  ][k  ] = 0.0;
			W[i  ][j  ][k-1] = -(W[i-1][j  ][k-1] + W[i  ][j-1][k-1]) * 0.5;
			break;
		case B_SOD:
			U[i  ][j  ][k  ] = 0.0;
			U[i-1][j  ][k  ] = -(U[i-1][j-1][k  ] + U[i-1][j  ][k-1]) * 0.5;
			V[i  ][j-1][k  ] = 0.0;
			V[i  ][j  ][k  ] = -(V[i+1][j  ][k  ] + V[i  ][j  ][k-1]) * 0.5;
			W[i  ][j  ][k-1] = 0.0;
			W[i  ][j  ][k  ] = -(W[i+1][j  ][k  ] + W[i  ][j-1][k  ]) * 0.5;
			break;
		case B_SWD:
			U[i-1][j  ][k  ] = 0.0;
			U[i  ][j  ][k  ] = -(U[i  ][j-1][k  ] + U[i  ][j  ][k-1]) * 0.5;
			V[i  ][j-1][k  ] = 0.0;
			V[i  ][j  ][k  ] = -(V[i-1][j  ][k  ] + V[i  ][j  ][k-1]) * 0.5;
			W[i  ][j  ][k-1] = 0.0;
			W[i  ][j  ][k  ] = -(W[i-1][j  ][k  ] + W[i  ][j-1][k  ]) * 0.5;
			break;
		default: //case 0
			//printf("case 0\n");
			break;
	}

}

void boundaryvalues_free_slip(
			int i,
			int j,
			int k,
			double ***U,
			double ***V,
			double ***W,
			int ***Flag
			){
        // Free slip boundary conditions for U, V and W.
	// 26 (6 + 12 + 8) cases in total.
	switch(getboundarytype(Flag[i  ][j  ][k  ])){
		case B_O:
			U[i  ][j  ][k  ] = 0.0;
			V[i  ][j-1][k  ] = V[i+1][j-1][k ];
			V[i  ][j  ][k  ] = V[i+1][j  ][k ];
			W[i  ][j  ][k-1] = W[i+1][j  ][k-1];
			W[i  ][j  ][k  ] = W[i+1][j  ][k ];
			break;
		case B_W:
			U[i-1][j  ][k  ] = 0.0;
			V[i  ][j-1][k  ] = V[i-1][j-1][k  ];
			V[i  ][j  ][k  ] = V[i-1][j  ][k  ];
			W[i  ][j  ][k-1] = W[i-1][j  ][k-1];
			W[i  ][j  ][k  ] = W[i-1][j  ][k  ];
			break;
		case B_N:
			V[i  ][j  ][k  ]   = 0.0;
			U[i-1][j  ][k  ]   = U[i-1][j+1][k  ];
			U[i  ][j  ][k  ]   = U[i  ][j+1][k  ];
			W[i  ][j  ][k-1]   = W[i  ][j+1][k-1];
			W[i  ][j  ][k  ]   = W[i  ][j+1][k  ];
			break;
		case B_S:
			V[i  ][j-1][k  ] = 0.0;
			U[i-1][j  ][k  ] = U[i-1][j-1][k  ];
			U[i  ][j  ][k  ] = U[i  ][j-1][k  ];
			W[i  ][j  ][k-1] = W[i  ][j-1][k-1];
			W[i  ][j  ][k  ] = W[i  ][j-1][k  ];
			break;
		case B_U:
			W[i  ][j  ][k  ]   = 0.0;
			U[i-1][j  ][k  ] = U[i-1][j  ][k+1];
			U[i  ][j  ][k  ]   = U[i  ][j  ][k+1];
			V[i  ][j-1][k  ] = V[i  ][j-1][k+1];
			V[i  ][j  ][k  ]   = V[i  ][j  ][k+1];
			break;
		case B_D:
			W[i  ][j  ][k-1] = 0.0;
			U[i-1][j  ][k  ] = U[i-1][j  ][k-1];
			U[i  ][j  ][k  ]   = U[i  ][j  ][k-1];
			V[i  ][j-1][k  ] = V[i  ][j-1][k-1];
			V[i  ][j  ][k  ]   = V[i  ][j  ][k-1];
			break;

		case B_NO:
			U[i  ][j  ][k  ]   = 0.0;
			U[i-1][j  ][k  ] = U[i-1][j+1][k  ];
			V[i  ][j  ][k  ]   = 0.0;
			V[i  ][j-1][k  ] = V[i+1][j-1][k  ];
			W[i  ][j  ][k  ]   = (W[i  ][j+1][k  ] + W[i+1][j  ][k  ]) * 0.5;
			W[i  ][j  ][k-1] = (W[i  ][j+1][k-1] + W[i+1][j  ][k-1]) * 0.5;
			break;
		case B_NW:
			U[i-1][j  ][k  ] = 0.0;
			U[i  ][j  ][k  ]   = U[i  ][j+1][k  ];
			V[i  ][j  ][k  ]   = 0.0;
			V[i  ][j-1][k  ] = V[i-1][j-1][k  ];
			W[i  ][j  ][k  ]   = (W[i  ][j+1][k  ] + W[i-1][j  ][k  ]) * 0.5;
			W[i  ][j  ][k-1] = (W[i  ][j+1][k-1] + W[i-1][j  ][k-1]) * 0.5;
			break;
		case B_NU:
			V[i  ][j  ][k  ]   = 0.0;
			V[i  ][j-1][k  ] = V[i  ][j-1][k+1];
			W[i  ][j  ][k  ]   = 0.0;
			W[i  ][j  ][k-1] = W[i  ][j+1][k-1];
			U[i  ][j  ][k  ]   = (U[i  ][j+1][k  ] + U[i  ][j  ][k+1]) * 0.5;
			U[i-1][j  ][k  ] = (U[i-1][j+1][k  ] + U[i-1][j  ][k+1]) * 0.5;
			break;
		case B_ND:
			V[i  ][j  ][k  ] = 0.0;
			V[i  ][j-1][k  ] = V[i  ][j-1][k-1];
			W[i  ][j  ][k-1] = 0.0;
			W[i  ][j  ][k  ] = W[i  ][j+1][k  ];
			U[i  ][j  ][k  ] = (U[i  ][j+1][k  ] + U[i  ][j  ][k-1]) * 0.5;
			U[i-1][j  ][k  ] = (U[i-1][j+1][k  ] + U[i-1][j  ][k-1]) * 0.5;
			break;
		case B_SO:
			U[i  ][j  ][k  ]   = 0.0;
			U[i-1][j  ][k  ] = U[i-1][j-1][k  ];
			V[i  ][j-1][k  ] = 0.0;
			V[i  ][j  ][k  ]   = V[i+1][j  ][k  ];
			W[i  ][j  ][k  ]   = (W[i  ][j-1][k  ] + W[i+1][j  ][k  ]) * 0.5;
			W[i  ][j  ][k-1] = (W[i  ][j-1][k-1] + W[i+1][j  ][k-1]) * 0.5;
			break;
		case B_SW:
			U[i-1][j  ][k  ] = 0.0;
			U[i  ][j  ][k  ]   = U[i  ][j-1][k  ];
			V[i  ][j-1][k  ] = 0.0;
			V[i  ][j  ][k  ]   = V[i-1][j  ][k  ];
			W[i  ][j  ][k  ]   = (W[i  ][j-1][k  ] + W[i-1][j  ][k  ]) * 0.5;
			W[i  ][j  ][k-1] = (W[i  ][j-1][k-1] + W[i-1][j  ][k-1]) * 0.5;
			break;
		case B_SU:
			V[i  ][j-1][k  ] = 0.0;
			V[i  ][j  ][k  ]   = V[i  ][j  ][k+1];
			W[i  ][j  ][k  ]   = 0.0;
			W[i  ][j  ][k-1] = W[i  ][j-1][k-1];
			U[i  ][j  ][k  ]   = (U[i  ][j-1][k  ] + U[i  ][j  ][k+1]) * 0.5;
			U[i-1][j  ][k  ] = (U[i-1][j-1][k  ] + U[i-1][j  ][k+1]) * 0.5;
			break;
		case B_SD:
			V[i  ][j-1][k  ] = 0.0;
			V[i  ][j  ][k  ]   = V[i  ][j  ][k-1];
			W[i  ][j  ][k-1] = 0.0;
			W[i  ][j  ][k  ]   = W[i  ][j-1][k  ];
			U[i  ][j  ][k  ]   = (U[i  ][j-1][k  ] + U[i  ][j  ][k-1]) * 0.5;
  		U[i-1][j  ][k  ]     = (U[i-1][j-1][k  ] + U[i-1][j  ][k-1]) * 0.5;
			break;
		case B_OU:
			U[i  ][j  ][k  ] = 0.0;
			U[i-1][j  ][k  ] = U[i-1][j  ][k+1];
			W[i  ][j  ][k  ] = 0.0;
			W[i  ][j  ][k-1] = W[i+1][j  ][k-1];
			V[i  ][j  ][k  ] = (V[i+1][j  ][k  ] + V[i  ][j  ][k+1]) * 0.5;
			V[i  ][j-1][k  ] = (V[i+1][j-1][k  ] + V[i  ][j-1][k+1]) * 0.5;
			break;
		case B_WU:
			U[i-1][j  ][k  ] = 0.0;
			U[i  ][j  ][k  ] = U[i  ][j  ][k+1];
			W[i  ][j  ][k  ] = 0.0;
			W[i  ][j  ][k-1] = W[i-1][j  ][k-1];
			V[i  ][j  ][k  ] = (V[i-1][j  ][k  ] + V[i  ][j  ][k+1]) * 0.5;
			V[i  ][j-1][k  ] = (V[i-1][j-1][k  ] + V[i  ][j-1][k+1]) * 0.5;
			break;
		case B_OD:
			U[i  ][j  ][k  ] = 0.0;
			U[i-1][j  ][k  ] = U[i-1][j  ][k-1];
			W[i  ][j  ][k-1] = 0.0;
			W[i  ][j  ][k  ] = W[i+1][j  ][k  ];
			V[i  ][j  ][k  ] = (V[i+1][j  ][k  ] + V[i  ][j  ][k-1]) * 0.5;
			V[i  ][j-1][k  ] = (V[i+1][j-1][k  ] + V[i  ][j-1][k-1]) * 0.5;
			break;
		case B_WD:
			U[i-1][j  ][k  ] = 0.0;
			U[i  ][j  ][k  ] = U[i  ][j  ][k-1];
			W[i  ][j  ][k-1] = 0.0;
			W[i  ][j  ][k  ] = W[i-1][j  ][k  ];
			V[i  ][j  ][k  ] = (V[i-1][j  ][k  ] + V[i  ][j  ][k-1]) * 0.5;
			V[i  ][j-1][k  ] = (V[i-1][j-1][k  ] + V[i  ][j-1][k-1]) * 0.5;
			break;

		case B_NOU:
			U[i  ][j  ][k  ] = 0.0;
			U[i-1][j  ][k  ] = (U[i-1][j+1][k  ] + U[i-1][j  ][k+1]) * 0.5;
			V[i  ][j  ][k  ] = 0.0;
			V[i  ][j-1][k  ] = (V[i+1][j-1][k  ] + V[i  ][j-1][k+1]) * 0.5;
			W[i  ][j  ][k  ] = 0.0;
			W[i  ][j  ][k-1] = (W[i+1][j  ][k-1] + W[i  ][j+1][k-1]) * 0.5;
			break;
		case B_NWU:
			U[i-1][j  ][k  ] = 0.0;
			U[i  ][j  ][k  ] = (U[i  ][j+1][k  ] + U[i  ][j  ][k+1]) * 0.5;
			V[i  ][j  ][k  ] = 0.0;
			V[i  ][j-1][k  ] = (V[i-1][j-1][k  ] + V[i  ][j-1][k+1]) * 0.5;
			W[i  ][j  ][k  ] = 0.0;
			W[i  ][j  ][k-1] = (W[i-1][j  ][k-1] + W[i  ][j+1][k-1]) * 0.5;
			break;
		case B_NOD:
			U[i  ][j  ][k  ] = 0.0;
			U[i-1][j  ][k  ] = (U[i-1][j+1][k  ] + U[i-1][j  ][k-1]) * 0.5;
			V[i  ][j  ][k  ] = 0.0;
			V[i  ][j-1][k  ] = (V[i+1][j-1][k  ] + V[i  ][j-1][k-1]) * 0.5;
			W[i  ][j  ][k-1] = 0.0;
			W[i  ][j  ][k  ] = (W[i+1][j  ][k  ] + W[i  ][j+1][k  ]) * 0.5;
			break;
		case B_NWD:
			U[i-1][j  ][k  ] = 0.0;
			U[i  ][j  ][k  ] = (U[i  ][j+1][k  ] + U[i  ][j  ][k-1]) * 0.5;
			V[i  ][j  ][k  ] = 0.0;
			V[i  ][j-1][k  ] = (V[i-1][j-1][k  ] + V[i  ][j-1][k-1]) * 0.5;
			W[i  ][j  ][k-1] = 0.0;
			W[i  ][j  ][k  ] = (W[i-1][j  ][k  ] + W[i  ][j+1][k  ]) * 0.5;
			break;
		case B_SOU:
			U[i  ][j  ][k  ] = 0.0;
			U[i-1][j  ][k  ] = (U[i-1][j-1][k  ] + U[i-1][j  ][k+1]) * 0.5;
			V[i  ][j-1][k  ] = 0.0;
			V[i  ][j  ][k  ] = (V[i+1][j  ][k  ] + V[i  ][j  ][k+1]) * 0.5;
			W[i  ][j  ][k  ] = 0.0;
			W[i  ][j  ][k-1] = (W[i+1][j  ][k-1] + W[i  ][j-1][k-1]) * 0.5;
			break;
		case B_SWU:
			U[i-1][j  ][k  ] = 0.0;
			U[i  ][j  ][k  ] = (U[i  ][j-1][k  ] + U[i  ][j  ][k+1]) * 0.5;
			V[i  ][j-1][k  ] = 0.0;
			V[i  ][j  ][k  ] = (V[i-1][j  ][k  ] + V[i  ][j  ][k+1]) * 0.5;
			W[i  ][j  ][k  ] = 0.0;
			W[i  ][j  ][k-1] = (W[i-1][j  ][k-1] + W[i  ][j-1][k-1]) * 0.5;
			break;
		case B_SOD:
			U[i  ][j  ][k  ] = 0.0;
			U[i-1][j  ][k  ] = (U[i-1][j-1][k  ] + U[i-1][j  ][k-1]) * 0.5;
			V[i  ][j-1][k  ] = 0.0;
			V[i  ][j  ][k  ] = (V[i+1][j  ][k  ] + V[i  ][j  ][k-1]) * 0.5;
			W[i  ][j  ][k-1] = 0.0;
			W[i  ][j  ][k  ] = (W[i+1][j  ][k  ] + W[i  ][j-1][k  ]) * 0.5;
			break;
		case B_SWD:
			U[i-1][j  ][k  ] = 0.0;
			U[i  ][j  ][k  ] = (U[i  ][j-1][k  ] + U[i  ][j  ][k-1]) * 0.5;
			V[i  ][j-1][k  ] = 0.0;
			V[i  ][j  ][k  ] = (V[i-1][j  ][k  ] + V[i  ][j  ][k-1]) * 0.5;
			W[i  ][j  ][k-1] = 0.0;
			W[i  ][j  ][k  ] = (W[i-1][j  ][k  ] + W[i  ][j-1][k  ]) * 0.5;
			break;

		default: //case 0
			break;
	}
}

void boundaryvalues_moving_wall(
			int i,
			int j,
			int k,
			double ***U,
			double ***V,
			double ***W,
			int ***Flag,
			double *velMW
			){
        // No-slip boundary conditions for U, V and W.
	// 18 (6 + 12) cases of moving wall in total.
	// For the remaining 8 cases of B_NOU, B_NWU, B_NOD, B_NWD, B_SOU, B_SWU, B_SOD, B_SWD,
	// moving wall boundary condition is not allowed, and no slip boundary condition applies.
	switch(getboundarytype(Flag[i  ][j  ][k  ])){

		case B_O: // wall is O. its moving direction is +/-y
			U[i  ][j  ][k  ] = 0.0;
			V[i  ][j-1][k  ] = 2.0*velMW[1] - V[i+1][j-1][k  ];
			V[i  ][j  ][k  ]   = 2.0*velMW[1] - V[i+1][j  ][k  ];
			W[i  ][j  ][k-1] = 2.0*velMW[2] - W[i+1][j  ][k-1];
			W[i  ][j  ][k  ]   = 2.0*velMW[2] - W[i+1][j  ][k  ];
			break;
		case B_W: // wall is W. its moving direction is +/-y
			U[i-1][j  ][k  ] = 0.0;
			V[i  ][j-1][k  ] = 2.0*velMW[1] - V[i-1][j-1][k  ];
			V[i  ][j  ][k  ]   = 2.0*velMW[1] - V[i-1][j  ][k  ];
			W[i  ][j  ][k-1] = 2.0*velMW[2] - W[i-1][j  ][k-1];
			W[i  ][j  ][k  ]   = 2.0*velMW[2] - W[i-1][j  ][k  ];
			break;
		case B_N: // wall is N. its moving direction is +/-z
			V[i  ][j  ][k  ] = 0.0;
			U[i-1][j  ][k  ] = 2.0*velMW[0] - U[i-1][j+1][k  ];
			U[i  ][j  ][k  ]   = 2.0*velMW[0] - U[i  ][j+1][k  ];
			W[i  ][j  ][k-1] = 2.0*velMW[2] - W[i  ][j+1][k-1];
			W[i  ][j  ][k  ]   = 2.0*velMW[2] - W[i  ][j+1][k  ];
			break;
		case B_S: // wall is S. its moving direction is +/-z
			V[i  ][j-1][k  ] = 0.0;
			U[i-1][j  ][k  ] = 2.0*velMW[0] - U[i-1][j-1][k  ];
			U[i  ][j  ][k  ]   = 2.0*velMW[0] - U[i  ][j-1][k  ];
			W[i  ][j  ][k-1] = 2.0*velMW[2] - W[i  ][j-1][k-1];
			W[i  ][j  ][k  ]   = 2.0*velMW[2] - W[i  ][j-1][k  ];
			break;
		case B_U: // wall is U. its moving direction is +/-x
			W[i  ][j  ][k  ]   = 0.0;
			U[i-1][j  ][k  ] = 2.0*velMW[0] - U[i-1][j  ][k+1];
			U[i  ][j  ][k  ]   = 2.0*velMW[0] - U[i  ][j  ][k+1];
			V[i  ][j-1][k  ] = 2.0*velMW[2] - V[i  ][j-1][k+1];
			V[i  ][j  ][k  ]   = 2.0*velMW[2] - V[i  ][j  ][k+1];
			break;
		case B_D: // wall is D. its moving direction is +/-x
			W[i  ][j  ][k-1] = 0.0;
			U[i-1][j  ][k  ] = 2.0*velMW[0] - U[i-1][j  ][k-1];
			U[i  ][j  ][k  ]   = 2.0*velMW[0] - U[i  ][j  ][k-1];
			V[i  ][j-1][k  ] = 2.0*velMW[1] - V[i  ][j-1][k-1];
			V[i  ][j  ][k  ]   = 2.0*velMW[1] - V[i  ][j  ][k-1];
			break;

		case B_NO: // circular moving direction of (O/W, N/S, U/D) is (+y,-x,0) or (-y,+x,0)
			U[i  ][j  ][k  ]   = 2.0*velMW[0] - U[i  ][j+1][k  ];
			U[i-1][j  ][k  ] = 2.0*velMW[0] - U[i-1][j+1][k  ];
			V[i  ][j  ][k  ]   = 2.0*velMW[1] - V[i+1][j  ][k  ];
			V[i  ][j-1][k  ] = 2.0*velMW[1] - V[i+1][j-1][k  ];
			W[i  ][j  ][k  ]   = - (W[i  ][j+1][k  ] + W[i+1][j  ][k  ]) * 0.5;
			W[i  ][j  ][k-1] = - (W[i  ][j+1][k-1] + W[i+1][j  ][k-1]) * 0.5;
			break;
		case B_NW: // circular moving direction of (O/W, N/S, U/D) is (+y,+x,0) or (-y,-x,0)
			U[i-1][j  ][k  ] = 2.0*velMW[0] - U[i-1][j+1][k  ];
			U[i  ][j  ][k  ]   = 2.0*velMW[0] - U[i  ][j+1][k  ];
			V[i  ][j  ][k  ]   = 2.0*velMW[1] - V[i-1][j  ][k  ];
			V[i  ][j-1][k  ] = 2.0*velMW[1] - V[i-1][j-1][k  ];
			W[i  ][j  ][k  ] = - (W[i  ][j+1][k  ] + W[i-1][j  ][k  ]) * 0.5;
			W[i  ][j  ][k-1] = - (W[i  ][j+1][k-1] + W[i-1][j  ][k-1]) * 0.5;
			break;
		case B_NU: // circular moving direction of (O/W, N/S, U/D) is (0,-z,+y) or (0,+z,-y)
			V[i  ][j  ][k  ] = 2.0*velMW[1] - V[i  ][j  ][k+1];
			V[i  ][j-1][k  ] = 2.0*velMW[1] - V[i  ][j-1][k+1];
			W[i  ][j  ][k  ] = 2.0*velMW[2] - W[i  ][j+1][k  ];
			W[i  ][j  ][k-1] = 2.0*velMW[2] - W[i  ][j+1][k-1];
			U[i  ][j  ][k  ] = - (U[i  ][j+1][k  ] + U[i  ][j  ][k+1]) * 0.5;
			U[i-1][j  ][k  ] = - (U[i-1][j+1][k  ] + U[i-1][j  ][k+1]) * 0.5;
			break;
		case B_ND: // circular moving direction of (O/W, N/S, U/D) is (0,+z,+y) or (0,-z,-y)
			V[i  ][j  ][k  ]   = 2.0*velMW[1] - V[i  ][j  ][k-1];
			V[i  ][j-1][k  ] = 2.0*velMW[1] - V[i  ][j-1][k-1];
			W[i  ][j  ][k-1] = 2.0*velMW[2] - W[i  ][j+1][k-1];
			W[i  ][j  ][k  ]   = 2.0*velMW[2] - W[i  ][j+1][k  ];
			U[i  ][j  ][k  ]   = -(U[i  ][j+1][k  ] + U[i  ][j  ][k-1]) * 0.5;
			U[i-1][j  ][k  ] = -(U[i-1][j+1][k  ] + U[i-1][j  ][k-1]) * 0.5;
			break;
		case B_SO: // circular moving direction of (O/W, N/S, U/D) is (+y,+x,0) or (-y,-x,0)
			U[i  ][j  ][k  ]   = 2.0*velMW[0] - U[i  ][j-1][k  ];
			U[i-1][j  ][k  ] = 2.0*velMW[0] - U[i-1][j-1][k  ];
			V[i  ][j-1][k  ] = 2.0*velMW[1] - V[i+1][j-1][k  ];
			V[i  ][j  ][k  ]   = 2.0*velMW[1] - V[i+1][j  ][k  ];
			W[i  ][j  ][k  ]   = - (W[i  ][j-1][k  ] + W[i+1][j  ][k  ]) * 0.5;
			W[i  ][j  ][k-1] = - (W[i  ][j-1][k-1] + W[i+1][j  ][k-1]) * 0.5;
			break;
		case B_SW: // circular moving direction of (O/W, N/S, U/D) is (+y,-x,0) or (-y,+x,0)
			U[i-1][j  ][k  ] = 2.0*velMW[0] - U[i-1][j-1][k  ];
			U[i  ][j  ][k  ]   = 2.0*velMW[0] - U[i  ][j-1][k  ];
			V[i  ][j-1][k  ] = 2.0*velMW[1] - V[i-1][j-1][k  ];
			V[i  ][j  ][k  ]   = 2.0*velMW[1] - V[i-1][j  ][k  ];
			W[i  ][j  ][k  ]   = - (W[i  ][j-1][k  ] + W[i-1][j  ][k  ]) * 0.5;
			W[i  ][j  ][k-1] = - (W[i  ][j-1][k-1] + W[i-1][j  ][k-1]) * 0.5;
			break;
		case B_SU: // circular moving direction of (O/W, N/S, U/D) is (0,+z,+y) or (0,-z,-y)
			V[i  ][j-1][k  ] = 2.0*velMW[1] - V[i  ][j-1][k+1];
			V[i  ][j  ][k  ]   = 2.0*velMW[1] - V[i  ][j  ][k+1];
			W[i  ][j  ][k  ]   = 2.0*velMW[2] - W[i  ][j-1][k  ];
			W[i  ][j  ][k-1] = 2.0*velMW[2] - W[i  ][j-1][k-1];
			U[i  ][j  ][k  ]   = - (U[i  ][j-1][k  ] + U[i  ][j  ][k+1]) * 0.5;
			U[i-1][j  ][k  ] = - (U[i-1][j-1][k  ] + U[i-1][j  ][k+1]) * 0.5;
			break;
		case B_SD: // circular moving direction of (O/W, N/S, U/D) is (0,-z,+y) or (0,+z,-y)
			V[i  ][j-1][k  ] = 2.0*velMW[1] - V[i  ][j-1][k-1];
			V[i  ][j  ][k  ]   = 2.0*velMW[1] - V[i  ][j  ][k-1];
			W[i  ][j  ][k-1] = 2.0*velMW[2] - W[i  ][j-1][k-1];
			W[i  ][j  ][k  ]   = 2.0*velMW[2] - W[i  ][j-1][k  ];
			U[i  ][j  ][k  ]   = - (U[i  ][j-1][k  ] + U[i  ][j  ][k-1]) * 0.5;
			U[i-1][j  ][k  ] = - (U[i-1][j-1][k  ] + U[i-1][j  ][k-1]) * 0.5;
			break;
		case B_OU: // circular moving direction of (O/W, N/S, U/D) is (+z,0,-x) or (-z,0,+x)
			U[i  ][j  ][k  ]   = 2.0*velMW[0] - U[i  ][j  ][k+1];
			U[i-1][j  ][k  ] = 2.0*velMW[0] - U[i-1][j  ][k+1];
			W[i  ][j  ][k  ]   = 2.0*velMW[2] - W[i+1][j  ][k  ];
			W[i  ][j  ][k-1] = 2.0*velMW[2] - W[i+1][j  ][k-1];
			V[i  ][j  ][k  ]   = - (V[i+1][j  ][k  ] + V[i  ][j  ][k+1]) * 0.5;
			V[i  ][j-1][k  ] = - (V[i+1][j-1][k  ] + V[i  ][j-1][k+1]) * 0.5;
			break;
		case B_WU: // circular moving direction of (O/W, N/S, U/D) is (+z,0,+x) or (-z,0,-x)
			U[i-1][j  ][k  ] = 2.0*velMW[0] - U[i-1][j  ][k+1];
			U[i  ][j  ][k  ]   = 2.0*velMW[0] - U[i  ][j  ][k+1];
			W[i  ][j  ][k  ]   = 2.0*velMW[2] - W[i-1][j  ][k  ];
			W[i  ][j  ][k-1] = 2.0*velMW[2] - W[i-1][j  ][k-1];
			V[i  ][j  ][k  ]   = - (V[i-1][j  ][k  ] + V[i  ][j  ][k+1]) * 0.5;
			V[i  ][j-1][k  ] = - (V[i-1][j-1][k  ] + V[i  ][j-1][k+1]) * 0.5;
			break;
		case B_OD: // circular moving direction of (O/W, N/S, U/D) is (+z,0,+x) or (-z,0,-x)
			U[i  ][j  ][k  ]   = 2.0*velMW[0] - U[i  ][j  ][k-1];
			U[i-1][j  ][k  ] = 2.0*velMW[0] - U[i-1][j  ][k-1];
			W[i  ][j  ][k-1] = 2.0*velMW[2] - W[i+1][j  ][k-1];
			W[i  ][j  ][k  ]   = 2.0*velMW[2] - W[i+1][j  ][k  ];
			V[i  ][j  ][k  ]   = - (V[i+1][j  ][k  ] + V[i  ][j  ][k-1]) * 0.5;
			V[i  ][j-1][k  ] = - (V[i+1][j-1][k  ] + V[i  ][j-1][k-1]) * 0.5;
			break;
		case B_WD: // circular moving direction of (O/W, N/S, U/D) is (+z,0,-x) or (-z,0,+x)
			U[i-1][j  ][k  ] = 2.0*velMW[0] - U[i-1][j  ][k-1];
			U[i  ][j  ][k  ]   = 2.0*velMW[0] - U[i  ][j  ][k-1];
			W[i  ][j  ][k-1] = 2.0*velMW[2] - W[i-1][j  ][k-1];
			W[i  ][j  ][k  ]   = 2.0*velMW[2] - W[i-1][j  ][k  ];
			V[i  ][j  ][k  ]   = - (V[i-1][j  ][k  ] + V[i  ][j  ][k-1]) * 0.5;
			V[i  ][j-1][k  ] = - (V[i-1][j-1][k  ] + V[i  ][j-1][k-1]) * 0.5;
			break;

		//before we find a better solution, the following 8 cases will be respectively same with no-slip.
		/*case B_NOU:
			printf("Warning: for now the moving wall condition is set same as no-slip, when the flag is B_NOU, B_NWU etc.");
			U[i  ][j  ][k  ] = 0.0;
			U[i-1][j  ][k  ] = -(U[i-1][j+1][k  ] + U[i-1][j  ][k+1]) * 0.5;
			V[i  ][j  ][k  ] = 0.0;
			V[i  ][j-1][k  ] = -(V[i+1][j-1][k  ] + V[i  ][j-1][k+1]) * 0.5;
			W[i  ][j  ][k  ] = 0.0;
			W[i  ][j  ][k-1] = -(W[i+1][j  ][k-1] + W[i  ][j+1][k-1]) * 0.5;
			break;
		case B_NWU:
			printf("Warning: for now the moving wall condition is set same as no-slip, when the flag is B_NOU, B_NWU etc.");
			U[i-1][j  ][k  ] = 0.0;
			U[i  ][j  ][k  ] = -(U[i  ][j+1][k  ] + U[i  ][j  ][k+1]) * 0.5;
			V[i  ][j  ][k  ] = 0.0;
			V[i  ][j-1][k  ] = -(V[i-1][j-1][k  ] + V[i  ][j-1][k+1]) * 0.5;
			W[i  ][j  ][k  ] = 0.0;
			W[i  ][j  ][k-1] = -(W[i-1][j  ][k-1] + W[i  ][j+1][k-1]) * 0.5;
			break;
		case B_NOD:
			printf("Warning: for now the moving wall condition is set same as no-slip, when the flag is B_NOU, B_NWU etc.");
			U[i  ][j  ][k  ] = 0.0;
			U[i-1][j  ][k  ] = -(U[i-1][j+1][k  ] + U[i-1][j  ][k-1]) * 0.5;
			V[i  ][j  ][k  ] = 0.0;
			V[i  ][j-1][k  ] = -(V[i+1][j-1][k  ] + V[i  ][j-1][k-1]) * 0.5;
			W[i  ][j  ][k-1] = 0.0;
			W[i  ][j  ][k  ] = -(W[i+1][j  ][k  ] + W[i  ][j+1][k  ]) * 0.5;
			break;
		case B_NWD:
			printf("Warning: for now the moving wall condition is set same as no-slip, when the flag is B_NOU, B_NWU etc.");
			U[i-1][j  ][k  ] = 0.0;
			U[i  ][j  ][k  ] = -(U[i  ][j+1][k  ] + U[i  ][j  ][k-1]) * 0.5;
			V[i  ][j  ][k  ] = 0.0;
			V[i  ][j-1][k  ] = -(V[i-1][j-1][k  ] + V[i  ][j-1][k-1]) * 0.5;
			W[i  ][j  ][k-1] = 0.0;
			W[i  ][j  ][k  ] = -(W[i-1][j  ][k  ] + W[i  ][j+1][k  ]) * 0.5;
			break;
		case B_SOU:
			printf("Warning: for now the moving wall condition is set same as no-slip, when the flag is B_NOU, B_NWU etc.");
			U[i  ][j  ][k  ] = 0.0;
			U[i-1][j  ][k  ] = -(U[i-1][j-1][k  ] + U[i-1][j  ][k+1]) * 0.5;
			V[i  ][j-1][k  ] = 0.0;
			V[i  ][j  ][k  ] = -(V[i+1][j  ][k  ] + V[i  ][j  ][k+1]) * 0.5;
			W[i  ][j  ][k  ] = 0.0;
			W[i  ][j  ][k-1] = -(W[i+1][j  ][k-1] + W[i  ][j-1][k-1]) * 0.5;
			break;
		case B_SWU:
			printf("Warning: for now the moving wall condition is set same as no-slip, when the flag is B_NOU, B_NWU etc.");
			U[i-1][j  ][k  ] = 0.0;
			U[i  ][j  ][k  ] = -(U[i  ][j-1][k  ] + U[i  ][j  ][k+1]) * 0.5;
			V[i  ][j-1][k  ] = 0.0;
			V[i  ][j  ][k  ] = -(V[i-1][j  ][k  ] + V[i  ][j  ][k+1]) * 0.5;
			W[i  ][j  ][k  ] = 0.0;
			W[i  ][j  ][k-1] = -(W[i-1][j  ][k-1] + W[i  ][j-1][k-1]) * 0.5;
			break;
		case B_SOD:
			printf("Warning: for now the moving wall condition is set same as no-slip, when the flag is B_NOU, B_NWU etc.");
			U[i  ][j  ][k  ] = 0.0;
			U[i-1][j  ][k  ] = -(U[i-1][j-1][k  ] + U[i-1][j  ][k-1]) * 0.5;
			V[i  ][j-1][k  ] = 0.0;
			V[i  ][j  ][k  ] = -(V[i+1][j  ][k  ] + V[i  ][j  ][k-1]) * 0.5;
			W[i  ][j  ][k-1] = 0.0;
			W[i  ][j  ][k  ] = -(W[i+1][j  ][k  ] + W[i  ][j-1][k  ]) * 0.5;
			break;
		case B_SWD:
			printf("Warning: for now the moving wall condition is set same as no-slip, when the flag is B_NOU, B_NWU etc.");
			U[i-1][j  ][k  ] = 0.0;
			U[i  ][j  ][k  ] = -(U[i  ][j-1][k  ] + U[i  ][j  ][k-1]) * 0.5;
			V[i  ][j-1][k  ] = 0.0;
			V[i  ][j  ][k  ] = -(V[i-1][j  ][k  ] + V[i  ][j  ][k-1]) * 0.5;
			W[i  ][j  ][k-1] = 0.0;
			W[i  ][j  ][k  ] = -(W[i-1][j  ][k  ] + W[i  ][j-1][k  ]) * 0.5;
			break;
		*/

		//before we find a better solution, the following 8 cases will be forbidden.
		/*case B_NOU:
			printf("Warning: It is forbidden for moving wall, when the flag is B_NOU, B_NWU etc.");
			break;
		case B_NWU:
			printf("Warning: It is forbidden for moving wall, when the flag is B_NOU, B_NWU etc.");
			break;
		case B_NOD:
			printf("Warning: It is forbidden for moving wall, when the flag is B_NOU, B_NWU etc.");
			break;
		case B_NWD:
			printf("Warning: It is forbidden for moving wall, when the flag is B_NOU, B_NWU etc.");
			break;
		case B_SOU:
			printf("Warning: It is forbidden for moving wall, when the flag is B_NOU, B_NWU etc.");
			break;
		case B_SWU:
			printf("Warning: It is forbidden for moving wall, when the flag is B_NOU, B_NWU etc.");
			break;
		case B_SOD:
			printf("Warning: It is forbidden for moving wall, when the flag is B_NOU, B_NWU etc.");
			break;
		case B_SWD:
			printf("Warning: It is forbidden for moving wall, when the flag is B_NOU, B_NWU etc.");
			break;
    */
		case 0: //case 0: inner cell -> do nothing
			break;
		default:
			printf("Warning: It is forbidden for moving wall, when the flag is B_NOU, B_NWU etc.");  //if we come in here, then were delaing with B_??? cell.
			break;
	}
}

void boundaryvalues_outflow(
			int i,
			int j,
			int k,
			double ***U,
			double ***V,
			double ***W,
			int ***Flag
			){
        // No-slip boundary conditions for U, V and W.
	// 26 (6 + 12 + 8) cases in total.
	switch(getboundarytype(Flag[i  ][j  ][k  ])){
		case B_O:
			U[i  ][j  ][k  ]   = U[i+1][j  ][k  ];

			V[i  ][j-1][k  ] = V[i+1][j-1][k  ];
			V[i  ][j  ][k  ]   = V[i+1][j  ][k  ];

			W[i  ][j  ][k-1] = W[i+1][j  ][k-1];
			W[i  ][j  ][k  ]   = W[i+1][j  ][k  ];
			break;
		case B_W: //step outflow
			U[i-1][j  ][k  ] = U[i-2][j  ][k  ];

			V[i  ][j-1][k  ] = V[i-1][j-1][k  ];
			V[i  ][j  ][k  ]   = V[i-1][j  ][k  ];

			W[i  ][j  ][k-1] = W[i-1][j  ][k-1];
			W[i  ][j  ][k  ]   = W[i-1][j  ][k  ];
			break;
		case B_N:
			V[i  ][j  ][k  ] = V[i  ][j+1][k  ];

			U[i-1][j  ][k  ] = U[i-1][j+1][k  ];
			U[i  ][j  ][k  ] = U[i  ][j+1][k  ];

			W[i  ][j  ][k-1] = W[i  ][j+1][k-1];
			W[i  ][j  ][k  ] = W[i  ][j+1][k  ];
			break;
		case B_S:
			V[i  ][j-1][k  ] = V[i  ][j-2][k  ];
			U[i-1][j  ][k  ] = U[i-1][j-1][k  ];
			U[i  ][j  ][k  ] = U[i  ][j-1][k  ];
			W[i  ][j  ][k-1] = W[i  ][j-1][k-1];
			W[i  ][j  ][k  ] = W[i  ][j-1][k  ];

			break;
		case B_U:
			W[i  ][j  ][k  ] = W[i  ][j  ][k+1];
			U[i-1][j  ][k  ] = U[i-1][j  ][k+1];
			U[i  ][j  ][k  ] = U[i  ][j  ][k+1];
			V[i  ][j-1][k  ] = V[i  ][j-1][k+1];
			V[i  ][j  ][k  ] = V[i  ][j  ][k+1];
			break;
		case B_D:
			W[i  ][j  ][k-1] = W[i  ][j  ][k-2];
			U[i-1][j  ][k  ] = U[i-1][j  ][k-1];
			U[i  ][j  ][k  ] = U[i  ][j  ][k-1];
			V[i  ][j-1][k  ] = V[i  ][j-1][k-1];
			V[i  ][j  ][k  ] = V[i  ][j  ][k-1];
			break;

		case B_NO:
			U[i  ][j  ][k  ] = U[i+1][j  ][k  ];
			U[i-1][j  ][k  ] = U[i-1][j+1][k  ];
			V[i  ][j  ][k  ] = V[i  ][j+1][k  ];
			V[i  ][j-1][k  ] = V[i+1][j-1][k  ];
			W[i  ][j  ][k  ] = (W[i  ][j+1][k  ] + W[i+1][j  ][k  ]) * 0.5;
			W[i  ][j  ][k-1] = (W[i  ][j+1][k-1] + W[i+1][j  ][k-1]) * 0.5;
			break;
		case B_NW:
			U[i-1][j  ][k  ] = U[i-2][j  ][k  ];
			U[i  ][j  ][k  ] = U[i  ][j+1][k  ];
			V[i  ][j  ][k  ] = V[i  ][j+1][k  ];
			V[i  ][j-1][k  ] = V[i-1][j-1][k  ];
			W[i  ][j  ][k  ] = (W[i  ][j+1][k  ] + W[i-1][j  ][k  ]) * 0.5;
			W[i  ][j  ][k-1] = (W[i  ][j+1][k-1] + W[i-1][j  ][k-1]) * 0.5;
			break;
		case B_NU:
			V[i  ][j  ][k  ] = V[i  ][j+1][k  ];
			V[i  ][j-1][k  ] = V[i  ][j-1][k+1];
			W[i  ][j  ][k  ] = W[i  ][j  ][k+1];
			W[i  ][j  ][k-1] = W[i  ][j+1][k-1];
			U[i  ][j  ][k  ] = (U[i  ][j+1][k  ] + U[i  ][j  ][k+1]) * 0.5;
			U[i-1][j  ][k  ] = (U[i-1][j+1][k  ] + U[i-1][j  ][k+1]) * 0.5;
			break;
		case B_ND:
			V[i  ][j  ][k  ] = V[i  ][j+1][k  ];
			V[i  ][j-1][k  ] = V[i  ][j-1][k-1];
			W[i  ][j  ][k-1] = W[i  ][j  ][k-2];
			W[i  ][j  ][k  ] = W[i  ][j+1][k  ];
			U[i  ][j  ][k  ] = (U[i  ][j+1][k  ] + U[i  ][j  ][k-1]) * 0.5;
			U[i-1][j  ][k  ] = (U[i-1][j+1][k  ] + U[i-1][j  ][k-1]) * 0.5;
			break;
		case B_SO:
			U[i  ][j  ][k  ] = U[i+1][j  ][k  ];
			U[i-1][j  ][k  ] = U[i-1][j-1][k  ];
			V[i  ][j-1][k  ] = V[i  ][j-2][k  ];
			V[i  ][j  ][k  ] = V[i+1][j  ][k  ];
			W[i  ][j  ][k  ] = (W[i  ][j-1][k  ] + W[i+1][j  ][k  ]) * 0.5;
			W[i  ][j  ][k-1] = (W[i  ][j-1][k-1] + W[i+1][j  ][k-1]) * 0.5;
			break;
		case B_SW:
			U[i-1][j  ][k  ] = U[i-2][j  ][k  ];
			U[i  ][j  ][k  ] = U[i  ][j-1][k  ];
			V[i  ][j-1][k  ] = V[i  ][j-2][k  ];
			V[i  ][j  ][k  ] = V[i-1][j  ][k  ];
			W[i  ][j  ][k  ] = (W[i  ][j-1][k  ] + W[i-1][j  ][k  ]) * 0.5;
			W[i  ][j  ][k-1] = (W[i  ][j-1][k-1] + W[i-1][j  ][k-1]) * 0.5;
			break;
		case B_SU:
			V[i  ][j-1][k  ] = V[i  ][j-2][k  ];
			V[i  ][j  ][k  ] = V[i  ][j  ][k+1];
			W[i  ][j  ][k  ] = W[i  ][j  ][k+1];
			W[i  ][j  ][k-1] = W[i  ][j-1][k-1];
			U[i  ][j  ][k  ] = (U[i  ][j-1][k  ] + U[i  ][j  ][k+1]) * 0.5;
			U[i-1][j  ][k  ] = (U[i-1][j-1][k  ] + U[i-1][j  ][k+1]) * 0.5;
			break;
		case B_SD:
			V[i  ][j-1][k  ] = V[i  ][j-2][k  ];
			V[i  ][j  ][k  ] = V[i  ][j  ][k-1];
			W[i  ][j  ][k-1] = W[i  ][j  ][k-2];
			W[i  ][j  ][k  ] = W[i  ][j-1][k  ];
			U[i  ][j  ][k  ] = (U[i  ][j-1][k  ] + U[i  ][j  ][k-1]) * 0.5;
			U[i-1][j  ][k  ] = (U[i-1][j-1][k  ] + U[i-1][j  ][k-1]) * 0.5;
			break;
		case B_OU:
			U[i  ][j  ][k  ] = U[i+1][j  ][k  ];
			U[i-1][j  ][k  ] = U[i-1][j  ][k+1];
			W[i  ][j  ][k  ] = W[i  ][j  ][k+1];
			W[i  ][j  ][k-1] = W[i+1][j  ][k-1];
			V[i  ][j  ][k  ] = (V[i+1][j  ][k  ] + V[i  ][j  ][k+1]) * 0.5;
			V[i  ][j-1][k  ] = (V[i+1][j-1][k  ] + V[i  ][j-1][k+1]) * 0.5;
			break;
		case B_WU:
			U[i-1][j  ][k  ] = U[i-2][j  ][k  ];
			U[i  ][j  ][k  ] = U[i  ][j  ][k+1];
			W[i  ][j  ][k  ] = W[i  ][j  ][k+1];
			W[i  ][j  ][k-1] = W[i-1][j  ][k-1];
			V[i  ][j  ][k  ] = (V[i-1][j  ][k  ] + V[i  ][j  ][k+1]) * 0.5;
			V[i  ][j-1][k  ] = (V[i-1][j-1][k  ] + V[i  ][j-1][k+1]) * 0.5;
			break;
		case B_OD:
			U[i  ][j  ][k  ] = U[i+1][j  ][k  ];
			U[i-1][j  ][k  ] = U[i-1][j  ][k-1];
			W[i  ][j  ][k-1] = W[i  ][j  ][k-2];
			W[i  ][j  ][k  ] = W[i+1][j  ][k  ];
			V[i  ][j  ][k  ] = (V[i+1][j  ][k  ] + V[i  ][j  ][k-1]) * 0.5;
			V[i  ][j-1][k  ] = (V[i+1][j-1][k  ] + V[i  ][j-1][k-1]) * 0.5;
			break;
		case B_WD:
			U[i-1][j  ][k  ] = U[i-2][j  ][k  ];
			U[i  ][j  ][k  ] = U[i  ][j  ][k-1];
			W[i  ][j  ][k-1] = W[i  ][j  ][k-2];
			W[i  ][j  ][k  ] = W[i-1][j  ][k  ];
			V[i  ][j  ][k  ] = (V[i-1][j  ][k  ] + V[i  ][j  ][k-1]) * 0.5;
			V[i  ][j-1][k  ] = (V[i-1][j-1][k  ] + V[i  ][j-1][k-1]) * 0.5;
			break;

		case B_NOU:
			U[i  ][j  ][k  ] = U[i+1][j  ][k  ];
			U[i-1][j  ][k  ] = (U[i-1][j+1][k  ] + U[i-1][j  ][k+1]) * 0.5;
			V[i  ][j  ][k  ] = V[i  ][j+1][k  ];
			V[i  ][j-1][k  ] = (V[i+1][j-1][k  ] + V[i  ][j-1][k+1]) * 0.5;
			W[i  ][j  ][k  ] = W[i  ][j  ][k+1];
			W[i  ][j  ][k-1] = (W[i+1][j  ][k-1] + W[i  ][j+1][k-1]) * 0.5;
			break;
		case B_NWU:
			U[i-1][j  ][k  ] = U[i-2][j  ][k  ];
			U[i  ][j  ][k  ] = (U[i  ][j+1][k  ] + U[i  ][j  ][k+1]) * 0.5;
			V[i  ][j  ][k  ] = V[i  ][j+1][k  ];
			V[i  ][j-1][k  ] = (V[i-1][j-1][k  ] + V[i  ][j-1][k+1]) * 0.5;
			W[i  ][j  ][k  ] = W[i  ][j  ][k+1];
			W[i  ][j  ][k-1] = (W[i-1][j  ][k-1] + W[i  ][j+1][k-1]) * 0.5;
			break;
		case B_NOD:
			U[i  ][j  ][k  ] = U[i+1][j  ][k  ];
			U[i-1][j  ][k  ] = (U[i-1][j+1][k  ] + U[i-1][j  ][k-1]) * 0.5;
			V[i  ][j  ][k  ] = V[i  ][j+1][k  ];
			V[i  ][j-1][k  ] = (V[i+1][j-1][k  ] + V[i  ][j-1][k-1]) * 0.5;
			W[i  ][j  ][k-1] = W[i  ][j  ][k-2];
			W[i  ][j  ][k  ] = (W[i+1][j  ][k  ] + W[i  ][j+1][k  ]) * 0.5;
			break;
		case B_NWD:
			U[i-1][j  ][k  ] = U[i-2][j  ][k  ];
			U[i  ][j  ][k  ] = (U[i  ][j+1][k  ] + U[i  ][j  ][k-1]) * 0.5;
			V[i  ][j  ][k  ] = V[i  ][j+1][k  ];
			V[i  ][j-1][k  ] = (V[i-1][j-1][k  ] + V[i  ][j-1][k-1]) * 0.5;
			W[i  ][j  ][k-1] = W[i  ][j  ][k-2];
			W[i  ][j  ][k  ] = (W[i-1][j  ][k  ] + W[i  ][j+1][k  ]) * 0.5;
			break;
		case B_SOU:
			U[i  ][j  ][k  ] = U[i+1][j  ][k  ];
			U[i-1][j  ][k  ] = (U[i-1][j-1][k  ] + U[i-1][j  ][k+1]) * 0.5;
			V[i  ][j-1][k  ] = V[i  ][j-2][k  ];
			V[i  ][j  ][k  ] = (V[i+1][j  ][k  ] + V[i  ][j  ][k+1]) * 0.5;
			W[i  ][j  ][k  ] = W[i  ][j  ][k+1];
			W[i  ][j  ][k-1] = (W[i+1][j  ][k-1] + W[i  ][j-1][k-1]) * 0.5;
			break;
		case B_SWU:
			U[i-1][j  ][k  ] = U[i-2][j  ][k  ];
			U[i  ][j  ][k  ] = (U[i  ][j-1][k  ] + U[i  ][j  ][k+1]) * 0.5;
			V[i  ][j-1][k  ] = V[i  ][j-2][k  ];
			V[i  ][j  ][k  ] = (V[i-1][j  ][k  ] + V[i  ][j  ][k+1]) * 0.5;
			W[i  ][j  ][k  ] = W[i  ][j  ][k+1];
			W[i  ][j  ][k-1] = (W[i-1][j  ][k-1] + W[i  ][j-1][k-1]) * 0.5;
			break;
		case B_SOD:
			U[i  ][j  ][k  ] = U[i+1][j  ][k  ];
			U[i-1][j  ][k  ] = (U[i-1][j-1][k  ] + U[i-1][j  ][k-1]) * 0.5;
			V[i  ][j-1][k  ] = V[i  ][j-2][k  ];
			V[i  ][j  ][k  ] = (V[i+1][j  ][k  ] + V[i  ][j  ][k-1]) * 0.5;
			W[i  ][j  ][k-1] = W[i  ][j  ][k-2];
			W[i  ][j  ][k  ] = (W[i+1][j  ][k  ] + W[i  ][j-1][k  ]) * 0.5;
			break;
		case B_SWD:
			U[i-1][j  ][k  ] = U[i-2][j  ][k  ];
			U[i  ][j  ][k  ] = (U[i  ][j-1][k  ] + U[i  ][j  ][k-1]) * 0.5;
			V[i  ][j-1][k  ] = V[i  ][j-2][k  ];
			V[i  ][j  ][k  ] = (V[i-1][j  ][k  ] + V[i  ][j  ][k-1]) * 0.5;
			W[i  ][j  ][k-1] = W[i  ][j  ][k-2];
			W[i  ][j  ][k  ] = (W[i-1][j  ][k  ] + W[i  ][j-1][k  ]) * 0.5;
			break;

		default: //case 0
			//printf("case 0\n");
			break;
	}

}

void boundaryvalues_inflow(
	int i,
	int j,
	int k,
	double ***U,
	double ***V,
	double ***W,
	int ***Flag,
	double velIN) { //carefull when setting up geometry file: here velIN is assumed to be velocity, perpendicular on the cell-fluid border(s)
		switch(getboundarytype(Flag[i  ][j  ][k  ])){
			case B_O:	U[i  ][j  ][k  ] = velIN;break;
			case B_W: U[i-1][j  ][k  ] = -velIN; break;
			case B_N: V[i  ][j  ][k  ] = velIN;	break;
			case B_S: V[i  ][j-1][k  ] = -velIN; break;
			case B_U: W[i  ][j  ][k  ] = velIN; break;
			case B_D: W[i  ][j  ][k-1] = -velIN; break;

			case 0: //when bound. cell is inner. this is only to be set if well change velIN to a vector.
				break;
			default: printf("Trying to set inflow at edge or corner cells. Not allowed!\n"); break; //when we have B_?? or B_??? we cant set it.
		}

}

void boundaryvalues_pressure(double ***P,int ***Flag,int imax,int jmax,int kmax){
	int i,j,k;

	/* set boundary values, here just for the 'real' boundaries - no air included yet (if even needed?) */
		for(i = 0; i <= imax+1; i++) {
			for(j = 0; j <= jmax+1; j++) {
				for(k = 0; k <= kmax+1; k++) {

					switch(getboundarytype(Flag[i  ][j  ][k  ])){
					case B_O: P[i  ][j  ][k  ]  = P[i+1][j  ][k  ]; break;

					case B_W: P[i  ][j  ][k  ]  = P[i-1][j  ][k  ]; break;

					case B_N: P[i  ][j  ][k  ]  = P[i  ][j+1][k  ]; break;

					case B_S: P[i  ][j  ][k  ]  = P[i  ][j-1][k  ]; break;

					case B_U: P[i  ][j  ][k  ]  = P[i  ][j  ][k+1]; break;

					case B_D: P[i  ][j  ][k  ]  = P[i  ][j  ][k-1]; break;

					case B_NO: P[i  ][j  ][k  ] = (P[i+1][j  ][k  ] + P[i  ][j+1][k  ])*0.5; break;
					case B_NW: P[i  ][j  ][k  ] = (P[i-1][j  ][k  ] + P[i  ][j+1][k  ])*0.5; break;

					case B_SO: P[i  ][j  ][k  ] = (P[i+1][j  ][k  ] + P[i  ][j-1][k  ])*0.5; break;
					case B_SW: P[i  ][j  ][k  ] = (P[i-1][j  ][k  ] + P[i  ][j-1][k  ])*0.5; break;

					case B_NU: P[i  ][j  ][k  ] = (P[i  ][j  ][k+1] + P[i  ][j+1][k  ])*0.5; break;
					case B_ND: P[i  ][j  ][k  ] = (P[i  ][j  ][k-1] + P[i  ][j+1][k  ])*0.5; break;

					case B_SU: P[i  ][j  ][k  ] = (P[i  ][j  ][k+1] + P[i  ][j-1][k  ])*0.5; break;
					case B_SD: P[i  ][j  ][k  ] = (P[i  ][j  ][k-1] + P[i  ][j-1][k  ])*0.5; break;

					case B_OU: P[i  ][j  ][k  ] = (P[i+1][j  ][k  ] + P[i  ][j  ][k+1])*0.5; break;
					case B_WU: P[i  ][j  ][k  ] = (P[i-1][j  ][k  ] + P[i  ][j  ][k+1])*0.5; break;

					case B_OD: P[i  ][j  ][k  ] = (P[i+1][j  ][k  ] + P[i  ][j  ][k-1])*0.5; break;
					case B_WD: P[i  ][j  ][k  ] = (P[i-1][j  ][k  ] + P[i  ][j  ][k-1])*0.5; break;

					case B_NOU: P[i  ][j  ][k  ] = (P[i+1][j  ][k  ] + P[i  ][j+1][k  ] + P[i  ][j  ][k+1])*1.0/3.0; break;
					case B_NWU: P[i  ][j  ][k  ] = (P[i-1][j  ][k  ] + P[i  ][j+1][k  ] + P[i  ][j  ][k+1])*1.0/3.0; break;

					case B_SOU: P[i  ][j  ][k  ] = (P[i+1][j  ][k  ] + P[i  ][j-1][k  ] + P[i  ][j  ][k+1])*1.0/3.0; break;
					case B_SWU: P[i  ][j  ][k  ] = (P[i-1][j  ][k  ] + P[i  ][j-1][k  ] + P[i  ][j  ][k+1])*1.0/3.0; break;

					case B_NOD: P[i  ][j  ][k  ] = (P[i+1][j  ][k  ] + P[i  ][j+1][k  ] + P[i  ][j  ][k-1])*1.0/3.0; break;
					case B_NWD: P[i  ][j  ][k  ] = (P[i-1][j  ][k  ] + P[i  ][j+1][k  ] + P[i  ][j  ][k-1])*1.0/3.0; break;

					case B_SOD: P[i  ][j  ][k  ] = (P[i+1][j  ][k  ] + P[i  ][j-1][k  ] + P[i  ][j  ][k-1])*1.0/3.0; break;
					case B_SWD: P[i  ][j  ][k  ] = (P[i-1][j  ][k  ] + P[i  ][j-1][k  ] + P[i  ][j  ][k-1])*1.0/3.0; break;


					default: break;
					}
				}
			}
		}

}

//the boundary values are set cell by cell
void boundaryvalues(
		int imax,
		int jmax,
		int kmax,
		double ***U,
		double ***V,
		double ***W,
		double ***P,
		double ***F,
		double ***G,
		double ***H,
		char *problem,  //should comment out? probably not needed.
		int ***Flag,
		double velIN,
		double *velMW
) {
	int i, j, k;
	for (i=0; i<imax+2; i++) {
		for (j=0; j<jmax+2; j++){
			for (k=0; k<kmax+2; k++){
				switch (getcelltype(Flag[i  ][j  ][k  ])) {
				case NO_SLIP:
					boundaryvalues_no_slip(i, j, k, U, V, W, Flag);
					/*printf("noslip\t");*/
					break;
				case FREE_SLIP:
					boundaryvalues_free_slip(i, j, k, U, V, W, Flag);
					/*printf("free slip\t");*/
					break;
				case INFLOW:
					boundaryvalues_no_slip(i, j, k, U, V, W, Flag);
					boundaryvalues_inflow(i, j, k, U, V, W, Flag, velIN);
					/*printf("inflow\t");*/
					break;
				case OUTFLOW:
					boundaryvalues_outflow(i, j, k, U, V, W, Flag);
					/*printf("outflow\n");*/
					break;
				case MOVING_WALL:
					boundaryvalues_moving_wall(i, j, k, U, V, W, Flag, velMW);
					/*printf("movingwall\t");*/
					break;
				default: //if we get to here, our cell is air or water. (temp>6) Maybe need to add something here when we do free surfaces.

					break;
				}

			}
		}
	}
}


/// Solve Poisson's equation for a rectangular box with Dirichlet
/// boundary conditions on each face.
/// \param source is a pointer to a flattened 3-D array for the source function
/// \param potential is a pointer to a flattened 3-D array for the calculated potential
/// \param xsize is the number of elements in the x-direction
/// \param ysize is the number of elements in the y-direction
/// \param zsize is the number of elements in the z-direction
/// \param delta is the voxel spacing in all directions
/// \param numiters is the number of iterations to perform
/// \param numcores is the number of CPU cores to use.  If 0, an optimal number is chosen
/// \return 0 on success.
#include <stdio.h>

// V_pot needs to be an array of i,j,k so each value can be saved. n is working kind of like time?? 
// Average the neighbours of the current point, and save it in V_pot, then loop again n times.

//TODO in each for loop branch, set boundry checks (not the n one cause it's time lololol)

int poisson_dirichlet (double * __restrict__ source, double *__restrict__ potential,
                       unsigned int xsize,  unsigned int ysize, unsigned int zsize,
                       double delta, unsigned int numiters, unsigned int numcores)
{
    // source[i, j, k] is accessed with source[((k * ysize) + j) * xsize + i]
  
	double *previous = potential;
	
	for(unsigned int n = 1;n < numiters;n++){
		//printf("-----next_iter-----\n");
			for(unsigned int i = 1; i < xsize - 1; i++){			
				for(unsigned int j = 1; j < ysize - 1; j++){	
							for(unsigned int k = 1; k < zsize - 1; k++){			
										
					potential[((k * ysize) + j) * xsize + i] =  previous[((k * ysize) + j) * xsize + (i+1)];
					potential[((k * ysize) + j) * xsize + i] += previous[((k * ysize) + j) * xsize + (i-1)];
					potential[((k * ysize) + j) * xsize + i] += previous[((k * ysize) + (j+1)) * xsize + i];
					potential[((k * ysize) + j) * xsize + i] += previous[((k * ysize) + (j-1)) * xsize + i];
					potential[((k * ysize) + j) * xsize + i] += previous[(((k+1) * ysize) + j) * xsize + i];
					potential[((k * ysize) + j) * xsize + i] += previous[(((k-1) * ysize) + j) * xsize + i];
					potential[((k * ysize) + j) * xsize + i] -= ((delta * delta) * source[((k * ysize) + j) * xsize + i]);
					potential[((k * ysize) + j) * xsize + i] = (potential[((k * ysize) + j) * xsize + i]) / 6.0;
					//printf("%d %d %d %f\n", i, j ,k, potential[((k * ysize) + j) * xsize + i]);
				}	
			}
		}
		previous = potential;
	}   
    return 0;
}



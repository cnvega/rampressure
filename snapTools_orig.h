/**
  \mainpage SnapTools
  Code developed by Cristian Vega cnvega(at)gmail.com
  under the GNU-GPL licence for working with snapshots.
*/

#ifndef _SNAPTOOLS_H
#define _SNAPTOOLS_H

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/**
  \file snapTools.h
  \author Cristian Vega
  \brief Several tools for working with snapshots.

  Available compilation flags:
  - ST_NOTIF_ON: Enables notification while loading files and allocating memory.
  - ST_SKIP_VEL: Doesn't load velocities.
 */

/**
  \brief Header-type struct. 
 */
typedef struct {
	int      npart[6];
	double   mass[6];
	double   time;
	double   redshift;
	int      flag_sfr;
	int      flag_feedback;
	int      npartTotal[6];
	int      flag_cooling;
	int      num_files;
	double   BoxSize;
	double   Omega0;
	double   OmegaLambda;
	double   HubbleParam;
	char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} HEADER;


/**
  \brief This prints in stdout the header loaded.
  \param *header Pointer to the header to print in stdout.
 */
void printHeader(HEADER *header);


/**
  \brief Returns offset of particles in the arrays.
  \param *headers Pointer to the list of headers to calculate offset.
  \param k particle type: 0..5
  \param i number of file: 0..N
  \return First index for finding particles of \e k in \ i.
 */
int calculate_offset(HEADER *headers, int k, int i);


/**
  \brief Returns offset of particles with masses in the arrays.
  \param *headers Pointer to the lit of headers to calculate offset.
  \param k particle type: 0..5
  \param i number of file: 0..N
  \return First index for finding particles of \e k in \ i.
 */
int calculate_offset_withmasses(HEADER *headers, int k, int i);


/**
  \brief Fragment a header into fraNum fragments.
  \param *header Pointer to header to fragment.
  \param fraNum Number of parts in which the header will be fragment.
  \return Pointer to an array of \e fracNum headers.
 */
HEADER * fragmentHeader(HEADER *header, int fraNum);


/**
  \brief Struct to storage particle positions and velocities.
 */
typedef struct {
	float  Pos[3];
	#ifndef ST_SKIP_VEL
	float  Vel[3];
	#endif
} PARTICLE;


/** 
  \brief Set position and velocity in a partcile.
  \param *part Pointer to particle for setting.
  \param pos Array with position values.
  \param vel Array with velocities values.
 */
void setParticle(PARTICLE *part, float pos[3], float vel[3]);


/**
  \brief Get position of a particle 
  \param *part Pointer to particle to extract positions array.
  \return Pointer to a new array with the positions of the particle.
 */
float * getPos(PARTICLE *part);


#ifndef ST_SKIP_VEL
/**
  \brief Get velocity of a particle.
  \param *part Pointer to particle to extract velocities array.
  \return Pointer to a new array with the velocities of the particle.
 */
float * getVel(PARTICLE *part);
#endif


/**
  \brief Copy particle source into destiny.
  \param destiny Pointer to the destination particle.
  \param source Pointer to the source data particle.
 */
void copyParticle(PARTICLE *destiny, PARTICLE *source);


/**
  \brief Struct for loading one snapshot file.
 */
typedef struct {
	/** Number of files in the snapshot. */
	int numFiles;
	/** List of file headers in with the snapshot is divided. */
	HEADER *headers; 
	/** Array with all the particles of the snapshot. */
	PARTICLE *particles;
	/** Array with all the IDs of the snapshot. */
	unsigned int *IDs;
   /** Array with particle masses (only if variable). */
	float *Mass;
   /** Array of energies (only gas particles). */
	float *Energy;  
   /** Array of densities (only gas particles). */
	float *Density; 
   /** Array of hsml (only gas particles). */
	float *HSML;
	/** IDs order tag. Is set 1 if yes, 0 if not. */
	int order;
} SNAPSHOT;


/** 
  \brief Constructor: Sets pointers to NULL.
  \return Pointer to a new snapshot.
 */
SNAPSHOT * createSnapshot();


/** 
  \brief Destructor: Frees memory used by the snapshot.
 */
void freeSnapshot(SNAPSHOT *snapshot);


/**
  \brief Load a fragmented file and returns a new snapshot
  \param *filename String containing the path and mean name of the snapshot file.
  \param files Number of files in with the snapshot is divided.
  \return Pointer to a new snapshot containing all the data loaded.
 */
SNAPSHOT * loadSnapshot(char *filename, int files);

/**
  \brief Load a fragment of a file and returns a new snapshot from it.
  \param *filename String containing the path and mean name of the snapshot file.
  \param file Number of fragment for loading as unique snapshot.
  \return Pointer to a new snapshot containing all the data loaded.
 */
SNAPSHOT * loadSnapshotFragment(char *filename, int file);

/**
  \brief Read the header of a snapshot file and returns a HEARDER type struct.
  \param *filename String containing the whole name of the snapshot file.
  \return Header data.
  */
HEADER * loadHeaderOfFile(char *filename);


#ifndef ST_SKIP_VEL
/**
  \brief Save the snapshot into filename.
  \param *snapshot Pointer to snapshot with data.
  \param *filename File name in with snapshot will be saved.
  \return 1 if an error ocurred, 0 if it was done.
 */
int saveSnapshot(SNAPSHOT *snapshot, char *filename);
#endif


/**
  \brief Rerturns memory used by snapshot.
  \param *snapshot Pointer to snapshot with data.
  \return Size in bytes used by snapshot in memory.
 */
double memStat(SNAPSHOT *snapshot);


/**
  \brief Combines all fragmented files into one virtual snapshot.
  \param *snapshot Pointer to snapshot with data.
 */
void unifySnapshot(SNAPSHOT *snapshot);


/**
  \brief Fragments a snapshot into fracNum number of files.
  \param *snapshot Pointer to snapshot with data.
  \param fracNum Number of files in which snapshot will be fragmented.
 */
void fragmentSnapshot(SNAPSHOT *snapshot, int fracNum);


/**
  \brief Prints in stdout the header of a snapshot file.
  \param *snapshot Pointer to snapshot with data.
  \param fileNum Number of the snapshot part that it will be printed. 
 */
void printSnapshotHeader(SNAPSHOT *snapshot, int fileNum);


/**
  \brief Total number of particles.
  \param *snapshot Pointer to snapshot with data
  \return Total number of particles calculated from header.
 */
int nParts(SNAPSHOT *snapshot);


/**
  \brief Total number of particles in file.
  \param *snapshot Pointer to snapshot with data
  \param file File for calculating number of particles.
  \return Total number of particles calculated from header.
 */
int fParts(SNAPSHOT *snapshot, int file);


/**
  \brief Number of particles of a particular type.
  \param *snapshot Pointer to snapshot with data.
  \param *type Particle type.
  \return Number of particles of particular type calculated from header.
 */
int tParts(SNAPSHOT *snapshot, int type);


/**
  \brief Number of particles of a particular type in a file.
  \param *snapshot Pointer to snapshot with data.
  \param file File for calculating number of particles.
  \param *type Particle type.
  \return Number of particles of particular type calculated from header.
 */
int ftParts(SNAPSHOT *snapshot, int file, int type);


/**
  \brief Mass of a particular particle type.
  \param *snapshot Pointer to snapshot with data.
  \param *type Particle type.
  \return Particle mass of a particular type.
 */
double massOfType(SNAPSHOT *snapshot, int type);


/**
  \brief Snapshot time.
  \param *snapshot Pointer to snapshot with data.
  \return Time of the snapshot saved in header.
 */
double get_time(SNAPSHOT *snapshot);


/**
  \brief Snapshot redshift.
  \param *snapshot Pointer to snapshot with data.
  \return Redshift of the snapshot saved in header.
 */
double redshift(SNAPSHOT *snapshot);


/**
  \brief Number of parts in which snapshot is divided.
  \param *snapshot Pointer to snapshot with data.
  \return Number of parts.
 */
int nFiles(SNAPSHOT *snapshot);


/**
  \brief Snapshot box size.
  \param *snapshot Pointer to snapshot with data.
  \return Box size.
 */
double boxSize(SNAPSHOT *snapshot);


/**
  \brief Snapshot Omega Matter.
  \param *snapshot Pointer to snapshot with data.
  \return Omega Matter.
 */
double omegaM(SNAPSHOT *snapshot);


/**
  \brief Snapshot Omega Lambda.
  \param *snapshot Pointer to snapshot with data.
  \return OmegaLambda.
 */
double omegaL(SNAPSHOT *snapshot);


/**
  \brief Snapshot Hubble constant.
  \param *snapshot Pointer to snapshot with data.
  \return Hubble constant.
 */
double hubble(SNAPSHOT *snapshot);


/**
  \brief Searchs the particle ID in the array and returns the index.
  \param *snapshot Pointer to snapshot with data.
  \param ID Tag of the desired particle.
  \return Index of particle in the array if it exists or nPart if don't.

  This function calculates the index using the minimum ID if the orderIDs()
  function has been executed before or perform an element-by-element searching
  of the particle if hasn't.
 */
unsigned int getIndex(SNAPSHOT *snapshot, unsigned int ID);


/**
  \brief Returns a pointer in the array to the desired particle.
  \param *snapshot Pointer to snapshot with data.
  \param ID Tag of the desired particle.
  \return Pointer to the particle.
  */
PARTICLE * pointerToParticle(SNAPSHOT *snapshot, unsigned int ID);


/**
  \brief Returns the desired particle.
  \param *snapshot Pointer to snapshot with data.
  \param ID Tag of the desired particle.
  \return A copy of the desired particle.
  */
PARTICLE getParticle(SNAPSHOT *snapshot, unsigned int ID);


/**
  \brief Set snapshot time.
  \param *snapshot Pointer to snapshot with data.
  \param Time to save in the snapshot header.

  This function also save redshift with z=1/a-1.
 */
void setTime(SNAPSHOT *snapshot, double a);


/**
  \brief Set snapshot redshift.
  \param *snapshot Pointer to snapshot with data.
  \param Redshift to save in the snapshot header.

  This function also save redshift with a=1/(1+z).
 */
void setRedshift(SNAPSHOT *snapshot, double z);


/**
  \brief Set snapshot box size.
  \param *snapshot Pointer to snapshot with data.
  \param Box size.
 */
void setBoxSize(SNAPSHOT *snapshot, double box);


/**
  \brief Set snapshot cosmology
  \param omegaM Matter density at z=0 in units of the critical density.
  \param omegaL Vacuum energy density at z=0 in units of the critical density.
  \param h The Hubble constant in units of 100 km/s/Mpc.
  */
void setCosmology(SNAPSHOT *snapshot, double omegaM, double omegaL, double h);


/**
  \brief Creates and initializes a new empty snapshot.
  \return Pointer to the new created snapashot.
  */
SNAPSHOT * newSnapshot();


/**
  \brief Flip particles A and B.
  \param *snapshot Pointer to snapshot with data.
  \param A Index of the particle to flip with B.
  \param B Index of the particle to flip with A.

  This function only works if the snapshot has just one type of particles and there
  are not gas particles.
 */
void flipParticles(SNAPSHOT *snapshot, int A, int B);


/**
  \brief Order the particles array according to IDs.
  \param *snapshot Pointer to snapshot with data.
  \return 1 if an error ocurred, 0 if it was done.
  
  This function only works if the snapshot has just one type of particles and there
  are not gas particles.
 */
int orderIDs(SNAPSHOT *snapshot);


/**
  \brief Adds one type particles to the snapashot.
  \param *snapshot Pointer to snapshot with data.
  \param *vector Array of particles for adding to snapshot.
  \param Npart Number of particles in the array.
  \param type Type of particles for adding particles.
  \param mass Mass of each particle.
  \return Pointer to the new array of particles added. After the first called of this
          function the original array is moved.
  
  This function only works if these are not gas particles and there is not mass
  array.
 */
PARTICLE * addParticles(SNAPSHOT *snapshot, PARTICLE *vector, 
                        int Npart, int type, float mass);


/**
  \brief Adds one all type of particles to the snapashot.
  \param *snapshot Pointer to snapshot with data.
  \param *vector Array of particles for adding to snapshot.
  \param Npart Number of particles of each type in the array.
  \param mass Mass of particles in each type.
  \return Pointer to the new array of particles.
  
  This function only works if these are not gas particles and there is not mass
  array.
 */
PARTICLE * addAllParticles(SNAPSHOT *snapshot, PARTICLE *vector, 
									int Npart[6], float mass[6]);


/**
  \brief Read file and write summary into file.
  \param *input String with the input file name.
  \param *output String with the output text file name.
  \return 1 if an error ocurred, 0 if it was done.
 */
int read_debug(char *input, char *output);

	
#endif


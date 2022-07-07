/*
 * Code developed by Cristian Vega-Martinez
 * under the GNU-GPL licence.
 */

#ifndef _SNAPTOOLS_CPP
#define _SNAPTOOLS_CPP

/* Include section */
#include "snapTools.h"

/* This prints in stdout the header loaded */
void Header::printHeader(Header *header)
{
	printf("\nSnapshot Header\n");
	printf("Npart:\n");
	printf("     Gas     = %d\n", header->npart[0]);
	printf("     Halo    = %d\n", header->npart[1]);
	printf("     Disk    = %d\n", header->npart[2]);
	printf("     Bulge   = %d\n", header->npart[3]);
	printf("     Stars   = %d\n", header->npart[4]);
	printf("     Bndry   = %d\n", header->npart[5]);
	printf("Mass: \n");
	printf("     Gas     = %g\n", header->mass[0]); 
	printf("     Halo    = %g\n", header->mass[1]); 
	printf("     Disk    = %g\n", header->mass[2]); 
	printf("     Bulge   = %g\n", header->mass[3]); 
	printf("     Stars   = %g\n", header->mass[4]); 
	printf("     Bndry   = %g\n", header->mass[5]); 	
	printf("Time            = %g\n", header->time);
	printf("Redshift        = %g\n", header->redshift);
	printf("flag_sfr        = %d\n", header->flag_sfr);
	printf("flag_feedback   = %d\n", header->flag_feedback);
	printf("NpartTotal:\n");
	printf("     Gas     = %d\n", header->npartTotal[0]);
	printf("     Halo    = %d\n", header->npartTotal[1]);
	printf("     Disk    = %d\n", header->npartTotal[2]);
	printf("     Bulge   = %d\n", header->npartTotal[3]);
	printf("     Stars   = %d\n", header->npartTotal[4]);
	printf("     Bndry   = %d\n", header->npartTotal[5]);
	printf("flag_cooling    = %d\n", header->flag_cooling);
	printf("numfiles        = %d\n", header->num_files);
	printf("boxsize         = %g\n", header->BoxSize);
	printf("Omega0          = %g\n", header->Omega0);
	printf("OmegaLambda     = %g\n", header->OmegaLambda);
	printf("Hubbleparam     = %g\n", header->HubbleParam);
}

/* Returns offset of particles in the arrays*/
int Header::calculate_offset(Header *headers, int pnum, int fnum)
{
	int offset = 0;
	for(int i = 0; i < pnum; i++)
		offset += headers[0].npartTotal[i];
	
	for(int i = 0; i < fnum; i++)
		offset += headers[i].npart[fnum];
	
	return offset;
}

/* Returns offset of particles with masses in the arrays*/
int Header::calculate_offset_withmasses(Header *headers, int pnum, int fnum)
{
	int offset = 0;
	for(int i = 0; i < pnum; i++)
		if(headers[0].mass[i] == 0)
			offset += headers[0].npartTotal[i];
	for(int i = 0; i < fnum; i++)
		offset += headers[i].npart[pnum];
	
	return offset;
}


/* Copy particle source into destiny */
void copyParticle(Particle *dst, Particle *src)
{
	int i;
	for(i=0; i<3; i++)
	{
		dst->Pos[i] = src->Pos[i];
		#ifndef PARTICLE_SKIP_VEL
		dst->Vel[i] = src->Vel[i];
		#endif
	}
}


/* Constructor: Sets pointers to NULL: */
Snapshot::Snapshot()
{
	this->numFiles = 0;
	this->headers = NULL;
	this->particles = NULL;
	this->IDs = NULL;
	this->Mass = NULL;
	this->Energy = NULL;
	this->Density = NULL;
	this->HSML = NULL;
	this->order = 0;
}

Snapshot::~Snapshot()
{
	free(this->headers);
	free(this->particles);
	free(this->IDs);
	free(this->Mass);
	free(this->Energy);
	free(this->Density);
	free(this->HSML);
}


/* Load a fragmented file and returns a pointer to snapshot*/
void Snapshot::loadSnapshot(char *fname, int nfiles)
{
	this->numFiles = nfiles;
	this->headers = (Header**) malloc(nfiles*sizeof(Header));
	for (int i = 0; i < nfiles; i++)
	{
		this->headers[i] = new Header();	
		// Reads over all the fragmented files:
		
		FILE *fd;
		char buf[255];
		if(nfiles > 1)
			sprintf(buf,"%s.%d", fname, i);
		else
			sprintf(buf,"%s",fname);

		if(!(fd=fopen(buf,"r")))
		{
			fprintf(stderr, "Can't open file `%s`\n",buf);
			return;
		}
		
		#ifdef ST_NOTIF_ON
		printf("Reading file: '%s'\n", buf);
		#endif

		fread(&skip, sizeof(skip), 1, fd);
		fread(&hds[i], sizeof(HEADER), 1, fd);
		fread(&skip, sizeof(skip), 1, fd);

	}



	
	int NumPart=0, Ngas=0, Nmasses=0;
	#ifdef PARTICLE_SKIP_VEL
	float skip_vel[3];
	#endif


	for(i=0; i<nfiles; i++)
   { 

		if(i == 0)
		{
		// Counting particles:
			for(k=0; k<5; k++)
			{
				NumPart += hds[0].npartTotal[k]; // all
				if(hds[0].mass[k] == 0)
					Nmasses+= hds[0].npartTotal[k]; // with masses
			}
			Ngas = hds[0].npartTotal[0];        // gas
		// Allocating memory:
			#ifdef ST_NOTIF_ON
			printf("\tAllocating memory...   "); fflush(stdout);
			#endif
			snapshot->headers = hds;
			if(!(snapshot->particles = malloc(NumPart*sizeof(PARTICLE)))
			|| !(snapshot->IDs = malloc(NumPart*sizeof(unsigned int))))
			{
				fprintf(stderr, "\nFailed to allocate memory: Particles Array\n");
				return NULL;
			}
			if(Nmasses != 0)
				if(!(snapshot->Mass =	malloc(Nmasses*sizeof(float))))
				{
					fprintf(stderr, "\nFailed to allocate memory: Mass Array\n");
					return NULL;
				}
			if(Ngas != 0)
				if(!(snapshot->Energy  = malloc(Ngas*sizeof(float)))
				|| !(snapshot->Density = malloc(Ngas*sizeof(float)))
				|| !(snapshot->HSML    = malloc(Ngas*sizeof(float))))
				{
					fprintf(stderr, "\nFailed to allocate memory: Gas Properties Arrays\n");
					return NULL;
				}
			#ifdef ST_NOTIF_ON
			printf("Done!\n");
			#endif
		}

		fread(&skip, sizeof(skip), 1, fd);

	// Reading the positions:
		#ifdef ST_NOTIF_ON
		printf("\tLoading positions...   "); fflush(stdout);
		#endif
		for(k=0; k<6; k++)
		{
			offset = calculate_offset(hds, k, i);
			for(n=0; n<hds[i].npart[k]; n++)
				fread(&snapshot->particles[offset+n].Pos[0],sizeof(float),3,fd);
		}
		#ifdef ST_NOTIF_ON
		printf("Done!\n");
		#endif
		fread(&skip, sizeof(skip), 1, fd);
		fread(&skip, sizeof(skip), 1, fd);
	// Reading the velocities:
		#ifdef ST_NOTIF_ON
		printf("\tLoading velocities...  "); fflush(stdout);
		#endif
		for(k=0; k<6; k++)
		{
			offset = calculate_offset(hds, k, i);
			for(n=0; n<hds[i].npart[k]; n++)
			{
				#ifdef PARTICLE_SKIP_VEL
				fread(&skip_vel[0],sizeof(float),3,fd);
				#else
				fread(&snapshot->particles[offset+n].Vel[0],sizeof(float),3,fd);
				#endif
			}
		}
		#ifdef ST_NOTIF_ON
		printf("Done!\n");
		#endif
		fread(&skip, sizeof(skip), 1, fd);
		fread(&skip, sizeof(skip), 1, fd);
	// Reading the IDs:
		#ifdef ST_NOTIF_ON
		printf("\tLoading IDs...         "); fflush(stdout);
		#endif
		for(k=0; k<6; k++)
		{
			offset = calculate_offset(hds, k, i);
			for(n=0; n<hds[i].npart[k]; n++)
				fread(&snapshot->IDs[offset+n], sizeof(unsigned int), 1, fd);
		}
		#ifdef ST_NOTIF_ON
		printf("Done!\n");
		#endif
		fread(&skip, sizeof(skip), 1, fd);
	// Reading the masses:
		if(Nmasses != 0)
		{
			#ifdef ST_NOTIF_ON
			printf("\tLoading masses...      "); fflush(stdout);
			#endif
			fread(&skip, sizeof(skip), 1, fd);
			for(k=0; k<6; k++)
			{
				offset = calculate_offset_withmasses(hds, k, i);
				for(n=0; n<hds[i].npart[k]; n++)
					fread(&snapshot->Mass[offset+n], sizeof(float), 1, fd);
			}
			#ifdef ST_NOTIF_ON
			printf("Done!\n");
			#endif
			fread(&skip, sizeof(skip), 1, fd);
		}
	// Reading gas properties:
		if(Ngas != 0)
		{
			fread(&skip, sizeof(skip), 1, fd);
		// Energy:
			#ifdef ST_NOTIF_ON
			printf("\tLoading gas energy...  "); fflush(stdout);
			#endif
			offset = calculate_offset(hds, 0, i);
			for(n=0; n<hds[i].npart[0]; n++)
				fread(&snapshot->Energy[offset+n], sizeof(float), 1, fd);
			#ifdef ST_NOTIF_ON
			printf("Done!\n");
			#endif
			fread(&skip, sizeof(skip), 1, fd);
			fread(&skip, sizeof(skip), 1, fd);
		// Density:
			#ifdef ST_NOTIF_ON
			printf("\tLoading gas density... "); fflush(stdout);
			#endif
			offset = calculate_offset(hds, 0, i);
			for(n=0; n<hds[i].npart[0]; n++)
				fread(&snapshot->Density[offset+n], sizeof(float), 1, fd);
			#ifdef ST_NOTIF_ON
			printf("Done!\n");
			#endif
			fread(&skip, sizeof(skip), 1, fd);
			fread(&skip, sizeof(skip), 1, fd);
		// HSML:
			#ifdef ST_NOTIF_ON
			printf("\tLoading gas hsml...    "); fflush(stdout);
			#endif
			offset = calculate_offset(hds, 0, i);
			for(n=0; n<hds[i].npart[0]; n++)
				fread(&snapshot->HSML[offset+n], sizeof(float), 1, fd);
			#ifdef ST_NOTIF_ON
			printf("Done!\n");
			#endif
			fread(&skip, sizeof(skip), 1, fd);
		}
		fclose(fd);

	}

	return snapshot;
}

/* Load a fragment of a file and returns a new snapshot from it*/
SNAPSHOT * loadSnapshotFragment(char *filename, int file)
{
	FILE *fd;
	char buf[200];
	HEADER *hds;
	SNAPSHOT *snapshot;
	int NumPart=0, Ngas=0, Nmasses=0;
	int k, n; //counters
	int offset=0;
	int skip;
	#ifdef PARTICLE_SKIP_VEL
	float skip_vel[3];
	#endif

	snapshot = createSnapshot();
	snapshot->numFiles = 1;

	hds = malloc(sizeof(HEADER));

// Reads over all the fragmented files:
	sprintf(buf,"%s.%d",filename,file);

	if(!(fd=fopen(buf,"r")))
	{
		fprintf(stderr, "Can't open file `%s`\n",buf);
		return NULL;
	}

	#ifdef ST_NOTIF_ON
	printf("Reading file: '%s'\n", buf);
	#endif

	fread(&skip, sizeof(skip), 1, fd);
	fread(&hds[0], sizeof(HEADER), 1, fd);
	fread(&skip, sizeof(skip), 1, fd);

// Counting particles:
	for(k=0; k<5; k++)
	{
		NumPart += hds[0].npart[k]; // all
		if(hds[0].mass[k] == 0)
			Nmasses+= hds[0].npart[k]; // with masses
	}
	Ngas = hds[0].npart[0];        // gas
// Allocating memory:
	#ifdef ST_NOTIF_ON
	printf("\tAllocating memory...   "); fflush(stdout);
	#endif
	snapshot->headers = hds;
	if(!(snapshot->particles = malloc(NumPart*sizeof(PARTICLE)))
		|| !(snapshot->IDs = malloc(NumPart*sizeof(unsigned int))))
	{
		fprintf(stderr, "\nFailed to allocate memory: Particles Array\n");
		return NULL;
	}
	if(Nmasses != 0)
		if(!(snapshot->Mass =	malloc(Nmasses*sizeof(float))))
		{
			fprintf(stderr, "\nFailed to allocate memory: Mass Array\n");
			return NULL;
		}
	if(Ngas != 0)
		if(!(snapshot->Energy  = malloc(Ngas*sizeof(float)))
			|| !(snapshot->Density = malloc(Ngas*sizeof(float)))
			|| !(snapshot->HSML    = malloc(Ngas*sizeof(float))))
		{
			fprintf(stderr, "\nFailed to allocate memory: Gas Properties Arrays\n");
			return NULL;
		}
	#ifdef ST_NOTIF_ON
	printf("Done!\n");
	#endif

	fread(&skip, sizeof(skip), 1, fd);

// Reading the positions:
	#ifdef ST_NOTIF_ON
	printf("\tLoading positions...   "); fflush(stdout);
	#endif
	for(k=0, offset=0; k<6; k++)
	{
		for(n=0; n<hds[0].npart[k]; n++)
			fread(&snapshot->particles[offset+n].Pos[0],sizeof(float),3,fd);
		offset += hds[0].npart[k];
	}
	#ifdef ST_NOTIF_ON
	printf("Done!\n");
	#endif
	fread(&skip, sizeof(skip), 1, fd);
	fread(&skip, sizeof(skip), 1, fd);
// Reading the velocities:
	#ifdef ST_NOTIF_ON
	printf("\tLoading velocities...  "); fflush(stdout);
	#endif
	for(k=0, offset=0; k<6; k++)
	{
		for(n=0; n<hds[0].npart[k]; n++)
		{
			#ifdef PARTICLE_SKIP_VEL
			fread(&skip_vel[0],sizeof(float),3,fd);
			#else
			fread(&snapshot->particles[offset+n].Vel[0],sizeof(float),3,fd);
			#endif
		}
		offset += hds[0].npart[k];
	}
	#ifdef ST_NOTIF_ON
	printf("Done!\n");
	#endif
	fread(&skip, sizeof(skip), 1, fd);
	fread(&skip, sizeof(skip), 1, fd);
// Reading the IDs:
	#ifdef ST_NOTIF_ON
	printf("\tLoading IDs...         "); fflush(stdout);
	#endif
	for(k=0, offset=0; k<6; k++)
	{
		for(n=0; n<hds[0].npart[k]; n++)
			fread(&snapshot->IDs[offset+n], sizeof(unsigned int), 1, fd);
		offset += hds[0].npart[k];
	}
	#ifdef ST_NOTIF_ON
	printf("Done!\n");
	#endif
	fread(&skip, sizeof(skip), 1, fd);
// Reading the masses:
	if(Nmasses != 0)
	{
		#ifdef ST_NOTIF_ON
		printf("\tLoading masses...      "); fflush(stdout);
		#endif
		fread(&skip, sizeof(skip), 1, fd);
		for(k=0, offset=0; k<6; k++)
		{
			for(n=0; n<hds[0].npart[k]; n++)
				fread(&snapshot->Mass[offset+n], sizeof(float), 1, fd);
			offset += hds[0].npart[k];
		}
		#ifdef ST_NOTIF_ON
		printf("Done!\n");
		#endif
		fread(&skip, sizeof(skip), 1, fd);
	}
// Reading gas properties:
	if(Ngas != 0)
	{
		fread(&skip, sizeof(skip), 1, fd);
	// Energy:
		#ifdef ST_NOTIF_ON
		printf("\tLoading gas energy...  "); fflush(stdout);
		#endif
		for(n=0; n<hds[0].npart[0]; n++)
			fread(&snapshot->Energy[n], sizeof(float), 1, fd);
		#ifdef ST_NOTIF_ON
		printf("Done!\n");
		#endif
		fread(&skip, sizeof(skip), 1, fd);
		fread(&skip, sizeof(skip), 1, fd);
	// Density:
		#ifdef ST_NOTIF_ON
		printf("\tLoading gas density... "); fflush(stdout);
		#endif
		for(n=0; n<hds[0].npart[0]; n++)
			fread(&snapshot->Density[n], sizeof(float), 1, fd);
		#ifdef ST_NOTIF_ON
		printf("Done!\n");
		#endif
		fread(&skip, sizeof(skip), 1, fd);
		fread(&skip, sizeof(skip), 1, fd);
	// HSML:
		#ifdef ST_NOTIF_ON
		printf("\tLoading gas hsml...    "); fflush(stdout);
		#endif
		for(n=0; n<hds[0].npart[0]; n++)
			fread(&snapshot->HSML[n], sizeof(float), 1, fd);
		#ifdef ST_NOTIF_ON
		printf("Done!\n");
		#endif
		fread(&skip, sizeof(skip), 1, fd);
	}
	fclose(fd);

	return snapshot;
}


/* Read the header of a snapshot file and returns a HEARDER type struct. */
HEADER * loadHeaderOfFile(char *filename)
{
	HEADER *header;
	FILE *fd;
	int skip;

	header = malloc(sizeof(HEADER));

	if(!(fd=fopen(filename,"r")))
	{
		fprintf(stderr, "Can't open file `%s`\n",filename);
		return NULL;
	}

	fread(&skip, sizeof(skip), 1, fd);
	fread(header, sizeof(HEADER), 1, fd);

	fclose(fd);

	return header;
}


#ifndef PARTICLE_SKIP_VEL
/* Save the snapshot into filename */
int saveSnapshot(SNAPSHOT *snapshot, char *filename)
{
	int i, k, n;
	char fname[255];
	FILE *output;
	int block=256;
	int NumPart=0, Ngas=0, Nmasses=0;
	int offset;

// Writing the new files:
	for(i=0; i<snapshot->headers[0].num_files; i++)
   { 
		if(snapshot->headers[0].num_files > 1)
			sprintf(fname,"%s.%d", filename, i);
		else
			sprintf(fname,"%s", filename);

		if(!(output=fopen(fname,"w")))
		{
			fprintf(stderr, "Can't write file `%s`\n",fname);
			return 1;
		}
		#ifdef ST_NOTIF_ON
		printf("Writing file: '%s'\n", fname);
		#endif
		// Counting particles:
		NumPart = 0;
		Nmasses = 0;
		Ngas = 0;
		for(k=0; k<5; k++)
		{
			NumPart += snapshot->headers[i].npart[k]; // all
			if(snapshot->headers[i].mass[k] == 0)
				Nmasses+= snapshot->headers[i].npart[k]; // with masses
		}
		Ngas = snapshot->headers[i].npart[0];        // gas

		// writing...
		fwrite(&block, sizeof(block), 1, output);
		fwrite(&snapshot->headers[i], sizeof(HEADER), 1, output);
		fwrite(&block, sizeof(block), 1, output);

    	block=NumPart*sizeof(float)*3; // pos block size

		fwrite(&block, sizeof(block), 1, output);
		for(k=0; k<6; k++)
		{
			offset = calculate_offset(snapshot->headers, k, i);
			for(n=0; n<snapshot->headers[i].npart[k]; n++)
				fwrite(snapshot->particles[offset+n].Pos,sizeof(float),3,output);
		}
		fwrite(&block, sizeof(block), 1, output); // pos block size

		fwrite(&block, sizeof(block), 1, output); // vel block size
		for(k=0; k<6; k++)
		{
			offset = calculate_offset(snapshot->headers, k, i);
			for(n=0; n<snapshot->headers[i].npart[k]; n++)
				fwrite(snapshot->particles[offset+n].Vel,sizeof(float),3,output);
		}
		fwrite(&block, sizeof(block), 1, output); // vel block size

		block = NumPart*sizeof(unsigned int); 
		fwrite(&block, sizeof(block), 1, output); // Ids block size
		for(k=0; k<6; k++)
		{
			offset = calculate_offset(snapshot->headers, k, i);
			for(n=0; n<snapshot->headers[i].npart[k]; n++)
				fwrite(&snapshot->IDs[offset+n], sizeof(unsigned int), 1, output);
		}
		fwrite(&block, sizeof(block), 1, output); // Ids block size

		if(Nmasses != 0)
		{
    		block = Nmasses*sizeof(float); 
			fwrite(&block, sizeof(block), 1, output); // masses block size
			for(k=0; k<6; k++)
			{
				offset = calculate_offset_withmasses(snapshot->headers, k, i);
				for(n=0; n<snapshot->headers[i].npart[k]; n++)
					fwrite(&snapshot->Mass[offset+n], sizeof(float), 1, output);
			}
			fwrite(&block, sizeof(block), 1, output); // masses block size
		}

		if(Ngas != 0)
		{
    		block = Ngas*sizeof(float);
			fwrite(&block, sizeof(block), 1, output); // Energies block size
			offset = calculate_offset(snapshot->headers, 0, i);
			for(n=0; n<snapshot->headers[i].npart[0]; n++)
				fwrite(&snapshot->Energy[offset+n], sizeof(float), 1, output);
			fwrite(&block, sizeof(block), 1, output); // Energies block size
			fwrite(&block, sizeof(block), 1, output); // Densities block size
			offset = calculate_offset(snapshot->headers, 0, i);
			for(n=0; n<snapshot->headers[i].npart[0]; n++)
				fwrite(&snapshot->Density[offset+n], sizeof(float), 1, output);
			fwrite(&block, sizeof(block), 1, output); // Densities block size
			fwrite(&block, sizeof(block), 1, output); // HSML block size
			offset = calculate_offset(snapshot->headers, 0, i);
			for(n=0; n<snapshot->headers[i].npart[0]; n++)
				fread(&snapshot->HSML[offset+n], sizeof(float), 1, output);
			fwrite(&block, sizeof(block), 1, output); // HSML block size
		}
		fclose(output);

	}
	return 0;
}
#endif


/* Rerturns memory used by a snapshot */
double memStat(SNAPSHOT *snapshot)
{
	double total=0.;
	int k;
	int NumPart=0, Ngas=0, Nmasses=0;

	for(k=0; k<5; k++)
	{
		NumPart += snapshot->headers[0].npartTotal[k]; // all
		if(snapshot->headers[0].mass[k] == 0)
			Nmasses+= snapshot->headers[0].npartTotal[k]; // with masses
	}
	Ngas = snapshot->headers[0].npartTotal[0];        // gas

	total += (snapshot->numFiles)*sizeof(HEADER);
	total += NumPart*sizeof(PARTICLE);
	total += NumPart*sizeof(unsigned int);
	total += Nmasses*sizeof(float);
	total += 3*Ngas*sizeof(float);

	#ifdef ST_NOTIF_ON
	printf("Memory used: %.1f Mb\n", total/1024.0/1024.0 );
	#endif

	return total;
}


/* Combines all fragmented files into one virtual snapshot */
void unifySnapshot(SNAPSHOT *snapshot)
{
	HEADER *newheader;
	int k;

	newheader = malloc(sizeof(HEADER));

	newheader[0] = snapshot->headers[0];

	for(k=0; k<6; k++)
		newheader[0].npart[k] = snapshot->headers[0].npartTotal[k];

	newheader[0].num_files = 1;

	free(snapshot->headers);
	snapshot->headers = newheader;
	
	snapshot->numFiles = 1;
}


/* Fragments a snapshot into fracNum number of files */
void fragmentSnapshot(SNAPSHOT *snapshot, int fracNum)
{
	HEADER *aux;

	snapshot->numFiles = fracNum;
	aux = fragmentHeader(snapshot->headers, fracNum);
	free(snapshot->headers);
	snapshot->headers = aux;
}

/* Prints snapshot header of a certain file */
void printSnapshotHeader(SNAPSHOT *snapshot, int fileNum)
{
	printHeader(&snapshot->headers[fileNum]);
}


/* Total number of particles. */
int nParts(SNAPSHOT *snapshot)
{
	int i, Nparts=0;
	for (i=0; i<6; i++)
		Nparts += snapshot->headers[0].npartTotal[i];
	return Nparts;
}


/* Total number of particles in file. */
int fParts(SNAPSHOT *snapshot, int file)
{
	int i, Nparts=0;
	for (i=0; i<6; i++)
		Nparts += snapshot->headers[file].npart[i];
	return Nparts;
}


/* Number of particles of a particular type. */
int tParts(SNAPSHOT *snapshot, int type)
{
	return snapshot->headers[0].npart[type];
}


/* Number of particles of a particular type in a file. */
int ftParts(SNAPSHOT *snapshot, int file, int type)
{
	return snapshot->headers[file].npart[type];
}


/* Mass of a particular particle type. */
double massOfType(SNAPSHOT *snapshot, int type)
{
	return snapshot->headers[0].mass[type];
}


/* Snapshot time. */
double get_time(SNAPSHOT *snapshot)
{
	return snapshot->headers[0].time;
}


/* Snapshot redshift. */
double redshift(SNAPSHOT *snapshot)
{
	return snapshot->headers[0].redshift;
}


/* Number of parts in which snapshot is divided. */
int nFiles(SNAPSHOT *snapshot)
{
	return snapshot->numFiles;
}


/* Snapshot box size. */
double boxSize(SNAPSHOT *snapshot)
{
	return snapshot->headers[0].BoxSize;
}


/* Snapshot Omega Matter. */
double omegaM(SNAPSHOT *snapshot)
{
	return snapshot->headers[0].Omega0;
}


/* Snapshot Omega Lambda. */
double omegaL(SNAPSHOT *snapshot)
{
	return snapshot->headers[0].OmegaLambda;
}


/* Snapshot Hubble constant. */
double hubble(SNAPSHOT *snapshot)
{
	return snapshot->headers[0].HubbleParam;
}


/* Searchs the particle ID in the array and returns the index or nParts. */
unsigned int getIndex(SNAPSHOT *snapshot, unsigned int ID)
{
	int i;
	if (snapshot->order == 1)
		return (ID - (snapshot->IDs[0]));
	else
		for (i=0; i<nParts(snapshot); i++)
			if(snapshot->IDs[i] == ID)
				return i;
	return i;
}


PARTICLE * pointerToParticle(SNAPSHOT *snapshot, unsigned int ID)
{
	PARTICLE *part;
	unsigned int index;

	index = getIndex(snapshot, ID);

	if (index == nParts(snapshot))
		return NULL;
	
	part = &snapshot->particles[ID];

	return part;
}


/* Returns the desired particle. */
PARTICLE getParticle(SNAPSHOT *snapshot, unsigned int ID)
{
	PARTICLE part;
	unsigned int index;
	float temp[3] = {0,0,0};

	index = getIndex(snapshot, ID);

	if (index == nParts(snapshot))
	{
		setParticle(&part, temp, temp);
		return part;
	}

	copyParticle(&part, &snapshot->particles[index]);

	return part;
}


/* Set snapshot time. */
void setTime(SNAPSHOT *snapshot, double a)
{
	int i;
	for (i=0; i<nFiles(snapshot); i++)
	{
		snapshot->headers[i].time = a;
		snapshot->headers[i].redshift = 1./a - 1.;
	}
}


/** Set snapshot redshift. */
void setRedshift(SNAPSHOT *snapshot, double z)
{
	int i;
	for (i=0; i<nFiles(snapshot); i++)
	{
		snapshot->headers[i].redshift = z;
		snapshot->headers[i].time = 1./(1. + z);
	}
}


/** Set snapshot box size. */
void setBoxSize(SNAPSHOT *snapshot, double box)
{
	int i;
	for (i=0; i<nFiles(snapshot); i++)
		snapshot->headers[i].BoxSize = box;
}


/* Set snapshot cosmology. */
void setCosmology(SNAPSHOT *snapshot, double omegaM, double omegaL, double h)
{
	int i;
	for (i=0; i<nFiles(snapshot); i++)
	{
		snapshot->headers[i].Omega0 = omegaM;
		snapshot->headers[i].OmegaLambda = omegaL;
		snapshot->headers[i].HubbleParam = h;
	}
}


/**
  \brief Creates and initializes a new empty snapshot.
  \return Pointer to the new created snapashot.
  */
SNAPSHOT * newSnapshot()
{
	SNAPSHOT *snapshot;
	HEADER *header;
	int i;

	snapshot = createSnapshot();
	header = malloc(sizeof(HEADER));

	for(i=0; i<6; i++)
	{
		header[0].npart[i] = 0;
		header[0].mass[i] = 0;
		header[0].npartTotal[i] = 0;
	}
	header[0].flag_sfr = 0;
	header[0].flag_feedback = 0;
	header[0].flag_cooling = 0;
	header[0].num_files = 1;

	snapshot->headers = header;
	snapshot->numFiles = 1;

	setRedshift(snapshot, 0.);
	setBoxSize(snapshot, 0.);
	setCosmology(snapshot, 0.3, 0.7, 0.7);

	return snapshot;
}


/* Flip particles A and B */
//TODO: this don't work with masses
void flipParticles(SNAPSHOT *snapshot, int A, int B)
{
	PARTICLE temp;
	unsigned int index;
	// data:
	copyParticle(&temp, &snapshot->particles[A]);
	copyParticle(&snapshot->particles[A], &snapshot->particles[B]);
	copyParticle(&snapshot->particles[B], &temp);
	// IDs:
	index = snapshot->IDs[A];
	snapshot->IDs[A] = snapshot->IDs[B];
	snapshot->IDs[B] = index;
}


/* Returns 1 if done, 0 if not */
//TODO: this don't work with masses
int orderIDs(SNAPSHOT *snapshot)
{
	int k, min, max;
	unsigned long p, np = 0;

// Has this function been executed before?
	if (snapshot->order == 1)
		return 0;
	if (snapshot->order == -1)
		return 1;

	// Counting particles:
	for(k=0; k<6; k++)
		np += snapshot->headers[0].npartTotal[k];
	
	// Checking limits:
	min = snapshot->IDs[0];
	max = snapshot->IDs[0];
	for(p=1; p<np; p++)
	{
		if(snapshot->IDs[p] < min)
			min = snapshot->IDs[p];
		if(snapshot->IDs[p] > max)
			max = snapshot->IDs[p];
	}

	// If is not possible to do it:
	if(np != (max - min + 1 ))
	{
		snapshot->order = -1;
		printf("WARNING: Particles can NOT be ordered.\n");
		return 1;
	}

	// Last option: order particles
	if (snapshot->headers[0].num_files != 1)
	{
		unifySnapshot(snapshot);
		printf("WARNING: Header has been unified.\n");
	}
	printf("\tOrdering particles...  "); fflush(stdout);
	p=0;
	while(p<np)
	{
		if(snapshot->IDs[p] != (p+min))
			flipParticles(snapshot, p, snapshot->IDs[p]-min);

		else p++;
	}
	snapshot->order = 1;
	printf("Done!\n");

	return 0;
}


/** Adds one type particles to the snapashot. */
PARTICLE * addParticles(SNAPSHOT *snapshot, PARTICLE *vector, 
                        int Npart, int type, float mass)
{
	int i;
	int NBefore = 0, NAfter = 0;
	int NOld = 0, NNew = 0;
	PARTICLE *newArray, *newType;

	if(nFiles(snapshot) != 1)
	{
		fprintf(stderr, "Cannot add particles with nFiles different than 1!\n");
		return NULL;
	}

	if (nParts(snapshot) == 0)
	{
		snapshot->particles = vector;
		snapshot->headers[0].npart[type] = Npart;
		snapshot->headers[0].npartTotal[type] = Npart;
		snapshot->headers[0].mass[type] = Npart;
		return vector;
	}
	else
	{
		for (i=0; i<=type; i++)
			NBefore += snapshot->headers[0].npartTotal[i];
		for (i=type+1; i<6; i++)
			NAfter += snapshot->headers[0].npartTotal[i];

		NOld = nParts(snapshot);
		NNew = NOld + Npart;

		newArray = malloc(NNew*sizeof(PARTICLE));
		
		memmove(newArray, snapshot->particles, NBefore*sizeof(PARTICLE));
		newType = memmove(&newArray[NBefore], vector, Npart*sizeof(PARTICLE));
		memmove(&newArray[NBefore+Npart], &snapshot->particles[NBefore], 
		        NAfter*sizeof(PARTICLE));

		snapshot->headers[0].npart[type] += Npart;
		snapshot->headers[0].npartTotal[type] += Npart;
		snapshot->headers[0].mass[type] = mass;

		snapshot->particles = newArray;

		return newType;
	}
}


/** Adds one all type of particles to the snapashot. */
PARTICLE * addAllParticles(SNAPSHOT *snapshot, PARTICLE *vector, 
									int Npart[6], float mass[6])
{
	int i;

	if(nFiles(snapshot) != 1)
	{
		fprintf(stderr, "Cannot add particles with nFiles different than 1!\n");
		return NULL;
	}

	if (nParts(snapshot) != 0)
		printf("WARNING: Overwriting %d particles in snapshot.\n", nParts(snapshot));

	snapshot->particles = vector;
	for (i=0; i<6; i++)
	{
		snapshot->headers[0].npart[i] = Npart[i];
		snapshot->headers[0].npartTotal[i] = Npart[i];
		snapshot->headers[0].mass[i] = mass[i];
	}

	return vector;
}


/* Read snapshot and write summary into file */
int read_debug(char *input, char *output)
{
	FILE *fd;
	FILE *data;
	int dummy;
	HEADER *header;
	float aux, aux3[3];
	int i, k, n; //counters
	int NumPart=0, Ngas=0, Nmasses=0;
	unsigned int auxi, auxc;

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);
#define PRINT fprintf(data, "DATA = %d\n", dummy);

	header = malloc(sizeof(HEADER));

	if(!(fd=fopen(input,"r")))
	{
		fprintf(stderr, "Can't open file `%s`\n", input);
		return 1;
	}
	#ifdef ST_NOTIF_ON
	printf("Reading file: '%s'\n", input);
	#endif

	if(!(data=fopen(output,"w")))
	{
		fprintf(stderr, "Can't open file `%s`\n", output);
		return 1;
	}
	#ifdef ST_NOTIF_ON
	printf("Writing file: '%s'\n", output);
	#endif

	fprintf(data, "#<<<<< %s >>>>>\n", input);

	SKIP; PRINT;
	fread(header, sizeof(HEADER), 1, fd);
	fprintf(data, "#<<< HEADER >>>\n");

	fprintf(data, "Npart:\n");
	fprintf(data, "     Gas     = %d\n",header->npart[0]);
	fprintf(data, "     Halo    = %d\n",header->npart[1]);
	fprintf(data, "     Disk    = %d\n",header->npart[2]);
	fprintf(data, "     Bulge   = %d\n",header->npart[3]);
	fprintf(data, "     Stars   = %d\n",header->npart[4]);
	fprintf(data, "     Bndry   = %d\n",header->npart[5]);
	fprintf(data, "Mass: \n");
	fprintf(data, "     Gas     = %g\n",header->mass[0]); 
	fprintf(data, "     Halo    = %g\n",header->mass[1]); 
	fprintf(data, "     Disk    = %g\n",header->mass[2]); 
	fprintf(data, "     Bulge   = %g\n",header->mass[3]); 
	fprintf(data, "     Stars   = %g\n",header->mass[4]); 
	fprintf(data, "     Bndry   = %g\n",header->mass[5]); 	
	fprintf(data, "Time            = %g\n",header->time);
	fprintf(data, "Redshift        = %g\n",header->redshift);
	fprintf(data, "flag_sfr        = %d\n",header->flag_sfr);
	fprintf(data, "flag_feedback   = %d\n",header->flag_feedback);
	fprintf(data, "NpartTotal:\n");
	fprintf(data, "     Gas     = %d\n",header->npartTotal[0]);
	fprintf(data, "     Halo    = %d\n",header->npartTotal[1]);
	fprintf(data, "     Disk    = %d\n",header->npartTotal[2]);
	fprintf(data, "     Bulge   = %d\n",header->npartTotal[3]);
	fprintf(data, "     Stars   = %d\n",header->npartTotal[4]);
	fprintf(data, "     Bndry   = %d\n",header->npartTotal[5]);
	fprintf(data, "flag_cooling    = %d\n",header->flag_cooling);
	fprintf(data, "numfiles        = %d\n",header->num_files);
	fprintf(data, "boxsize         = %g\n",header->BoxSize);
	fprintf(data, "Omega0          = %g\n",header->Omega0);
	fprintf(data, "OmegaLambda     = %g\n",header->OmegaLambda);
	fprintf(data, "Hubbleparam     = %g\n",header->HubbleParam);

	SKIP; PRINT;
	// Counting particles:
	for(k=0; k<5; k++)
	{
		NumPart += header->npartTotal[k]; // all
		if(header->mass[k] == 0)
			Nmasses+= header->npartTotal[k]; // with masses
	}
	Ngas = header->npartTotal[0];        // gas
	
	SKIP; PRINT;

	// Reading the positions:
	for(i=0, k=0; k<6; k++)
		for(n=0; n<header->npart[k]; n++)
		{
			fread(&aux3[0],sizeof(float),3,fd);
			i++;
		}
	fprintf(data, "Read %d particle positions\n", i);

	SKIP; PRINT;
	SKIP; PRINT;

	// Reading the velocities:
	for(i=0, k=0; k<6; k++)
		for(n=0; n<header->npart[k]; n++)
		{
			fread(&aux3[0],sizeof(float),3,fd);
			i++;
		}
	fprintf(data, "Read %d particle velocities\n", i);
		
	SKIP; PRINT;
	SKIP; PRINT;

	// Reading the IDs:
	for(auxc=0, k=0; k<6; k++)
		for(n=0; n<header->npart[k]; n++)
		{
			fread(&auxi, sizeof(unsigned int), 1, fd);
			auxc++;
		}
	fprintf(data, "Read %u particle IDs\n", auxc);
		
	SKIP; PRINT;

	// Reading the masses:
	if(Nmasses != 0)
	{
		SKIP; PRINT;
		for(i=0, k=0; k<6; k++)
			for(n=0; n<header->npart[k]; n++)
			{
				fread(&aux, sizeof(float), 1, fd);
				i++;
			}
		fprintf(data, "Read %d particle masses\n ", i);
		SKIP; PRINT;
	}

	// Reading gas properties:
	if(Ngas != 0)
	{
		SKIP; PRINT;
	// Energy:
		for(n=0; n<header->npart[0]; n++)
			fread(&aux, sizeof(float), 1, fd);
		fprintf(data, "Read %d gas particle energy values", n);
		SKIP; PRINT;
		SKIP; PRINT;
	// Density:
		for(n=0; n<header->npart[0]; n++)
			fread(&aux, sizeof(float), 1, fd);
		fprintf(data, "Read %d gas particle density values", n);
		SKIP; PRINT;
		SKIP; PRINT;
	// HSML:
		for(n=0; n<header->npart[0]; n++)
			fread(&aux, sizeof(float), 1, fd);
		fprintf(data, "Read %d gas particle hslm values", n);
		SKIP; PRINT;
	}
	fclose(fd);
	fclose(data);

	return 0;
}

#endif


#include "snapTools.h"

int main(int argc, char **argv)
{
	// Several declarations:

	SNAPSHOT *a;
	PARTICLE list[10];

	// reading from file:
	char filename[255];
	int nfiles;

	sprintf(filename, "/path/to/file/snapshot");
	nfiles = 16; // n fragments of snapshot.

	a = loadSnapshot(filename, nfiles);

	// Unify fragments:
	unifySnapshot(a);

	// Order particles:
	orderIDs(a);

	// Copy first 10 particles:
	for (int i=0; i<10; i++)
		list[i] = getParticle(a, i);
	
	// Free snapshot:
	freeSnapshot(a);

	// Create new snapshot:
	a = newSnapshot();

	// Set values:
	setRedshift(a, 0);
	setCosmology(a, 0.3, 0.7, 0.7);

	// Add particles:
	addParticles(a, list, 10, 1, 0.003);

	// Save it in 2 fragments:
	fragmentSnapshot(a, 2);
	sprintf(filename, "/path/to/output/snapshot");
	saveSnapshot(a, filename);

	freeSnapshot(a);

	return 0;
}

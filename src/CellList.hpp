#ifndef CELL_LIST_H
#define CELL_LIST_H

#include <vector>
#include "Particle.hpp"

class CellList {
public:
    std::vector<int> BucketFirst;   // Index of the first particle in each bucket
    std::vector<int> BucketLast;    // Index of the last particle in each bucket
    std::vector<int> Nextof;        // Linked list to traverse particles in the same bucket
    std::vector<int> BucketIndex;   // The bucket index for each particle

    double DB ; //the one side of the bucket
    int nBx ; //the number of the buckets in x axis
    int nBy ; //the number of the buckets in y axis
    int nBz ; //the number of the buckets in z axis
    int nBxy ;
    int nBxyz ;

    double MIN_X ; //the minimum side of x in analysis area
    double MIN_Y ;
    double MIN_Z ;
    double MAX_X ; //the maximum side of x in analysis area
    double MAX_Y ;
    double MAX_Z ;

    double Lx;
    double Ly;
    double Lz;

    // Constructor
    CellList(int numParticles, const CellList& celllist);

    // Methods
    void initializeBuckets();
    void assignParticlesToBuckets(const std::vector<Particle>& particles);
    // void updateParticlePosition(int particleIndex, const std::vector<double>& newPosition, const CellList& celllist);
    // std::vector<int> getNeighboringParticles(int particleIndex, const CellList& celllist) const;
};

#endif // CELL_LIST_H

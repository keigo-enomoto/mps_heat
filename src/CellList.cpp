
/*

まだ、何もやっていない。後回し
 */

#include <iostream>
#include "CellList.hpp"

CellList::CellList(int numParticles, const CellList& celllist)
    : numParticles(numParticles) {

    // Initialize bucket dimensions and cell sizes
    numBucketsX = celllist.grid_size_x;
    numBucketsY = celllist.grid_size_y;
    numBucketsZ = celllist.grid_size_z;

    domainSizeX = celllist.Lx;
    domainSizeY = celllist.Ly;
    domainSizeZ = celllist.Lz;

    cellSizeX = domainSizeX / numBucketsX;
    cellSizeY = domainSizeY / numBucketsY;
    cellSizeZ = domainSizeZ / numBucketsZ;

    // Resize vectors for buckets and particle lists
    BucketFirst.resize(numBucketsX * numBucketsY * numBucketsZ, -1);
    BucketLast.resize(numBucketsX * numBucketsY * numBucketsZ, -1);
    Nextof.resize(numParticles, -1);
    BucketIndex.resize(numParticles, -1);
}

void CellList::initializeBuckets() {
    // Reset the bucket indices
    std::fill(BucketFirst.begin(), BucketFirst.end(), -1);
    std::fill(BucketLast.begin(), BucketLast.end(), -1);
    std::fill(Nextof.begin(), Nextof.end(), -1);
}

void CellList::assignParticlesToBuckets(const std::vector<Particle>& particles) {
    initializeBuckets();

    for (int i = 0; i < numParticles; ++i) {
        const auto& pos = particles[i].position;

        // Determine which bucket the particle belongs to
        int ix = static_cast<int>(pos[0] / cellSizeX);
        int iy = static_cast<int>(pos[1] / cellSizeY);
        int iz = static_cast<int>(pos[2] / cellSizeZ);

        // Compute the bucket index
        int bucketIndex = ix + numBucketsX * (iy + numBucketsY * iz);
        BucketIndex[i] = bucketIndex;

        // Insert particle into the bucket's linked list
        if (BucketFirst[bucketIndex] == -1) {
            BucketFirst[bucketIndex] = i;
        } else {
            Nextof[BucketLast[bucketIndex]] = i;
        }
        BucketLast[bucketIndex] = i;
    }
}

void CellList::updateParticlePosition(int particleIndex, const std::vector<double>& newPosition, const CellList& celllist) {
    int oldBucket = BucketIndex[particleIndex];

    // Determine the new bucket based on the updated position
    int ix = static_cast<int>(newPosition[0] / cellSizeX);
    int iy = static_cast<int>(newPosition[1] / cellSizeY);
    int iz = static_cast<int>(newPosition[2] / cellSizeZ);
    int newBucket = ix + numBucketsX * (iy + numBucketsY * iz);

    // If the particle has moved to a new bucket, update the linked list
    if (newBucket != oldBucket) {
        int prevParticle = -1;
        int currentParticle = BucketFirst[oldBucket];

        // Traverse the linked list to find the particle
        while (currentParticle != -1) {
            if (currentParticle == particleIndex) {
                if (prevParticle == -1) {
                    BucketFirst[oldBucket] = Nextof[currentParticle];
                } else {
                    Nextof[prevParticle] = Nextof[currentParticle];
                }
                break;
            }
            prevParticle = currentParticle;
            currentParticle = Nextof[currentParticle];
        }

        // Add the particle to the new bucket's linked list
        if (BucketFirst[newBucket] == -1) {
            BucketFirst[newBucket] = particleIndex;
        } else {
            Nextof[BucketLast[newBucket]] = particleIndex;
        }
        BucketLast[newBucket] = particleIndex;

        BucketIndex[particleIndex] = newBucket;
    }
}

std::vector<int> CellList::getNeighboringParticles(int particleIndex, const CellList& celllist) const {
    std::vector<int> neighbors;

    // Get the bucket of the particle
    int bucketIndex = BucketIndex[particleIndex];
    int ix = bucketIndex % numBucketsX;
    int iy = (bucketIndex / numBucketsX) % numBucketsY;
    int iz = bucketIndex / (numBucketsX * numBucketsY);

    // Loop over neighboring buckets (including itself)
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                int neighborIx = (ix + dx + numBucketsX) % numBucketsX;
                int neighborIy = (iy + dy + numBucketsY) % numBucketsY;
                int neighborIz = (iz + dz + numBucketsZ) % numBucketsZ;

                int neighborBucketIndex = neighborIx + numBucketsX * (neighborIy + numBucketsY * neighborIz);
                int currentParticle = BucketFirst[neighborBucketIndex];

                // Traverse the particles in the neighboring bucket
                while (currentParticle != -1) {
                    neighbors.push_back(currentParticle);
                    currentParticle = Nextof[currentParticle];
                }
            }
        }
    }

    return neighbors;
}

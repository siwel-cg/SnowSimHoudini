#pragma once

#include <vector>
#include <Eigen/Dense>
#include "particle.h"
#include "grid.h"
#pragma once

// #ifndef MPMSOLVER_H
// #define MPMSOLVER_H

class MPMSolver
{
public:
    MPMSolver();
    MPMSolver(Eigen::Vector3f gridDim, float spacing, Eigen::Vector3f gridOrigin, float groundPlane, float dt,
        float critCompression, float critStretch,
        float hardeningCoeff, float initialDensity, float youngsMod,
        float poissonRatio);
    void addParticle(const MPMParticle& particle);
    void computeForcesAndIntegrate();

	const std::vector<MPMParticle>& getParticles() const;
    
    float stepSize;
    float critCompression;
    float critStretch;
    float hardeningCoeff;
    float initialDensity;
    float youngsMod;
    float poissonRatio;

    float groundPlaneY;

    float mu0;
    float lambda0;

    void computeInitialDensity();

    // THIS IS THE MAIN FUNCTION THAT SHOULD UPDATE THE SIM EACH FRAME
    void step();

private:
    std::vector<MPMParticle> particles;
    mpmgrid grid;


    Eigen::Vector3f computeGravity(const MPMParticle& p);
    Eigen::Vector3f computeCohesion(const MPMParticle& p);
    void integrate(MPMParticle& p, Eigen::Vector3f force);

    // PARTICLE FUNCTIONS
    void computeSigma();
    void updateParticleDefGrad();

    // GRID FUNCTIONS
    void particleToGridTransfer();
    void computeForce();
    void updateGridVel();

};

//#endif // MPMSOLVER_H

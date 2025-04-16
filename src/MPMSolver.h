#pragma once

// #ifndef MPMSOLVER_H
// #define MPMSOLVER_H

#include "particle.h"
#include "grid.h"
#include <Eigen/Dense>


class MPMSolver
{
public:
    MPMSolver(Eigen::Vector3f gridDim, float spacing, Eigen::Vector3f gridOrigin, float dt,
        float critCompression, float critStretch,
        float hardeningCoeff, float initialDensity, float youngsMod,
        float poissonRatio);
    MPMSolver();
    void addParticle(const particle& particle);
    void computeForcesAndIntegrate();

    const vector<particle>& getParticles() const;

    float stepSize;
    float critCompression;
    float critStretch;
    float hardeningCoeff;
    float initialDensity;
    float youngsMod;
    float poissonRatio;

    float mu0;
    float lambda0;

    void computeInitialDensity();

    // THIS IS THE MAIN FUNCTION THAT SHOULD UPDATE THE SIM EACH FRAME
    void step();

private:
    vector<particle> particles;
    mpmgrid grid;


    Eigen::Vector3f computeGravity(const particle& p);
    Eigen::Vector3f computeCohesion(const particle& p);
    void integrate(particle& p, Eigen::Vector3f force);

    // PARTICLE FUNCTIONS
    void computeSigma();
    void updateParticleDefGrad();

    // GRID FUNCTIONS
    void particleToGridTransfer();
    void computeForce();
    void updateGridVel();

};

//#endif // MPMSOLVER_H

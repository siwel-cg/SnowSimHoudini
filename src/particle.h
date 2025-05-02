#pragma once
#include <vector>
#include <Eigen/Dense>

class MPMParticle
{
public:
    Eigen::Vector3f position;
    Eigen::Vector3f velocity;
    float mass;
    float density; 
    float volume;

    Eigen::Matrix3f FE;
    Eigen::Matrix3f FP;
    Eigen::Matrix3f sigma;

    MPMParticle(Eigen::Vector3f pos, Eigen::Vector3f vel, float m);

    void computeFE();
    void computeFP();
};




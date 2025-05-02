#include "particle.h"

MPMParticle::MPMParticle(Eigen::Vector3f pos, Eigen::Vector3f vel, float m)
    : position(pos), velocity(vel), mass(m), FE(Eigen::Matrix3f::Identity()), FP(Eigen::Matrix3f::Identity()), sigma(Eigen::Matrix3f::Identity()) {
}


#ifndef MPMGRID_H
#define MPMGRID_H

#include <Eigen/Dense>

struct GridNode {
    Eigen::Vector3f velocity;
    Eigen::Vector3f prevVelocity;
    float mass;
    float density;
    Eigen::Vector3f force;
    float velocityMass;

    Eigen::Vector3f worldPos;
    //Eigen::Vector3i idx;
    GridNode();
};

class mpmgrid
{
private:

public:
    mpmgrid();
    mpmgrid(Eigen::Vector3f dim, float spacing, Eigen::Vector3f center);
    Eigen::Vector3f dimension;
    float spacing;
    Eigen::Vector3f center;
    std::vector<GridNode> gridNodes;

    // NUMBER OF CELLS IN EACH DIM
    int nx;
    int ny;
    int nz;

    GridNode* getGridNode(float x, float y, float z);
    void clearGrid();
    void divideMass();
};

#endif // MPMGRID_H

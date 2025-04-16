#include "grid.h"
#include <cmath>

GridNode::GridNode()
    : velocity(0.f, 0.f, 0.f), prevVelocity(0.f, 0.f, 0.f), mass(0.f), density(0.f), force(0.f, 0.f, 0.f),
    worldPos(0.f, 0.f, 0.f), velocityMass(0.f)
{
}

mpmgrid::mpmgrid()
    : dimension(Eigen::Vector3f(3.0f, 3.0f, 3.0f)),
    spacing(1.f),
    center(Eigen::Vector3f(0.f, 0.f, 0.f))
{
}

mpmgrid::mpmgrid(Eigen::Vector3f dim, float spc, Eigen::Vector3f cent)
    : dimension(dim),
    spacing(spc),
    center(cent)
{
    nx = static_cast<int>(std::floor(dimension[0] / spacing));
    ny = static_cast<int>(std::floor(dimension[1] / spacing));
    nz = static_cast<int>(std::floor(dimension[2] / spacing));

    if (nx < 1) nx = 1;
    if (ny < 1) ny = 1;
    if (nz < 1) nz = 1;

    gridNodes.resize(nx * ny * nz);

    // INITIALIZE EACH GRID NODES POSITION AND INDEX
    Eigen::Vector3f minCorner = center - 0.5f * dimension;
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int idx = i + nx * (j + ny * k);

                //gridNodes[idx].idx = Eigen::Vector3i(i, j, k);

                // CENTER POS IN WORLD SPACE
                float xNode = minCorner[0] + (i + 0.5f) * spacing;
                float yNode = minCorner[1] + (j + 0.5f) * spacing;
                float zNode = minCorner[2] + (k + 0.5f) * spacing;
                gridNodes[idx].worldPos = Eigen::Vector3f(xNode, yNode, zNode);
            }
        }
    }
}


GridNode* mpmgrid::getGridNode(float x, float y, float z)
{
    float xMin = center[0] - 0.5f * dimension[0];
    float yMin = center[1] - 0.5f * dimension[1];
    float zMin = center[2] - 0.5f * dimension[2];

    int i = static_cast<int>(std::floor((x - xMin) / spacing));
    int j = static_cast<int>(std::floor((y - yMin) / spacing));
    int k = static_cast<int>(std::floor((z - zMin) / spacing));

    if (i < 0)    i = 0;
    if (i >= nx)  i = nx - 1;
    if (j < 0)    j = 0;
    if (j >= ny)  j = ny - 1;
    if (k < 0)    k = 0;
    if (k >= nz)  k = nz - 1;

    int idx = i + nx * (j + ny * k);

    return &gridNodes[idx];
}

void mpmgrid::clearGrid() {
    // CLEAR GRID FOR NEW VALUES
    for (GridNode& node : gridNodes) {
        node.mass = 0.0f;
        node.velocity = Eigen::Vector3f(0.f, 0.f, 0.f);
    }
}

// THIS FUNCTION IS USED WHEN TRANSFERING PARTICLES TO THE GRID
// THIS DIVIDES THE VELOCITY OF EACH NODE BY ITS MASS
//
// ONLY CALL THIS WHEN YOU ARE SURE THE STORED VELOCITY IS MOMENTUM
void mpmgrid::divideMass() {
    for (GridNode& n : gridNodes) {
        if (n.mass != 0.0) {
            n.velocity /= n.mass;
        }
    }
}




#include "grid.h"
#include <cmath>

GridNode::GridNode()
    : velocity(Eigen::Vector3f::Zero()), prevVelocity(Eigen::Vector3f::Zero()),
    mass(0.f), density(0.f),
    force(Eigen::Vector3f::Zero()), velocityMass(0.f),
    worldPos(Eigen::Vector3f::Zero()),
    idx(Eigen::Vector3i::Zero())
{
}

mpmgrid::mpmgrid()
    : dimension(3.0f, 3.0f, 3.0f),
    spacing(1.f),
    center(0, 0, 0)
{
}

mpmgrid::mpmgrid(const Eigen::Vector3f& dim,
    float spc,
    const Eigen::Vector3f& cent)
    : dimension(dim),
    spacing(spc),
    center(cent)
{
    nx = std::max(1, int(std::floor(dimension.x() / spacing)));
    ny = std::max(1, int(std::floor(dimension.y() / spacing)));
    nz = std::max(1, int(std::floor(dimension.z() / spacing)));

    gridNodes.resize(nx * ny * nz);

    Eigen::Vector3f minCorner = center - 0.5f * dimension;
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int idx = i + nx * (j + ny * k);
                gridNodes[idx].idx = Eigen::Vector3i(i, j, k);
                float xNode = minCorner.x() + (i + 0.5f) * spacing;
                float yNode = minCorner.y() + (j + 0.5f) * spacing;
                float zNode = minCorner.z() + (k + 0.5f) * spacing;
                gridNodes[idx].worldPos = { xNode,yNode,zNode };
            }
        }
    }
}

GridNode* mpmgrid::getGridNode(float x, float y, float z) {
    float xMin = center.x() - 0.5f * dimension.x();
    float yMin = center.y() - 0.5f * dimension.y();
    float zMin = center.z() - 0.5f * dimension.z();

    int i = std::clamp(int(std::floor((x - xMin) / spacing)), 0, nx - 1);
    int j = std::clamp(int(std::floor((y - yMin) / spacing)), 0, ny - 1);
    int k = std::clamp(int(std::floor((z - zMin) / spacing)), 0, nz - 1);

    int idx = i + nx * (j + ny * k);
    return &gridNodes[idx];
}

void mpmgrid::clearGrid() {
    for (auto& n : gridNodes) {
        n.mass = 0.0f;
        n.density = 0.0f;
        n.velocity = Eigen::Vector3f::Zero();
        n.prevVelocity = Eigen::Vector3f::Zero();
        n.force = Eigen::Vector3f::Zero();
        n.velocityMass = 0.0f;
    }
}

void mpmgrid::divideMass() {
    for (auto& n : gridNodes) {
        if (n.mass != 0.0f)
            n.velocity /= n.mass;
    }
}

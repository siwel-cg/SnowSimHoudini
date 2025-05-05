#include "mpmsolver.h"
#include <iostream>

#include <Eigen/Dense>

const float dt = 0.01f; // Time step
const float gravity = 30.81f; // Gravity acceleration
const float cohesionStrength = 0.1f; // Cohesion force factor

MPMSolver::MPMSolver() :
     stepSize(0.00001f)
    , grid(Eigen::Vector3f(1.0f, 1.0f, 1.0f), 0.05f, Eigen::Vector3f(0.0f, 0.0f, 0.0f))
    , critCompression(0.025f)
    , critStretch(0.0075f)
    , hardeningCoeff(10.0f)
    , initialDensity(400.0f)
    , youngsMod(140000.0f)
    , poissonRatio(0.2f)
{
    groundPlaneY = 0.0;
    mu0 = youngsMod / (2.f * (1.f + poissonRatio));
    lambda0 = (youngsMod * poissonRatio) / ((1.f + poissonRatio) * (1.f - 2.f * poissonRatio));
    const GU_PrimVDB* vdbPrimSDF = nullptr;
} 


MPMSolver::MPMSolver(Eigen::Vector3f gridDim, float spacing, Eigen::Vector3f gridOrigin, float groundPlane, float dt,
    float critCompression, float critStretch, float hardeningCoeff,
    float initialDensity, float youngsMod, float poissonRatio, const GU_PrimVDB* vdbPrimSDF)
    : stepSize(dt), grid(Eigen::Vector3f(gridDim), spacing, Eigen::Vector3f(gridOrigin)), groundPlaneY(groundPlane),
    critCompression(critCompression), critStretch(critStretch), hardeningCoeff(hardeningCoeff),
    initialDensity(initialDensity), youngsMod(youngsMod), poissonRatio(poissonRatio), vdbPrimSDF(vdbPrimSDF)
{
    mu0 = youngsMod / (2.f * (1.f + poissonRatio));
    lambda0 = (youngsMod * poissonRatio) / ((1.f + poissonRatio) * (1.f - 2.f * poissonRatio));
}

void MPMSolver::addParticle(const MPMParticle& particle) {
    particles.push_back(particle);
}

// THIS IS HELPER FUNCTION FOR COMPUTING WEIGHTING
// TODO : ENSURE ALL WEIGHTS IN KERNEL ADD UP TO 1
static float weightFun(float x) {
    x = abs(x);
    if (x >= 0.f && x < 1.f) {
        return 0.5f * (x * x * x) - (x * x) + 2.f / 3.f;
    }
    else if (x >= 1.f && x < 2.f) {
        return -(1.f / 6.f) * (x * x * x) + (x * x) - 2.f * x + (4.f / 3.f);
    }
    else {
        return 0.0f;
    }
}

// static float weightFunGradient (float x) {
//     x = abs(x);
//     if (x > 0.f && x < 1.f) {
//         return 3.f*abs(x)/(2.f*x) - 2.f*x;
//     } else if (x >= 1.f && x < 2.f) {
//         return -abs(x)/(2.f*x)*x*x + 2.f*x - 2.f*(abs(x)/x);
//     } else {
//         return 0.0f;
//     }
// }

static float weightFunGradient(float x)
{
    const float ax = std::abs(x);
    const float s = (x < 0.f) ? -1.f : 1.f;
    if (ax >= 2.f) {
        return 0.f;
    }
    if (ax < 1.f) {      // |x| < 1
        return (1.5f * ax - 2.f) * s;
    }
    /* 1 <= |x| < 2 */
    return (-0.5f * ax + 2.f - (2.f / 3.f)) * s;
}

// [======] PARTICLE FUNCTIONS [======]

void MPMSolver::computeSigma() {
    for (MPMParticle& p : particles) {
        // CONVERT GLM MATRIX TO EIGEN
        Eigen::Matrix3f Fe = p.FE;
        Eigen::Matrix3f Fp = p.FP;

        // COMPUTE POLAR DECOMP TO GET ROTATIONAL PART OF F
        Eigen::JacobiSVD<Eigen::Matrix3f> svd(Fe, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Matrix3f U = svd.matrixU();
        Eigen::Matrix3f V = svd.matrixV();
        Eigen::Matrix3f R = U * V.transpose();

        float Je = Fe.determinant();  // Elastic volume
        float Jp = Fp.determinant();  // Plastic volume
        float Jtot = (Fe * Fp).determinant(); // TOTAL DEFORMATION

        // Check if deformation determinant is negative (Flips particles inside out)
        // Reset to keep sim stable
        auto bad = [](float x) {
            return (x <= 0.0f) || !std::isfinite(x);
            };

        if (bad(Je) || bad(Jp) || bad(Jtot))
        {
            //std::cout << "WAHOOO" << std::endl;
            p.FE = Eigen::Matrix3f::Identity();
            p.FP = Eigen::Matrix3f::Identity();
            p.sigma = Eigen::Matrix3f::Zero();
            continue;
        }

        float xi = hardeningCoeff;//10.0f; // default value, tweak as needed

        // Plastic-hardening-modified Lame parameters
        float mu = mu0 * std::exp(xi * (1.0f - Jp));
        float lambda = lambda0 * std::exp(xi * (1.0f - Jp));

        Eigen::Matrix3f strain = Fe - R;
        // Eigen::Matrix3f sigma = 2.f * mu * strain + lambda * (Je - 1.0f) * Eigen::Matrix3f::Identity();
        // sigma *= (1.f / Je); // Cauchy stress

        Eigen::Matrix3f sigma = 2.f * mu * (Fe - R) + lambda * std::log(Je) * Eigen::Matrix3f::Identity();
        sigma *= 1.f / Je;

        p.sigma = sigma;
    }
}

void MPMSolver::updateParticleDefGrad() {
    //Fn+1 = (I + ∆t∇v^n+1_p)Fn

    float xMin = grid.center[0] - 0.5f * grid.dimension[0];
    float yMin = grid.center[1] - 0.5f * grid.dimension[1];
    float zMin = grid.center[2] - 0.5f * grid.dimension[2];

    for (MPMParticle& p : particles) {
        Eigen::Matrix3f velGrad = Eigen::Matrix3f::Zero();

        // WORLD SPACE POSITIONS
        float x = p.position[0];
        float y = p.position[1];
        float z = p.position[2];

        int i = static_cast<int>(std::floor((x - xMin) / grid.spacing));
        int j = static_cast<int>(std::floor((y - yMin) / grid.spacing));
        int k = static_cast<int>(std::floor((z - zMin) / grid.spacing));

        Eigen::Vector3f vPic = Eigen::Vector3f::Zero();
        Eigen::Vector3f vFlip = Eigen::Vector3f::Zero();
        // YOU NEED TO DO THIS FOR ALL CELLS WITHIN SOME RADIUS
        for (int di = -2; di <= 2; di++) {
            for (int dj = -2; dj <= 2; dj++) {
                for (int dk = -2; dk <= 2; dk++) {
                    // INDEX OF CURRENT NODE WE ARE LOOKING AT
                    int iNode = i + di;
                    int jNode = j + dj;
                    int kNode = k + dk;

                    // CLAMP TO BE INSIDE THE GRID DOMAIN
                    if (iNode < 0 || iNode >= grid.nx) continue;
                    if (jNode < 0 || jNode >= grid.ny) continue;
                    if (kNode < 0 || kNode >= grid.nz) continue;

                    int idx = iNode + grid.nx * (jNode + grid.ny * kNode);
                    GridNode& curNode = grid.gridNodes[idx];

                    // POSITION IN GRID SPACE
                    float xGrid = (x - curNode.worldPos[0]) / grid.spacing;
                    float yGrid = (y - curNode.worldPos[1]) / grid.spacing;
                    float zGrid = (z - curNode.worldPos[2]) / grid.spacing;

                    // WEIGHT OF GRID CELL RELATIVE TO PARTICLE GRID
                    float weight = weightFun(xGrid) * weightFun(yGrid) * weightFun(zGrid);
                    if (weight == 0.0) continue;

                    // TRANSFER WEIGHTED VELOCITY FROM GRID TO PARTICLE
                    //p.velocity += curNode.velocity * weight;
                    //if (glm::any(glm::isnan(curNode.velocity))) continue; // THIS IS WHERE NANS ARE COMING FROM


                    Eigen::Vector3f gradWeight;
                    gradWeight[0] = 1.f / grid.spacing * weightFunGradient(xGrid) * weightFun(yGrid) * weightFun(zGrid);
                    gradWeight[1] = 1.f / grid.spacing * weightFun(xGrid) * weightFunGradient(yGrid) * weightFun(zGrid);
                    gradWeight[2] = 1.f / grid.spacing * weightFun(xGrid) * weightFun(yGrid) * weightFunGradient(zGrid);

                    // THIS IS v * gradW^T
                    velGrad += curNode.velocity * gradWeight.transpose();

                    vPic += curNode.velocity * weight;
                    vFlip += (curNode.velocity - curNode.prevVelocity) * weight;

                    // THIS IS THE OLD INCCORECT WAY OF UPDATING VELOCITY
                    // IT GIVES MOTION, BUT IS WRONG
                    //p.velocity += curNode.velocity * weight;
                }
            }
        }

        vFlip += p.velocity;
        float alpha = 0.80f;

        p.velocity = (1.f - alpha) * vPic + alpha * vFlip;
        const float damping = 0.01f;     // 1% velocity loss each step
        p.velocity *= (1.0f - damping);

        // UPDATE DEFORMATION GRADIENT
        p.FE = (Eigen::Matrix3f::Identity() + stepSize * velGrad) * p.FE;


        // Predict the new total deformation gradient (FE * FP)
        //glm::mat3 F_total = (glm::mat3(1.0f) + stepSize * velGrad) * p.FE * p.FP;

        // Predict the new elastic deformation gradient
        Eigen::Matrix3f FE_hat = p.FE;


        Eigen::Matrix3f FE_hat_eigen = FE_hat;


        // SVD: FE_hat = U * Σ * V^T
        Eigen::JacobiSVD<Eigen::Matrix3f> svd(FE_hat_eigen, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Matrix3f U = svd.matrixU();
        Eigen::Matrix3f V = svd.matrixV();
        Eigen::Vector3f sigma_hat = svd.singularValues(); // Σ̂

        // Clamp singular values to [1 - θc, 1 + θs]
        Eigen::Vector3f sigma_clamped = sigma_hat;
        for (int i = 0; i < 3; ++i)
            sigma_clamped[i] = std::clamp(sigma_hat[i], 1.0f - critCompression, 1.0f + critStretch);

        // Reconstruct clamped FE
        Eigen::Matrix3f Sigma_clamped = Eigen::Matrix3f::Zero();
        for (int i = 0; i < 3; ++i) {
            Sigma_clamped(i, i) = sigma_clamped[i];
        }
    
        Eigen::Matrix3f FE_new = U * Sigma_clamped * V.transpose();
        Eigen::Matrix3f Fp_old = p.FP;
        /*for (int c = 0; c < 3; ++c) {
            for (int r = 0; r < 3; ++r) {
                Fp_old(r, c) = p.FP[r][c];
            }
        }*/

        // Eigen::Matrix3f F_total_eigen;
        // for (int col = 0; col < 3; ++col)
        //     for (int row = 0; row < 3; ++row)
        //         F_total_eigen(row, col) = F_total[col][row];

        // Update FP using: FP = V * Σ⁻¹ * Uᵀ * F_total
        Eigen::Matrix3f Sigma_inv = Eigen::Matrix3f::Zero();
        for (int i = 0; i < 3; ++i)
            Sigma_inv(i, i) = 1.0f / sigma_clamped[i];

        Eigen::Matrix3f F_tot_new = FE_new * Fp_old;
        Eigen::Matrix3f FP_new = V * Sigma_inv * U.transpose() * F_tot_new;

        p.FE = FE_new;
        p.FP = FP_new;


        // UPDATE POINT POSITIONS
        p.position += stepSize * p.velocity;
    }


    // glm::vec3 minCorner = grid.center - 0.5f * grid.dimension;
    // glm::vec3 maxCorner = grid.center + 0.5f * grid.dimension;


    // //THIS IS DOING NOTHING : :
    // float damping = 0.0f; // or try 0.01f, 0.1f for bounciness
    // for (MPMParticle &p : particles) {
    //     for (int axis = 0; axis < 3; ++axis) {
    //         if (p.position[axis] < minCorner[axis]) {
    //             p.position[axis] = minCorner[axis];
    //             if (p.velocity[axis] < 0.f) {
    //                 p.velocity[axis] *= -damping;
    //             }
    //         }

    //         if (p.position[axis] > maxCorner[axis]) {
    //             p.position[axis] = maxCorner[axis];
    //             if (p.velocity[axis] > 0.f) {
    //                 p.velocity[axis] *= -damping;
    //             }
    //         }
    //     }
    // }
}

// [======] GRID FUNCTIONS [======]
void MPMSolver::particleToGridTransfer() {
    float xMin = grid.center[0] - 0.5f * grid.dimension[0];
    float yMin = grid.center[1] - 0.5f * grid.dimension[1];
    float zMin = grid.center[2] - 0.5f * grid.dimension[2];

    //grid.clearGrid();
    for (MPMParticle& p : particles) {

        // World space positions
        float x = p.position[0];
        float y = p.position[1];
        float z = p.position[2];

        int i = static_cast<int>(std::floor((x - xMin) / grid.spacing));
        int j = static_cast<int>(std::floor((y - yMin) / grid.spacing));
        int k = static_cast<int>(std::floor((z - zMin) / grid.spacing));

        // YOU NEED TO DO THIS FOR ALL CELLS WITHIN SOME RADIUS
        for (int di = -2; di <= 2; di++) {
            for (int dj = -2; dj <= 2; dj++) {
                for (int dk = -2; dk <= 2; dk++) {
                    // INDEX OF CURRENT NODE WE ARE LOOKING AT
                    int iNode = i + di;
                    int jNode = j + dj;
                    int kNode = k + dk;

                    // CLAMP TO BE INSIDE THE GRID DOMAIN
                    if (iNode < 0 || iNode >= grid.nx) continue;
                    if (jNode < 0 || jNode >= grid.ny) continue;
                    if (kNode < 0 || kNode >= grid.nz) continue;


                    int idx = iNode + grid.nx * (jNode + grid.ny * kNode);
                    GridNode& curNode = grid.gridNodes[idx];

                    // POSITION IN GRID SPACE
                    float xGrid = (x - curNode.worldPos[0]) / grid.spacing;
                    float yGrid = (y - curNode.worldPos[1]) / grid.spacing;
                    float zGrid = (z - curNode.worldPos[2]) / grid.spacing;

                    // WEIGHT OF GRID CELL RELATICE TO PARTICLE GRID
                    float weight = weightFun(xGrid) * weightFun(yGrid) * weightFun(zGrid);
                    if (weight == 0.0) continue;


                    // UPDATE CURRENT NODE WE ARE LOOKING AT
                    curNode.mass += p.mass * weight;
                    curNode.velocity += p.velocity * p.mass * weight;

                    // if (glm::length(p.velocity) > 1e-6f) { // Only if velocity is non-zero
                    //     curNode.velocityMass += curNode.mass;
                    // }
                }
            }
        }

    }

    // // Normalize grid node velocities using only effective mass
    // for (GridNode& node : grid.gridNodes) {
    //     if (node.velocityMass > 0.f) {
    //         node.velocity /= node.velocityMass;
    //     } else {
    //         node.velocity = glm::vec3(0.f);
    //     }
    // }

    // Currently the velocity stored is actually the total weighted momentum
    // To convert it to actual vel we divide each gridCell by its mass
    grid.divideMass();
}

// THIS SHOULD ONLY BE CALLED ONCE AT t=0
void MPMSolver::computeInitialDensity() {
    float xMin = grid.center[0] - 0.5f * grid.dimension[0];
    float yMin = grid.center[1] - 0.5f * grid.dimension[1];
    float zMin = grid.center[2] - 0.5f * grid.dimension[2];

    //grid.clearGrid();
    for (MPMParticle& p : particles) {
        // World space positions
        float x = p.position[0];
        float y = p.position[1];
        float z = p.position[2];

        int i = static_cast<int>(std::floor((x - xMin) / grid.spacing));
        int j = static_cast<int>(std::floor((y - yMin) / grid.spacing));
        int k = static_cast<int>(std::floor((z - zMin) / grid.spacing));

        // YOU NEED TO DO THIS FOR ALL CELLS WITHIN SOME RADIUS
        for (int di = -2; di <= 2; di++) {
            for (int dj = -2; dj <= 2; dj++) {
                for (int dk = -2; dk <= 2; dk++) {
                    // INDEX OF CURRENT NODE WE ARE LOOKING AT
                    int iNode = i + di;
                    int jNode = j + dj;
                    int kNode = k + dk;

                    // CLAMP TO BE INSIDE THE GRID DOMAIN
                    if (iNode < 0 || iNode >= grid.nx) continue;
                    if (jNode < 0 || jNode >= grid.ny) continue;
                    if (kNode < 0 || kNode >= grid.nz) continue;

                    int idx = iNode + grid.nx * (jNode + grid.ny * kNode);
                    GridNode& curNode = grid.gridNodes[idx];

                    // POSITION IN GRID SPACE
                    float xGrid = (x - curNode.worldPos[0]) / grid.spacing;
                    float yGrid = (y - curNode.worldPos[1]) / grid.spacing;
                    float zGrid = (z - curNode.worldPos[2]) / grid.spacing;

                    // WEIGHT OF GRID CELL RELATICE TO PARTICLE GRID
                    float weight = weightFun(xGrid) * weightFun(yGrid) * weightFun(zGrid);
                    if (weight == 0.0) continue;

                    // UPDATE CURRENT NODE WE ARE LOOKING AT
                    curNode.mass += p.mass * weight;
                }
            }
        }
    }

    // COMPUTE GRID DENSITY
    float volume = grid.spacing * grid.spacing * grid.spacing;
    for (GridNode& g : grid.gridNodes) {
        if (volume != 0.0) {
            g.density = (g.mass / volume);
        }
    }

    // TRANSFER GRID DENSITY TO PARTICLES
    for (MPMParticle& p : particles) {
        // World space positions
        float x = p.position[0];
        float y = p.position[1];
        float z = p.position[2];

        int i = static_cast<int>(std::floor((x - xMin) / grid.spacing));
        int j = static_cast<int>(std::floor((y - yMin) / grid.spacing));
        int k = static_cast<int>(std::floor((z - zMin) / grid.spacing));


        // YOU NEED TO DO THIS FOR ALL CELLS WITHIN SOME RADIUS
        for (int di = -2; di <= 2; di++) {
            for (int dj = -2; dj <= 2; dj++) {
                for (int dk = -2; dk <= 2; dk++) {
                    // INDEX OF CURRENT NODE WE ARE LOOKING AT
                    int iNode = i + di;
                    int jNode = j + dj;
                    int kNode = k + dk;

                    // CLAM TO BE INSIDE THE GRID DOMAIN
                    if (iNode < 0 || iNode >= grid.nx) continue;
                    if (jNode < 0 || jNode >= grid.ny) continue;
                    if (kNode < 0 || kNode >= grid.nz) continue;

                    int idx = iNode + grid.nx * (jNode + grid.ny * kNode);
                    GridNode& curNode = grid.gridNodes[idx];

                    // POSITION IN GRID SPACE
                    float xGrid = (x - curNode.worldPos[0]) / grid.spacing;
                    float yGrid = (y - curNode.worldPos[1]) / grid.spacing;
                    float zGrid = (z - curNode.worldPos[2]) / grid.spacing;

                    // WEIGHT OF GRID CELL RELATICE TO PARTICLE GRID
                    float weight = weightFun(xGrid) * weightFun(yGrid) * weightFun(zGrid);
                    if (weight == 0.0) continue;

                    // SUM UP DENSITY CONTRIBUTIONS
                    p.density += curNode.density * weight;
                }
            }
        }

        // DIVIDE TO GET INITIAL VOLUME OF EACH PARTICLE
        p.volume = p.mass / p.density;
    }
}

void MPMSolver::computeForce() {
    // ZERO OUT FORCE FOR EACH GRID CELL
    for (GridNode& g : grid.gridNodes) {
        g.force = Eigen::Vector3f(0.f, 0.f, 0.f);
    }

    UT_Vector3 gravity = GRAVITY(context.getTime());

    float xMin = grid.center[0] - 0.5f * grid.dimension[0];
    float yMin = grid.center[1] - 0.5f * grid.dimension[1];
    float zMin = grid.center[2] - 0.5f * grid.dimension[2];

    // UPDATE CAUCHY STRESS FOR EACH PARTICLE
    computeSigma();

    for (MPMParticle& p : particles) {
        // World space positions
        float x = p.position[0];
        float y = p.position[1];
        float z = p.position[2];

        int i = static_cast<int>(std::floor((x - xMin) / grid.spacing));
        int j = static_cast<int>(std::floor((y - yMin) / grid.spacing));
        int k = static_cast<int>(std::floor((z - zMin) / grid.spacing));


        // YOU NEED TO DO THIS FOR ALL CELLS WITHIN SOME RADIUS
        for (int di = -1; di < 2; di++) {
            for (int dj = -1; dj < 2; dj++) {
                for (int dk = -1; dk < 2; dk++) {
                    // INDEX OF CURRENT NODE WE ARE LOOKING AT
                    int iNode = i + di;
                    int jNode = j + dj;
                    int kNode = k + dk;

                    // CLAM TO BE INSIDE THE GRID DOMAIN
                    if (iNode < 0 || iNode >= grid.nx) continue;
                    if (jNode < 0 || jNode >= grid.ny) continue;
                    if (kNode < 0 || kNode >= grid.nz) continue;

                    int idx = iNode + grid.nx * (jNode + grid.ny * kNode);
                    GridNode& curNode = grid.gridNodes[idx];

                    // POSITION IN GRID SPACE
                    float xGrid = (x - curNode.worldPos[0]) / grid.spacing;
                    float yGrid = (y - curNode.worldPos[1]) / grid.spacing;
                    float zGrid = (z - curNode.worldPos[2]) / grid.spacing;

                    Eigen::Vector3f gradWeight;
                    gradWeight[0] = 1.f / grid.spacing * weightFunGradient(xGrid) * weightFun(yGrid) * weightFun(zGrid);
                    gradWeight[1] = 1.f / grid.spacing * weightFun(xGrid) * weightFunGradient(yGrid) * weightFun(zGrid);
                    gradWeight[2] = 1.f / grid.spacing * weightFun(xGrid) * weightFun(yGrid) * weightFunGradient(zGrid);

                    float weight = weightFun(xGrid) * weightFun(yGrid) * weightFun(zGrid);
                    //
                    curNode.force -= (p.volume * p.sigma * gradWeight);
                    curNode.force += weight * Eigen::Vector3f(gravity.x(), gravity.y(), gravity.z()) * p.mass;
                }
            }
        }
    }
}

float sphereSDF(const Eigen::Vector3f& pos, const Eigen::Vector3f& center, float radius) {
    //Eigen::Vector3f posVec(pos.x, pos.y, pos.z);
    return (pos - center).norm() - radius;
    /*Eigen::Vector3f extent(0.5f, 0.5f, 0.5f);
    Eigen::Vector3f q = (pos - center).cwiseAbs() - extent;
    return std::max(q.maxCoeff(), 0.0f) + std::min(std::max(q.x(), std::max(q.y(), q.z())), 0.0f);*/
}


float MPMSolver::querySdf(Eigen::Vector3f point)
{
    if (!vdbPrimSDF) {
        return std::numeric_limits<float>::max();
    }

    const UT_Vector3& pointVDB = UT_Vector3(point.x(), point.y(), point.z());
    return vdbPrimSDF->getValueF(pointVDB);
}

Eigen::Vector3f MPMSolver::sdfNormal(Eigen::Vector3f point)
{
    if (!vdbPrimSDF) {
        return Eigen::Vector3f(0, 1, 0);
    }

    const UT_Vector3& pointVDB = UT_Vector3(point.x(), point.y(), point.z());

    UT_Vector3 normal = vdbPrimSDF->getGradient(pointVDB);
    return Eigen::Vector3f(normal.x(), normal.y(), normal.z()).normalized();
}

// QUICK NOTE ABOUT THIS: The full implementation described in the paper uses a semi-implicit method for solving
//                        which results in a more stable and complex system. Since this is a push to get something out
//                        I am just doing an explicit solve with the most basic collision handling possible.
//                        For future implementation, fix this so that it fully implements the paper
void MPMSolver::updateGridVel() {
    Eigen::Vector3f minCorner = grid.center - 0.5f * grid.dimension;
    Eigen::Vector3f maxCorner = grid.center + 0.5f * grid.dimension;

    float spacing = 0.14f;
    Eigen::Vector3f dim(12.f, 12.f, 12.f);
    Eigen::Vector3f origin = dim * (-0.5f * spacing);
    float sphereRadius = 0.25f * std::min({ dim.x(), dim.y(), dim.z() }) * spacing;
    Eigen::Vector3f sphereCenter = origin + 0.5f * dim * spacing;
    sphereCenter.y() -= 1.f;

    // EXPLICIT UPDATE JUST TO TEST
    for (GridNode& g : grid.gridNodes) {
        if (g.mass > 0.f) {
            // COMPUTE VEL FROM GRID FORCES
            g.prevVelocity = g.velocity;
            g.velocity += stepSize * (1.f / g.mass) * g.force;

            // // --- SDF COLLISIONS ---
            if (g.mass > 0.f) {
                g.prevVelocity = g.velocity;
                g.velocity += stepSize * (1.f / g.mass) * g.force;

                float sdfVal = querySdf(g.worldPos);
                if (sdfVal < 0.f) {
                    Eigen::Vector3f normal = sdfNormal(g.worldPos);

                    float vn = g.velocity.dot(normal);
                    if (vn < 0.f) {
                        float restitution = 1.f;
                        g.velocity -= (1.f + restitution) * vn * normal;
                    }
                }
            }
            // CLAMP VELOCITY AT BOUNDS
            for (int i = 0; i < 3; i++) {
                if (i == 1 && g.worldPos[1] <= groundPlaneY) {
                    g.velocity[i] = 0.f;
                }
                if (g.worldPos[i] <= minCorner[i] || g.worldPos[i] >= maxCorner[i]) {
                    g.velocity[i] = 0.f;
                }
            }
        }
    }
}



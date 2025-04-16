#include "MPMSolver.h"
#include <Eigen/Dense>

const float dt = 0.01f; // Time step
const float gravity = -9.81f; // Gravity acceleration
const float cohesionStrength = 0.1f; // Cohesion force factor

MPMSolver::MPMSolver(Eigen::Vector3f gridDim, float spacing, Eigen::Vector3f gridOrigin, float dt,
    float critCompression, float critStretch, float hardeningCoeff,
    float initialDensity, float youngsMod, float poissonRatio)
    : stepSize(dt), grid(Eigen::Vector3f(gridDim), spacing, Eigen::Vector3f(gridOrigin)),
    critCompression(critCompression), critStretch(critStretch), hardeningCoeff(hardeningCoeff),
    initialDensity(initialDensity), youngsMod(youngsMod), poissonRatio(poissonRatio)
{
 
    mu0 = youngsMod / (2.f * (1.f + poissonRatio));
    lambda0 = (youngsMod * poissonRatio) / ((1.f + poissonRatio) + (1.f - 2.f * poissonRatio));
}

void MPMSolver::addParticle(const particle& particle) {
    particles.push_back(particle);
}


// THIS IS ORIGINAL PARTICLE STUFF ?? IDK IT WAS HERE WHEN I GOT HERE
void MPMSolver::computeForcesAndIntegrate() {
    for (particle& p : particles) {
        Eigen::Vector3f force = computeGravity(p) + computeCohesion(p);
        integrate(p, force);
    }
}

Eigen::Vector3f MPMSolver::computeGravity(const particle& p) {
    return Eigen::Vector3f(0.0, p.mass * gravity, 0.0);
}

Eigen::Vector3f MPMSolver::computeCohesion(const particle& p) {
    Eigen::Vector3f cohesionForce(0.0, 0.0, 0.0);
    for (const particle& neighbor : particles) {
        if (&p == &neighbor) continue;
        Eigen::Vector3f diff = neighbor.position - p.position;
        float distance = diff.size();
        if (distance > 0 && distance < grid.spacing) {
            cohesionForce += (diff * (1.f / diff.size())) * cohesionStrength / distance;
        }
    }
    return cohesionForce;
}

void MPMSolver::integrate(particle& p, Eigen::Vector3f force) {
    Eigen::Vector3f acceleration = force / p.mass;
    p.velocity += acceleration * dt;
    p.position += p.velocity * dt;
}

const vector<particle>& MPMSolver::getParticles() const {
    return particles;
}

void MPMSolver::step() {
    grid.clearGrid();
    particleToGridTransfer();
    computeForce();
    updateGridVel();
    updateParticleDefGrad();
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

static float weightFunGradient(float x) {
    x = abs(x);
    if (x > 0.f && x < 1.f) {
        return 3.f * abs(x) / (2.f * x) - 2.f * x;
    }
    else if (x >= 1.f && x < 2.f) {
        return -abs(x) / (2.f * x) * x * x + 2.f * x - 2.f * (abs(x) / x);
    }
    else {
        return 0.0f;
    }
}

// [======] PARTICLE FUNCTIONS [======]

void MPMSolver::computeSigma() {
    for (particle& p : particles) {
        Eigen::Matrix3f Fe = p.FE;// Eigen::Matrix3f::Identity();
        Eigen::Matrix3f Fp = p.FP;// Eigen::Matrix3f::Identity();

        // COMPUTE POLAR DECOMP TO GET ROTATIONAL PART OF F
        Eigen::JacobiSVD<Eigen::Matrix3f> svd(Fe, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Matrix3f U = svd.matrixU();
        Eigen::Matrix3f V = svd.matrixV();

        Eigen::Matrix3f R = U * V.transpose();
        float Je = Fe.determinant();
        float Jp = Fp.determinant();


        float xi = 10.0f; // default value, tweak as needed
 
        // Plastic-hardening-modified Lame parameters
        float mu = mu0 * std::exp(xi * (1.0f - Jp));
        float lambda = lambda0 * std::exp(xi * (1.0f - Jp));
 
        Eigen::Matrix3f strain = Fe - R;
        Eigen::Matrix3f sigma = 2.0f * mu * strain + lambda * (Je - 1.0f) * Eigen::Matrix3f::Identity();
        sigma *= (1.0f / Je); // Cauchy stress

        p.sigma = sigma;
    }
}

void MPMSolver::updateParticleDefGrad() {
    //Fn+1 = (I + ∆t∇v^n+1_p)Fn

    float xMin = grid.center[0] - 0.5f * grid.dimension[0];
    float yMin = grid.center[1] - 0.5f * grid.dimension[1];
    float zMin = grid.center[2] - 0.5f * grid.dimension[2];

    for (particle& p : particles) {
        Eigen::Matrix3f velGrad = Eigen::Matrix3f(0.f);

        // WORLD SPACE POSITIONS
        float x = p.position[0];
        float y = p.position[1];
        float z = p.position[2];

        int i = static_cast<int>(std::floor((x - xMin) / grid.spacing));
        int j = static_cast<int>(std::floor((y - yMin) / grid.spacing));
        int k = static_cast<int>(std::floor((z - zMin) / grid.spacing));

        Eigen::Vector3f vPic = Eigen::Vector3f(0.f);
        Eigen::Vector3f vFlip = Eigen::Vector3f(0.f);

        // YOU NEED TO DO THIS FOR ALL CELLS WITHIN SOME RADIUS
        for (int di = -2; di < 2; di++) {
            for (int dj = -2; dj < 2; dj++) {
                for (int dk = -2; dk < 2; dk++) {
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

                    Eigen::Vector3f gradWeight;
                    gradWeight[0] = 1.f / grid.spacing * weightFunGradient(xGrid) * weightFun(yGrid) * weightFun(zGrid);
                    gradWeight[1] = 1.f / grid.spacing * weightFun(xGrid) * weightFunGradient(yGrid) * weightFun(zGrid);
                    gradWeight[2] = 1.f / grid.spacing * weightFun(xGrid) * weightFun(yGrid) * weightFunGradient(zGrid);

                    // THIS IS v * gradW^T aka outer product
                    velGrad += curNode.velocity * gradWeight.transpose();

                    vPic += curNode.velocity * weight;
                    vFlip += (curNode.velocity - curNode.prevVelocity) * weight;
                }
            }
        }

        vFlip += p.velocity;
        float alpha = 0.95f; // Recommended 0.95
 
        p.velocity = (1.0f - alpha) * vPic + alpha * vFlip;


        // UPDATE DEFORMATION GRADIENT
        p.FE = (Eigen::Matrix3f::Identity() + stepSize * velGrad) * p.FE;

        // Predict the new total deformation gradient (FE * FP)
        Eigen::Matrix3f F_total = (Eigen::Matrix3f::Identity() + stepSize * velGrad) * p.FE * p.FP;
        
        // Predict the new elastic deformation gradient
        Eigen::Matrix3f FE_hat = (Eigen::Matrix3f::Identity() + stepSize * velGrad) * p.FE;

        // Convert to Eigen for SVD
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
        for (int i = 0; i < 3; ++i)
            Sigma_clamped(i, i) = sigma_clamped[i];

        Eigen::Matrix3f FE_new = U * Sigma_clamped * V.transpose();
        Eigen::Matrix3f F_total_eigen = F_total;


        // Update FP using: FP = V * Σ⁻¹ * Uᵀ * F_total
        Eigen::Matrix3f Sigma_inv = Eigen::Matrix3f::Zero();
        for (int i = 0; i < 3; ++i)
            Sigma_inv(i, i) = 1.0f / sigma_clamped[i];

        Eigen::Matrix3f FP_new = V * Sigma_inv * U.transpose() * F_total_eigen;

        p.FE = FE_new;
        p.FP = FP_new;

        // UPDATE POINT POSITIONS
        p.position += stepSize * p.velocity;
    }

    Eigen::Vector3f minCorner = grid.center - 0.5f * grid.dimension;
    Eigen::Vector3f maxCorner = grid.center + 0.5f * grid.dimension;

    for (particle& p : particles) {
        for (int axis = 0; axis < 3; ++axis) {
            if (p.position[axis] < minCorner[axis]) {
                p.position[axis] = minCorner[axis];
                if (p.velocity[axis] < 0.f) {
                    p.velocity[axis] *= -damping;
                }
            }

            if (p.position[axis] > maxCorner[axis]) {
                p.position[axis] = maxCorner[axis];
                if (p.velocity[axis] > 0.f) {
                    p.velocity[axis] *= -damping;
                }
            }
        }
    }

}



// [======] GRID FUNCTIONS [======]

void MPMSolver::particleToGridTransfer() {
    float xMin = grid.center[0] - 0.5f * grid.dimension[0];
    float yMin = grid.center[1] - 0.5f * grid.dimension[1];
    float zMin = grid.center[2]- 0.5f * grid.dimension[2];

    grid.clearGrid();
    for (particle& p : particles) {
        // World space positions
        float x = p.position[0];
        float y = p.position[1];
        float z = p.position[2];

        int i = static_cast<int>(std::floor((x - xMin) / grid.spacing));
        int j = static_cast<int>(std::floor((y - yMin) / grid.spacing));
        int k = static_cast<int>(std::floor((z - zMin) / grid.spacing));

        // YOU NEED TO DO THIS FOR ALL CELLS WITHIN SOME RADIUS
        for (int di = -2; di < 2; di++) {
            for (int dj = -2; dj < 2; dj++) {
                for (int dk = -2; dk < 2; dk++) {
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

                    if (p.velocity.size() > 1e-6f) { // Only if velocity is non-zero
                        curNode.velocityMass += curNode.mass;
                    }
                }
            }
        }
    }

    // Normalize grid node velocities using only effective mass
    for (GridNode& node : grid.gridNodes) {
        if (node.velocityMass > 0.f) {
            node.velocity /= node.velocityMass;
        }
        else {
            node.velocity = Eigen::Vector3f(0.f);
        }
    }

    // Currently the velocity stored is actually the total weighted momentum
    // To convert it to actual vel we divide each gridCell by its mass
    //grid.divideMass();
}

// THIS SHOULD ONLY BE CALLED ONCE AT t=0
void MPMSolver::computeInitialDensity() {
    float xMin = grid.center[0] - 0.5f * grid.dimension[0];
    float yMin = grid.center[1] - 0.5f * grid.dimension[1];
    float zMin = grid.center[2] - 0.5f * grid.dimension[2];

    grid.clearGrid();
    for (particle& p : particles) {
        // World space positions
        float x = p.position[0];
        float y = p.position[1];
        float z = p.position[2];

        int i = static_cast<int>(std::floor((x - xMin) / grid.spacing));
        int j = static_cast<int>(std::floor((y - yMin) / grid.spacing));
        int k = static_cast<int>(std::floor((z - zMin) / grid.spacing));

        // YOU NEED TO DO THIS FOR ALL CELLS WITHIN SOME RADIUS
        for (int di = -2; di < 2; di++) {
            for (int dj = -2; dj < 2; dj++) {
                for (int dk = -2; dk < 2; dk++) {
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
        g.density = (g.mass / volume);
    }

    // TRANSFER GRID DENSITY TO PARTICLES
    for (particle& p : particles) {
        // World space positions
        float x = p.position[0];
        float y = p.position[1];
        float z = p.position[2];

        int i = static_cast<int>(std::floor((x - xMin) / grid.spacing));
        int j = static_cast<int>(std::floor((y - yMin) / grid.spacing));
        int k = static_cast<int>(std::floor((z - zMin) / grid.spacing));


        // YOU NEED TO DO THIS FOR ALL CELLS WITHIN SOME RADIUS
        for (int di = -2; di < 2; di++) {
            for (int dj = -2; dj < 2; dj++) {
                for (int dk = -2; dk < 2; dk++) {
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

    float xMin = grid.center[0] - 0.5f * grid.dimension[0];
    float yMin = grid.center[1] - 0.5f * grid.dimension[1];
    float zMin = grid.center[2] - 0.5f * grid.dimension[2];

    // UPDATE CAUCHY STRESS FOR EACH PARTICLE
    computeSigma();

    for (particle& p : particles) {
        // World space positions
        float x = p.position[0];
        float y = p.position[1];
        float z = p.position[2];

        int i = static_cast<int>(std::floor((x - xMin) / grid.spacing));
        int j = static_cast<int>(std::floor((y - yMin) / grid.spacing));
        int k = static_cast<int>(std::floor((z - zMin) / grid.spacing));


        // YOU NEED TO DO THIS FOR ALL CELLS WITHIN SOME RADIUS
        for (int di = -2; di < 2; di++) {
            for (int dj = -2; dj < 2; dj++) {
                for (int dk = -2; dk < 2; dk++) {
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

                    curNode.force -= p.volume * p.sigma * gradWeight + (Eigen::Vector3f(0.f, gravity, 0.f) * p.mass);
                }
            }
        }
    }
}


// QUICK NOTE ABOUT THIS: The full implementation described in the paper uses a semi-implicit method for solving
//                        which results in a more stable and complex system. Since this is a push to get something out
//                        I am just doing an explicit solve with the most basic collision handling possible.
//                        For future implementation, fix this so that it fully implements the paper
void MPMSolver::updateGridVel() {
    Eigen::Vector3f minCorner = grid.center - 0.5f * grid.dimension;
    Eigen::Vector3f maxCorner = grid.center + 0.5f * grid.dimension;
    // EXPLICIT UPDATE JUST TO TEST
    for (GridNode& g : grid.gridNodes) {
        if (g.mass > 0.f) {
            // COMPUTE VEL FROM GRID FORCES
            g.prevVelocity = g.velocity;
            g.velocity += stepSize * (1.f / g.mass) * g.force;

            // VERY SIMPLE BOUNDING COLLISION
            // CLAMP VELOCITY AT BOUNDS
            if (g.worldPos[0] <= minCorner[0] || g.worldPos[0] >= maxCorner[0]) {
                g.velocity[0] = 0.f;
            }

            if (g.worldPos[1] <= minCorner[1] || g.worldPos[1] >= maxCorner[1]) {
                g.velocity[1] = 0.f;
            }

            if (g.worldPos[2] <= minCorner[2] || g.worldPos[2] >= maxCorner[2]) {
                g.velocity[2] = 0.f;
            }
        }
    }
}


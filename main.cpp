#include "surfaceFitting.h"

int main(){
    std::unique_ptr<FileReader> pc1 = std::make_unique<FileReader>("pc1.csv");
    std::unique_ptr<FileReader> pc2 = std::make_unique<FileReader>("pc2.csv");

    std::unique_ptr<SurfaceFitting> surface1 = std::make_unique<SurfaceFitting>("pc1", pc1->x, pc1->y, pc1->z);
    Eigen::Matrix3d CovarianceMatrix = surface1->covarianceMatrix();
    auto [directionNormal, magnitudeNormal] = surface1->computeNormal(CovarianceMatrix);

    std::cout << "directionNormal: " << directionNormal << "\n";
    std::cout << "magnitudeNormal: " << magnitudeNormal << "\n";
    
    surface1->leastSquaredFitting();

    std::unique_ptr<SurfaceFitting> surface2 = std::make_unique<SurfaceFitting>("pc2", pc2->x, pc2->y, pc2->z);
    surface2->leastSquaredFitting();

    surface1->totalLeastSquared();
    surface2->totalLeastSquared();

    surface1->ransac();
    surface2->ransac();
}
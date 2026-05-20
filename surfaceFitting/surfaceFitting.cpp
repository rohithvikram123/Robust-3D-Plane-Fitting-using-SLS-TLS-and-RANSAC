#include "surfaceFitting.h"

void FileReader::readFile(){
    std::ifstream file(filename);
    if (!file){
        std::cout << "Error Opening file" << "\n";
        return;
    }

    std::string line;

    while (std::getline(file, line)){
        std::stringstream ss(line);
        std::string value;
        std::vector<double> row;

        while (std::getline(ss, value, ',')){
            row.push_back(std::stod(value));
        }
        
        if (row.size() == 3){
            x_vals.push_back(row[0]);
            y_vals.push_back(row[1]);
            z_vals.push_back(row[2]);
        }
    }
    file.close();

    x = Eigen::Map<Eigen::VectorXd>(x_vals.data(), x_vals.size());
    y = Eigen::Map<Eigen::VectorXd>(y_vals.data(), y_vals.size());
    z = Eigen::Map<Eigen::VectorXd>(z_vals.data(), z_vals.size());

}

Eigen::Matrix3d SurfaceFitting::covarianceMatrix(){
    double x_variance = (x_.array() - x_mean_).square().sum() / x_.size();
    double y_variance = (y_.array() - y_mean_).square().sum() / y_.size();
    double z_variance = (z_.array() - z_mean_).square().sum() / z_.size();

    // double xy_covariance = ((x_ - Eigen::VectorXd::Constant(x_.size(), x_mean_)).transpose() * (y_ - Eigen::VectorXd::Constant(y_.size(), y_mean_)))(0,0)/x_.size();
    Eigen::VectorXd dx = x_.array() - x_mean_;
    Eigen::VectorXd dy = y_.array() - y_mean_;
    Eigen::VectorXd dz = z_.array() - z_mean_;

    double xy_covariance = dx.dot(dy) / x_.size();
    double xz_covariance = dx.dot(dz) / x_.size();
    double yz_covariance = dy.dot(dz) / x_.size();

    Eigen::Matrix3d covarianceMatrix;
    covarianceMatrix << x_variance, xy_covariance, xz_covariance,
                        xy_covariance, y_variance, yz_covariance,
                        xz_covariance, yz_covariance, z_variance;

    return covarianceMatrix;
}

std::pair<Eigen::Vector3d, double> SurfaceFitting::computeNormal(Eigen::Matrix3d covarianceMatrix){
// void SurfaceFitting::computeNormal(Eigen::Matrix3d covarianceMatrix){
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(covarianceMatrix);
    Eigen::Matrix3d eigenVectors = solver.eigenvectors();
    Eigen::Vector3d eigenValues = solver.eigenvalues();

    int minIndex;
    double magnitudeNormal = eigenValues.minCoeff(&minIndex);
    Eigen::Vector3d directionNormal = eigenVectors.col(minIndex);

    if (directionNormal(2) < 0){
        directionNormal = -directionNormal;
    }
    
    return {directionNormal, magnitudeNormal};
}

Eigen::MatrixXd SurfaceFitting::pseudoInverse(Eigen::MatrixXd X_matrix){
    Eigen::MatrixXd X_TX = X_matrix.transpose() * X_matrix;
    X_TX = 0.5 * (X_TX + X_TX.transpose());  // enforce geometry
     
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver_1(X_TX);
    Eigen::MatrixXd V = solver_1.eigenvectors();
    Eigen::VectorXd singular_values = solver_1.eigenvalues();
    
    // We are going to reverse the singular, U and V matrices since the the singular values should be from greater to smaller and not in asscending order
    singular_values = singular_values.reverse().eval();
    V = V.rowwise().reverse().eval();

    singular_values = singular_values.cwiseSqrt();
    // std::cout << "Singular values: " << singular_values << "\n";
    
    Eigen::MatrixXd singular_inverseMatrix = Eigen::MatrixXd::Zero(V.cols(), V.cols());

    // std::cout << "U_matrix: " << U << "\n";
    // std::cout << "V_matrix: " << V << "\n";

    for (int i{0}; i < singular_values.size(); i++){
        singular_inverseMatrix(i, i) = (1/singular_values(i));
    }
    // std::cout << "Singular_inverseMatrix: \n" << singular_inverseMatrix << "\n";
    
    // Computing thin U
    Eigen::MatrixXd U = X_matrix * V * singular_inverseMatrix;
    // std::cout << "U_matrix: " << U << "\n";
    Eigen::MatrixXd X_pseudoInverse = V * singular_inverseMatrix * U.transpose();

    // std::cout << "X_pseudoInverse" << X_pseudoInverse << "\n";
    
    return X_pseudoInverse;
}

void SurfaceFitting::plot(std::string fittingType, std::vector<double>& x_ovec, std::vector<double>& y_ovec, std::vector<double>& z_ovec, std::vector<std::vector<double>>& x_vec, std::vector<std::vector<double>>& y_vec, std::vector<std::vector<double>>& z_vec){
    plt::figure();
    plt::scatter(x_ovec, y_ovec, z_ovec, 10.0);
    plt::title(fittingType + "- Scatter");
    plt::save(data_name_ + "/" + fittingType + "_scatter.png");
    plt::close();

    plt::figure();
    plt::plot_surface(x_vec, y_vec, z_vec, {{"alpha", "1.0"}});
    plt::title(fittingType + "- Surface");
    plt::save(data_name_ + "/" + fittingType + "_surface.png");
    plt::close();
}

void SurfaceFitting::leastSquaredFitting(){
    Eigen::MatrixXd X_matrix(x_.size(), 3);
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(x_.rows());

    X_matrix.col(0) = x_;
    X_matrix.col(1) = y_;
    X_matrix.col(2) = ones;

    Eigen::MatrixXd X_inverseMatrix = pseudoInverse(X_matrix);

    std::cout << "X_inverseMatrix: " << "col: " << X_inverseMatrix.cols() << "rows: " << X_inverseMatrix.rows() << "\n";

    Eigen::VectorXd B_matrix = X_inverseMatrix * z_;

    std::cout << "B_matrix: " << B_matrix << "\n";

    int N{300};

    double x_min = x_.minCoeff();
    double x_max = x_.maxCoeff();
    double y_min = y_.minCoeff();
    double y_max = y_.maxCoeff();

    std::vector<double> x_ovec(x_.data(), x_.data() + x_.size());
    std::vector<double> y_ovec(y_.data(), y_.data() + y_.size());
    std::vector<double> z_ovec(z_.data(), z_.data() + z_.size());

    Eigen::VectorXd x_plot = Eigen::VectorXd::LinSpaced(N, x_min, x_max);
    Eigen::VectorXd y_plot = Eigen::VectorXd::LinSpaced(N, y_min, y_max);

    std::vector<double> x_vec(x_plot.data(), x_plot.data() + x_plot.size());
    std::vector<double> y_vec(y_plot.data(), y_plot.data() + y_plot.size());

    std::vector<std::vector<double>> X, Y, Z;

    for (int i{0}; i < x_vec.size(); i++){
        std::vector<double> x_row, y_row, z_row;

        for (int j{0}; j < y_vec.size(); j++){
            double xv = x_vec[i];
            double yv = y_vec[j];

            x_row.push_back(xv);
            y_row.push_back(yv);

            double zv = (xv * B_matrix(0)) + (yv * B_matrix(1)) + B_matrix(2);
            z_row.push_back(zv);
        }
        X.push_back(x_row);
        Y.push_back(y_row);
        Z.push_back(z_row);
    }

    plot("Least_Squared", x_ovec, y_ovec, z_ovec, X, Y, Z);
    
}

void SurfaceFitting::totalLeastSquared(){
    Eigen::MatrixXd X_matrix(x_.size(), 3);
    X_matrix.col(0) = x_.array() - x_mean_;
    X_matrix.col(1) = y_.array() - y_mean_;
    X_matrix.col(2) = z_.array() - z_mean_;

    Eigen::Matrix3d C_matrix = X_matrix.transpose() * X_matrix;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(C_matrix);
    Eigen::MatrixXd eigen_vectors = solver.eigenvectors();
    tls_coeffs_ = eigen_vectors.col(0);
    
    int N{x_.size()};

    double x_min = x_.minCoeff();
    double x_max = x_.maxCoeff();
    double y_min = y_.minCoeff();
    double y_max = y_.maxCoeff();

    std::vector<double> x_ovec(x_.data(), x_.data() + x_.size());
    std::vector<double> y_ovec(y_.data(), y_.data() + y_.size());
    std::vector<double> z_ovec(z_.data(), z_.data() + z_.size());

    Eigen::VectorXd x_plot = Eigen::VectorXd::LinSpaced(N, x_min, x_max);
    Eigen::VectorXd y_plot = Eigen::VectorXd::LinSpaced(N, y_min, y_max);

    std::vector<double> x_vec(x_plot.data(), x_plot.data() + x_plot.size());
    std::vector<double> y_vec(y_plot.data(), y_plot.data() + y_plot.size());

    std::vector<std::vector<double>> X, Y, Z;

    for (int i{0}; i < x_vec.size(); i++){
        std::vector<double> x_row, y_row, z_row;

        for (int j{0}; j < y_vec.size(); j++){
            double xv = x_vec[i];
            double yv = y_vec[j];

            x_row.push_back(xv);
            y_row.push_back(yv);

            double zv = z_mean_ -(tls_coeffs_(0)/tls_coeffs_(2))*(xv - x_mean_) - (tls_coeffs_(1)/tls_coeffs_(2))*(yv - y_mean_);
            z_row.push_back(zv);
        }
        X.push_back(x_row);
        Y.push_back(y_row);
        Z.push_back(z_row);
    }
    plot("Total_Least_Squared", x_ovec, y_ovec, z_ovec, X, Y, Z);
}

void SurfaceFitting::ransac(){
    Eigen::Vector3d unit_normal = tls_coeffs_.normalized();
    Eigen::Vector3d centroid(x_mean_, y_mean_, z_mean_);

    double d = -unit_normal.dot(centroid);
    std::vector<Eigen::Vector3d> points;

    for (int i{0}; i < x_.size(); i++){
        points.emplace_back(x_(i), y_(i), z_(i));
    }

    std::vector<double> residuals;

    for (const auto &p: points){
        double r = unit_normal.dot(p) + d;
        residuals.push_back(r);
    }

    // Using MAD(Median Absolute Deviation) instead of mean, because it is robust to outliers and extream values.

    std::sort(residuals.begin(), residuals.end());

    double median = residuals[residuals.size()/2];

    std::vector<double> absolute_deviation;

    for(double r: residuals){
        absolute_deviation.push_back(std::abs(r - median));
    }

    std::sort(absolute_deviation.begin(), absolute_deviation.end());

    double mad = absolute_deviation[absolute_deviation.size()/2];

    double std_dev = 1.4826 * mad;

    double threshold = 1.96 * std_dev;

    int max_iterations{1000}; // setting a very large iterations value (Adaptive ransac procedure)

    int sample_count{0};
    int max_inlier_count{0};
    double e_outlier_prob;
    while (max_iterations > sample_count){
        // Using random integer generator to select 3 unique indices from the point vector
        std::random_device rd; //seed
        std::mt19937 gen(rd()); // random engine
        
        std::uniform_int_distribution<> dist(0, points.size() - 1);

        int i1 = dist(gen);
        int i2 = dist(gen);
        int i3 = dist(gen);

        if(i1 == i2 || i2 == i3 || i1 == i3){
            continue;
        }

        Eigen::Vector3d p1 = points[i1];
        Eigen::Vector3d p2 = points[i2];
        Eigen::Vector3d p3 = points[i3];

        Eigen::Vector3d s_normal = (p3 - p1).cross(p2 - p1);

        if (s_normal.norm() < 1e-9){      // checking if the points are colinear
            continue;
        }

        s_normal.normalize();

        double s_d = -s_normal.dot(p1);

        int inlier_count{0};

        for(const auto &p: points){
            double dist = std::abs(s_normal.dot(p) + s_d);

            if(dist < threshold){
                inlier_count++;
            }
        }
        std::cout << "inlier_count: " << inlier_count <<"\n"; 
        std::cout << "points_size: " << points.size();

        if (inlier_count > max_inlier_count){
            max_inlier_count = inlier_count;
            e_outlier_prob = 1 - (static_cast<double>(inlier_count)/points.size());

            std::cout << "oulier prob: " << e_outlier_prob << "\n";

            max_iterations = std::log(1-desired_prob)/std::log(1 - std::pow((1-e_outlier_prob),s_min));
            ransac_normal = s_normal;
            ransac_d = s_d;
            
        }
        sample_count++;
        std::cout << "max_iterations: " << max_iterations << ", sample_count: " << sample_count << "\n";
    }

    int N{x_.size()};

    double x_min = x_.minCoeff();
    double x_max = x_.maxCoeff();
    double y_min = y_.minCoeff();
    double y_max = y_.maxCoeff();

    std::vector<double> x_ovec(x_.data(), x_.data() + x_.size());
    std::vector<double> y_ovec(y_.data(), y_.data() + y_.size());
    std::vector<double> z_ovec(z_.data(), z_.data() + z_.size());

    Eigen::VectorXd x_plot = Eigen::VectorXd::LinSpaced(N, x_min, x_max);
    Eigen::VectorXd y_plot = Eigen::VectorXd::LinSpaced(N, y_min, y_max);

    std::vector<double> x_vec(x_plot.data(), x_plot.data() + x_plot.size());
    std::vector<double> y_vec(y_plot.data(), y_plot.data() + y_plot.size());

    std::vector<std::vector<double>> X, Y, Z;

    for (int i{0}; i < x_vec.size(); i++){
        std::vector<double> x_row, y_row, z_row;

        for (int j{0}; j < y_vec.size(); j++){
            double xv = x_vec[i];
            double yv = y_vec[j];

            x_row.push_back(xv);
            y_row.push_back(yv);

            double zv = ((-ransac_normal(0) * xv) - (ransac_normal(1) * yv) - ransac_d)/ransac_normal(2);
            z_row.push_back(zv);
        }
        X.push_back(x_row);
        Y.push_back(y_row);
        Z.push_back(z_row);
    }

    plot("RANSAC", x_ovec, y_ovec, z_ovec, X, Y, Z);

}
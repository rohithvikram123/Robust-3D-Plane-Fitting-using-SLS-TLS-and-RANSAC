#ifndef __SURFACEFITTING__H__
#define __SURFACEFITTING__H__

#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <vector>
#include <cmath>
#include <random>
#include <Eigen/Dense>
#include <opencv2/opencv.hpp>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

class FileReader{
    public:
        // Constructor
        FileReader(const std::string& fname):filename{fname}{
            readFile();
        }

        Eigen::VectorXd x;
        Eigen::VectorXd y;
        Eigen::VectorXd z;
        std::string filename;
            
    private:
        std::vector<double> x_vals;
        std::vector<double> y_vals;
        std::vector<double> z_vals;
        
        void readFile();
};

class SurfaceFitting{
    public:
        SurfaceFitting(const std::string data_name, const Eigen::VectorXd& x, const Eigen::VectorXd& y, const Eigen::VectorXd& z): data_name_{data_name}, x_{x}, y_{y}, z_{z}, x_mean_{x_.mean()}, y_mean_{y_.mean()}, z_mean_{z_.mean()}{
        }

        Eigen::Matrix3d covarianceMatrix();
        std::pair<Eigen::Vector3d, double> computeNormal(Eigen::Matrix3d covarianceMatrix);
        void leastSquaredFitting();
        void totalLeastSquared();
        void ransac();
        std::string data_name_;

    private:
        Eigen::VectorXd x_;
        Eigen::VectorXd y_;
        Eigen::VectorXd z_;
        double x_mean_;
        double y_mean_;
        double z_mean_;
        Eigen::Vector3d tls_coeffs_;
        Eigen::Vector3d ransac_normal;
        double ransac_d;
        Eigen::MatrixXd pseudoInverse(Eigen::MatrixXd X_matrix);

        // RANSAC parameters
        int s_min{3}; // minimum number of points required to fit a plane
        double desired_prob{0.99}; // desired probability that we get a good sample
        
        
        void plot(std::string fittingType, std::vector<double>& x_ovec, std::vector<double>& y_ovec, std::vector<double>& z_ovec, std::vector<std::vector<double>>& x_vec, std::vector<std::vector<double>>& y_vec, std::vector<std::vector<double>>& z_vec);
};


#endif
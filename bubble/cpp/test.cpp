#include <iostream>
#include <armadillo>

int main() {
    arma::mat A = arma::randu<arma::mat>(4, 5);
    std::cout << "Armadillo version: " << arma::arma_version::as_string() << std::endl;
    return 0;
}

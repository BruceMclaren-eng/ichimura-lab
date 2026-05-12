#include <iostream>
#include <vector>

double dot_product(const std::vector<double>& x, const std::vector<double>& y) {
    double result = 0.0;
    for (int i = 0; i < x.size(); i++) {
        result += x[i] * y[i];
    }
    return result;
}

void scale(std::vector<double>& x, double alpha) {
    for (int i = 0; i < x.size(); i++) {
        x[i] *= alpha;
    }
}

void print_vector(const std::vector<double>& x) {
    for (int i = 0; i < x.size(); i++) {
        std::cout << "x[" << i << "] =" << x[i] << "\n";
    }
}

int main() {
    std::vector<double> x = {1.0, 2.0, 3.0};
    std::vector<double> y = {4.0, 5.0, 6.0};

    std::cout << "dot product: " << dot_product(x,y) << "\n";

    scale(x, 2.0);
    std::cout << "scaled x:\n";
    print_vector(x);

    return 0;
}
#include <fstream>
#include <iostream>

int main(int argc, char *argv[])
{
    if (argc < 2) return 1;
    
    std::ifstream binary_float_input(argv[1], std::ios::binary);

    std::cout.precision(14);
    std::cout.setf(std::ios::scientific, std::ios::floatfield);
    
    while (binary_float_input) {
        double temp;
        binary_float_input.read(reinterpret_cast<char*>(&temp), sizeof(double));
        if (binary_float_input) {
            //no failure
            std::cout << temp << "\n";
        }
    }
    
    return 0;
}

#include "NJL_2flavor.h"
#include <fstream>
#include <sstream>
#include <string>

int main()
{
    std::ifstream infile("input.txt");
    std::string line;
    int case_id;

    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        std::string key;
        std::string value;

        std::getline(iss, key, '=');
        std::getline(iss, value);

        if (key == "case")
        {
            case_id = std::stoi(value);
        }
    }

    NJL_2flavors::Constants A;
    A.compute_NJL2(case_id);

    return 0;
}

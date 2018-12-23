#ifndef ANA_NDD_UTILS_H
#define ANA_NDD_UTILS_H

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

using std::size_t;
namespace ANA {
namespace NDD {

    class Modes {
    public:
        Modes(std::string const &filename);

        void get_modes_from_raw(std::string_view const texto);

        std::vector<std::vector<float>> evectors;
        std::vector<float> evals;
        size_t i, j;
    };

} // namespace NDD
} // namespace ANA

#endif // _H
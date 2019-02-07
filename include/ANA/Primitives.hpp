#ifndef ANA_PRIMITIVES_H
#define ANA_PRIMITIVES_H
#include <ANA/CGALUtils.hpp>
#include <ANA/Includes.hpp>
#include <ANA/Options.hpp>

namespace ANA {

struct Points {
public:
    Points() = default;

    Points(int const i) : i(_cnt) {
        _x.reserve(_cnt);
        _y.reserve(_cnt);
        _z.reserve(_cnt);

        _index.reserve(_cnt);
        _radius.reserve(_cnt);
        _resn.reserve(_cnt);
        _resi.reserve(_cnt);
    }

    int const _cnt;
    std::vector<double> _x, _y, _z;

    std::vector<int> _index;
    std::vector<double> _radius;
    // atom's residue number
    std::vector<int> _resn;
    // atom's residue name in 3 letter format
    std::vector<std::string> _resi;
};

struct Tetrahedron {
public:
    Tetrahedron() noexcept = default;

    Tetrahedron(
        Point const &p0, Point const &p1, Point const &p2, Point const &p3) :
        _data({p0, p1, p2, p3}) {}

    Tetrahedron(Point &&p0, Point &&p1, Point &&p2, Point &&p3) :
        _data({p0, p1, p2, p3}) {}

    std::array<Point, 4> _data;
};

struct TriangularPrism {
public:
    TriangularPrism() noexcept = default;

    TriangularPrism(Point const &p0, Point const &p1, Point const &p2,
        Point const &p3, Point const &p4, Point const &p5) :
        _data({p0, p1, p2, p3, p4, p5}) {}

    TriangularPrism(Point &&p0, Point &&p1, Point &&p2, Point &&p3, Point &&p4,
        Point &&p5) :
        _data({p0, p1, p2, p3, p4, p5}) {}

    std::array<Point, 6> _data;
};
}

#endif // _H
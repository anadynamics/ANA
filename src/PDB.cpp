#include <ANA/PDB.hpp>

namespace ANA {

void draw(ConvexHull const &CH, std::string const &filename) {

    FILE *out_file = std::fopen(filename.c_str(), "w");
    if (out_file) {
        int idx = 1, resid = 1;
        for (auto const &triangle : CH._triangles) {
            draw(triangle, out_file, idx, resid);
        }
        connect_triangle(out_file, 1, resid);
    } else {
        printf("Could not open %s.\n", filename.c_str());
    }
    std::fclose(out_file);

    return;
}

void draw(Triangle const &t, FILE *out_file, int &idx, int &resid) {

    auto const i = idx++;
    auto const j = idx++;
    auto const k = idx++;

    draw(t.vertex(0), out_file, i, resid);
    draw(t.vertex(1), out_file, j, resid);
    draw(t.vertex(2), out_file, k, resid);
    ++resid;
    return;
}

void connect_triangle(FILE *out_file, int const first_t, int const last_t) {

    for (auto r = first_t; r < last_t; ++r) {
        auto const i = (r - 1) * 3 + 1;
        auto const j = i + 1;
        auto const k = i + 2;

        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", i, j, k);
        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", j, i, k);
        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4}\n", k, i, j);
    }
    return;
}

void draw(Cavity const &hueco, std::string const &filename) {
    FILE *out_file = std::fopen(filename.c_str(), "w");
    if (out_file) {
        int idx = 1, resid = 1;
        for (auto const &cell : hueco._included_cells) {
            draw(cell, out_file, idx, resid);
        }
        connect_cell(out_file, 1, resid);
        // for (auto const &polyhedron : hueco._border) {
        //     draw(polyhedron, out_file, idx, resid);
        // }
    } else {
        printf("Could not open %s.\n", filename.c_str());
    }
    std::fclose(out_file);
    return;
}

void draw(
    Finite_cells_iterator const &cell, FILE *out_file, int &idx, int &resid) {

    auto const i = idx++;
    auto const j = idx++;
    auto const k = idx++;
    auto const l = idx++;
    draw(cell->vertex(0)->point(), out_file, i, resid);
    draw(cell->vertex(1)->point(), out_file, j, resid);
    draw(cell->vertex(2)->point(), out_file, k, resid);
    draw(cell->vertex(3)->point(), out_file, l, resid);
    ++resid;
    return;
}

void draw(Polyhedron const &polyhedron, FILE *out_file, int &idx, int &resid) {

    auto v_end = polyhedron.vertices_end();
    for (auto v_ite = polyhedron.vertices_begin(); v_ite != v_end; ++v_ite) {
        draw(v_ite->point(), out_file, idx, resid);
    }
    ++resid;
    return;
}

void draw(Point const &punto, FILE *out_file, int const idx, int const resid) {
    fmt::print(out_file,
        "{: <6}{: >5} {: <4s} {:3} {:1}{: >4}    "
        "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {: >2s}\n",
        "HETATM", idx, "H", "ANA", "A", resid, punto.x(), punto.y(), punto.z(),
        1.0, 0.0, "H");
    return;
}

void connect_cell(FILE *out_file, int const first_cell, int const last_cell) {
    std::cout << first_cell << "  " << last_cell << '\n';
    assert(first_cell < last_cell);

    for (auto r = first_cell; r < last_cell; ++r) {
        auto const i = (r - 1) * 4 + 1;
        auto const j = i + 1;
        auto const k = i + 2;
        auto const l = i + 3;

        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", i, j, k, l);
        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", j, k, l, i);
        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", k, l, i, j);
        fmt::print(out_file, "CONECT {:>4} {:>4} {:>4} {:>4}\n", l, i, j, k);
    }
    return;
}

} // namespace ANA
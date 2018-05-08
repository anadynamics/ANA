// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include "catch.hpp"
#include "helpers.hpp"
#include "chemfiles.hpp"
using namespace chemfiles;

TEST_CASE("Read files in LAMMPS .lammpstrj format using Molfile"){
    Trajectory file("data/lammps/polymer.lammpstrj");
    Frame frame = file.read();
    double eps = 1e-3;

    CHECK(frame.natoms() == 1714);
    auto positions = frame.positions();
    CHECK(approx_eq(positions[0], Vector3D(51.8474, 100.348, 116.516), eps));
    // this one has a non zero image index (1 0 0)
    CHECK(approx_eq(positions[1189], Vector3D(116.829, 91.2404, 79.8858), eps));
}

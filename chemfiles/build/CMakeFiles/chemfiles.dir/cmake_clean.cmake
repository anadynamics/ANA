file(REMOVE_RECURSE
  "libchemfiles.pdb"
  "libchemfiles.a"
)

# Per-language clean rules from dependency scanning.
foreach(lang C CXX)
  include(CMakeFiles/chemfiles.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()

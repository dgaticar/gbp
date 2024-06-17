Compact Representations of Spatial Hierarchical Structures with Support for Topological Queries
=========

Requirements
------------

The implementations provided requires:

* Installing the following version of SDSL: https://github.com/jfuentess/sdsl-lite

------------

This project provides three succinct data structures to represent a multigranular model based on the use of rules of inference, plus two variant for the second and third approach.

The structures can be used as follows:

Structure based on Gloud

```cpp
#include "GloudO.hpp"


using namespace std;
using namespace sdsl;

gloud<> *Granularities;

int main(int argc, char **argv){
  // argv[1] is the path to a file with the subsumption graph.
  // argv[2] is the path to a file with the disjoint graph.
  // argv[3] is the path to a file with the not subsumption graph.
  // argv[4] is the path to a file with the not disjoint graph.
Granularities = new gloud<>(argv[1],argv[2],argv[3],argv[4]);

if(Granularities->reglasub(2,8))cout << "Granule 2 Inside in granule 8"\n;

if(Granularities->reglanotdis(2,8))cout << "Granule 2 not disjoint with granule 8"\n;

return 0;
}
```

Structure based on GBP
```cpp
#include "GloudBP.hpp"


using namespace std;
using namespace sdsl;

gbp<> *Granularities;

int main(int argc, char **argv){
  // argv[1] is the path to a file with the subsumption graph.
  // argv[2] is the path to a file with the disjoint graph.
  // argv[3] is the path to a file with the not subsumption graph.
  // argv[4] is the path to a file with the not disjoint graph.
Granularities = new gbp<>(argv[1],argv[2],argv[3],argv[4]);

if(Granularities->reglasub(2,8))cout << "Granule 2 Inside in granule 8"\n;

if(Granularities->reglanotdis(2,8))cout << "Granule 2 not disjoint with granule 8"\n;

return 0;
}
```

Structure based on GBP + matrices
```cpp
#include "GloudMatrix.hpp"


using namespace std;
using namespace sdsl;

gbp<> *Granularities;

int main(int argc, char **argv){
  // argv[1] is the path to a file with the subsumption graph.
  // argv[2] is the path to a file with the disjoint graph.
  // argv[3] is the path to a file with the not subsumption graph.
  // argv[4] is the path to a file with the not disjoint graph.
Granularities = new gbp<>(argv[1],argv[2],argv[3],argv[4]);

if(Granularities->reglasub(2,8))cout << "Granule 2 Inside in granule 8"\n;

if(Granularities->reglanotdis(2,8))cout << "Granule 2 not disjoint with granule 8"\n;

return 0;
}
```



to compile, just run:


```sh
g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib program.cpp -o program -lsdsl
```

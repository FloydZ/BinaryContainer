IMPORTANT NOTE
----
This is only a pre alpha version.

Usage
---
```C++
#include "binary_container"
using BC = BinaryContainer<1024, uint64_t>;
BC a,b,c;
    
```

Flags
----
The following compiler flags are supported:
```bash
USE_BRANCH_PREDICTION
USE_PREFETCH
USE_LOOP_UNROLL
```

Additionally if you include [M4RI](https://bitbucket.org/malb/m4ri) (alternatively you can define `M4RI_M4RI_H` ) before add this project supporting functions to convert vectors from or to M4RI.

TODO
----



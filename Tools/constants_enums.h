//#####################################################################
// Copyright 2015, Lakshman Anumolu.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef __CONSTANTS_ENUMS__
#define __CONSTANTS_ENUMS__

namespace PhysBAM {

// enums
enum TemporalDiscretization {
    FORWARD_EULER
};

enum SpatialDiscretization {
    GAUSS
};

enum SurfaceInterpolationScheme {
    INTERFACE_COMPRESSION,
    LINEAR,
    LIMITED_LINEAR,
    MINMOD,
    MUSCL,
    SUPERBEE,
    UPWIND,
    VANALBADA,
    VANLEER
};

}

#endif // __CONSTANTS_ENUMS__

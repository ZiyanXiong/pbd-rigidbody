#pragma once
#pragma once
#include "Utils.h"

namespace _2psp {

    class Contact {
    public:
        Vector3 _xi1;     // contact point in body1 frame.
        Vector3 _xi2;     // contact point in body2 frame.
        //Vector3 _xw;     // contact point in world frame.
        dtype _d;        // penetration depth.
        Vector3 _normal; // contact normal

        Contact(Vector3 xi1, Vector3 xi2, dtype d, Vector3 normal) {
            _xi1 = xi1;
            _xi2 = xi2;
            //_xw = xw;
            _d = d;
            _normal = normal;
        }

    };
}
#include <iostream>
#include <set>

#include <doctest/doctest.h>

#include "contour.h"

TEST_CASE("contour.basic")
{
    #define COLLIDE_POINT(c_, point_, result_, .../*allowed ray dirs*/) \
        do \
        { \
            [[maybe_unused]] auto point = point_; \
            [[maybe_unused]] auto col = c_.MakeCollider(point_); \
            REQUIRE(std::set<int>{__VA_ARGS__}.contains(col.DebugRayDirection())); \
            REQUIRE_EQ(col.CollidePoint(point), result_); \
        } \
        while (false)

    ContourShape<int> c;

    // An octagon with no axis-aligned edges.
    constexpr int halfsize = 100, gap = 20;
    c.AddLoop(std::array{
        ivec2(halfsize        ,                0),
        ivec2(halfsize*2 - gap,              gap),
        ivec2(halfsize*2      , halfsize        ),
        ivec2(halfsize*2 - gap, halfsize*2 - gap),
        ivec2(halfsize        , halfsize*2      ),
        ivec2(             gap, halfsize*2 - gap),
        ivec2(               0, halfsize        ),
        ivec2(             gap,              gap),
    });

    REQUIRE_EQ(c.Bounds(), ivec2(0).rect_size(next_value(ivec2(halfsize * 2))));

    COLLIDE_POINT(c, ivec2(191, 60), false, 0);
    COLLIDE_POINT(c, ivec2(190, 60), true , 0);
    COLLIDE_POINT(c, ivec2(189, 60), true , 0);
    COLLIDE_POINT(c, ivec2(191,140), false, 0);
    COLLIDE_POINT(c, ivec2(190,140), true , 0);
    COLLIDE_POINT(c, ivec2(189,140), true , 0);

    COLLIDE_POINT(c, ivec2( 60,191), false, 1);
    COLLIDE_POINT(c, ivec2( 60,190), true , 1);
    COLLIDE_POINT(c, ivec2( 60,189), true , 1);
    COLLIDE_POINT(c, ivec2(140,191), false, 1);
    COLLIDE_POINT(c, ivec2(140,190), true , 1);
    COLLIDE_POINT(c, ivec2(140,189), true , 1);

    COLLIDE_POINT(c, ivec2(  9, 60), false, 2);
    COLLIDE_POINT(c, ivec2( 10, 60), true , 2);
    COLLIDE_POINT(c, ivec2( 11, 60), true , 2);
    COLLIDE_POINT(c, ivec2(  9,140), false, 2);
    COLLIDE_POINT(c, ivec2( 10,140), true , 2);
    COLLIDE_POINT(c, ivec2( 11,140), true , 2);

    COLLIDE_POINT(c, ivec2( 60,  9), false, 3);
    COLLIDE_POINT(c, ivec2( 60, 10), true , 3);
    COLLIDE_POINT(c, ivec2( 60, 11), true , 3);
    COLLIDE_POINT(c, ivec2(140,  9), false, 3);
    COLLIDE_POINT(c, ivec2(140, 10), true , 3);
    COLLIDE_POINT(c, ivec2(140, 11), true , 3);

    COLLIDE_POINT(c, ivec2(190, 100), true, 0);
    COLLIDE_POINT(c, ivec2(100, 190), true, 1);
    COLLIDE_POINT(c, ivec2( 10, 100), true, 2);
    COLLIDE_POINT(c, ivec2(100,  10), true, 3);

    #undef COLLIDE_POINT
}

TEST_CASE("contour.pixelperfect")
{
    ContourShape<int> c;
    c.AddLoop(std::array{
        ivec2( 4, 0),
        ivec2( 8, 0),
        ivec2(12, 4),
        ivec2(12, 8),
        ivec2( 8,12),
        ivec2( 8,12),
        // #error finish this, and render to bitmap
    });
}

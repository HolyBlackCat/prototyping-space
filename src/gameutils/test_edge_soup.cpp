#include <iostream>
#include <set>

#include <doctest/doctest.h>

#include "edge_soup.h"

// Some basic point collision tests.
TEST_CASE("edge_soup.basic")
{
    #define COLLIDE_POINT(c_, point_, result_, .../*allowed ray dirs*/) \
        do \
        { \
            [[maybe_unused]] auto point = point_; \
            [[maybe_unused]] auto col = c_.MakePointCollider(point_); \
            REQUIRE(std::set<int>{__VA_ARGS__}.contains(col.DebugRayDirection())); \
            REQUIRE_EQ(col.CollidePoint(point), result_); \
        } \
        while (false)

    EdgeSoup<int> c;

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

// Rasterize using point collisions and compare against reference.
TEST_CASE("edge_soup.pixelperfect.basic")
{
    EdgeSoup<int> c;
    c.AddLoop(std::array{
        ivec2( 6, 0),
        ivec2(12, 6),
        ivec2( 6,12),
        ivec2( 0, 6),
    });
    // Hole:
    c.AddLoop(std::array{
        ivec2( 2, 6),
        ivec2( 6,10),
        ivec2(10, 6),
        ivec2( 6, 2),
    });

    constexpr ivec2 corner_a(-1), corner_b(13);

    for (int ray_dir = 0; ray_dir < 4; ray_dir++)
    {
        CAPTURE(ray_dir);

        std::string str;
        for (ivec2 pos : corner_a <= vector_range <= corner_b)
        {
            str += ".#"[c.MakePointColliderWithRayDir(ray_dir).CollidePoint(pos)];
            if (pos.x == corner_b.x)
                str += '\n';
        }

        REQUIRE_EQ(str, &1[R"(
...............
.......#.......
......###......
.....#####.....
....###.###....
...###...###...
..###.....###..
.###.......###.
..###.....###..
...###...###...
....###.###....
.....#####.....
......###......
.......#.......
...............
)"]);
    }
}

// Rasterize using point collisions and compare against reference.
// This additionally includes axis-aligned edges.
TEST_CASE("edge_soup.pixelperfect.axis_aligned")
{
    EdgeSoup<int> c;
    c.AddLoop(std::array{
        ivec2( 4, 0),
        ivec2( 8, 0),
        ivec2(12, 4),
        ivec2(12, 8),
        ivec2( 8,12),
        ivec2( 4,12),
        ivec2( 0, 8),
        ivec2( 0, 4),
    });
    // Hole:
    c.AddLoop(std::array{
        ivec2( 2, 5),
        ivec2( 2, 7),
        ivec2( 5,10),
        ivec2( 7,10),
        ivec2(10, 7),
        ivec2(10, 5),
        ivec2( 7, 2),
        ivec2( 5, 2),
    });

    constexpr ivec2 corner_a(-1), corner_b(13);

    for (int ray_dir = 0; ray_dir < 4; ray_dir++)
    {
        CAPTURE(ray_dir);

        std::string str;
        for (ivec2 pos : corner_a <= vector_range <= corner_b)
        {
            str += ".#"[c.MakePointColliderWithRayDir(ray_dir).CollidePoint(pos)];
            if (pos.x == corner_b.x)
                str += '\n';
        }

        REQUIRE_EQ(str, &1[R"(
...............
.....#####.....
....#######....
...#########...
..####...####..
.####.....####.
.###.......###.
.###.......###.
.###.......###.
.####.....####.
..####...####..
...#########...
....#######....
.....#####.....
...............
)"]);
    }
}

// Different kinds of overhangs.
TEST_CASE("edge_soup.overhands")
{
    EdgeSoup<int> c;
    c.AddLoop(std::array{
        ivec2(  0, 0),
        ivec2(  5,20),
        ivec2(  0,15), // Sharp edge, angles backwards.
        ivec2(  5,25),
        ivec2(  0,25), // This edge is parallel to the ray.
        ivec2(  5,30),
        ivec2(  0,35), // Dull edge, angles forwards.
        ivec2(  5,45),
        ivec2(  3,45), // An edge parallel to the ray, followed by a backward edge.
        ivec2(  0,42),
        ivec2(  0,60),
        ivec2(-10,60),
        ivec2(-10, 0),
    });

    constexpr ivec2 corner_a(-11,-1), corner_b(6,61);

    for (int ray_dir = 0; ray_dir < 4; ray_dir++)
    {
        CAPTURE(ray_dir);

        std::string str;
        for (ivec2 pos : corner_a <= vector_range <= corner_b)
        {
            str += ".#"[c.MakePointColliderWithRayDir(ray_dir).CollidePoint(pos)];
            if (pos.x == corner_b.x)
                str += '\n';
        }

        // Note that in this string, all teeth must reach their desired X position.
        // A bad algorithm can make them fall short of that.
        REQUIRE_EQ(str, &1[R"(
..................
.###########......
.###########......
.###########......
.###########......
.############.....
.############.....
.############.....
.############.....
.#############....
.#############....
.#############....
.#############....
.##############...
.##############...
.##############...
.##############...
.###############..
.###############..
.############.##..
.#############.#..
.#############..#.
.##############...
.##############...
.###############..
.###############..
.################.
.############.....
.#############....
.##############...
.###############..
.################.
.###############..
.##############...
.#############....
.############.....
.###########......
.###########......
.############.....
.############.....
.#############....
.#############....
.##############...
.##############...
.###############..
.###########.###..
.###########..###.
.###########......
.###########......
.###########......
.###########......
.###########......
.###########......
.###########......
.###########......
.###########......
.###########......
.###########......
.###########......
.###########......
.###########......
.###########......
..................
)"]);
    }
}

// Single edge to edge collisions.
TEST_CASE("edge_soup.edge_to_edge")
{
    using T = EdgeSoup<int>;

    auto sameInAllModes = [](T::Edge e1, T::Edge e2) -> bool
    {
        [[maybe_unused]] bool r1 = e1.CollideWithEdgeInclusive(e2, T::Edge::EdgeCollisionMode::parallelRejected);
        [[maybe_unused]] bool r2 = e1.CollideWithEdgeInclusive(e2, T::Edge::EdgeCollisionMode::parallelUnspecified);
        REQUIRE_EQ(r1, r2);
        return r1;
    };

    // True intersection.
    REQUIRE(sameInAllModes({ivec2(9, 20), ivec2(11,20)}, {ivec2(10,19), ivec2(10,21)}));

    // Point 1 touches.
    REQUIRE(sameInAllModes({ivec2(10, 20), ivec2(12,20)}, {ivec2(10,19), ivec2(10,21)}));
    REQUIRE_FALSE(sameInAllModes({ivec2(11, 20), ivec2(12,20)}, {ivec2(10,19), ivec2(10,21)}));
    // Point 2 touches.
    REQUIRE(sameInAllModes({ivec2(8, 20), ivec2(10,20)}, {ivec2(10,19), ivec2(10,21)}));
    REQUIRE_FALSE(sameInAllModes({ivec2(8, 20), ivec2(9,20)}, {ivec2(10,19), ivec2(10,21)}));
    // Point 3 touches.
    REQUIRE(sameInAllModes({ivec2(9, 20), ivec2(11,20)}, {ivec2(10,18), ivec2(10,20)}));
    REQUIRE_FALSE(sameInAllModes({ivec2(9, 20), ivec2(11,20)}, {ivec2(10,18), ivec2(10,19)}));
    // Point 4 touches.
    REQUIRE(sameInAllModes({ivec2(9, 20), ivec2(11,20)}, {ivec2(10,20), ivec2(10,22)}));
    REQUIRE_FALSE(sameInAllModes({ivec2(9, 20), ivec2(11,20)}, {ivec2(10,21), ivec2(10,22)}));

    // Parallel overlap.
    REQUIRE(T::Edge{ivec2(10, 20), ivec2(14,20)}.CollideWithEdgeInclusive({ivec2(12,20), ivec2(16,20)}, T::Edge::EdgeCollisionMode::parallelUnspecified));
    // Parallel overlap rejected.
    REQUIRE_FALSE(T::Edge{ivec2(10, 20), ivec2(14,20)}.CollideWithEdgeInclusive({ivec2(12,20), ivec2(16,20)}, T::Edge::EdgeCollisionMode::parallelRejected));
    // Parallel side by side.
    REQUIRE_FALSE(sameInAllModes({ivec2(10, 20), ivec2(11,20)}, {ivec2(10,21), ivec2(11,21)}));
    // Parallel, same line but no overlap.
    REQUIRE_FALSE(T::Edge{ivec2(10, 20), ivec2(11,20)}.CollideWithEdgeInclusive({ivec2(12,20), ivec2(13,20)}, T::Edge::EdgeCollisionMode::parallelRejected));
    REQUIRE_FALSE(T::Edge{ivec2(11, 20), ivec2(10,20)}.CollideWithEdgeInclusive({ivec2(12,20), ivec2(13,20)}, T::Edge::EdgeCollisionMode::parallelRejected));
    REQUIRE_FALSE(T::Edge{ivec2(11, 20), ivec2(10,20)}.CollideWithEdgeInclusive({ivec2(13,20), ivec2(12,20)}, T::Edge::EdgeCollisionMode::parallelRejected));
}

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

    auto SameInAllModes = [](T::Edge e1, T::Edge e2) -> bool
    {
        [[maybe_unused]] bool r1 = e1.CollideWithEdgeInclusive(e2, T::Edge::EdgeCollisionMode::parallel_rejected);
        [[maybe_unused]] bool r2 = e1.CollideWithEdgeInclusive(e2, T::Edge::EdgeCollisionMode::parallel_unspecified);
        REQUIRE_EQ(r1, r2);
        return r1;
    };

    // True intersection.
    REQUIRE(SameInAllModes({ivec2(9, 20), ivec2(11,20)}, {ivec2(10,19), ivec2(10,21)}));

    // Point 1 touches.
    REQUIRE(SameInAllModes({ivec2(10, 20), ivec2(12,20)}, {ivec2(10,19), ivec2(10,21)}));
    REQUIRE_FALSE(SameInAllModes({ivec2(11, 20), ivec2(12,20)}, {ivec2(10,19), ivec2(10,21)}));
    // Point 2 touches.
    REQUIRE(SameInAllModes({ivec2(8, 20), ivec2(10,20)}, {ivec2(10,19), ivec2(10,21)}));
    REQUIRE_FALSE(SameInAllModes({ivec2(8, 20), ivec2(9,20)}, {ivec2(10,19), ivec2(10,21)}));
    // Point 3 touches.
    REQUIRE(SameInAllModes({ivec2(9, 20), ivec2(11,20)}, {ivec2(10,18), ivec2(10,20)}));
    REQUIRE_FALSE(SameInAllModes({ivec2(9, 20), ivec2(11,20)}, {ivec2(10,18), ivec2(10,19)}));
    // Point 4 touches.
    REQUIRE(SameInAllModes({ivec2(9, 20), ivec2(11,20)}, {ivec2(10,20), ivec2(10,22)}));
    REQUIRE_FALSE(SameInAllModes({ivec2(9, 20), ivec2(11,20)}, {ivec2(10,21), ivec2(10,22)}));

    // Parallel overlap.
    REQUIRE(T::Edge{ivec2(10, 20), ivec2(14,20)}.CollideWithEdgeInclusive({ivec2(12,20), ivec2(16,20)}, T::Edge::EdgeCollisionMode::parallel_unspecified));
    // Parallel overlap rejected.
    REQUIRE_FALSE(T::Edge{ivec2(10, 20), ivec2(14,20)}.CollideWithEdgeInclusive({ivec2(12,20), ivec2(16,20)}, T::Edge::EdgeCollisionMode::parallel_rejected));
    // Parallel side by side.
    REQUIRE_FALSE(SameInAllModes({ivec2(10, 20), ivec2(11,20)}, {ivec2(10,21), ivec2(11,21)}));
    // Parallel, same line but no overlap.
    REQUIRE_FALSE(T::Edge{ivec2(10, 20), ivec2(11,20)}.CollideWithEdgeInclusive({ivec2(12,20), ivec2(13,20)}, T::Edge::EdgeCollisionMode::parallel_rejected));
    REQUIRE_FALSE(T::Edge{ivec2(11, 20), ivec2(10,20)}.CollideWithEdgeInclusive({ivec2(12,20), ivec2(13,20)}, T::Edge::EdgeCollisionMode::parallel_rejected));
    REQUIRE_FALSE(T::Edge{ivec2(11, 20), ivec2(10,20)}.CollideWithEdgeInclusive({ivec2(13,20), ivec2(12,20)}, T::Edge::EdgeCollisionMode::parallel_rejected));
}

// Collisions between edge soups.
TEST_CASE("edge_soup.soup_to_soup")
{
    using T = EdgeSoup<float>;

    static constexpr auto PointsNear = [](fvec2 a, fvec2 b) -> bool
    {
        return (a - b).len_sq() < 0.0001f;
    };

    struct CollidingPoint
    {
        fvec2 pos;
        std::optional<int> self_edge;
        std::optional<int> other_edge;
    };

    auto Collide = [](const T &shape_a, const T &shape_b, float t, fvec2 pos_a, float angle_a, fvec2 offset_a, float rot_a, fvec2 pos_b, float angle_b, fvec2 offset_b, float rot_b, const std::vector<CollidingPoint> &points)
    {
        auto HalfCollide = [](bool inv, float t, const T &shape_a, const T &shape_b, fvec2 pos_a, float angle_a, fvec2 offset_a, float rot_a, fvec2 pos_b, float angle_b, fvec2 offset_b, float rot_b, std::vector<CollidingPoint> points)
        {
            std::vector<char> memory;
            std::size_t memory_pos = 0;
            T::EdgeSoupCollider::Params params{
                .memory_pool = &memory,
                .memory_pool_offset = &memory_pos,
                .self_pos = pos_a,
                .other_pos = pos_b,
                .self_vel = offset_a,
                .other_vel = offset_b,
                .self_rot = fmat2::rotate(angle_a),
                .other_rot = fmat2::rotate(angle_b),
                .self_angular_vel_abs_upper_bound = abs(rot_a),
                .other_angular_vel_abs_upper_bound = abs(rot_b),
            };
            auto collider = shape_a.MakeEdgeSoupCollider(shape_b, params);
            fvec2 col_pos_a = pos_a + offset_a * t;
            fvec2 col_pos_b = pos_b + offset_b * t;
            fmat2 col_mat_a = fmat2::rotate(angle_a + rot_a * t);
            fmat2 col_mat_b = fmat2::rotate(angle_b + rot_b * t);
            collider.Collide(col_pos_a, col_mat_a, col_pos_b, col_mat_b,
                [&](
                    T::EdgeIndex self_edge_index, const T::Edge &self_edge, const T::Edge &world_self_edge,
                    T::EdgeIndex other_edge_index, const T::Edge &other_edge, const T::Edge &world_other_edge,
                    T::scalar num_self, T::scalar num_other, T::scalar den
                ) -> bool
                {
                    (void)self_edge_index;
                    (void)self_edge;
                    (void)other_edge_index;
                    (void)other_edge;
                    (void)world_other_edge;

                    REQUIRE_MESSAGE(!points.empty(), "More collisions than expected.");

                    if (inv)
                        std::swap(num_self, num_other);

                    fvec2 point = world_self_edge.Point(num_self, den);
                    REQUIRE(PointsNear(point, world_other_edge.Point(num_other, den)));
                    REQUIRE(PointsNear(point, col_mat_a * self_edge.Point(num_self, den) + col_pos_a));
                    REQUIRE(PointsNear(point, col_mat_b * other_edge.Point(num_other, den) + col_pos_b));

                    auto it = std::find_if(points.begin(), points.end(), [&](const CollidingPoint &p)
                    {
                        return PointsNear(p.pos, point);
                    });
                    REQUIRE(it != points.end());

                    if (it->self_edge)
                        REQUIRE(shape_a.GetEdge(self_edge_index).dense_index == *it->self_edge);
                    if (it->other_edge)
                        REQUIRE(shape_a.GetEdge(other_edge_index).dense_index == *it->other_edge);

                    points.erase(it);
                    return false;
                }
            );
            REQUIRE_MESSAGE(points.empty(), "Less collisions than expected.");
        };

        HalfCollide(false, t, shape_a, shape_b, pos_a, angle_a, offset_a, rot_a, pos_b, angle_b, offset_b, rot_b, points);
        HalfCollide(true, t, shape_b, shape_a, pos_b, angle_b, offset_b, rot_b, pos_a, angle_a, offset_a, rot_a, points);
    };

    T shape_a;
    shape_a.AddLoop({fvec2(4,4),fvec2(-4,4),fvec2(-4,-4),fvec2(4,-4)});

    T shape_b;
    shape_b.AddLoop({fvec2(3,0),fvec2(0,3),fvec2(-3,0),fvec2(0,-3)});

    Collide(shape_a, shape_b, 1, fvec2(), 0, fvec2(), 0, fvec2(), 0, fvec2(), 0, {});
    Collide(shape_a, shape_b, 1, fvec2(), 0, fvec2(), 0, fvec2(-5,0), 0, fvec2(), 0, {
        {.pos = fvec2(-4,-2), .self_edge = 1, .other_edge = 3},
        {.pos = fvec2(-4, 2), .self_edge = 1, .other_edge = 0},
    });
}

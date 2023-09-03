#include <iostream>
#include <set>

#include <doctest/doctest.h>

#include "edge_soup.h"

// Tests the edge editing interface.
TEST_CASE("edge_soup.editor")
{
    using T = EdgeSoup<int>;

    static constexpr auto CheckDenseIndices = [](const T &soup)
    {
        for (int i = 0; i < soup.NumEdges(); i++)
            REQUIRE(soup.GetEdge(soup.GetEdgeIndex(i)).dense_index == i);
    };

    static constexpr auto AssertCorrectness = [](const T &soup)
    {
        CheckDenseIndices(soup);

        for (int i = 0; i < soup.NumEdges(); i++)
        {
            CAPTURE(i);
            const auto &edge = soup.GetEdge(soup.GetEdgeIndex(i));
            REQUIRE(edge.prev != T::null_edge);
            REQUIRE(edge.next != T::null_edge);
            [[maybe_unused]] const auto &prev_edge = soup.GetEdge(edge.prev);
            [[maybe_unused]] const auto &next_edge = soup.GetEdge(edge.next);
            REQUIRE_EQ(edge.a, prev_edge.b);
            REQUIRE_EQ(edge.b, next_edge.a);
        }
    };

    T c;
    c.Edit([&]([[maybe_unused]] T::Editor &e)
    {
        REQUIRE((e.IsComplete() && e.NumDetachedPoints() == 0));
    });

    c.AddClosedLoop({{10,10}, {20,10}, {20,20}, {10, 20}});
    AssertCorrectness(c);
    ASSERT(c.NumEdges() == 4);
    ASSERT(c.GetEdge(c.GetEdgeIndex(0)).a == ivec2(10,10));
    ASSERT(c.GetEdge(c.GetEdgeIndex(1)).a == ivec2(20,10));
    ASSERT(c.GetEdge(c.GetEdgeIndex(2)).a == ivec2(20,20));
    ASSERT(c.GetEdge(c.GetEdgeIndex(3)).a == ivec2(10,20));

    c.Edit([&](T::Editor &e)
    {
        [[maybe_unused]] T::EdgeIndex edge{};

        REQUIRE((e.IsComplete() && e.NumDetachedPoints() == 0));
        e.RemoveEdge(c.GetEdgeIndex(1));
        CheckDenseIndices(c);
        REQUIRE((!e.IsComplete() && e.NumDetachedPoints() == 2));

        edge = e.AddEdge({{22,12},{22,18}});
        REQUIRE(c.GetEdge(edge).dense_index == 3);
        REQUIRE((!e.IsComplete() && e.NumDetachedPoints() == 4));

        edge = e.AddEdge({{20,10},{22,12}});
        REQUIRE(c.GetEdge(edge).dense_index == 4);
        REQUIRE((!e.IsComplete() && e.NumDetachedPoints() == 2));

        edge = e.AddEdge({{22,18},{20,20}});
        REQUIRE(c.GetEdge(edge).dense_index == 5);
        REQUIRE((e.IsComplete() && e.NumDetachedPoints() == 0));
    });
}

// Some basic point collision tests.
TEST_CASE("edge_soup.point_to_soup.basic")
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
    c.AddClosedLoop(std::array{
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
TEST_CASE("edge_soup.point_to_soup.pixelperfect.basic")
{
    EdgeSoup<int> c;
    c.AddClosedLoop(std::array{
        ivec2( 6, 0),
        ivec2(12, 6),
        ivec2( 6,12),
        ivec2( 0, 6),
    });
    // Hole:
    c.AddClosedLoop(std::array{
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
TEST_CASE("edge_soup.point_to_soup.pixelperfect.axis_aligned")
{
    EdgeSoup<int> c;
    c.AddClosedLoop(std::array{
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
    c.AddClosedLoop(std::array{
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
    c.AddClosedLoop(std::array{
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
        Storage::MonotonicPool persistent_pool;
        Storage::MonotonicPool temp_pool;

        auto HalfCollide = [&persistent_pool, &temp_pool](bool inv, float t, const T &shape_a, const T &shape_b, fvec2 pos_a, float angle_a, fvec2 offset_a, float rot_a, fvec2 pos_b, float angle_b, fvec2 offset_b, float rot_b, std::vector<CollidingPoint> points)
        {
            CAPTURE(inv);


            T::EdgeSoupCollider::Params params{
                .persistent_pool = &persistent_pool,
                .temp_pool = temp_pool,
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
                [&](T::EdgeSoupCollider::CallbackData data) -> bool
                {
                    (void)data;

                    REQUIRE_MESSAGE(!points.empty(), "More collisions than expected.");

                    fvec2 point = data.world_self_edge.Point(data.num_self, data.den);
                    REQUIRE(PointsNear(point, data.world_other_edge.Point(data.num_other, data.den)));
                    REQUIRE(PointsNear(point, col_mat_a * data.self_edge.Point(data.num_self, data.den) + col_pos_a));
                    REQUIRE(PointsNear(point, col_mat_b * data.other_edge.Point(data.num_other, data.den) + col_pos_b));

                    [[maybe_unused]] bool has_near_points = false;
                    auto it = std::find_if(points.begin(), points.end(), [&](const CollidingPoint &p)
                    {
                        if (PointsNear(p.pos, point))
                            has_near_points = true;
                        else
                            return false;

                        std::optional<int> expected_self_edge = p.self_edge;
                        std::optional<int> expected_other_edge = p.other_edge;
                        if (inv)
                            std::swap(expected_self_edge, expected_other_edge);

                        if (expected_self_edge && shape_a.GetEdge(data.self_edge_index).dense_index != *expected_self_edge)
                            return false;
                        if (expected_other_edge && shape_b.GetEdge(data.other_edge_index).dense_index != *expected_other_edge)
                            return false;
                        return true;
                    });
                    REQUIRE_MESSAGE(has_near_points, "Unexpected collision point position.");
                    REQUIRE_MESSAGE(it != points.end(), "Unexpected collision point metadata.");

                    points.erase(it);
                    return false;
                }
            );
            REQUIRE_MESSAGE(points.empty(), "Less collisions than expected.");
        };

        HalfCollide(false, t, shape_a, shape_b, pos_a, angle_a, offset_a, rot_a, pos_b, angle_b, offset_b, rot_b, points);
        HalfCollide(true, t, shape_b, shape_a, pos_b, angle_b, offset_b, rot_b, pos_a, angle_a, offset_a, rot_a, points);
    };

    // Static objects:

    { // Simple collision.
        T shape_a;
        shape_a.AddClosedLoop({fvec2(4,4),fvec2(-4,4),fvec2(-4,-4),fvec2(4,-4)});

        T shape_b;
        shape_b.AddClosedLoop({fvec2(3,0),fvec2(0,3),fvec2(0,-3)});

        // No collision, we're fully inside.
        Collide(shape_a, shape_b, 1, fvec2(), 0, fvec2(), 0, fvec2(), 0, fvec2(), 0, {});
        // No collision, we're outside.
        Collide(shape_a, shape_b, 1, fvec2(), 0, fvec2(), 0, fvec2(-8,0), 0, fvec2(), 0, {});
        // Penetrating collision.
        Collide(shape_a, shape_b, 1, fvec2(), 0, fvec2(), 0, fvec2(-5,0), 0, fvec2(), 0, {
            {.pos = fvec2(-4,-2), .self_edge = 1, .other_edge = 2},
            {.pos = fvec2(-4, 2), .self_edge = 1, .other_edge = 0},
        });
        // Touching collision.
        Collide(shape_a, shape_b, 1, fvec2(), 0, fvec2(), 0, fvec2(-7,0), 0, fvec2(), 0, {
            {.pos = fvec2(-4, 0), .self_edge = 1, .other_edge = 2},
            {.pos = fvec2(-4, 0), .self_edge = 1, .other_edge = 0},
        });
        // Touching corner collision.
        Collide(shape_a, shape_b, 1, fvec2(), 0, fvec2(), 0, fvec2(-7,4), 0, fvec2(), 0, {
            {.pos = fvec2(-4, 4), .self_edge = 1, .other_edge = 2},
            {.pos = fvec2(-4, 4), .self_edge = 1, .other_edge = 0},
            {.pos = fvec2(-4, 4), .self_edge = 0, .other_edge = 2},
            {.pos = fvec2(-4, 4), .self_edge = 0, .other_edge = 0},
        });
    }

    { // Box to box with parallel lines.
        T shape_a;
        shape_a.AddClosedLoop({fvec2(4,4),fvec2(-4,4),fvec2(-4,-4),fvec2(4,-4)});
        T shape_b;
        shape_b.AddClosedLoop({fvec2(4,4),fvec2(-4,4),fvec2(-4,-4),fvec2(4,-4)});

        // Exact match, the corners collide.
        Collide(shape_a, shape_b, 1, fvec2(), 0, fvec2(), 0, fvec2(), 0, fvec2(), 0, {
            {.pos = fvec2( 4, 4), .self_edge = 3, .other_edge = 0},
            {.pos = fvec2( 4, 4), .self_edge = 0, .other_edge = 3},
            {.pos = fvec2(-4, 4), .self_edge = 0, .other_edge = 1},
            {.pos = fvec2(-4, 4), .self_edge = 1, .other_edge = 0},
            {.pos = fvec2(-4,-4), .self_edge = 1, .other_edge = 2},
            {.pos = fvec2(-4,-4), .self_edge = 2, .other_edge = 1},
            {.pos = fvec2( 4,-4), .self_edge = 2, .other_edge = 3},
            {.pos = fvec2( 4,-4), .self_edge = 3, .other_edge = 2},
        });
        // Single axis offset. The corners of the overlapping area collide.
        Collide(shape_a, shape_b, 1, fvec2(), 0, fvec2(), 0, fvec2(-4,0), 0, fvec2(), 0, {
            {.pos = fvec2( 0, 4), .self_edge = 0, .other_edge = 3},
            {.pos = fvec2( 0,-4), .self_edge = 2, .other_edge = 3},
            {.pos = fvec2(-4, 4), .self_edge = 1, .other_edge = 0},
            {.pos = fvec2(-4,-4), .self_edge = 1, .other_edge = 2},
        });
        // Two axis offset. Two collision points.
        Collide(shape_a, shape_b, 1, fvec2(), 0, fvec2(), 0, fvec2(4,4), 0, fvec2(), 0, {
            {.pos = fvec2( 4, 0), .self_edge = 3, .other_edge = 2},
            {.pos = fvec2( 0, 4), .self_edge = 0, .other_edge = 1},
        });
    }

    { // Movement.
        T shape_a;
        shape_a.AddClosedLoop({fvec2(4,4),fvec2(-4,4),fvec2(-4,-4),fvec2(4,-4)});

        T shape_a_offset;
        shape_a_offset.AddClosedLoop({fvec2(-30+4,4),fvec2(-30-4,4),fvec2(-30-4,-4),fvec2(-30+4,-4)});

        T shape_b;
        shape_b.AddClosedLoop({fvec2(3,0),fvec2(0,3),fvec2(-3,0),fvec2(0,-3)});

        T shape_b_offset;
        shape_b_offset.AddClosedLoop({fvec2(20+3,0),fvec2(20,3),fvec2(20-3,0),fvec2(20,-3)});

        for (int rot_index_a = 0; rot_index_a < 4; rot_index_a++)
        {
            CAPTURE(rot_index_a);
            float angle_a = rot_index_a * f_pi / 2;

            [[maybe_unused]] int
                a0 = (4 - rot_index_a) % 4,
                a1 = (5 - rot_index_a) % 4,
                a2 = (6 - rot_index_a) % 4,
                a3 = (7 - rot_index_a) % 4;

            for (int rot_index_b = 0; rot_index_b < 4; rot_index_b++)
            {
                CAPTURE(rot_index_b);
                float angle_b = rot_index_b * f_pi / 2;

                [[maybe_unused]] int
                    b0 = (4 - rot_index_b) % 4,
                    b1 = (5 - rot_index_b) % 4,
                    b2 = (6 - rot_index_b) % 4,
                    b3 = (7 - rot_index_b) % 4;

                // No movement.
                Collide(shape_a, shape_b, 1, fvec2(), angle_a, fvec2(), 0, fvec2(-5,0), angle_b, fvec2(), 0, {
                    {.pos = fvec2(-4,-2), .self_edge = a1, .other_edge = b3},
                    {.pos = fvec2(-4, 2), .self_edge = a1, .other_edge = b0},
                });
                Collide(shape_a, shape_b, 1, fvec2(), 0, fvec2(), angle_a, fvec2(-5,0), 0, fvec2(), angle_b, {
                    {.pos = fvec2(-4,-2), .self_edge = a1, .other_edge = b3},
                    {.pos = fvec2(-4, 2), .self_edge = a1, .other_edge = b0},
                });

                // Approaching static body.
                Collide(shape_a, shape_b, 1, fvec2(), angle_a, fvec2(), 0, fvec2(-25,0), angle_b, fvec2(20,0), 0, {
                    {.pos = fvec2(-4,-2), .self_edge = a1, .other_edge = b3},
                    {.pos = fvec2(-4, 2), .self_edge = a1, .other_edge = b0},
                });
                Collide(shape_a, shape_b, 1, fvec2(), 0, fvec2(), angle_a, fvec2(-25,0), 0, fvec2(20,0), angle_b, {
                    {.pos = fvec2(-4,-2), .self_edge = a1, .other_edge = b3},
                    {.pos = fvec2(-4, 2), .self_edge = a1, .other_edge = b0},
                });

                // Overshooting static body.
                Collide(shape_a, shape_b, 0.5f, fvec2(), angle_a, fvec2(), 0, fvec2(-25,0), angle_b, fvec2(40,0), 0, {
                    {.pos = fvec2(-4,-2), .self_edge = a1, .other_edge = b3},
                    {.pos = fvec2(-4, 2), .self_edge = a1, .other_edge = b0},
                });
                Collide(shape_a, shape_b, 0.5f, fvec2(), 0, fvec2(), angle_a*2, fvec2(-25,0), 0, fvec2(40,0), angle_b*2, {
                    {.pos = fvec2(-4,-2), .self_edge = a1, .other_edge = b3},
                    {.pos = fvec2(-4, 2), .self_edge = a1, .other_edge = b0},
                });

                // Approaching moving body from the front.
                Collide(shape_a, shape_b, 1, fvec2(10,0), angle_a, fvec2(-10,0), 0, fvec2(-35,0), angle_b, fvec2(30,0), 0, {
                    {.pos = fvec2(-4,-2), .self_edge = a1, .other_edge = b3},
                    {.pos = fvec2(-4, 2), .self_edge = a1, .other_edge = b0},
                });
                Collide(shape_a, shape_b, 0.5f, fvec2(10,0), angle_a, fvec2(-20,0), 0, fvec2(-35,0), angle_b, fvec2(60,0), 0, {
                    {.pos = fvec2(-4,-2), .self_edge = a1, .other_edge = b3},
                    {.pos = fvec2(-4, 2), .self_edge = a1, .other_edge = b0},
                });

                // Approaching moving body from the side.
                Collide(shape_a, shape_b, 1, fvec2(0,-20), angle_a, fvec2(0,20), 0, fvec2(-5,10), angle_b, fvec2(0,-10), 0, {
                    {.pos = fvec2(-4,-2), .self_edge = a1, .other_edge = b3},
                    {.pos = fvec2(-4, 2), .self_edge = a1, .other_edge = b0},
                });
                Collide(shape_a, shape_b, 0.5f, fvec2(0,-20), angle_a, fvec2(0,40), 0, fvec2(-5,10), angle_b, fvec2(0,-20), 0, {
                    {.pos = fvec2(-4,-2), .self_edge = a1, .other_edge = b3},
                    {.pos = fvec2(-4, 2), .self_edge = a1, .other_edge = b0},
                });

                // Approaching moving body diagonally.
                Collide(shape_a, shape_b, 1, fvec2(10,-20), angle_a, fvec2(-10,20), 0, fvec2(-35,10), angle_b, fvec2(30,-10), 0, {
                    {.pos = fvec2(-4,-2), .self_edge = a1, .other_edge = b3},
                    {.pos = fvec2(-4, 2), .self_edge = a1, .other_edge = b0},
                });
                Collide(shape_a, shape_b, 0.5f, fvec2(10,-20), angle_a, fvec2(-20,40), 0, fvec2(-35,10), angle_b, fvec2(60,-20), 0, {
                    {.pos = fvec2(-4,-2), .self_edge = a1, .other_edge = b3},
                    {.pos = fvec2(-4, 2), .self_edge = a1, .other_edge = b0},
                });
            }

            // Offset body, and the other is fixed.

            // No collision.
            Collide(shape_a, shape_b_offset, 1, fvec2(), angle_a, fvec2(), 0, fvec2(), 0, fvec2(), 0, {});

            // Collision without movement.
            Collide(shape_a, shape_b_offset, 1, fvec2(), angle_a, fvec2(), 0, fvec2(-25,0), 0, fvec2(), 0, {
                {.pos = fvec2(-4,-2), .self_edge = a1, .other_edge = 3},
                {.pos = fvec2(-4, 2), .self_edge = a1, .other_edge = 0},
            });

            // Collision with straight movement.
            Collide(shape_a, shape_b_offset, 1, fvec2(), angle_a, fvec2(), 0, fvec2(-45,0), 0, fvec2(20,0), 0, {
                {.pos = fvec2(-4,-2), .self_edge = a1, .other_edge = 3},
                {.pos = fvec2(-4, 2), .self_edge = a1, .other_edge = 0},
            });
            // Collision with sideways movement.
            Collide(shape_a, shape_b_offset, 1, fvec2(), angle_a, fvec2(), 0, fvec2(-25,-20), 0, fvec2(0,20), 0, {
                {.pos = fvec2(-4,-2), .self_edge = a1, .other_edge = 3},
                {.pos = fvec2(-4, 2), .self_edge = a1, .other_edge = 0},
            });
            Collide(shape_a, shape_b_offset, 0.5f, fvec2(), angle_a, fvec2(), 0, fvec2(-25,-20), 0, fvec2(0,40), 0, {
                {.pos = fvec2(-4,-2), .self_edge = a1, .other_edge = 3},
                {.pos = fvec2(-4, 2), .self_edge = a1, .other_edge = 0},
            });
            // Collision with diagonal movement.
            Collide(shape_a, shape_b_offset, 1, fvec2(), angle_a, fvec2(), 0, fvec2(-55,-20), 0, fvec2(30,20), 0, {
                {.pos = fvec2(-4,-2), .self_edge = a1, .other_edge = 3},
                {.pos = fvec2(-4, 2), .self_edge = a1, .other_edge = 0},
            });

            { // Only rotation.
                // Collision at t=0.
                Collide(shape_a, shape_b_offset, 0, fvec2(), angle_a, fvec2(), 0, fvec2(-25,0), 0, fvec2(), f_pi, {
                    {.pos = fvec2(-4,-2), .self_edge = a1, .other_edge = 3},
                    {.pos = fvec2(-4, 2), .self_edge = a1, .other_edge = 0},
                });

                // Collision at t=0.5, quarter circle.
                Collide(shape_a, shape_b_offset, 0.5f, fvec2(), angle_a, fvec2(), 0, fvec2(0,-25), 0, fvec2(), f_pi, {
                    {.pos = fvec2( 2,-4), .self_edge = a2, .other_edge = 3},
                    {.pos = fvec2(-2,-4), .self_edge = a2, .other_edge = 0},
                });
                // Collision at t=0.5, quarter circle, starting from a rotated position.
                Collide(shape_a, shape_b_offset, 0.5f, fvec2(), angle_a, fvec2(), 0, fvec2(-25,0), -f_pi/2, fvec2(), f_pi, {
                    {.pos = fvec2(-4,-2), .self_edge = a1, .other_edge = 3},
                    {.pos = fvec2(-4, 2), .self_edge = a1, .other_edge = 0},
                });

                // Collision at t=1, quarter circle.
                Collide(shape_a, shape_b_offset, 1, fvec2(), angle_a, fvec2(), 0, fvec2(0,-25), 0, fvec2(), f_pi/2, {
                    {.pos = fvec2( 2,-4), .self_edge = a2, .other_edge = 3},
                    {.pos = fvec2(-2,-4), .self_edge = a2, .other_edge = 0},
                });
                // Collision at t=1, quarter circle, starting from a rotated position.
                Collide(shape_a, shape_b_offset, 1, fvec2(), angle_a, fvec2(), 0, fvec2(-25,0), -f_pi/2, fvec2(), f_pi/2, {
                    {.pos = fvec2(-4,-2), .self_edge = a1, .other_edge = 3},
                    {.pos = fvec2(-4, 2), .self_edge = a1, .other_edge = 0},
                });

                // Collision at t=1, half circle.
                Collide(shape_a, shape_b_offset, 1, fvec2(), angle_a, fvec2(), 0, fvec2(25,0), 0, fvec2(), f_pi, {
                    {.pos = fvec2( 4, 2), .self_edge = a3, .other_edge = 3},
                    {.pos = fvec2( 4,-2), .self_edge = a3, .other_edge = 0},
                });
                // Collision at t=1, quarter circle, starting from a rotated position.
                Collide(shape_a, shape_b_offset, 1, fvec2(), angle_a, fvec2(), 0, fvec2(-25,0), -f_pi, fvec2(), f_pi, {
                    {.pos = fvec2(-4,-2), .self_edge = a1, .other_edge = 3},
                    {.pos = fvec2(-4, 2), .self_edge = a1, .other_edge = 0},
                });
            }

            { // Rotation with movement.
                // Sideways slap.
                Collide(shape_a, shape_b_offset, 0.5f, fvec2(), angle_a, fvec2(), 0, fvec2(-25,-20), -f_pi/2, fvec2(0,40), f_pi, {
                    {.pos = fvec2(-4,-2), .self_edge = a1, .other_edge = 3},
                    {.pos = fvec2(-4, 2), .self_edge = a1, .other_edge = 0},
                });
                Collide(shape_a, shape_b_offset, 1, fvec2(), angle_a, fvec2(), 0, fvec2(-25,-20), -f_pi/2, fvec2(0,20), f_pi/2, {
                    {.pos = fvec2(-4,-2), .self_edge = a1, .other_edge = 3},
                    {.pos = fvec2(-4, 2), .self_edge = a1, .other_edge = 0},
                });

                // Forward slap.
                Collide(shape_a, shape_b_offset, 0.5f, fvec2(), angle_a, fvec2(), 0, fvec2(-45,0), -f_pi/2, fvec2(40,0), f_pi, {
                    {.pos = fvec2(-4,-2), .self_edge = a1, .other_edge = 3},
                    {.pos = fvec2(-4, 2), .self_edge = a1, .other_edge = 0},
                });

                // Forward touch, starting a quarter circle back.
                Collide(shape_a, shape_b_offset, 1, fvec2(), angle_a, fvec2(), 0, fvec2(-45,0), -f_pi/2, fvec2(20,0), f_pi/2, {
                    {.pos = fvec2(-4,-2), .self_edge = a1, .other_edge = 3},
                    {.pos = fvec2(-4, 2), .self_edge = a1, .other_edge = 0},
                });
                // Forward touch, starting a half circle back.
                Collide(shape_a, shape_b_offset, 1, fvec2(), angle_a, fvec2(), 0, fvec2(-45,0), -f_pi, fvec2(20,0), f_pi, {
                    {.pos = fvec2(-4,-2), .self_edge = a1, .other_edge = 3},
                    {.pos = fvec2(-4, 2), .self_edge = a1, .other_edge = 0},
                });
            }

            // Offset body, and the other is moving.
            Collide(shape_a, shape_b_offset, 1, fvec2(30,10), angle_a, fvec2(-30,-10), 0, fvec2(-25,-20), -f_pi/2, fvec2(0,20), f_pi/2, {
                {.pos = fvec2(-4,-2), .self_edge = a1, .other_edge = 3},
                {.pos = fvec2(-4, 2), .self_edge = a1, .other_edge = 0},
            });
        }

        // Both bodies have offset.
        Collide(shape_a_offset, shape_b_offset, 1, fvec2(0,40), 0, fvec2(0,-10), f_pi/2, fvec2(-5,-40), 0, fvec2(0,20), f_pi/2, {
            {.pos = fvec2(-4,-2), .self_edge = 0, .other_edge = 2},
            {.pos = fvec2(-4, 2), .self_edge = 0, .other_edge = 3},
        });
    }
}

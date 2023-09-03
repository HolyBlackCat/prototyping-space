#pragma once

#include <algorithm>
#include <iterator>
#include <ranges>

#include <program/errors.h>
#include <utils/aabb_tree.h>
#include <utils/mat.h>

#define IMP_CONTOUR_DEBUG(...) // (std::cout << FMT(__VA_ARGS__) << '\n')

// Represents a 2D contour.
// Internally stores independent edges, but the collisions only work correctly if the contour is watertight.
// The contour is sensitive to the winding direction. The correct direction is CW if the Y points down.
template <Math::scalar T>
class ContourShape
{
  public:
    using scalar = T;
    using vector = vec2<T>;
    using rect = rect2<T>;

    struct Edge
    {
        vector a;
        vector b;

        // The AABB for this edge.
        [[nodiscard]] rect Bounds() const
        {
            rect ret = a.rect_to(b).fix();
            ret.b = next_value(ret.b);
            return ret;
        }
    };

    using AabbTree = AabbTree<vec2<T>, Edge>;
    using EdgeIndex = typename AabbTree::NodeIndex;

    static constexpr EdgeIndex null_index = AabbTree::null_index;

  private:
    AabbTree aabb_tree;

  public:
    constexpr ContourShape() : aabb_tree(typename AabbTree::Params(vector(0))) {}

    [[nodiscard]] bool IsEmpty() const
    {
        return aabb_tree.IsEmpty();
    }

    // Returns the shape bounds. Those may or may not be slightly larger than the actual shape.
    // Returns a default-constructed rect if the shape is empty.
    [[nodiscard]] rect Bounds() const
    {
        return aabb_tree.Bounds();
    }

    // Inserts a new edge.
    // Save the index to be able to remove it later.
    EdgeIndex AddEdge(Edge edge)
    {
        return aabb_tree.AddNode(edge.Bounds(), edge);
    }

    // Adds an edge loop from the vertices.
    template <std::input_iterator I, std::sentinel_for<I> J = I>
    void AddLoop(I begin, J end)
    {
        vector prev;
        vector first_point = *begin;
        bool first = true;
        while (begin != end)
        {
            if (!first)
            {
                AddEdge({.a = prev, .b = *begin});
            }
            first = false;
            prev = *begin;
            ++begin;
        }
        if (!first)
            AddEdge({.a = prev, .b = first_point});
    }
    void AddLoop(std::ranges::range auto &&range)
    {
        AddLoop(std::ranges::begin(range), std::ranges::end(range));
    }

    // Removes an edge by index.
    // Returns false if no such edge.
    bool RemoveEdge(EdgeIndex index)
    {
        return aabb_tree.RemoveNode(index);
    }

    // Performs collision tests.
    class Collider
    {
        const ContourShape *target = nullptr;

        // The ray parameters for point collisions.
        // The ray always points in one of the 4 cardinal directions.
        bool ray_vertical = false; // X = false, Y = true.
        bool ray_negative = false; // false = positive, true = negative.

      public:
        // See `MakeCollider()` below.
        Collider(const ContourShape &new_target, vector reference_point)
            : target(&new_target)
        {
            rect bounds = target->Bounds();

            int distances[4] = {
                bounds.b.x - reference_point.x,
                bounds.b.y - reference_point.y,
                reference_point.x - bounds.a.x,
                reference_point.y - bounds.a.y,
            };

            int ray_dir = std::min_element(std::begin(distances), std::end(distances)) - std::begin(distances);
            ray_vertical = ray_dir % 2;
            ray_negative = ray_dir >= 2;
        }

        // Returns true if the point is inside the shape, or is exactly on the boundary.
        [[nodiscard]] bool CollidePoint(vector point) const
        {
            IMP_CONTOUR_DEBUG("------------------------------- shape to point collision test:");

            rect bounds = target->Bounds();
            if (!bounds.contains(point))
            {
                IMP_CONTOUR_DEBUG("-----> not in the bounding box");
                return false; // We're outside of the bounds, abort early.
            }

            // We emit a ray in one of the 4 cardinal directions.
            // Compute a bounding box of the ray for the AABB tree query.
            rect ray_aabb = point.tiny_rect();
            if (ray_negative)
                ray_aabb.a[ray_vertical] = bounds.a[ray_vertical];
            else
                ray_aabb.b[ray_vertical] = bounds.b[ray_vertical];

            int counter = 0;
            target->aabb_tree.CollideAabb(ray_aabb, [&](AabbTree::NodeIndex edge_index) -> bool
            {
                const Edge &edge = target->aabb_tree.GetNodeUserData(edge_index);

                // If the edge is entirely on one side of the ray line, OR exactly on the ray line, there's no collision.
                int a_side_sign = sign(edge.a[!ray_vertical] - point[!ray_vertical]);
                int b_side_sign = sign(edge.b[!ray_vertical] - point[!ray_vertical]);
                if (a_side_sign == b_side_sign)
                {
                    IMP_CONTOUR_DEBUG("{} doesn't cross the ray line", edge_index);
                    return false;
                }

                // 1 = we're exiting the edge from behind.
                // -1 = we're entering the edge from the front.
                int ray_edge_sign = (edge.b[!ray_vertical] > edge.a[!ray_vertical]) != ray_vertical != ray_negative ? 1 : -1;

                // Run a precise collision test.
                // We increment the counter, with the sign depending on the ray-to-edge direction.
                // If one of the edge points exactly overlaps the ray, this counts as a half-collision.
                if ((edge.b - edge.a) /cross/ (point - edge.a) * ray_edge_sign >= 0)
                {
                    int delta = (bool(a_side_sign) + bool(b_side_sign)) * ray_edge_sign;
                    counter += delta;
                    IMP_CONTOUR_DEBUG("edge {} PASSES precise collision: {:+}", edge_index, delta);
                    return false;
                }

                IMP_CONTOUR_DEBUG("edge {} fails precise collision (sign = {:+})", edge_index, ray_edge_sign);

                return false;
            });

            // 0 is no collision. 1 tangent collision with the edge.
            // 2 is a collision. Everything else is a geometry issue or an altorithm issue.
            // This holds even if the geometry has holes.
            ASSERT(counter == 0 || counter == 1 || counter == 2, "Bad contour geometry, or broken algorithm.");

            IMP_CONTOUR_DEBUG("// counter = {}", counter);
            IMP_CONTOUR_DEBUG("-----> {}", counter == 0 ? "no collision" : counter == 1 ? "TANGENT collision" : counter == 2 ? "COLLISION" : "(error?)");

            return counter != 0;
        }

        // Which way we throw rays. This is a 4-direction.
        [[nodiscard]] int DebugRayDirection() const
        {
            return ray_vertical + ray_negative * 2;
        }
    };

    // Makes a collider object that can be used to test collisions against this shape.
    // The collider should be a short-lived temporary.
    // The `reference_point` should be near the things you're going to test against, this improves performance.
    // If you're testing against several objects in different locations, create a separate collider for each one.
    [[nodiscard]] Collider MakeCollider(vector reference_point) const
    {
        return Collider(*this, reference_point);
    }

    // Get the AABB tree of edges, mostly for debug purposes.
    [[nodiscard]] const AabbTree &EdgeTree() const
    {
        return aabb_tree;
    }
};

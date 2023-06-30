#pragma once

#include <algorithm>
#include <iterator>
#include <new>
#include <ranges>
#include <span>

#include "program/errors.h"
#include "utils/aabb_tree.h"
#include "utils/alignment.h"
#include "utils/mat.h"

#define IMP_CONTOUR_DEBUG(...) // (std::cout << FMT(__VA_ARGS__) << '\n')

// Represents a 2D contour.
// The edges are stored independently* of each other, in an AABB tree.
// The collisions only work correctly if the contour is watertight.
// The contour is sensitive to the winding direction. The correct direction is CW if the Y points down.
// The edge coorinates are inclusive.
// * (hence "soup", compare with a common phrase "triangle soup")
template <Math::scalar T>
class EdgeSoup
{
  public:
    using scalar = T;
    using vector = vec2<T>;
    using matrix = mat2<T>;
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

        enum class EdgeCollisionMode
        {
            parallelUnspecified, // Parallel edges are handled in an unspecified way. (In practice, should be true if they are on the same line.)
            parallelRejected, // Parallel edges are always considered not colliding.
        };

        [[nodiscard]] bool CollideWithEdgeInclusive(Edge other, EdgeCollisionMode mode) const
        {
            // This is essentially `all_of(0 <= Math::point_dir_intersection_factor_two_way(...)[i] <= 1)`, but without division.
            vector delta_a = other.a - a;
            vector delta_self = b - a;
            vector delta_other = other.b - other.a;
            scalar d = delta_self /cross/ delta_other;

            if (mode == EdgeCollisionMode::parallelRejected && d == 0)
                return false;

            scalar abs_d = abs(d);

            scalar u = delta_a /cross/ delta_other;
            scalar v = delta_a /cross/ delta_self;
            return abs(u) <= abs_d && abs(v) <= abs_d && u * d >= 0 && v * d >= 0;
        }
    };

    struct EdgeWithData : Edge
    {
        int dense_index = 0; // 0 <= ... < NumEdges()

        // How far the edge is from the local origin. The maximum of the two points' distances is used.
        // For integral `T`s, an upper bound approximation is used.
        scalar max_distance_to_origin = 0;
    };

    using AabbTree = AabbTree<vec2<T>, EdgeWithData>;
    enum class EdgeIndex : typename AabbTree::NodeIndex {};

    static constexpr EdgeIndex null_index = EdgeIndex(AabbTree::null_index);

    // For integral `T`s, an upper bound approximation is used.
    [[nodiscard]] static scalar VertexDistanceToOrigin(vector point)
    {
        if constexpr (Math::floating_point_scalar<T>)
        {
            return point.len();
        }
        else
        {
            // *looks at butterfly* Is this a distance approximation?
            return point.abs().sum();
        }
    }

    // This is precomputed and stored in edges as `.max_distance_to_origin`.
    [[nodiscard]] static scalar MaxEdgeDistanceToOrigin(Edge edge)
    {
        return max(VertexDistanceToOrigin(edge.a), VertexDistanceToOrigin(edge.b));
    }

    // "Inverts" the rotation matrix `m`.
    // In practice, merely transposes it, for speed.
    [[nodiscard]] static matrix InvertRotationMatrix(matrix m)
    {
        // Good enough for rotation and/or mirror matrices.
        // Also does work for integral `T`s, unlike `.inverse()`.
        return m.transpose();
    }

  private:
    AabbTree aabb_tree;

    // NOTE: The elements here are `EdgeIndex`es, and their indices is what's stored in `EdgeWithData::dense_index`.
    // This is not exactly intuitive.
    SparseSet<int> dense_edge_indices;

  public:
    constexpr EdgeSoup() : aabb_tree(typename AabbTree::Params(vector(0))) {}

    [[nodiscard]] bool IsEmpty() const
    {
        return aabb_tree.IsEmpty();
    }

    // The exact number of edges.
    [[nodiscard]] int NumEdges() const
    {
        return dense_edge_indices.ElemCount();
    }

    // Returns the shape bounds. Those may or may not be slightly larger than the actual shape.
    // Returns a default-constructed rect if the shape is empty.
    [[nodiscard]] rect Bounds() const
    {
        return aabb_tree.Bounds();
    }

    // Returns the edge by index.
    [[nodiscard]] const EdgeWithData &GetEdge(EdgeIndex edge_index) const
    {
        return aabb_tree.GetNodeUserData(typename AabbTree::NodeIndex(edge_index));
    }

    // Returns `EdgeIndex` (which are not consecutive) from a consecutive index (which are `0 <= ... < NumEdges()`).
    // Raises a debug assertion if the index is invalid.
    [[nodiscard]] EdgeIndex GetEdgeIndex(int dense_edge_index) const
    {
        return EdgeIndex(dense_edge_indices.GetElem(dense_edge_index));
    }

    // Inserts a new edge.
    // Save the index to be able to remove it later.
    // NOTE: Those indices are not consecutive, that's why `EdgeIndex` is a `enum class`.
    // For consecutive indices, use `GetEdge(edge_index).dense_index`.
    EdgeIndex AddEdge(Edge edge)
    {
        EdgeWithData edge_with_data;
        edge_with_data.Edge::operator=(edge);
        edge_with_data.max_distance_to_origin = MaxEdgeDistanceToOrigin(edge_with_data);
        dense_edge_indices.Insert(dense_edge_indices.ElemCount()); // Sic.
        auto ret = aabb_tree.AddNode(edge.Bounds(), std::move(edge_with_data));

        // `SparseSet` doesn't auto-increase capacity.
        // Since `AabbTree::AddNode()` already calculates a suitable capacity, we just reuse it here.
        dense_edge_indices.Reserve(aabb_tree.Nodes().Capacity());
        dense_edge_indices.Insert(ret);

        return EdgeIndex(ret);
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
    // Raises a debug assertion if no such edge.
    void RemoveEdge(EdgeIndex edge_index)
    {
        [[maybe_unused]] bool ok = dense_edge_indices.EraseUnordered(AabbTree::NodeIndex(edge_index)); // Sic.
        ASSERT(ok);
        ok = aabb_tree.RemoveNode(typename AabbTree::NodeIndex(edge_index));
        ASSERT(ok);
    }

    // Get the AABB tree of edges, mostly for debug purposes.
    [[nodiscard]] const AabbTree &EdgeTree() const
    {
        return aabb_tree;
    }

    // Performs collision tests on points.
    class PointCollider
    {
        const EdgeSoup *target = nullptr;

        // The ray parameters for point collisions.
        // The ray always points in one of the 4 cardinal directions.
        bool ray_vertical = false; // X = false, Y = true.
        bool ray_negative = false; // false = positive, true = negative.

      public:
        // See (and use) `MakePointCollider()` below.
        PointCollider(const EdgeSoup &new_target, int ray_dir)
            : target(&new_target)
        {
            ray_dir &= 3;
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
            target->aabb_tree.CollideAabb(ray_aabb, [&](typename AabbTree::NodeIndex edge_index_raw) -> bool
            {
                const Edge &edge = target->GetEdge(EdgeIndex(edge_index_raw));

                // If the edge is entirely on one side of the ray line, OR exactly on the ray line, there's no collision.
                int a_side_sign = sign(edge.a[!ray_vertical] - point[!ray_vertical]);
                int b_side_sign = sign(edge.b[!ray_vertical] - point[!ray_vertical]);
                if (a_side_sign == b_side_sign)
                {
                    IMP_CONTOUR_DEBUG("edge {} doesn't cross the ray line", edge_index_raw);
                    return false;
                }

                // true = we're exiting the edge from behind.
                // false = we're entering the edge from the front.
                bool ray_exits_shape = (edge.b[!ray_vertical] > edge.a[!ray_vertical]) != ray_vertical != ray_negative;

                // Run a precise collision test.
                // We increment the counter, with the sign depending on the ray-to-edge direction.
                // If one of the edge points exactly overlaps the ray, this counts as a half-collision.
                if (((edge.b - edge.a) /cross/ (point - edge.a) >= 0) == ray_exits_shape)
                {
                    int delta = (bool(a_side_sign) + bool(b_side_sign)) * (ray_exits_shape ? 1 : -1);
                    counter += delta;
                    IMP_CONTOUR_DEBUG("edge {} PASSES precise collision: {:+}", edge_index_raw, delta);
                    return false;
                }

                IMP_CONTOUR_DEBUG("edge {} fails precise collision (sign = {:+})", edge_index_raw, ray_exits_shape ? 1 : -1);

                return false;
            });

            // 0 = no collision
            // 1 = tangent collision with an outside edge (the point is on an outside edge, and the ray is parallel to it).
            // 2 = collision
            // 3 = convex corner collision (the point is exactly in a convex corner, and the adjacent edge directions point to different sides of the ray)
            // Everything else is a geometry issue or an altorithm issue.
            // This holds even if the geometry has holes.
            ASSERT(counter == 0 || counter == 1 || counter == 2 || counter == 3, "Bad contour geometry, or broken algorithm.");

            IMP_CONTOUR_DEBUG("// counter = {}", counter);
            IMP_CONTOUR_DEBUG("-----> {}",
                counter == 0 ? "no collision" :
                counter == 1 ? "TANGENT collision" :
                counter == 2 ? "COLLISION" :
                counter == 3 ? "CONVEX CORNER collision" :
                "(error?)");

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
    [[nodiscard]] PointCollider MakePointCollider(vector reference_point) const
    {
        rect bounds = Bounds();

        scalar distances[4] = {
            bounds.b.x - reference_point.x,
            bounds.b.y - reference_point.y,
            reference_point.x - bounds.a.x,
            reference_point.y - bounds.a.y,
        };

        int ray_dir = std::min_element(std::begin(distances), std::end(distances)) - std::begin(distances);
        return MakePointColliderWithRayDir(ray_dir);
    }
    // Same, but will cast rays in a custom 4-direction `ray_dir`.
    [[nodiscard]] PointCollider MakePointColliderWithRayDir(int ray_dir) const
    {
        return PointCollider(*this, ray_dir);
    }

    // Collides with an other edge soup.
    // Only edges are tested. If you want to test for a complete overlap, a separate test is needed.
    template <typename AabbOrOtherToSelf, typename SelfToOther, typename EdgeBoundsInOther, typename F>
    bool CollideEdgeSoupLow(
        const EdgeSoup &other,
        // If `other` is in the same coorinate system, nullptr. Otherwise,
        // either a `rect` in self coorinate system that participates in the collision,
        // or `(vec point) -> vec` that maps from `other` coorinates to self;
        AabbOrOtherToSelf &&aabb_or_other_to_self,
        // If `other` is in the same coorinate system, nullptr. Otherwise `(vec point) -> vec` that maps from self coorinates to `other`.
        SelfToOther &&self_to_other,
        // If this is a single collision, nullptr. Otherwise `(EdgeIndex edge_index, const Edge &transformed_edge) -> rect`
        // that returns the AABB of the edge in `other` coordinates. `transformed_edge` is edge at `edge_index` already transformed to those coordinates.
        // Normally it should return `transformed_edge.Bounds()` extended in some direction, to then perform multiple collisions tests.
        EdgeBoundsInOther &&edge_bounds_in_other,
        // This receives the possibly overlapping edge pairs.
        // `(EdgeIndex self_edge_index, const Edge &self_transformed_edge, EdgeIndex other_edge_index) -> bool`,
        // where `{self,other}_edge_index` are the respective edge indices, and `self_transformed_edge` is the self edge transformed to `other` coorinates.
        // If it returns true, everything stops immediately and also returns true.
        F &&receive_edge_pairs) const
    {
        // Get AABB of other, transform with `map_point`,
        // then get AABB of the resulting polygon.
        rect other_bounds;
        if constexpr (std::is_null_pointer_v<AabbOrOtherToSelf>)
        {
            other_bounds = other.Bounds();
        }
        else if constexpr (std::is_same_v<std::remove_cvref_t<AabbOrOtherToSelf>, rect>)
        {
            other_bounds = aabb_or_other_to_self;
        }
        else
        {
            for (bool first = true; vec point : other.Bounds().to_contour())
            {
                if (first)
                {
                    other_bounds = aabb_or_other_to_self(point).tiny_rect();
                    first = false;
                }
                else
                {
                    other_bounds = other_bounds.combine(aabb_or_other_to_self(point));
                }
            }
        }

        return aabb_tree.CollideAabb(other_bounds, [&](typename AabbTree::NodeIndex edge_index_raw)
        {
            const EdgeIndex edge_index = EdgeIndex(edge_index_raw);
            Edge edge = GetEdge(edge_index);
            if constexpr (!std::is_null_pointer_v<SelfToOther>)
            {
                edge.a = self_to_other(edge.a);
                edge.b = self_to_other(edge.b);
            }

            rect edge_bounds;
            if constexpr (std::is_null_pointer_v<EdgeBoundsInOther>)
                edge_bounds = edge.Bounds();
            else
                edge_bounds = edge_bounds_in_other(edge_index, edge);

            return other.EdgeTree().CollideAabb(edge_bounds, [&](typename AabbTree::NodeIndex other_edge_index)
            {
                return receive_edge_pairs(edge_index, edge, EdgeIndex(other_edge_index));
            });
        });
    }

    // A simple one-time collision between two edge soups.
    // Only edges are tested. If you want to test for a complete overlap, a separate test is needed.
    // `self_matrix` can only contain rotation and mirroring.
    [[nodiscard]] bool CollideEdgeSoupSimple(const EdgeSoup &other, vector self_offset, matrix self_matrix) const
    {
        matrix inv_matrix = InvertRotationMatrix(self_matrix.transpose);
        return CollideEdgeSoupLow(other,
            [&](vector v){return inv_matrix * (v - self_offset);},
            [&](vector v){return self_matrix * v + self_offset;},
            nullptr,
            [&](EdgeIndex self_edge_index, const Edge &self_transformed_edge, EdgeIndex other_edge_index)
            {
                (void)self_edge_index;
                return self_transformed_edge.CollideWithEdgeInclusive(other.GetEdge(other_edge_index), Edge::EdgeCollisionMode::parallelRejected);
            }
        );
    }

    // Performs multiple near collision tests on two edge soups.
    class EdgeSoupCollider
    {
      public:
        struct Params
        {
            // We use this to store the temporary edge lists. The old value is ignored.
            // We increase the size of the vector if it's not large enough.
            std::vector<char> *memory_pool = nullptr;
            // Current byte offset in `memory_pool`.
            std::size_t *memory_pool_offset = nullptr;

            const EdgeSoup *self = nullptr;
            const EdgeSoup *other = nullptr;

            // Other position relative to self.
            vector other_delta_pos;
            // Other velocity relative to self.
            vector other_delta_vel;

            // Rotation matrixes. Can include mirroring.
            matrix self_rot;
            matrix other_rot;

            // Angular velocity in radians, absolute. Can be larger than the true value.
            scalar self_angular_vel_abs_upper_bound = 0;
            scalar other_angular_vel_abs_upper_bound = 0;
        };

        Params params;

        template <typename U>
        [[nodiscard]] U *AllocInPool(std::size_t n = 1)
        {
            static_assert(std::is_trivially_copyable_v<U> && alignof(U) <= __STDCPP_DEFAULT_NEW_ALIGNMENT__);
            std::size_t this_offset = Storage::Align<alignof(U)>(*params->memory_pool_offset);
            *params->memory_pool_offset = this_offset + sizeof(U) * n;
            if (*params->memory_pool_offset > params.memory_pool->size())
                params.memory_pool->resize(std::max(*params->memory_pool_offset, params.memory_pool->size() * 2));
            return ::new((void *)(params.memory_pool->data() + this_offset)) U[n];
        }

        struct EdgeEntry
        {
            Edge transformed_edge;
            std::span<const Edge> collision_candidates;
        };

        // One entry per those self edges that can participate in the collision.
        std::span<const EdgeEntry> edge_entries;

      public:
        EdgeSoupCollider(Params params) : params(std::move(params))
        {
            EdgeEntry *edge_entries_ptr = AllocInPool<EdgeEntry>(params.self->NumEdges());
            #error continue
        }


    };


};

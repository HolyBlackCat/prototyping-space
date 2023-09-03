#pragma once

#include <algorithm>
#include <iterator>
#include <new>
#include <span>
#include <utility>

IMP_PLATFORM_IF(clang)(
    IMP_DIAGNOSTICS_PUSH
    IMP_DIAGNOSTICS_IGNORE("-Wdeprecated-builtins")
    IMP_DIAGNOSTICS_IGNORE("-Wdeprecated-declarations")
)
#include <parallel_hashmap/phmap.h>
IMP_PLATFORM_IF(clang)(
    IMP_DIAGNOSTICS_POP
)

#include "program/errors.h"
#include "utils/aabb_tree.h"
#include "utils/alignment.h"
#include "utils/mat.h"
#include "utils/monotonic_pool.h"

#define IMP_EDGESOUP_DEBUG(...) // (std::cout << FMT(__VA_ARGS__) << '\n')

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

        // Returns a point between `a` (t=0) and `b` (t=1).
        [[nodiscard]] vector Point(scalar t) const requires Math::floating_point_scalar<scalar>
        {
            return a + (b - a) * t;
        }
        // Same but the position is represented by a fraction.
        [[nodiscard]] vector Point(scalar num, scalar den) const
        {
            return a + (b - a) * num / den;
        }

        enum class EdgeCollisionMode
        {
            parallel_unspecified, // Parallel edges are handled in an unspecified way. (In practice, should be true if they are on the same line.)
            parallel_rejected, // Parallel edges are always considered not colliding.
        };

        // If specified, `func` is called on collision, with three parameters: `(scalar num_self, scalar num_other, scalar den) -> void`.
        // `num_self / den` is the fractional position on self edge, and `num_other / den` is the fractional position on the other edge.
        template <typename F = std::nullptr_t>
        [[nodiscard]] bool CollideWithEdgeInclusive(Edge other, EdgeCollisionMode mode, F &&func = nullptr) const
        {
            // This is essentially `all_of(0 <= Math::point_dir_intersection_factor_two_way(...)[i] <= 1)`, but without division.
            vector delta_a = other.a - a;
            vector delta_self = b - a;
            vector delta_other = other.b - other.a;
            scalar d = delta_self /cross/ delta_other;

            if (mode == EdgeCollisionMode::parallel_rejected && d == 0)
                return false;

            scalar abs_d = abs(d);

            scalar u = delta_a /cross/ delta_other;
            scalar v = delta_a /cross/ delta_self;
            bool ret = abs(u) <= abs_d && abs(v) <= abs_d && u * d >= 0 && v * d >= 0;
            if constexpr (!std::is_null_pointer_v<F>)
            {
                if (ret)
                    std::forward<F>(func)(u, v, d);
            }
            return ret;
        }
    };

    struct EdgeWithData;

    using AabbTree = AabbTree<vec2<T>, EdgeWithData>;

    // NOTE: Those indices are not consecutive, that's why it's `enum class`.
    // For consecutive indices, use `GetEdge(edge_index).dense_index`.
    enum class EdgeIndex : typename AabbTree::NodeIndex {};

    // This ID is guaranteed to be unused.
    static constexpr EdgeIndex null_edge = EdgeIndex(AabbTree::null_index);

    struct EdgeWithData : Edge
    {
        int dense_index = 0; // 0 <= ... < NumEdges()

        // Previous and next edges in the loop, or `null_edge` if there are somehow none.
        EdgeIndex prev = null_edge;
        EdgeIndex next = null_edge;

        // How far the edge is from the local origin. The maximum of the two points' distances is used.
        // For integral `T`s, an upper bound approximation is used.
        scalar max_distance_to_origin = 0;
    };

    static constexpr EdgeIndex null_index = EdgeIndex(AabbTree::null_index);

    // For integral `T`s, an upper bound approximation is used.
    [[nodiscard]] static scalar VertexDistanceToOrigin(vector point)
    {
        if constexpr (Math::floating_point_scalar<T>)
        {
            return Math::next_value(point.len());
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

    // Returns a mutable reference to the edge by index.
    [[nodiscard]] EdgeWithData &GetEdgeMutable(EdgeIndex edge_index)
    {
        return aabb_tree.GetNodeUserData(typename AabbTree::NodeIndex(edge_index));
    }

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

    // `EdgeIndex` values are less than this number. Use for sparse arrays indexed by edge indices.
    // This probably should be `<= NumEdges() * 2`, but this is an implementation detail.
    [[nodiscard]] int NumSparseEdges() const
    {
        return aabb_tree.Nodes().ElemCount();
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

    // This is used to edit the edges.
    // Can't be created directly, instead call `Edit()` on the target object.
    class Editor
    {
        friend EdgeSoup;
        EdgeSoup *target = nullptr;

        // Maps point to edge that starts/ends at that point, and doesn't have an adjacent edge in that direction.
        phmap::flat_hash_map<vector, EdgeIndex> missing_prev, missing_next;

        constexpr Editor() {}

        void Finish()
        {
            ASSERT(IsComplete(), "The edge soup is not watertight after `Edit()`.");
        }

        // Inserts a single edge.
        // Doesn't update `missing_{prev,next}`.
        EdgeIndex AddEdgeLow(Edge edge, EdgeIndex prev_index = null_edge, EdgeIndex next_index = null_edge)
        {
            EdgeWithData edge_with_data;
            edge_with_data.Edge::operator=(edge);
            edge_with_data.dense_index = target->dense_edge_indices.ElemCount(); // Sic.
            edge_with_data.prev = prev_index;
            edge_with_data.next = next_index;
            edge_with_data.max_distance_to_origin = MaxEdgeDistanceToOrigin(edge_with_data);
            auto ret = target->aabb_tree.AddNode(edge.Bounds(), std::move(edge_with_data));

            // `SparseSet` doesn't auto-increase capacity.
            // Since `AabbTree::AddNode()` already calculates a suitable capacity, we just reuse it here.
            target->dense_edge_indices.Reserve(target->aabb_tree.Nodes().Capacity());
            target->dense_edge_indices.Insert(ret);

            return EdgeIndex(ret);
        }

        // Removes an edge by index, but doesn't touch `prev` and `next` and doesn't update `missing_{prev,next}`.
        // Raises a debug assertion if no such edge.
        void RemoveEdgeLow(EdgeIndex edge_index)
        {
            int dense_index = target->GetEdge(edge_index).dense_index;
            [[maybe_unused]] bool ok = target->dense_edge_indices.EraseUnordered(typename AabbTree::NodeIndex(edge_index)); // Sic.
            ASSERT(ok);

            // Deleting the edge means we need to reassign its dense index to a different edge.
            // That is, unless its index was equal to the number of edges.
            if (dense_index < target->NumEdges())
                target->GetEdgeMutable(target->GetEdgeIndex(dense_index)).dense_index = dense_index;

            ok = target->aabb_tree.RemoveNode(typename AabbTree::NodeIndex(edge_index));
            ASSERT(ok);
        }

        // Fills `prev` or `next` (`End == 0` and `1` respectively) in the specified newly created edge.
        // Updates the `missing_{prev,next}` maps.
        template <bool Detach, bool End>
        void AttachOrDetach(EdgeIndex edge_index)
        {
            // This edge, as specified in `edge_index`.
            EdgeWithData &edge = target->GetEdgeMutable(edge_index);
            // This `End` of the `edge`.
            const vector point = End ? edge.b : edge.a;

            // The reference to the index of the adjacent edge sharing this `point` (in direction `End`), or `null_edge` if none.
            EdgeIndex &other_edge_index = (End ? edge.next : edge.prev);

            // The map to add this edge to.
            auto &self_map = End ? missing_next : missing_prev;
            // THe map to add the other edge to.
            auto &other_map = End ? missing_prev : missing_next;

            if constexpr (Detach)
            {
                if (other_edge_index == null_edge)
                {
                    // No adjacent edge in this direction. Self should be in the map, remove it from there.
                    auto it = self_map.find(point);
                    ASSERT(it != self_map.end(), "Internal error: Expected the edge that's being deleted to be the map, since it doesn't have an adjacent one, but it's not there.");
                    ASSERT(it->second == edge_index, "Internal error: Expected the edge that's being deleted to be the map, since it doesn't have an adjacent one. The map has an entry for this point, but the edge index is wrong.");
                    self_map.erase(it);
                }
                else
                {
                    // Have adjacent edge in this direction. Add it to the map.
                    EdgeWithData &other_edge = target->GetEdgeMutable(other_edge_index);
                    EdgeIndex &self_edge_index_in_other = (End ? other_edge.prev : other_edge.next);
                    ASSERT(self_edge_index_in_other == edge_index, "Internal error: While deleting an edge, the adjacent edge doesn't contain the correct index of the edge being deleted.");
                    self_edge_index_in_other = null_edge;
                    [[maybe_unused]] bool ok = other_map.try_emplace(point, other_edge_index).second;
                    ASSERT(ok, "While deleting an edge, the resulting detached point immediately overlapped with an other detached point. This is weird and not supported.");
                    other_edge_index = null_edge;
                }
            }
            else
            {
                ASSERT(other_edge_index == null_edge, "Internal error: This edge is already attached in this direction.");

                auto it = other_map.find(point);
                if (it != other_map.end())
                {
                    // Have something to attach to. Remove other from the map.
                    other_edge_index = it->second;

                    EdgeWithData &other_edge = target->GetEdgeMutable(other_edge_index);
                    (End ? other_edge.prev : other_edge.next) = edge_index;

                    other_map.erase(it);
                }
                else
                {
                    // Don't have anything to attach to, add self to the map.
                    self_map.try_emplace(point, edge_index);
                }
            }
        }

        // Inserts a sequence of fully initialized edges.
        // If `Loop`, inserts a final edge to close the loop.
        template <bool Loop, typename I, typename J>
        void AddPathLow(I begin, J end)
        {
            if (begin == end)
                return;

            vector first_point = *begin++;
            vector prev_point = first_point;
            vector this_point;

            if (begin == end)
                return;

            bool first = true;
            EdgeIndex first_edge = null_edge;
            EdgeIndex last_edge = null_edge;

            do
            {
                this_point = *begin++;
                // The order of the assignment targets is important here.
                EdgeIndex this_edge = AddEdgeLow({prev_point, this_point}, last_edge, null_edge);
                if (!first)
                    target->GetEdgeMutable(last_edge).next = this_edge;
                last_edge = this_edge;
                prev_point = this_point;

                if (first)
                {
                    first_edge = this_edge;
                    first = false;
                }
            }
            while (begin != end);

            if constexpr (Loop)
            {
                EdgeIndex final_edge = AddEdgeLow({this_point, first_point}, last_edge, first_edge);
                target->GetEdgeMutable(first_edge).prev = final_edge;
                target->GetEdgeMutable(last_edge).next = final_edge;
            }
            else
            {
                AttachOrDetach<false, false>(first_edge);
                AttachOrDetach<false, true>(last_edge);
            }
        }

      public:
        Editor(const Editor &) = delete;
        Editor &operator=(const Editor &) = delete;

        // If true, the edges are currently watertight, and finishing now wouldn't produce any errors.
        [[nodiscard]] bool IsComplete() const
        {
            return missing_prev.empty() && missing_next.empty();
        }

        // The total number of edge ends that are not attached yet.
        // 0 means that we're watertight.
        [[nodiscard]] std::size_t NumDetachedPoints() const
        {
            return missing_prev.size() + missing_next.size();
        }

        EdgeIndex AddEdge(Edge edge)
        {
            EdgeIndex ret = AddEdgeLow(edge);
            AttachOrDetach<false, false>(ret);
            AttachOrDetach<false, true>(ret);
            return ret;
        }

        void RemoveEdge(EdgeIndex edge_index)
        {
            AttachOrDetach<true, false>(edge_index);
            AttachOrDetach<true, true>(edge_index);
            RemoveEdgeLow(edge_index);
        }

        // Add non-closed path.
        // This can handle closed loops too (if the start and end points overlap), but for those prefer `AddLoop`, which is faster.
        template <std::ranges::range R = std::initializer_list<vector>>
        requires std::same_as<std::ranges::range_value_t<R>, vector>
        void AddPath(R &&range)
        {
            AddPathLow<false>(range.begin(), range.end());
        }

        // Add a closed path.
        // The final closing edge is inserted automatically.
        template <std::ranges::range R = std::initializer_list<vector>>
        requires std::same_as<std::ranges::range_value_t<R>, vector>
        void AddClosedLoop(R &&range)
        {
            AddPathLow<true>(range.begin(), range.end());
        }
    };

    // Can be used to edit the shape.
    // When done, triggers an assertion if the mesh is not watertight.
    // `func` is `(Editor &e) -> eoid`.
    template <typename F>
    void Edit(F &&func)
    {
        Editor editor;
        editor.target = this;
        std::forward<F>(func)(editor);
        editor.Finish();
    }

    // Since this can't generate non-watertight edge soups, we can expose this directly outside of `Edit()`.
    template <std::ranges::range R = std::initializer_list<vector>>
    requires std::same_as<std::ranges::range_value_t<R>, vector>
    void AddClosedLoop(R &&range)
    {
        Edit([&](Editor &e){e.AddClosedLoop(std::forward<R>(range));});
    }

    // Get the AABB tree of edges, mostly for debug purposes.
    [[nodiscard]] const AabbTree &EdgeTree() const
    {
        return aabb_tree;
    }

    // Returns the radius enclosing all edges, of a circle centered at the local origin.
    [[nodiscard]] scalar EnclosingRadius() const
    {
        rect bounds = Bounds();
        scalar sq = bounds.corner(0).len_sq();
        for (int i = 1; i < 3; i++)
            clamp_var_min(sq, bounds.corner(i).len_sq());
        return sqrt(sq);
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
            IMP_EDGESOUP_DEBUG("------------------------------- shape to point collision test:");

            rect bounds = target->Bounds();
            if (!bounds.contains(point))
            {
                IMP_EDGESOUP_DEBUG("-----> not in the bounding box");
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
                    IMP_EDGESOUP_DEBUG("edge {} doesn't cross the ray line", edge_index_raw);
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
                    IMP_EDGESOUP_DEBUG("edge {} PASSES precise collision: {:+}", edge_index_raw, delta);
                    return false;
                }

                IMP_EDGESOUP_DEBUG("edge {} fails precise collision (sign = {:+})", edge_index_raw, ray_exits_shape ? 1 : -1);

                return false;
            });

            // 0 = no collision
            // 1 = tangent collision with an outside edge (the point is on an outside edge, and the ray is parallel to it).
            // 2 = collision
            // 3 = convex corner collision (the point is exactly in a convex corner, and the adjacent edge directions point to different sides of the ray)
            // Everything else is a geometry issue or an altorithm issue.
            // This holds even if the geometry has holes.
            ASSERT(counter == 0 || counter == 1 || counter == 2 || counter == 3, "Bad contour geometry, or broken algorithm.");

            IMP_EDGESOUP_DEBUG("// counter = {}", counter);
            IMP_EDGESOUP_DEBUG("-----> {}",
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
    template <typename AabbOrOtherToSelf, typename SelfToOther, typename EdgeBoundsInOther, typename BeginF, typename StepF, typename EndF>
    bool CollideEdgeSoupLow(
        const EdgeSoup &other,
        // If `other` is in the same coorinate system, nullptr. Otherwise,
        // either a `rect` in self coorinate system that participates in the collision (often, other bounds transformed to self coordinates),
        // or `(vector point) -> vector` that maps from `other` coorinates to self.
        AabbOrOtherToSelf &&aabb_or_other_to_self,
        // If `other` is in the same coorinate system, nullptr. Otherwise `(vector point) -> vector` that maps from self coorinates to `other`.
        SelfToOther &&self_to_other,
        // If this is a single collision, nullptr. Otherwise `(EdgeIndex edge_index, const EdgeWithData &edge, const Edge &transformed_edge) -> rect`
        // that returns the AABB of the edge in `other` coordinates. `transformed_edge` is edge at `edge_index` already transformed to those coordinates.
        // `edge` is redundant when we have `edge_index`, but it's still passed for convenience.
        // Normally it should return `transformed_edge.Bounds()` extended in some direction, to then perform multiple collisions tests.
        EdgeBoundsInOther &&edge_bounds_in_other,
        // This is called before processing each potentially colliding edge of self.
        // It's GUARANTEED to be followed by at least one `add_other_edge`.
        // Either nullptr or `(EdgeIndex self_edge_index, const EdgeWithData &edge, const Edge &self_transformed_edge) -> bool`,
        // where `self_transformed_edge` is the self edge transformed to other coordinates.
        // If this returns true, everything stops and the function also returns true.
        BeginF &&begin_self_edge,
        // This receives other edges possibly colliding with the edge previously passed to `begin_self_edge`.
        // This will be called at least once per `begin_self_edge`.
        // `(EdgeIndex self_edge_index, const EdgeWithData &edge, const Edge &self_transformed_edge, EdgeIndex other_edge_index) -> bool`.
        // If this returns true, everything stops and the function also returns true.
        StepF &&add_other_edge,
        // This is the counterpart of `begin_self_edge`.
        // Either nullptr or `(EdgeIndex self_edge_index, const EdgeWithData &edge, const Edge &self_transformed_edge) -> bool`.
        // If this returns true, everything stops and the function also returns true.
        EndF &&end_self_edge
    ) const
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
            for (bool first = true; vector point : other.Bounds().to_contour())
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
            const EdgeWithData &original_edge = GetEdge(edge_index);
            Edge edge = original_edge;
            if constexpr (!std::is_null_pointer_v<SelfToOther>)
            {
                edge.a = self_to_other(edge.a);
                edge.b = self_to_other(edge.b);
            }

            rect edge_bounds;
            if constexpr (std::is_null_pointer_v<EdgeBoundsInOther>)
                edge_bounds = edge.Bounds();
            else
                edge_bounds = edge_bounds_in_other(edge_index, original_edge, edge);

            bool first = true;

            bool stop = other.EdgeTree().CollideAabb(edge_bounds, [&](typename AabbTree::NodeIndex other_edge_index)
            {
                if (first)
                {
                    first = false;
                    if constexpr (!std::is_null_pointer_v<BeginF>)
                    {
                        if (begin_self_edge(edge_index, original_edge, edge))
                            return true;
                    }
                }
                return add_other_edge(edge_index, original_edge, edge, EdgeIndex(other_edge_index));
            });
            if (stop)
                return true;

            if constexpr (!std::is_null_pointer_v<EndF>)
            {
                if (first)
                    IMP_EDGESOUP_DEBUG("edge {} is in the box, but has no candidates", edge_index_raw);

                if (!first && end_self_edge(edge_index, original_edge, edge))
                    return true;
            }

            return false;
        });
    }

    // A simple one-time collision between two edge soups.
    // Only edges are tested. If you want to test for a complete overlap, a separate test is needed.
    // `self_matrix` can only contain rotation and mirroring.
    [[nodiscard]] bool CollideEdgeSoupSimple(const EdgeSoup &other, vector other_offset, matrix other_matrix) const
    {
        IMP_EDGESOUP_DEBUG("------------------------------- simple shape to shape collision:");

        matrix inv_matrix = InvertRotationMatrix(other_matrix);
        return CollideEdgeSoupLow(other,
            [&](vector v){return other_matrix * v + other_offset;},
            [&](vector v){return inv_matrix * (v - other_offset);},
            nullptr,
            nullptr,
            [&](EdgeIndex self_edge_index, const EdgeWithData &self_edge, const Edge &self_transformed_edge, EdgeIndex other_edge_index)
            {
                (void)self_edge_index;
                (void)self_edge;
                bool collides = self_transformed_edge.CollideWithEdgeInclusive(other.GetEdge(other_edge_index), Edge::EdgeCollisionMode::parallel_rejected);
                IMP_EDGESOUP_DEBUG("{} to {} -> {}", typename AabbTree::NodeIndex(self_edge_index), typename AabbTree::NodeIndex(other_edge_index), collides);
                return collides;
            },
            nullptr
        );
    }

    // Performs multiple near collision tests on two edge soups.
    class EdgeSoupCollider
    {
      public:
        struct Params
        {
            // This is used during construction, and then must be kept alive for all subsequent collision tests.
            Storage::MonotonicPool *persistent_pool = nullptr;
            // This is used only during construction, and then cleared.
            Storage::TemporaryPoolRef temp_pool;

            // Positions.
            vector self_pos;
            vector other_pos;
            // Velocities.
            vector self_vel;
            vector other_vel;

            // Rotation matrixes. Can include mirroring.
            matrix self_rot;
            matrix other_rot;

            // Angular velocity in radians, absolute. Can be larger than the true value.
            scalar self_angular_vel_abs_upper_bound = 0;
            scalar other_angular_vel_abs_upper_bound = 0;

            scalar hitbox_expansion_epsilon = []{
                if constexpr (Math::floating_point_scalar<scalar>)
                    return scalar(0.01); // Intended as a value measured in pixels.
                else
                    return 1;
            }();
        };

      private:
        const EdgeSoup *self = nullptr;
        const EdgeSoup *other = nullptr;

        struct CollisionCandidate
        {
            EdgeIndex edge_index;
            Edge edge;
        };

        struct EdgeEntry
        {
            EdgeIndex edge_index;
            Edge edge;
            std::span<CollisionCandidate> collision_candidates;
        };

        // One entry per those self edges that can participate in the collision.
        std::span<EdgeEntry> edge_entries;

      public:
        EdgeSoupCollider(const EdgeSoup &new_self, const EdgeSoup &new_other, const Params &params)
            : self(&new_self), other(&new_other)
        {
            ASSERT(params.self_angular_vel_abs_upper_bound >= 0);
            ASSERT(params.other_angular_vel_abs_upper_bound >= 0);

            edge_entries = params.persistent_pool->template AllocateArray<EdgeEntry>(self->NumEdges());
            // The rest of `edge_entries` is unused in this function, but we shrink it at the end.
            std::size_t num_edge_entries = 0;

            matrix inv_self_rot = InvertRotationMatrix(params.self_rot);
            matrix inv_other_rot = InvertRotationMatrix(params.other_rot);

            auto map_other_point_to_self = [&](vector point){return inv_self_rot * (params.other_rot * point + params.other_pos - params.self_pos);};
            auto map_self_point_to_other = [&](vector point){return inv_other_rot * (params.self_rot * point + params.self_pos - params.other_pos);};

            // Other velocity relative to self, in world coordinates.
            vector other_delta_vel = params.other_vel - params.self_vel;
            // Same, but in other coordinates.
            vector other_delta_vel_in_other = inv_other_rot * other_delta_vel;

            // Would need `next_value` here, if not for `hitbox_expansion_epsilon`.
            scalar max_distance_to_other = sqrt(max((params.other_pos - params.self_pos).len_sq(), ((params.other_pos + other_delta_vel) - (params.self_pos + params.self_vel)).len_sq()));

            // Collision candidates are temporarily stored here.
            auto temp_other_edges = params.temp_pool->template AllocateArray<EdgeIndex>(other->NumEdges());
            std::size_t num_temp_other_edges = 0;

            rect other_aabb_in_self_coords = [&]{
                rect other_bounds = other->Bounds();
                // Compute the outer radius of the other hitbox.
                scalar other_outer_radius = other->EnclosingRadius();
                // Expand according to other angular velocity, but not larger than a box enclosing the enclosing circle.
                // NOTE: This must be the first operation on the box, because we clamp the resulting box.
                other_bounds = vector().centered_rect_halfsize(other_outer_radius).intersect(other_bounds.expand(other_outer_radius * params.other_angular_vel_abs_upper_bound));
                // Expand by the relative velocity.
                other_bounds = other_bounds.expand_dir(other_delta_vel_in_other);
                // Expand according to self angular velocity.
                other_bounds = other_bounds.expand((other_outer_radius + max_distance_to_other) * params.self_angular_vel_abs_upper_bound);

                rect ret;
                for (bool first = true; vector point : other_bounds.to_contour())
                {
                    if (first)
                    {
                        ret = map_other_point_to_self(point).tiny_rect();
                        first = false;
                    }
                    else
                    {
                        ret = ret.combine(map_other_point_to_self(point));
                    }
                }
                return ret.expand(params.hitbox_expansion_epsilon);
            }();

            self->CollideEdgeSoupLow(
                *other,
                other_aabb_in_self_coords,
                map_self_point_to_other,
                [&](EdgeIndex edge_index, const EdgeWithData &edge, const Edge &transformed_edge)
                {
                    (void)edge_index;
                    return transformed_edge.Bounds().expand_dir(-other_delta_vel_in_other).expand(
                        params.hitbox_expansion_epsilon
                        + edge.max_distance_to_origin * params.self_angular_vel_abs_upper_bound
                        + (edge.max_distance_to_origin + max_distance_to_other) * params.other_angular_vel_abs_upper_bound
                    );
                },
                nullptr,
                [&](EdgeIndex self_edge_index, const EdgeWithData &self_edge, const Edge &self_transformed_edge, EdgeIndex other_edge_index)
                {
                    (void)self_edge_index;
                    (void)self_edge;
                    (void)self_transformed_edge;
                    temp_other_edges[num_temp_other_edges++] = other_edge_index;
                    return false;
                },
                [&](EdgeIndex self_edge_index, const EdgeWithData &self_edge, const Edge &self_transformed_edge)
                {
                    (void)self_transformed_edge;
                    EdgeEntry &edge_entry = edge_entries[num_edge_entries];
                    edge_entry.edge_index = self_edge_index;
                    edge_entry.edge = self_edge;
                    edge_entry.collision_candidates = params.persistent_pool->template AllocateArray<CollisionCandidate>(num_temp_other_edges);
                    for (std::size_t i = 0; i < num_temp_other_edges; i++)
                        edge_entry.collision_candidates[i] = {.edge_index = temp_other_edges[i], .edge = other->GetEdge(temp_other_edges[i])};
                    num_edge_entries++;
                    num_temp_other_edges = 0;
                    return false;
                }
            );

            // Shrink `edge_entries` to its proper size.
            edge_entries = {edge_entries.data(), num_edge_entries};
        }

        struct CallbackData
        {
            // Self edge.
            EdgeIndex self_edge_index;
            const Edge &self_edge; // In self coords.
            const Edge &world_self_edge; // In world coords.

            // Other edge.
            EdgeIndex other_edge_index;
            const Edge &other_edge; // In other coords.
            const Edge &world_other_edge; // In world coords.

            // The relative intersection position on the edges, like in `Edge::CollideWithEdgeInclusive()`.
            scalar num_self = 0;
            scalar num_other = 0;
            scalar den = 1;
        };

        // If `func` is not specified, returns true on collision.
        // If `func` is specified, it's called for every edge collision point.
        // If it returns true, the function stops immediately and also returns true.
        // `func` is `(EdgeSoupCollider::CallbackData data) -> bool`.
        template <typename F = std::nullptr_t>
        bool Collide(vector self_pos, matrix self_rot, vector other_pos, matrix other_rot, F &&func = nullptr)
        {
            IMP_EDGESOUP_DEBUG("------------------------------- shape to shape collision in a collider:");

            for (const EdgeEntry &entry : edge_entries)
            {
                IMP_EDGESOUP_DEBUG("edge {}", typename AabbTree::NodeIndex(entry.edge_index));

                Edge world_self_edge = entry.edge;
                world_self_edge.a = self_pos + self_rot * world_self_edge.a;
                world_self_edge.b = self_pos + self_rot * world_self_edge.b;
                IMP_EDGESOUP_DEBUG("{}-{} -> {}-{}", entry.edge.a, entry.edge.b, world_self_edge.a, world_self_edge.b);

                for (const CollisionCandidate &candidate : entry.collision_candidates)
                {
                    IMP_EDGESOUP_DEBUG("  candidate {}", typename AabbTree::NodeIndex(candidate.edge_index));

                    Edge world_other_edge = candidate.edge;
                    world_other_edge.a = other_pos + other_rot * world_other_edge.a;
                    world_other_edge.b = other_pos + other_rot * world_other_edge.b;
                    IMP_EDGESOUP_DEBUG("{}-{} -> {}-{}", candidate.edge.a, candidate.edge.b, world_other_edge.a, world_other_edge.b);
                    if constexpr (std::is_null_pointer_v<F>)
                    {
                        if (world_self_edge.CollideWithEdgeInclusive(world_other_edge, Edge::EdgeCollisionMode::parallel_rejected))
                            return true;
                    }
                    else
                    {
                        bool stop = false;
                        (void)world_self_edge.CollideWithEdgeInclusive(world_other_edge, Edge::EdgeCollisionMode::parallel_rejected, [&](scalar num_self, scalar num_other, scalar den)
                        {
                            if (func(CallbackData{
                                .self_edge_index  = entry.edge_index,
                                .self_edge        = entry.edge,
                                .world_self_edge  = world_self_edge,
                                .other_edge_index = candidate.edge_index,
                                .other_edge       = candidate.edge,
                                .world_other_edge = world_other_edge,
                                .num_self         = num_self,
                                .num_other        = num_other,
                                .den              = den,
                            }))
                            {
                                stop = true;
                            }
                        });
                        if (stop)
                            return true;
                    }
                }
            }

            return false;
        }
    };

    // Used to perform multiple near collision tests on two edge soups.
    [[nodiscard]] EdgeSoupCollider MakeEdgeSoupCollider(const EdgeSoup &other, const EdgeSoupCollider::Params &params) const
    {
        return EdgeSoupCollider(*this, other, params);
    }

    // Accumulates `CallbackData` objects, and produces a collision summary from them.
    class CollisionPointsAccumulator
    {
      public:
        struct Params
        {
            // This needs to stay alive as long as the object is used.
            Storage::MonotonicPool *persistent_pool = nullptr;

            // The starting capacity of point lists in edges.
            // This is used when the first point is inserted, before that the capacity is zero.
            std::size_t edge_points_capacity_initial = 4;
            // When a point list of an edge runs out of capacity (and it's already positive),
            // the capacity is multiplied by this number.
            std::size_t edge_points_capacity_factor = 2;
        };

        struct PointDesc
        {
            // If true, the other edge enters the self shape.
            // If false, the other edge exits the self shape.
            bool in = false;

            // The relative position of the point on the self edge, as a fraction.
            scalar self_num = 0;
            scalar den = 1;
        };

        struct EdgeDesc
        {
            std::span<PointDesc> points;
            // The number of objects we can store at `points.data()`.
            std::size_t points_capacity = 0;
        };

        struct State
        {
            const EdgeSoup *self = nullptr;
            Params params;

            // This is sparse, use `EdgeIndex` as indices.
            std::span<EdgeDesc> edge_entries;
        };
        State state;

      public:
        CollisionPointsAccumulator(const EdgeSoup &self, Params params)
        {
            state.self = self;
            state.params = std::move(params);
            state.edge_entries = state.params.persistent_pool->template AllocateArray<EdgeDesc>(state.self->NumSparseEdges());
        }

        CollisionPointsAccumulator(CollisionPointsAccumulator &&other) noexcept : state(std::exchange(other.state, {})) {}
        CollisionPointsAccumulator &operator=(CollisionPointsAccumulator other) noexcept
        {
            std::swap(state, other.state);
            return *this;
        }

        #error continue (hovering cat.jpt)
    };
};

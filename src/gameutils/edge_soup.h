#pragma once

#include <algorithm>
#include <iterator>
#include <new>
#include <optional>
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
#include "utils/sparse_set.h"

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
    using float_scalar = Math::floating_point_t<T>;
    using vector = vec2<scalar>;
    using float_vector = vec2<float_scalar>;
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

        enum class CollideMode
        {
            normal,
            retry_parallel_prev_self,
            retry_parallel_prev_other,
        };

        // Check if this edge collides with `other`. The edges are half-open, including `a` but excluding `b`.
        // Parallel edges are never considered colliding, even if they are on the same line.
        // If specified, `func` is called on collision, with three parameters: `(scalar num_self, scalar num_other, scalar den) -> void`.
        // `num_self / den` is the fractional position on self edge, and `num_other / den` is the fractional position on the other edge.
        // If `same_line_func` is specified, it's called as `() -> void` if the edges are on the same line.
        // In that case, you're recommended to retry the previous self edge against the same other edge (with `Mode == retry_parallel_prev_self`),
        // and also the same self edge againt the previous other edge (with `Mode == retry_parallel_prev_other`). In those modes, `same_line_func` is ignored.
        template <CollideMode Mode = CollideMode::normal, typename Func = std::nullptr_t, typename SameLineFunc = std::nullptr_t>
        [[nodiscard]] bool CollideWithEdge(Edge other, Func &&func = nullptr, SameLineFunc &&same_line_func = nullptr) const
        {
            if constexpr (Mode == CollideMode::retry_parallel_prev_self)
            {
                Edge self{.a = b, .b = a};
                auto lambda = [&]{
                    if constexpr (std::is_null_pointer_v<Func>)
                    {
                        return nullptr;
                    }
                    else
                    {
                        return [&func](scalar num_self, scalar num_other, scalar den)
                        {
                            std::forward<Func>(func)(den - num_self, num_other, den);
                        };
                    }
                }();
                return self.CollideWithEdge(other, std::move(lambda), std::forward<SameLineFunc>(same_line_func));
            }
            else if constexpr (Mode == CollideMode::retry_parallel_prev_other)
            {
                Edge new_other{.a = other.b, .b = other.a};

                auto lambda = [&]{
                    if constexpr (std::is_null_pointer_v<Func>)
                    {
                        return nullptr;
                    }
                    else
                    {
                        return [&func](scalar num_self, scalar num_other, scalar den)
                        {
                            std::forward<Func>(func)(num_self, den - num_other, den);
                        };
                    }
                }();
                return CollideWithEdge(new_other, std::move(lambda), std::forward<SameLineFunc>(same_line_func));
            }
            else
            {
                // This is essentially `all_of(0 <= Math::point_dir_intersection_factor_two_way(...)[i] <= 1)`, but without division.
                vector delta_a = other.a - a;
                vector delta_self = b - a;
                vector delta_other = other.b - other.a;
                scalar d = delta_self /cross/ delta_other;

                // If you want edges to be inclusive on both sides, you need to do `if (d == 0) return false;` here,
                //   otherwise the parallel edges lying on the same line would always return true.
                // This isn't needed when the second point is exclusive, because of the `< abs_d` test below (which becomes `<= abs_d` if the second point is inclusive).

                if (d == 0)
                {
                    // Handle parallel or overlapping lines.

                    if constexpr (!std::is_null_pointer_v<SameLineFunc> && Mode == CollideMode::normal)
                    {
                        // If the lines overlap, rather than just being parallell...
                        if (delta_a /cross/ delta_other == 0)
                            std::forward<SameLineFunc>(same_line_func)();
                    }

                    return false;
                }

                scalar abs_d = abs(d);

                scalar u = delta_a /cross/ delta_other;
                scalar v = delta_a /cross/ delta_self;

                bool ret = abs(u) < abs_d && abs(v) < abs_d && u * d >= 0 && v * d >= 0;
                if constexpr (!std::is_null_pointer_v<Func>)
                {
                    if (ret)
                    {
                        // Make sure the denominator is not negative, since this is highly unintuitive.
                        if (d < 0)
                        {
                            d = -d;
                            u = -u;
                            v = -v;
                        }

                        std::forward<Func>(func)(u, v, d);
                    }
                }
                return ret;
            }
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

        // The normalized normal.
        float_vector normal;

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

    [[nodiscard]] static std::underlying_type_t<EdgeIndex> EdgeIndexToUnderlying(EdgeIndex edge_index)
    {
        return std::underlying_type_t<EdgeIndex>(edge_index);
    }

    // Returns a mutable reference to the edge by index.
    [[nodiscard]] EdgeWithData &GetEdgeMutable(EdgeIndex edge_index)
    {
        return aabb_tree.GetNodeUserData(EdgeIndexToUnderlying(edge_index));
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
        return aabb_tree.GetNodeUserData(EdgeIndexToUnderlying(edge_index));
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
        // If `precalculated_normal` isn't null, uses it instead of computing the normal.
        EdgeIndex AddEdgeLow(Edge edge, EdgeIndex prev_index, EdgeIndex next_index, std::optional<float_vector> precalculated_normal)
        {
            EdgeWithData edge_with_data;
            edge_with_data.Edge::operator=(edge);
            edge_with_data.dense_index = target->dense_edge_indices.ElemCount(); // Sic.
            edge_with_data.prev = prev_index;
            edge_with_data.next = next_index;
            edge_with_data.normal = precalculated_normal ? *precalculated_normal : (edge.b - edge.a).rot90(-1).norm();
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
            [[maybe_unused]] bool ok = target->dense_edge_indices.EraseUnordered(EdgeIndexToUnderlying(edge_index)); // Sic.
            ASSERT(ok);

            // Deleting the edge means we need to reassign its dense index to a different edge.
            // That is, unless its index was equal to the number of edges.
            if (dense_index < target->NumEdges())
                target->GetEdgeMutable(target->GetEdgeIndex(dense_index)).dense_index = dense_index;

            ok = target->aabb_tree.RemoveNode(EdgeIndexToUnderlying(edge_index));
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
                EdgeIndex this_edge = AddEdgeLow({prev_point, this_point}, last_edge, null_edge, {});
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
                EdgeIndex final_edge = AddEdgeLow({this_point, first_point}, last_edge, first_edge, {});
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

        // Adds a single edge.
        // If `precalculated_normal` isn't null, it's used instead of calcualting the normal. You can use this as an optimization.
        EdgeIndex AddEdge(Edge edge, std::optional<float_vector> precalculated_normal = {})
        {
            EdgeIndex ret = AddEdgeLow(edge, null_edge, null_edge, precalculated_normal);
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

        // Forget unfinished contours, letting you exit the editor without triggering an assertion.
        // Reentering the editor will not let you attach more edges to those contours.
        void DebugForgetUnfinishedContours()
        {
            missing_prev.clear();
            missing_next.clear();
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
            target->aabb_tree.CollideAabb(ray_aabb, [&](AabbTree::NodeIndex edge_index_raw) -> bool
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

        return aabb_tree.CollideAabb(other_bounds, [&](AabbTree::NodeIndex edge_index_raw)
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

            bool stop = other.EdgeTree().CollideAabb(edge_bounds, [&](AabbTree::NodeIndex other_edge_index)
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
        auto other_to_self = [&](vector v){return other_matrix * v + other_offset;};
        auto self_to_other = [&](vector v){return inv_matrix * (v - other_offset);};
        return CollideEdgeSoupLow(other,
            other_to_self,
            self_to_other,
            nullptr,
            nullptr,
            [&](EdgeIndex self_edge_index, const EdgeWithData &self_edge, const Edge &self_transformed_edge, EdgeIndex other_edge_index)
            {
                (void)self_edge_index;
                (void)self_edge;
                bool collides = false;
                const EdgeWithData &other_edge = other.GetEdge(other_edge_index);
                if (self_transformed_edge.CollideWithEdge(other_edge, nullptr, [&]{
                    if (self_edge.prev != null_edge)
                    {
                        Edge prev_self_edge = GetEdge(self_edge.prev);
                        prev_self_edge.a = self_to_other(prev_self_edge.a);
                        prev_self_edge.b = self_to_other(prev_self_edge.b);
                        if (prev_self_edge.template CollideWithEdge<Edge::CollideMode::retry_parallel_prev_self>(other_edge))
                        {
                            collides = true;
                            return;
                        }
                    }
                    if (other_edge.prev != null_edge)
                    {
                        const Edge &prev_other_edge = other.GetEdge(other_edge.prev);
                        if (self_transformed_edge.template CollideWithEdge<Edge::CollideMode::retry_parallel_prev_other>(prev_other_edge))
                        {
                            collides = true;
                            return;
                        }
                    }
                }))
                    collides = true;
                IMP_EDGESOUP_DEBUG("{} to {} -> {}", EdgeIndexToUnderlying(self_edge_index), EdgeIndexToUnderlying(other_edge_index), collides);
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
            const EdgeWithData *edge = nullptr;
        };

        struct EdgeEntry
        {
            EdgeIndex edge_index;
            const EdgeWithData *edge = nullptr;
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
                    edge_entry.edge = &self_edge;
                    edge_entry.collision_candidates = params.persistent_pool->template AllocateArray<CollisionCandidate>(num_temp_other_edges);
                    for (std::size_t i = 0; i < num_temp_other_edges; i++)
                        edge_entry.collision_candidates[i] = {.edge_index = temp_other_edges[i], .edge = &other->GetEdge(temp_other_edges[i])};
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
            const EdgeWithData &self_edge; // In self coords.
            const Edge &world_self_edge; // In world coords.

            // Other edge.
            EdgeIndex other_edge_index;
            const EdgeWithData &other_edge; // In other coords.
            const Edge &world_other_edge; // In world coords.

            // The relative intersection position on the edges, like in `Edge::CollideWithEdge()`.
            scalar num_self = 0;
            scalar num_other = 0;
            scalar den = 1;

            // Returns the same data, but with self and other swapped.
            [[nodiscard]] CallbackData Mirror() const
            {
                return {
                    .self_edge_index  = other_edge_index,
                    .self_edge        = other_edge,
                    .world_self_edge  = world_other_edge,
                    .other_edge_index = self_edge_index,
                    .other_edge       = self_edge,
                    .world_other_edge = world_self_edge,
                    .num_self         = num_other,
                    .num_other        = num_self,
                    .den              = den,
                };
            }
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
                IMP_EDGESOUP_DEBUG("edge {}", EdgeIndexToUnderlying(entry.edge_index));

                Edge world_self_edge = *entry.edge;
                world_self_edge.a = self_pos + self_rot * world_self_edge.a;
                world_self_edge.b = self_pos + self_rot * world_self_edge.b;
                IMP_EDGESOUP_DEBUG("{}-{} -> {}-{}", entry.edge.a, entry.edge.b, world_self_edge.a, world_self_edge.b);

                for (const CollisionCandidate &candidate : entry.collision_candidates)
                {
                    IMP_EDGESOUP_DEBUG("  candidate {}", EdgeIndexToUnderlying(candidate.edge_index));

                    Edge world_other_edge = *candidate.edge;
                    world_other_edge.a = other_pos + other_rot * world_other_edge.a;
                    world_other_edge.b = other_pos + other_rot * world_other_edge.b;
                    IMP_EDGESOUP_DEBUG("{}-{} -> {}-{}", candidate.edge.a, candidate.edge.b, world_other_edge.a, world_other_edge.b);
                    if constexpr (std::is_null_pointer_v<F>)
                    {
                        if (world_self_edge.CollideWithEdge(world_other_edge))
                            return true;
                    }
                    else
                    {
                        bool stop = false;
                        (void)world_self_edge.CollideWithEdge(world_other_edge,
                            [&](scalar num_self, scalar num_other, scalar den)
                            {
                                if (func(CallbackData{
                                    .self_edge_index  = entry.edge_index,
                                    .self_edge        = *entry.edge,
                                    .world_self_edge  = world_self_edge,
                                    .other_edge_index = candidate.edge_index,
                                    .other_edge       = *candidate.edge,
                                    .world_other_edge = world_other_edge,
                                    .num_self         = num_self,
                                    .num_other        = num_other,
                                    .den              = den,
                                }))
                                {
                                    stop = true;
                                }
                            },
                            [&]
                            {
                                // Retry previous self edge.
                                if (const EdgeIndex prev_self_edge_index = entry.edge->prev; prev_self_edge_index != null_edge)
                                {
                                    const EdgeWithData &prev_self_edge = self->GetEdge(prev_self_edge_index);
                                    Edge prev_world_self_edge = prev_self_edge;
                                    prev_world_self_edge.a = self_pos + self_rot * prev_world_self_edge.a;
                                    prev_world_self_edge.b = self_pos + self_rot * prev_world_self_edge.b;

                                    (void)prev_world_self_edge.template CollideWithEdge<Edge::CollideMode::retry_parallel_prev_self>(world_other_edge,
                                        [&](scalar num_self, scalar num_other, scalar den)
                                        {
                                            if (func(CallbackData{
                                                .self_edge_index  = prev_self_edge_index,
                                                .self_edge        = prev_self_edge,
                                                .world_self_edge  = prev_world_self_edge,
                                                .other_edge_index = candidate.edge_index,
                                                .other_edge       = *candidate.edge,
                                                .world_other_edge = world_other_edge,
                                                .num_self         = num_self,
                                                .num_other        = num_other,
                                                .den              = den,
                                            }))
                                            {
                                                stop = true;
                                                return;
                                            }
                                        }
                                    );
                                }

                                // Retry previous other edge.
                                if (const EdgeIndex prev_other_edge_index = candidate.edge->prev; prev_other_edge_index != null_edge)
                                {
                                    const EdgeWithData &prev_other_edge = other->GetEdge(prev_other_edge_index);
                                    Edge prev_world_other_edge = prev_other_edge;
                                    prev_world_other_edge.a = other_pos + other_rot * prev_world_other_edge.a;
                                    prev_world_other_edge.b = other_pos + other_rot * prev_world_other_edge.b;

                                    (void)world_self_edge.template CollideWithEdge<Edge::CollideMode::retry_parallel_prev_other>(prev_world_other_edge,
                                        [&](scalar num_self, scalar num_other, scalar den)
                                        {
                                            if (func(CallbackData{
                                                .self_edge_index  = entry.edge_index,
                                                .self_edge        = *entry.edge,
                                                .world_self_edge  = world_self_edge,
                                                .other_edge_index = prev_other_edge_index,
                                                .other_edge       = prev_other_edge,
                                                .world_other_edge = prev_world_other_edge,
                                                .num_self         = num_self,
                                                .num_other        = num_other,
                                                .den              = den,
                                            }))
                                            {
                                                stop = true;
                                                return;
                                            }
                                        }
                                    );
                                }
                            }
                        );

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

    // Describes collisions between two edge soups.
    struct CollisionData
    {
        struct Segment
        {
            // The world positions of the start and end of the segment.
            vector world_a, world_b;

            // The true normal.
            float_vector world_normal;

            // The edges describing this point.
            EdgeIndex self_edge_index_a = null_edge;
            EdgeIndex other_edge_index_a = null_edge;
            // Those can be equal to their `_a` counterparts above, if the whole collision segment is entirely on one edge.
            // The equality has to be checked separately for self and other.
            EdgeIndex self_edge_index_b = null_edge;
            EdgeIndex other_edge_index_b = null_edge;
        };

        std::span<Segment> segments;

        // This should be zero. Larger numbers indicate an unlucky collision, or broken edge soups, or our algorithm being broken.
        int num_lone_points = 0;
    };

    // Accumulates `EdgeSoupCollider::CallbackData` objects, and produces a `CollisionData` from them.
    class CollisionPointsAccumulator
    {
      public:
        struct Params
        {
            // This needs to stay alive as long as the object is used (until you call `Finalize()`).
            Storage::MonotonicPool *temp_pool = nullptr;

            // Rotation matrices for self and other.
            matrix self_rot;
            matrix other_rot;

            // The starting capacity of point lists in edges.
            // This is used when the first point is inserted, before that the capacity is zero.
            std::size_t edge_points_capacity_initial = 4;
            // When a point list of an  edge runs out of capacity (and it's already positive),
            // the capacity is multiplied by this number.
            std::size_t edge_points_capacity_factor = 2;
        };

        struct PointDesc
        {
            // If true, the other edge enters the self shape, while the self edge exits the other shape.
            // If false, the other edge exits the self shape, while the self edge enters the other shape.
            bool other_enters_self = false;

            // The other edge. The self edge is implied by where the point is stored.
            EdgeIndex other_edge_index = null_edge;
            // The normal of the other edge.
            float_vector other_normal;

            // The fractional position on the self edge.
            scalar num_self = 0;
            scalar den = 1;

            // The world positition of this point.
            vector world_pos;
        };

        struct EdgeEntry
        {
            // Those are sorted by the position on the edge.
            std::span<PointDesc> points;
            // The number of objects we can store at `points.data()`.
            std::size_t points_capacity = 0;

            // The normalized normal, in the world coordinates.
            float_vector normal;
        };

      private:
        struct State
        {
            const EdgeSoup *self = nullptr;
            const EdgeSoup *other = nullptr;
            Params params;

            // Information about each edge.
            // This is sparse, use self `EdgeIndex` as indices.
            std::span<EdgeEntry> edge_entries;
            // This lists all edges with at least one point.
            std::span<EdgeIndex> edges_with_points;

            // The number of collision points.
            int total_num_points = 0;
        };
        State state;

      public:
        CollisionPointsAccumulator() {}

        CollisionPointsAccumulator(const EdgeSoup &self, const EdgeSoup &other, Params params)
        {
            state.self = &self;
            state.other = &other;
            state.params = std::move(params);
            state.edge_entries = state.params.temp_pool->template AllocateArray<EdgeEntry>(state.self->NumSparseEdges());
            state.edges_with_points = {state.params.temp_pool->template AllocateArray<EdgeIndex>(state.self->NumEdges()).data(), 0};
        }

        CollisionPointsAccumulator(CollisionPointsAccumulator &&other) noexcept : state(std::exchange(other.state, {})) {}
        CollisionPointsAccumulator &operator=(CollisionPointsAccumulator other) noexcept
        {
            std::swap(state, other.state);
            return *this;
        }

        // Returns all edges that have at least one point, out of order.
        [[nodiscard]] std::span<const EdgeIndex> GetEdgesWithPoints() const
        {
            return state.edges_with_points;
        }

        // Returns the information about the specified edge.
        [[nodiscard]] const EdgeEntry &GetEdgeEntry(EdgeIndex index) const
        {
            return state.edge_entries[EdgeIndexToUnderlying(index)];
        }

        // Note, the collider must use the same object as 'self' as the one passed to the constructor.
        void AddPoint(const EdgeSoupCollider::CallbackData &data)
        {
            state.total_num_points++;

            EdgeEntry &edge_entry = state.edge_entries[EdgeIndexToUnderlying(data.self_edge_index)];

            // Check that we have enough capacity in the point list.
            ASSERT(edge_entry.points.size() <= edge_entry.points_capacity);
            if (edge_entry.points.size() == edge_entry.points_capacity)
            {
                // Out of capacity.
                if (edge_entry.points_capacity == 0)
                {
                    edge_entry.points_capacity = state.params.edge_points_capacity_initial;
                    // If this is the first time we see this edge, add it to the list.
                    // This shouldn't overflow.
                    state.edges_with_points.data()[state.edges_with_points.size()] = data.self_edge_index;
                    state.edges_with_points = {state.edges_with_points.data(), state.edges_with_points.size() + 1};

                    edge_entry.normal = state.params.self_rot * data.self_edge.normal;
                }
                else
                {
                    edge_entry.points_capacity *= state.params.edge_points_capacity_factor;
                }

                PointDesc *new_points = state.params.temp_pool->template AllocateArray<PointDesc>(edge_entry.points_capacity).data();
                std::move(edge_entry.points.begin(), edge_entry.points.end(), new_points);
                edge_entry.points = {new_points, edge_entry.points.size()/*sic*/};
            }
            ASSERT(edge_entry.points.size() < edge_entry.points_capacity);

            PointDesc &new_point = edge_entry.points.data()[edge_entry.points.size()];
            edge_entry.points = {edge_entry.points.data(), edge_entry.points.size() + 1};

            new_point.num_self = data.num_self;
            new_point.den = data.den;
            new_point.world_pos = data.world_self_edge.Point(data.num_self, data.den);
            new_point.other_enters_self = (data.world_self_edge.b - data.world_self_edge.a) /cross/ (data.world_other_edge.b - data.world_other_edge.a) > 0;
            new_point.other_edge_index = data.other_edge_index;
            new_point.other_normal = data.other_edge.normal;
        }

        [[nodiscard]] CollisionData Finalize(Storage::MonotonicPool &persistent_pool)
        {
            CollisionData ret;
            const int max_num_collisions = state.total_num_points / 2;
            ret.segments = {persistent_pool.AllocateArray<typename CollisionData::Segment>(max_num_collisions).data(), 0};

            // Those contain `EdgeIndex`es of edges with lone points at the end or at the beginning respectively.
            SparseSetNonOwning<int> lone_begin_points(persistent_pool.AllocateArray<int>(state.self->NumSparseEdges() * 2)); // Note, need x2 storage for a sparse set.
            SparseSetNonOwning<int> lone_end_points  (persistent_pool.AllocateArray<int>(state.self->NumSparseEdges() * 2));

            // Creates a new collision and returns the index.
            // If `precalculated_normal` is not null, it's used as the normal (as an optimization). Otherwise it's calculated from the points.
            auto AddNewCollision = [&](
                EdgeIndex self_edge_index_a, EdgeIndex self_edge_index_b,
                const PointDesc &a, const PointDesc &b,
                const float_vector *precalculated_normal
            ) -> int
            {
                int index = int(ret.segments.size());
                ASSERT(index < max_num_collisions);
                ret.segments = {ret.segments.data(), ret.segments.size() + 1};

                typename CollisionData::Segment &seg = ret.segments.back();

                seg.world_a = a.world_pos;
                seg.world_b = b.world_pos;
                seg.world_normal = precalculated_normal ? *precalculated_normal : (seg.world_b - seg.world_a).rot90(-1).norm();

                seg.self_edge_index_a = self_edge_index_a;
                seg.self_edge_index_b = self_edge_index_b;
                seg.other_edge_index_a = a.other_edge_index;
                seg.other_edge_index_b = b.other_edge_index;

                return index;
            };

            // Process individual edges.
            for (EdgeIndex edge_index : state.edges_with_points)
            {
                EdgeEntry &edge_entry = state.edge_entries[EdgeIndexToUnderlying(edge_index)];

                if (!edge_entry.points.empty())
                {
                    // Sort the points by the position on the edge.
                    // Resolve the ties so that `.other_enters_self == true` points go first. (Not sure if strictly necessary, but doesn't hurt.)
                    std::sort(edge_entry.points.begin(), edge_entry.points.end(), [](const PointDesc &a, const PointDesc &b)
                    {
                        if (auto d = a.num_self * b.den - b.num_self * a.den)
                            return d < 0;
                        if (auto d = a.other_enters_self - b.other_enters_self)
                            return d > 0;
                        return false;
                    });

                    // The valid point indices on this edge. We remove the first and/or last indices from this if they form incomplete collisions.
                    std::size_t begin = 0;
                    std::size_t end = edge_entry.points.size();

                    // Handle unfinished collisions at the beginning and the end of the edge.
                    if (edge_entry.points.front().other_enters_self)
                    {
                        begin++;
                        lone_end_points.Insert(EdgeIndexToUnderlying(edge_index));
                    }
                    if (!edge_entry.points.back().other_enters_self)
                    {
                        end--;
                        lone_begin_points.Insert(EdgeIndexToUnderlying(edge_index));
                    }

                    // Handle collisions that are entirely on this edge.
                    for (std::size_t i = begin; i + 1 < end; i++)
                    {
                        if (!edge_entry.points[i].other_enters_self && edge_entry.points[i+1].other_enters_self)
                        {
                            // Create the collision.
                            AddNewCollision(edge_index, edge_index/*again*/, edge_entry.points[i], edge_entry.points[i+1], &edge_entry.normal);
                            i++;
                        }
                    }
                }
            }

            // Handle collisions that span multiple edges.
            // Sic, note `&&`. If exactly one of the sets is empty, we can't process points from the other one here anyway. This shouldn't happen often.
            while (lone_begin_points.ElemCount() > 0 && lone_end_points.ElemCount() > 0)
            {
                // Decide which set to start with.
                // We try to consume the smaller set first, to break out of the outer loop faster.
                bool search_back = lone_end_points.ElemCount() < lone_begin_points.ElemCount();

                auto &target_set = search_back ? lone_end_points : lone_begin_points;
                auto &other_set = search_back ? lone_begin_points : lone_end_points;

                // Grab one of the edges.
                const EdgeIndex target_edge_index = EdgeIndex(target_set.GetElem(target_set.ElemCount() - 1));
                target_set.EraseUnordered(EdgeIndexToUnderlying(target_edge_index));

                const EdgeEntry &target_edge_entry = state.edge_entries[EdgeIndexToUnderlying(target_edge_index)];

                EdgeIndex cur_edge_index = target_edge_index;

                // Walk the edge loop.
                while (true)
                {
                    const EdgeIndex next_edge_index = state.self->GetEdge(cur_edge_index).*(search_back ? &EdgeWithData::prev : &EdgeWithData::next);
                    if (next_edge_index == null_edge)
                    {
                        // Hmm.
                        ret.num_lone_points++;
                        break;
                    }

                    const EdgeEntry &next_edge_entry = state.edge_entries[EdgeIndexToUnderlying(next_edge_index)];

                    if (next_edge_entry.points.empty())
                    {
                        // No points on this edge, continue to the next one.
                        cur_edge_index = next_edge_index;
                        continue;
                    }

                    if (!other_set.Contains(EdgeIndexToUnderlying(next_edge_index)))
                    {
                        // Found an edge with points, but no matching unfinished collision here.
                        ret.num_lone_points++;
                        break;
                    }

                    // Found a matching point.
                    other_set.EraseUnordered(EdgeIndexToUnderlying(next_edge_index));

                    // Get the two points.
                    const PointDesc &point_a = search_back ? next_edge_entry.points.back() : target_edge_entry.points.back();
                    const PointDesc &point_b = search_back ? target_edge_entry.points.front() : next_edge_entry.points.front();

                    if (point_a.world_pos == point_b.world_pos)
                    {
                        // The points overlap exactly, so we can't compute the normal.
                        ret.num_lone_points += 2; // We already deleted two points, hence `2`.
                        break;
                    }

                    // Success.

                    EdgeIndex self_edge_index_a = target_edge_index;
                    EdgeIndex self_edge_index_b = next_edge_index;
                    if (search_back)
                        std::swap(self_edge_index_a, self_edge_index_b);

                    AddNewCollision(self_edge_index_a, self_edge_index_b, point_a, point_b, nullptr);
                    break;
                }
            }
            ret.num_lone_points += lone_begin_points.ElemCount() + lone_end_points.ElemCount();

            return ret;
        }

        // Find the edge normal closest to `world_normal` and return it.
        // This ignores the normal directions, and considers them as lines.
        [[nodiscard]] float_vector GetClosestEdgeNormal(const CollisionData::Segment &seg) const
        {
            float_vector ret;
            float_scalar max_value = 0;

            EdgeIndex edge_indices[4] = {
                // Don't reorder, we assume the first two to be 'self' below.
                seg.self_edge_index_a,
                seg.self_edge_index_b,
                seg.other_edge_index_a,
                seg.other_edge_index_b,
            };

            bool first = true;

            for (int i = 0; i < 4; i++)
            {
                EdgeIndex edge_index = edge_indices[i];
                bool is_self = i < 2;
                float_vector normal = (is_self ? state.params.self_rot : state.params.other_rot) * state.self->GetEdge(edge_index).normal;
                float_scalar value = abs(normal /dot/ seg.world_normal);
                if (first || value > max_value)
                {
                    ret = normal;
                    max_value = value;
                }

                first = false;
            }

            return ret;
        }
    };
};

//
// Created by user on 31/07/22.
//

#ifndef LIBALGEBRA_LITE_FREE_TENSOR_H
#define LIBALGEBRA_LITE_FREE_TENSOR_H

#include <libalgebra_lite/implementation_types.h>
#include "tensor_basis.h"
#include "basis_traits.h"
#include <libalgebra_lite/coefficients.h>
#include "vector_traits.h"
#include "dense_vector.h"
#include <libalgebra_lite/algebra.h>

#include <algorithm>
#include <utility>
#include <vector>

namespace lal {

template<typename,
         template<typename, typename, typename...> class,
         template <typename> class,
         typename...>
class free_tensor;


#define LAL_IS_TENSOR(V) std::is_same<typename V::basis_type, \
    tensor_basis>::value

#define LAL_SAME_COEFFS(V1, V2) std::is_same<typename V1::coefficient_ring, \
    typename V2::coefficient_ring>::value

#define LAL_TENSOR_COMPAT_RVV(R, V1, V2) \
    std::enable_if_t<LAL_IS_TENSOR(R) && LAL_IS_TENSOR(V1) && LAL_IS_TENSOR(V2) \
                      && LAL_SAME_COEFFS(R, V1) && LAL_SAME_COEFFS(R, V2)>

template<typename Coefficients>
class dense_multiplication_helper
{
    using traits = coefficient_trait<Coefficients>;
    using scalar_type = typename traits::scalar_type;

    std::vector<scalar_type> left_read_buffer;
    std::vector<scalar_type> right_read_buffer;
    std::vector<scalar_type> write_buffer;
    std::vector<scalar_type> reverse_buffer;
    const tensor_basis* p_basis;
    const scalar_type* left_ptr;
    const scalar_type* right_ptr;

    scalar_type* out_ptr;
    deg_t lhs_deg;
    deg_t rhs_deg;
    using key_type = typename tensor_basis::key_type;
public:

    dimn_t tile_width;
    dimn_t tile_size;
    deg_t tile_letters;

    template <typename Alloc>
    dense_multiplication_helper(
            free_tensor<Coefficients, dense_vector_view, Alloc>& out,
            const free_tensor<Coefficients, dense_vector_view, Alloc>& lhs,
            const free_tensor<Coefficients, dense_vector_view, Alloc>& rhs
            )
        : p_basis(&out.basis()), lhs_deg(lhs.degree()), rhs_deg(rhs.degree())
    {
        const auto& powers = p_basis->powers();

        // TODO: replace with logic
        tile_letters = 1;

        left_ptr = lhs.as_ptr();
        right_ptr = rhs.as_ptr();
        out_ptr = out.as_mut_ptr();

        tile_width = powers[tile_letters];
        tile_size = tile_width*tile_width;

        left_read_buffer.resize(tile_width);
        right_read_buffer.resize(tile_width);
        write_buffer.resize(tile_size);
    }

    const scalar_type& left_unit() const noexcept
    { return *left_ptr; }
    const scalar_type& right_unit() const noexcept
    { return *right_ptr; }

    const scalar_type* left_tile() const noexcept
    { return left_read_buffer.data(); }
    const scalar_type* right_tile() const noexcept
    { return right_read_buffer.data(); }
    scalar_type* write_tile() noexcept
    { return write_buffer.data(); }

    deg_t lhs_degree() const noexcept {  return lhs_deg; }
    deg_t rhs_degree() const noexcept {  return rhs_deg; }

    const scalar_type* left_fwd_read(key_type k) const noexcept
    {
        auto offset = p_basis->start_of_degree(static_cast<deg_t>(k.degree()));
        return left_ptr + k.index()*tile_width + offset;
    }
    const scalar_type* right_fwd_read(key_type k) const noexcept
    {
        auto offset = p_basis->start_of_degree(static_cast<deg_t>(k.degree()));
        return right_ptr + k.index()*tile_width + offset;
    }
    scalar_type* fwd_write(key_type k) const noexcept
    {
        auto offset = p_basis->start_of_degree(static_cast<deg_t>(k.degree()));
        return out_ptr + k.index()*tile_width + offset;
    }

    void read_left_tile(key_type k) noexcept
    {
        auto offset = p_basis->start_of_degree(static_cast<deg_t>(k.degree()));
        const auto* reverse_ptr = reverse_buffer.data() + k.index()*tile_width + offset;
        std::copy(reverse_ptr, reverse_ptr + tile_width, left_read_buffer.data());
    }
    void read_right_tile(key_type k) noexcept
    {
        auto offset = p_basis->start_of_degree(static_cast<deg_t>(k.degree()));
        const auto* fwd_ptr = right_ptr + k.index()*tile_width + offset;
        std::copy(fwd_ptr, fwd_ptr + tile_width, right_read_buffer.data());
    }
    void write_tile_in(key_type k, key_type kr) noexcept
    {
        const auto offset = p_basis->start_of_degree(k.degree() + 2*tile_letters);
        const auto* in_ptr = out_ptr + k.index()*tile_width + offset;
        auto* tile_ptr = write_tile();
        const auto stride = p_basis->powers()[k.degree()+tile_letters];

        for (dimn_t i=0; i<tile_width; ++i) {
            for (dimn_t j=0; j<tile_width; ++j) {
                tile_ptr[i*tile_width+j] = in_ptr[i*stride+j];
            }
        }
    }

    void write_tile_out(key_type k, key_type kr) noexcept
    {
        const auto deg = k.degree();
        const auto offset = p_basis->start_of_degree(deg + 2*tile_letters);
        const auto stride = p_basis->powers()[deg + tile_letters];

        auto* ptr = out_ptr + k.index()*tile_width + offset;
        auto* tile_ptr = write_tile();

        for (dimn_t i=0; i<tile_width; ++i) {
            for (dimn_t j=0; j<tile_width; ++j) {
                ptr[i*stride + j] = tile_ptr[i*tile_width + j];
            }
        }

        if (deg < p_basis->depth()) {
            // Write reverse data
        }


    }

    key_type reverse(key_type k) const noexcept
    {
        const auto width = p_basis->width();
        auto idx = k.index();

        typename key_type ::index_type result_idx = 0;
        while (idx) {
            result_idx *= width;
            result_idx += idx % width;
            idx /= width;
        }

        return key_type {k.degree(), result_idx};
    }
    std::pair<key_type, key_type> split_key(key_type k, deg_t lhs_size) const noexcept
    {
        auto rhs_size = k.degree() - lhs_size;
        auto split = p_basis->powers()[rhs_size];
        return { key_type(lhs_size, k.index() / split),
                 key_type(rhs_size, k.index() % split)};
    }

    dimn_t stride(deg_t deg) const noexcept {
        return p_basis->powers()[deg - tile_letters];
    }

    dimn_t combine(dimn_t lhs, dimn_t rhs, deg_t rh_deg)
    {
        const auto shift = p_basis->powers()[rh_deg];
        return lhs*shift + rhs;
    }
    key_type combine(key_type lhs, key_type rhs)
    {
        const auto rhs_deg = rhs.degree();
        const auto shift = p_basis->powers()[rhs_deg];
        return key_type{lhs.degree() + rhs_deg, lhs.index()*shift+rhs.index()};
    }
    dimn_t combine(dimn_t lhs, key_type rhs)
    {
        const auto rhs_deg = rhs.degree();
        const auto shift = p_basis->powers()[rhs_deg];
        return lhs*shift + rhs.index();
    }


    std::pair<dimn_t, dimn_t> range_size(deg_t lhs, deg_t rhs) const noexcept
    {
        const auto& powers = p_basis->powers();
        return {powers[lhs], powers[rhs]};
    }

    dimn_t range_size(deg_t deg) const noexcept
    { return p_basis->powers()[deg]; }


};


class free_tensor_multiplication
{

    template <typename Coefficients, typename Fn>
    void fma_dense_traditional(
            dense_multiplication_helper<Coefficients>& helper,
            Fn fn, deg_t out_degree
            ) noexcept
    {
        using key_type = tensor_basis::key_type;
        auto lhs_deg = helper.lhs_degree();
        auto rhs_deg = helper.rhs_degree();

        for (deg_t out_deg = out_degree; out_deg>=0; --out_deg) {
            auto lhs_deg_min = std::max(0, out_deg-rhs_deg);
            auto lhs_deg_max = std::min(out_deg, lhs_deg);

            auto* out_ptr = helper.fwd_write(key_type(out_deg, 0));

            for (deg_t lh_deg = lhs_deg_max; lh_deg>=0; --lh_deg) {
                auto rh_deg = out_deg-lh_deg;

                auto lhs_ptr = helper.left_fwd_read(key_type(lh_deg, 0));
                auto rhs_ptr = helper.right_fwd_read(key_type(rh_deg, 0));

                auto range_sizes = helper.range_size(lh_deg, rh_deg);

                auto* p = out_ptr;
                for (dimn_t i = 0; i<range_sizes.first; ++i) {
                    for (dimn_t j = 0; j<range_sizes.second; ++j) {
                        *(p++) += fn(lhs_ptr[i]*rhs_ptr[j]);
                    }
                }
            }
        }
    }


    template <typename Coefficients, typename Fn>
    void fma_dense_giles(dense_multiplication_helper<Coefficients>& helper,
            Fn fn, deg_t out_degree) noexcept
    {
        using key_type = tensor_basis::key_type;

        auto lhs_deg = helper.lhs_degree();
        auto rhs_deg = helper.rhs_degree();

        auto* tile = helper.write_tile();
        const auto* left_rtile = helper.left_tile();
        const auto* right_rtile = helper.right_tile();

        for (deg_t out_deg=out_degree; out_deg > 2*helper.tile_letters; --out_deg) {
            const auto stride = helper.stride(out_deg);
            const auto adj_deg = out_deg - 2*helper.tile_letters;

            // end is not actually a valid key, but it serves as a marker.
            key_type start{adj_deg, 0}, end{adj_deg, helper.range_size(adj_deg)};

            for (auto k = start; k<end; ++k) {
                auto k_reverse = helper.reverse(k);

                helper.write_tile_in(k, k_reverse);

                {
                    const auto& lhs_unit = helper.lhs_unit();
                    const auto* rhs_ptr = helper.right_fwd_read(k);

                    for (dimn_t i=0; i<helper.tile_width; ++i) {
                        for (dimn_t j=0; j<helper.tile_width; ++j) {
                            tile[i*helper.tile_width+j] += fn(lhs_unit*rhs_ptr[i*stride+j]);
                        }
                    }
                }

                {
                    const auto* lhs_ptr = helper.lhs_fwd_read(k);
                    const auto& rhs_unit = helper.rhs_unit();
                    for (dimn_t i = 0; i<helper.tile_width; ++i) {
                        for (dimn_t j = 0; j<helper.tile_width; ++j) {
                            tile[i*helper.tile_width+j] += fn(lhs_ptr[i*stride+j]*rhs_unit);
                        }
                    }
                }

                for (deg_t lh_deg=1; lh_deg < helper.tile_letters; ++lh_deg) {
                    auto rh_deg = adj_deg - lh_deg;
                    for (dimn_t i=0; i<helper.tile_width; ++i) {
                        const auto split = helper.split_key(k);
                        const auto& lhs_val = *helper.left_fwd_read(split.first);
                        helper.read_right_tile(helper.combine(key_type(helper.tile_letters-lh_deg, split.first), k));
                        for (dimn_t j=0; j<helper.tile_width; ++j) {
                            tile[i*helper.tile_width+j] += fn(lhs_val*right_rtile[j]);
                        }
                    }
                }

                for (deg_t lh_deg=0; lh_deg<adj_deg; ++lh_deg) {
                    const auto rh_deg = adj_deg - lh_deg;
                    auto split = helper.split_key(rh_deg, k);
                    helper.read_left_tile(helper.reverse(split.first));
                    helper.read_right_tile(split.second);

                    for (dimn_t i=0; i<helper.tile_width; ++i) {
                        for (dimn_t j=0; j<helper.tile_width; ++j) {
                            tile[i*helper.tile_width+j] += fn(left_rtile[i]*right_rtile[j]);
                        }
                    }
                }

                for (deg_t rh_deg=1; rh_deg<helper.tile_letters; ++rh_deg) {
                    const auto lh_deg = adj_deg - rh_deg;
                    for (dimn_t j=0; j<helper.tile_width; ++j) {
                        const auto split = helper.split_key(key_type(rh_deg, j));
                        const auto& rhs_val = *helper.right_fwd_read(helper.combine(k_reverse, helper.reverse(split.first)));
                        helper.read_left_tile(split.second);

                        for (dimn_t i=0; i<helper.tile_width; ++i) {
                            tile[i*helper.tile_width+j] += fn(left_rtile[i]*rhs_val);
                        }
                    }
                }

                helper.write_tile_out(k, k_reverse);
            }

        }

        fma_dense_traditional(helper, fn, 2*helper.tile_letters);
    }

    template <typename Coefficients, typename Alloc, typename Fn>
    void fma_dense(
            free_tensor<Coefficients, dense_vector_view, Alloc>& result,
            const free_tensor<Coefficients, dense_vector_view, Alloc>& lhs,
            const free_tensor<Coefficients, dense_vector_view, Alloc>& rhs,
            Fn fn,
            deg_t max_deg
            )
    {
        using key_type = typename tensor_basis::key_type;
        if (max_deg>result.depth()) {
            max_deg = result.depth();
        }

        auto lhs_deg = lhs.degree();
        auto rhs_deg = rhs.degree();

        deg_t out_degree = std::min(max_deg, lhs_deg+rhs_deg);
        const tensor_basis* basis = &result.basis();
        result.resize(basis->size(out_degree));

        dense_multiplication_helper<Coefficients> helper(result, lhs, rhs);
        if (out_degree > 2*helper.tile_letters) {
            fma_dense_giles(helper, fn, out_degree);
        } else {
            fma_dense_traditional(helper, fn, out_degree);
        }

    }


public:

    template <typename Result, typename Vector1, typename Vector2, typename Fn>
    LAL_TENSOR_COMPAT_RVV(Result, Vector1, Vector2)
    multiply_and_add(Result& result, const Vector1& lhs, const Vector2& rhs, Fn op)
    {

    }

    template <typename Coefficients, typename Alloc, typename Fn>
    void multiply_inplace(dense_vector_view<tensor_basis, Coefficients, Alloc>& lhs,
            const dense_vector_view<tensor_basis, Coefficients, Alloc>& rhs,
            Fn op) noexcept
    {

    }

};

#undef LAL_TENSOR_COMPAT_RVV
#undef LAL_SAME_COEFFS
#undef LAL_IS_TENSOR


template <typename Coefficients,
          template <typename, typename, typename...> class VectorType,
          template <typename> class StorageModel,
          typename... Args>
class free_tensor : public algebra<tensor_basis, Coefficients, free_tensor_multiplication, VectorType, StorageModel>
{
    using algebra_type = algebra<tensor_basis, Coefficients, free_tensor_multiplication, VectorType, StorageModel>;

    using algebra_type::algebra_type;




};

} // namespace lal


#endif //LIBALGEBRA_LITE_FREE_TENSOR_H

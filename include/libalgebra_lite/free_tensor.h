//
// Created by user on 31/07/22.
//

#ifndef LIBALGEBRA_LITE_FREE_TENSOR_H
#define LIBALGEBRA_LITE_FREE_TENSOR_H

#include "implementation_types.h"
#include "libalgebra_lite_export.h"

#include <algorithm>
#include <memory>
#include <mutex>
#include <unordered_map>
#include <utility>
#include <vector>

#include <boost/functional/hash.hpp>
#include <boost/mpl/vector.hpp>

#include "algebra.h"
#include "basis_traits.h"
#include "coefficients.h"
#include "dense_vector.h"
#include "registry.h"
#include "tensor_basis.h"
#include "vector_traits.h"

namespace lal {

template<typename, template<typename, typename> class,
        template<typename> class>
class free_tensor;

namespace dtl {

#define LAL_IS_TENSOR(V)                                                       \
  std::is_same<typename V::basis_type, tensor_basis>::value

#define LAL_SAME_COEFFS(V1, V2)                                                \
  std::is_same<typename V1::coefficient_ring,                                  \
               typename V2::coefficient_ring>::value

#define LAL_TENSOR_COMPAT_RVV(R, V1, V2)                                       \
  std::enable_if_t<LAL_IS_TENSOR(R) && LAL_IS_TENSOR(V1) &&                    \
                   LAL_IS_TENSOR(V2) && LAL_SAME_COEFFS(R, V1) &&              \
                   LAL_SAME_COEFFS(R, V2)>

template<typename Coefficients> class dense_multiplication_helper {
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

    using tensor_type = dense_vector<tensor_basis, Coefficients>;

public:
    dimn_t tile_width;
    dimn_t tile_size;
    deg_t tile_letters;

    dense_multiplication_helper(tensor_type& out, const tensor_type& lhs,
            const tensor_type& rhs)
            :p_basis(&out.basis()), lhs_deg(lhs.degree()), rhs_deg(rhs.degree())
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

        if (lhs.degree()>1) {
            reverse_buffer.resize(p_basis->size(lhs.degree()-1));
        }
    }

    const scalar_type& left_unit() const noexcept { return *left_ptr; }
    const scalar_type& right_unit() const noexcept { return *right_ptr; }

    const scalar_type* left_tile() const noexcept
    {
        return left_read_buffer.data();
    }
    const scalar_type* right_tile() const noexcept
    {
        return right_read_buffer.data();
    }
    scalar_type* write_tile() noexcept { return write_buffer.data(); }

    deg_t lhs_degree() const noexcept { return lhs_deg; }
    deg_t rhs_degree() const noexcept { return rhs_deg; }

    const scalar_type* left_fwd_read(key_type k) const noexcept
    {
        auto offset = p_basis->start_of_degree(static_cast<deg_t>(k.degree()));
        return left_ptr+k.index()*tile_width+offset;
    }
    const scalar_type* right_fwd_read(key_type k) const noexcept
    {
        auto offset = p_basis->start_of_degree(static_cast<deg_t>(k.degree()));
        return right_ptr+k.index()*tile_width+offset;
    }
    scalar_type* fwd_write(key_type k) const noexcept
    {
        auto offset = p_basis->start_of_degree(static_cast<deg_t>(k.degree()));
        return out_ptr+k.index()*tile_width+offset;
    }

    void read_left_tile(key_type k) noexcept
    {
        auto offset = p_basis->start_of_degree(static_cast<deg_t>(k.degree()));
        const auto* reverse_ptr =
                reverse_buffer.data()+k.index()*tile_width+offset;
        std::copy(reverse_ptr, reverse_ptr+tile_width, left_read_buffer.data());
    }
    void read_right_tile(key_type k)
    {
        auto offset = p_basis->start_of_degree(static_cast<deg_t>(k.degree()));
        const auto* fwd_ptr = right_ptr+k.index()*tile_width+offset;
        std::copy(fwd_ptr, fwd_ptr+tile_width, right_read_buffer.data());
    }
    void write_tile_in(key_type k, key_type kr)
    {
        const auto offset = p_basis->start_of_degree(k.degree()+2*tile_letters);
        const auto* in_ptr = out_ptr+k.index()*tile_width+offset;
        auto* tile_ptr = write_tile();
        const auto stride = p_basis->powers()[k.degree()+tile_letters];

        for (dimn_t i = 0; i<tile_width; ++i) {
            for (dimn_t j = 0; j<tile_width; ++j) {
                tile_ptr[i*tile_width+j] = in_ptr[i*stride+j];
            }
        }
    }

    void write_tile_out(key_type k, key_type kr)
    {
        const auto deg = k.degree();
        const auto offset = p_basis->start_of_degree(deg+2*tile_letters);
        const auto stride = p_basis->powers()[deg+tile_letters];

        auto* ptr = out_ptr+k.index()*tile_width+offset;
        auto* tile_ptr = write_tile();

        for (dimn_t i = 0; i<tile_width; ++i) {
            for (dimn_t j = 0; j<tile_width; ++j) {
                ptr[i*stride+j] = tile_ptr[i*tile_width+j];
            }
        }

        if (deg<p_basis->depth()) {
            // Write reverse data
        }
    }

    key_type reverse(key_type k) const noexcept
    {
        const auto width = p_basis->width();
        auto idx = k.index();

        typename key_type::index_type result_idx = 0;
        while (idx) {
            result_idx *= width;
            result_idx += idx%width;
            idx /= width;
        }

        return key_type{k.degree(), result_idx};
    }
    std::pair<key_type, key_type> split_key(key_type k,
            deg_t lhs_size) const noexcept
    {
        auto rhs_size = k.degree()-lhs_size;
        auto split = p_basis->powers()[rhs_size];
        return {key_type(lhs_size, k.index()/split),
                key_type(rhs_size, k.index()%split)};
    }

    dimn_t stride(deg_t deg) const noexcept
    {
        return p_basis->powers()[deg-tile_letters];
    }

    key_type combine(key_type lhs, key_type rhs)
    {
        const auto rhs_deg = rhs.degree();
        const auto shift = p_basis->powers()[rhs_deg];
        return key_type{lhs.degree()+rhs_deg, lhs.index()*shift+rhs.index()};
    }
    dimn_t combine(dimn_t lhs, key_type rhs)
    {
        const auto rhs_deg = rhs.degree();
        const auto shift = p_basis->powers()[rhs_deg];
        return lhs*shift+rhs.index();
    }

    std::pair<dimn_t, dimn_t> range_size(deg_t lhs, deg_t rhs) const noexcept
    {
        const auto& powers = p_basis->powers();
        return {powers[lhs], powers[rhs]};
    }

    dimn_t range_size(deg_t deg) const noexcept { return p_basis->powers()[deg]; }
};

} // namespace dtl


class LIBALGEBRA_LITE_EXPORT free_tensor_multiplier
    : public base_multiplier<free_tensor_multiplier, tensor_basis> {
    deg_t m_width;

public:
    using key_type = typename tensor_basis::key_type;
    using basis_type = tensor_basis;

    explicit free_tensor_multiplier(deg_t width)
        : m_width(width) {}

    static key_type concat_product(const tensor_basis &basis, key_type lhs,
                                   key_type rhs) {
        const auto lhs_deg = lhs.degree();
        const auto rhs_deg = rhs.degree();
        const auto shift = basis.powers()[rhs_deg];

        const auto idx = lhs.index() * shift + rhs.index();
        return key_type(lhs_deg + rhs_deg, idx);
    }

    using product_type =
        boost::container::small_vector<std::pair<key_type, int>, 1>;

    product_type operator()(const tensor_basis &basis, key_type lhs,
                            key_type rhs) const ;
};

class LIBALGEBRA_LITE_EXPORT free_tensor_multiplication
        : public base_multiplication<free_tensor_multiplier, free_tensor_multiplication> {
    using base_type = base_multiplication<free_tensor_multiplier, free_tensor_multiplication>;

    template<typename C> using ctraits = coefficient_trait<C>;

    template<typename C> using sca_ref = typename ctraits<C>::scalar_type&;
    template<typename C>
    using sca_cref = const typename ctraits<C>::scalar_type&;
    template<typename C> using sca_ptr = typename ctraits<C>::scalar_type*;
    template<typename C>
    using sca_rptr = typename ctraits<C>::scalar_type* LAL_RESTRICT;
    template<typename C>
    using sca_cptr = const typename ctraits<C>::scalar_type*;
    template<typename C>
    using sca_crptr = const typename ctraits<C>::scalar_type* LAL_RESTRICT;

    using key_type = typename tensor_basis::key_type;

    template<typename Coefficients, typename Fn>
    void
    fma_dense_traditional(dtl::dense_multiplication_helper<Coefficients>& helper,
            Fn fn, deg_t out_degree) const
    {
        auto lhs_deg = helper.lhs_degree();
        auto rhs_deg = helper.rhs_degree();

        for (deg_t out_deg = out_degree; out_deg>=0; --out_deg) {
            auto lhs_deg_min = std::max(0, out_deg-rhs_deg);
            auto lhs_deg_max = std::min(out_deg, lhs_deg);

            auto* out_ptr = helper.fwd_write(key_type(out_deg, 0));

            for (deg_t lh_deg = lhs_deg_max; lh_deg>=lhs_deg_min; --lh_deg) {
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

    template<typename C, typename Fn>
    LAL_INLINE_ALWAYS static void impl_db0(sca_rptr<C> tile, sca_crptr<C> lhs_ptr,
            sca_cref<C> rhs_unit, dimn_t stride,
            dimn_t tile_width, Fn op) noexcept { }
    template<typename C, typename Fn>
    LAL_INLINE_ALWAYS static void impl_0bd(sca_rptr<C> tile, sca_cref<C> lhs_unit,
            sca_crptr<C> rhs_ptr, dimn_t stride,
            dimn_t tile_width, Fn op) noexcept { }

    template<typename C, typename Fn>
    LAL_INLINE_ALWAYS static void
    impl_mid(sca_rptr<C> tile, sca_crptr<C> lhs_tile, sca_crptr<C> rhs_tile,
            dimn_t stride, dimn_t tile_width, Fn op) noexcept { }

    template<typename C, typename Fn>
    LAL_INLINE_ALWAYS static void
    impl_lb1(sca_rptr<C> tile, sca_cref<C> lhs_val, sca_crptr<C> rhs_tile,
            dimn_t lhs_index, dimn_t tile_width, Fn op) noexcept { }

    template<typename C, typename Fn>
    LAL_INLINE_ALWAYS static void
    impl_1br(sca_rptr<C> tile, sca_crptr<C> lhs_tile, sca_cref<C> rhs_val,
            dimn_t index, dimn_t tile_width, Fn op) noexcept { }

    template<typename Coefficients, typename Fn>
    void fma_dense_tiled(dtl::dense_multiplication_helper<Coefficients>& helper,
            Fn fn, deg_t out_degree) const
    {
        using key_type = tensor_basis::key_type;

        auto lhs_deg = helper.lhs_degree();
        auto rhs_deg = helper.rhs_degree();

        auto* tile = helper.write_tile();
        const auto* left_rtile = helper.left_tile();
        const auto* right_rtile = helper.right_tile();

        for (deg_t out_deg = out_degree; out_deg>2*helper.tile_letters;
             --out_deg) {
            const auto stride = helper.stride(out_deg);
            const auto adj_deg = out_deg-2*helper.tile_letters;

            // end is not actually a valid key, but it serves as a marker.
            key_type start{adj_deg, 0}, end{adj_deg, helper.range_size(adj_deg)};

            for (auto k = start; k<end; ++k) {
                auto k_reverse = helper.reverse(k);

                helper.write_tile_in(k, k_reverse);

                {
                    const auto& lhs_unit = helper.left_unit();
                    const auto* rhs_ptr = helper.right_fwd_read(k);

                    for (dimn_t i = 0; i<helper.tile_width; ++i) {
                        for (dimn_t j = 0; j<helper.tile_width; ++j) {
                            tile[i*helper.tile_width+j] +=
                                    fn(lhs_unit*rhs_ptr[i*stride+j]);
                        }
                    }
                }

                {
                    const auto* lhs_ptr = helper.left_fwd_read(k);
                    const auto& rhs_unit = helper.right_unit();
                    for (dimn_t i = 0; i<helper.tile_width; ++i) {
                        for (dimn_t j = 0; j<helper.tile_width; ++j) {
                            tile[i*helper.tile_width+j] +=
                                    fn(lhs_ptr[i*stride+j]*rhs_unit);
                        }
                    }
                }

                for (deg_t lh_deg = 1; lh_deg<helper.tile_letters; ++lh_deg) {
                    auto rh_deg = adj_deg-lh_deg;
                    for (dimn_t i = 0; i<helper.tile_width; ++i) {
                        const auto split = helper.split_key(k, lh_deg);
                        const auto& lhs_val = *helper.left_fwd_read(split.first);
                        helper.read_right_tile(helper.combine(split.first, k));
                        for (dimn_t j = 0; j<helper.tile_width; ++j) {
                            tile[i*helper.tile_width+j] += fn(lhs_val*right_rtile[j]);
                        }
                    }
                }

                for (deg_t lh_deg = 0; lh_deg<adj_deg; ++lh_deg) {
                    const auto rh_deg = adj_deg-lh_deg;
                    auto split = helper.split_key(k, lh_deg);
                    helper.read_left_tile(helper.reverse(split.first));
                    helper.read_right_tile(split.second);

                    for (dimn_t i = 0; i<helper.tile_width; ++i) {
                        for (dimn_t j = 0; j<helper.tile_width; ++j) {
                            tile[i*helper.tile_width+j] +=
                                    fn(left_rtile[i]*right_rtile[j]);
                        }
                    }
                }

                for (deg_t rh_deg = 1; rh_deg<helper.tile_letters; ++rh_deg) {
                    const auto lh_deg = adj_deg-rh_deg;
                    for (dimn_t j = 0; j<helper.tile_width; ++j) {
                        const auto split = helper.split_key(key_type(rh_deg, j), lh_deg);
                        const auto& rhs_val = *helper.right_fwd_read(
                                helper.combine(k_reverse, helper.reverse(split.first)));
                        helper.read_left_tile(split.second);

                        for (dimn_t i = 0; i<helper.tile_width; ++i) {
                            tile[i*helper.tile_width+j] += fn(left_rtile[i]*rhs_val);
                        }
                    }
                }

                helper.write_tile_out(k, k_reverse);
            }
        }

        fma_dense_traditional(helper, fn, 2*helper.tile_letters);
    }

public:
    template<typename Coeff>
    using dense_tensor_vec = dense_vector<tensor_basis, Coeff>;

    using base_type::base_type;

    using base_type::fma;
    using base_type::multiply_inplace;

    template<typename Coeff, typename Op>
    void fma(dense_tensor_vec<Coeff>& out, const dense_tensor_vec<Coeff>& lhs,
            const dense_tensor_vec<Coeff>& rhs, Op op) const
    {
        fma(out, lhs, rhs, op, out.basis().depth());
    }

    template<typename Coeff, typename Op>
    void fma(dense_tensor_vec<Coeff>& out, const dense_tensor_vec<Coeff>& lhs,
            const dense_tensor_vec<Coeff>& rhs, Op op, deg_t max_degree) const
    {
        const auto& basis = out.basis();
        if (max_degree>=basis.depth()) {
            max_degree = basis.depth();
        }

        deg_t out_degree = std::min(max_degree, lhs.degree()+rhs.degree());

        const auto out_size = basis.size(out_degree);
        if (out.size()<out_size) {
            /*
             * The resize function will look for the smallest dimension larger than
             * the requested dim, so if we give it size, it will look for the smallest
             * dimension greater than size, not simply size. Thus, we subtract 1 to make
             * sure the next smallest dimension is equal to size.
             */
            out.resize(out_size-1);
        }

        dtl::dense_multiplication_helper<Coeff> helper(out, lhs, rhs);
        //        if (out_degree > 2*helper.tile_letters) {
        //            fma_dense_tiled(helper, op, out_degree);
        //        } else {
        fma_dense_traditional(helper, op, out_degree);
        //        }
    }


};



#undef LAL_TENSOR_COMPAT_RVV
#undef LAL_SAME_COEFFS
#undef LAL_IS_TENSOR

template<typename Coefficients, template<typename, typename> class VectorType,
        template<typename> class StorageModel>
class free_tensor
        : public algebra<tensor_basis, Coefficients, free_tensor_multiplication,
                         VectorType, StorageModel> {
    using algebra_type =
            algebra<tensor_basis, Coefficients, free_tensor_multiplication,
                    VectorType, StorageModel>;


    static void resize_to_degree(free_tensor<Coefficients, dense_vector, StorageModel>& arg, deg_t degree) {
        auto size = arg.basis().size(degree);
        /*
         * The resize function will look for the smallest dimension larger than
         * the requested dim, so if we give it size, it will look for the smallest
         * dimension greater than size, not simply size. Thus, we subtract 1 to make
         * sure the next smallest dimension is equal to size.
         */
        arg.base_vector().resize(size-1);
    }

    template <template <typename, typename> class OVT>
    static void resize_to_degree(free_tensor<Coefficients, OVT, StorageModel>& arg, deg_t degree) {}

public:
    using typename algebra_type::basis_type;
    using typename algebra_type::coefficient_ring;
    using typename algebra_type::key_type;
    using typename algebra_type::rational_type;
    using typename algebra_type::scalar_type;

    using typename algebra_type::basis_pointer;
    using typename algebra_type::multiplication_pointer;

    using algebra_type::algebra_type;

    free_tensor(basis_pointer basis, multiplication_pointer mul, scalar_type arg)
            :algebra_type(basis, mul, key_type(0, 0), std::move(arg)) { }

    free_tensor(basis_pointer basis, scalar_type arg)
            :algebra_type(basis, key_type(0, 0), std::move(arg)) { }


    free_tensor& fmexp_inplace(const free_tensor& exp_arg) {
        free_tensor original(*this), x(exp_arg);

        x[key_type(0, 0)] = scalar_type(0);

        auto degree = this->basis().depth();
        resize_to_degree(*this, degree);
        for (deg_t i=degree; i >= 1; --i) {
            this->mul_scal_div(x, rational_type(i), degree - i + 1);
            *this += original;
        }

        return *this;
    }

    free_tensor fmexp(const free_tensor& exp_arg) const {
        free_tensor result(*this), x(exp_arg);

        x[key_type(0, 0)] = scalar_type(0);

        auto degree = this->basis().depth();
        resize_to_degree(result, degree);

        for (deg_t i=degree; i>= 1; --i) {
            result.mul_scal_div(x, rational_type(i), degree - i + 1);
            result += *this;
        }

        return result;
    }

    friend free_tensor exp(const free_tensor &arg) {
        free_tensor result(arg.get_basis(), arg.multiplication(), scalar_type(1));
        free_tensor one(arg.get_basis(), arg.multiplication(), scalar_type(1));

        const auto degree = arg.basis().depth();
        resize_to_degree(result, degree);
        for (deg_t i = degree; i >= 1; --i) {
            result.mul_scal_div(arg, rational_type(i));
            result += one;
        }

        return result;
    }

    friend free_tensor log(const free_tensor &arg) {

        auto x = arg;
        x[typename tensor_basis::key_type(0, 0)] = scalar_type(0);

        free_tensor result(arg.get_basis(), arg.multiplication());
        const auto degree = arg.basis().depth();
        resize_to_degree(result, degree);

        free_tensor one(arg.get_basis(), arg.multiplication(), scalar_type(1));
        for (deg_t i = degree; i >= 1; --i) {
            if (i % 2 == 0) {
                result.sub_scal_div(one, rational_type(i));
            } else {
                result.add_scal_div(one, rational_type(i));
            }
            result *= x;
        }

        return result;
    }

    friend free_tensor inverse(const free_tensor &arg) {

        const auto &a = arg[key_type(0, 0)];
        assert(a != Coefficients::zero());
        auto x = arg;
        x[key_type(0, 0)] = Coefficients::zero();

        const auto degree = arg.basis().depth();
        free_tensor a_inverse(arg.get_basis(), arg.multiplication(), scalar_type(1) / a);
        free_tensor result(a_inverse);
        resize_to_degree(result, degree);

        auto z = x / a;
        for (deg_t d=0; d<degree; ++d) {
            result = a_inverse + z*result;
        }

        return result;
    }

    friend free_tensor antipode(const free_tensor& arg) {
        // TODO: replace with implementation of antipode.
        return inverse(arg);
    }

};





LAL_EXPORT_TEMPLATE(multiplication_registry, free_tensor_multiplication)


} // namespace lal

#endif // LIBALGEBRA_LITE_FREE_TENSOR_H

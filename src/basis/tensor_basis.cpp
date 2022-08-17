//
// Created by user on 27/07/22.
//
#include <libalgebra_lite/basis/tensor_basis.h>


namespace alg {


alg::tensor_basis::tensor_basis(deg_t width, deg_t depth) :
        m_width (width), m_depth(depth)
{
    m_powers.reserve(depth+1);
    m_sizes.reserve(depth+2);
    m_powers.push_back(1);
    m_sizes.push_back(1);

    for (dimn_t d = 1; d<=depth; ++d) {
        m_powers.push_back(m_powers.back()*width);
        m_sizes.push_back(1+width*m_sizes.back());
    }
    m_sizes.push_back(1+width*m_sizes.back());
}

}

//
// Created by sam on 07/09/22.
//

#ifndef LIBALGEBRA_LITE_TENSOR_FIXTURE_H
#define LIBALGEBRA_LITE_TENSOR_FIXTURE_H


#include <libalgebra_lite/free_tensor.h>
#include <libalgebra_lite/shuffle_tensor.h>

#include <gtest/gtest.h>

struct TensorFixture : public ::testing::Test
{
    lal::deg_t width = 5;
    lal::deg_t depth = 5;
    std::shared_ptr<const lal::tensor_basis> basis;

    using key_type = typename lal::tensor_basis::key_type;


    TensorFixture() : basis(new lal::tensor_basis(width, depth))
    {}

    key_type key(std::initializer_list<lal::let_t> arg) {
        typename key_type::index_type idx = 0;
        for (auto let : arg) {
            idx *= width;
            idx += let - 1;
        }
        return key_type(arg.size(), idx);
    }

    template <typename... Ints>
    key_type key(Ints... ints) {
        return key({lal::let_t(ints)...});
    }

};




#endif //LIBALGEBRA_LITE_TENSOR_FIXTURE_H

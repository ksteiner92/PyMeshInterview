#ifndef PYMESH_QUANTITY_H
#define PYMESH_QUANTITY_H

//#include <SI/detail/unit.h>
#include <type_traits>

namespace {
//template <template<typename...> class T, typename... Ts>
//concept is_A = std::is_same_v<SI::detail::unit_t<Ts...>,T<Ts...>>;
}
namespace mesh {

template<typename T, typename Unit>
class Quantity {
private:
};

} // namespace mesh

#endif // PYMESH_QUANTITY_H

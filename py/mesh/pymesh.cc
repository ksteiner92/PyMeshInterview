#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "attribute.h"
#include "diffoperators.h"
#include "mesh.h"
#include "poisson.h"
#include "system.h"

namespace py = pybind11;
using rvp = py::return_value_policy;
using namespace pybind11::literals;

PYBIND11_MAKE_OPAQUE(std::vector<mesh::Segment<2, 2>>);
PYBIND11_MAKE_OPAQUE(std::vector<mesh::Segment<2, 1>>);
PYBIND11_MAKE_OPAQUE(std::vector<mesh::Segment<1, 1>>);
PYBIND11_MAKE_OPAQUE(std::vector<mesh::Segment<3, 1>>);
PYBIND11_MAKE_OPAQUE(std::vector<mesh::Segment<3, 2>>);

template <uint Dim>
static py::class_<mesh::Mesh<Dim, 0>, mesh::MeshBase>
declareMesh0D(const py::module &m) {
  using Class = mesh::Mesh<Dim, 0>;
  using PyClass = py::class_<Class, mesh::MeshBase>;
  const std::string name =
      Dim == 0 ? "Mesh0D" : "Mesh" + std::to_string(Dim) + "0D";

  py::class_<typename Class::VerticesProxy, mesh::MeshElementsProxy> proxy(
      m, ("VerticesProxy" + std::to_string(Dim) + "D").c_str());
  proxy.def("add_points",
            py::overload_cast<const EigenDRef<const Eigen::MatrixXd> &>(
                &Class::VerticesProxy::add),
            rvp::reference_internal, py::arg().noconvert());

  PyClass cls(m, name.c_str());
  cls.def(py::init<>())
      .def(py::init<mesh::Mesh<Dim, 0> *>())
      .def_property_readonly("vertices", &Class::vertices,
                             rvp::reference_internal);
  return cls;
}

template <uint Dim>
static py::class_<mesh::Mesh<Dim, 1>, mesh::Mesh<Dim, 0>>
declareMesh1D(const py::module &m) {
  using Class = mesh::Mesh<Dim, 1>;
  using PyClass = py::class_<Class, mesh::Mesh<Dim, 0>>;
  const std::string name =
      Dim == 1 ? "Mesh1D" : "Mesh" + std::to_string(Dim) + "1D";

  py::class_<typename Class::EdgesProxy, mesh::MeshElementsProxy> proxy(
      m, ("EdgesProxy" + std::to_string(Dim) + "D").c_str());

  PyClass cls(m, name.c_str());
  cls.def(py::init<>())
      .def(py::init<mesh::Mesh<Dim, 1> *>())
      .def_property_readonly("edges", &Class::edges, rvp::reference_internal);

  return cls;
}

template <uint Dim>
static py::class_<mesh::Mesh<Dim, 2>, mesh::Mesh<Dim, 1>>
declareMesh2D(const py::module &m) {
  using Class = mesh::Mesh<Dim, 2>;
  using PyClass = py::class_<Class, mesh::Mesh<Dim, 1>>;
  const std::string name =
      Dim == 2 ? "Mesh2D" : "Mesh" + std::to_string(Dim) + "2D";

  py::class_<typename Class::FacesProxy, mesh::MeshElementsProxy> proxy(
      m, ("FacesProxy" + std::to_string(Dim) + "D").c_str());

  PyClass cls(m, name.c_str());
  cls.def(py::init<>())
      .def(py::init<mesh::Mesh<Dim, 2> *>())
      .def_property_readonly("faces", &Class::faces, rvp::reference_internal);
  return cls;
}

template <uint Dim, uint TopDim>
static void declareSimplex(const py::module &m) {
  using Class = mesh::Simplex<Dim, TopDim>;
  using PyClass = py::class_<Class, mesh::MeshElement>;

  std::stringstream ss;
  ss << "Simplex";
  if (Dim == TopDim)
    ss << TopDim << 'D';
  else
    ss << Dim << TopDim << 'D';
  PyClass cls(m, ss.str().c_str());
}

template <uint Dim, uint TopDim>
static void declareSegment(const py::module &m) {
  using Class = mesh::Segment<Dim, TopDim>;
  using PyClass = py::class_<Class>;

  std::stringstream ss;
  ss << "Segment";
  if (Dim == TopDim)
    ss << TopDim << 'D';
  else
    ss << Dim << TopDim << 'D';
  PyClass cls(m, ss.str().c_str());
  cls.def_property_readonly("id", &Class::getID);
  cls.def("get_or_create_attribute",
          &Class::template getOrCreateAttribute<double>,
          rvp::reference_internal);
  cls.def_property_readonly("mesh", py::overload_cast<>(&Class::mesh),
                            rvp::reference_internal);
}

template <uint Dim, uint TopDim>
static void declareInterface(const py::module &m) {
  using Class = mesh::Interface<Dim, TopDim>;
  using PyClass = py::class_<Class, mesh::Segment<Dim, TopDim - 1>>;

  std::stringstream ss;
  ss << "Interface";
  if (Dim == TopDim)
    ss << TopDim << 'D';
  else
    ss << Dim << TopDim << 'D';
  PyClass cls(m, ss.str().c_str());
  cls.def_property_readonly("id", &Class::getID);
  cls.def_property_readonly("mesh", py::overload_cast<>(&Class::mesh),
                            rvp::reference_internal);
}

template <uint Dim, uint TopDim>
static void declareSurface(const py::module &m) {
  using Class = mesh::Surface<Dim, TopDim>;
  using PyClass = py::class_<Class, mesh::Segment<Dim, TopDim - 1>>;

  std::stringstream ss;
  ss << "Surface";
  if (Dim == TopDim)
    ss << TopDim << 'D';
  else
    ss << Dim << TopDim << 'D';
  PyClass cls(m, ss.str().c_str());
  cls.def_property_readonly("id", &Class::getID);
  cls.def_property_readonly("mesh", py::overload_cast<>(&Class::mesh),
                            rvp::reference_internal);
}

template <uint Dim, uint TopDim>
static void declarePoisson(const py::module &m) {
  using PoissonClass = math::Poisson<Dim, TopDim>;
  using PyClassPoisson = py::class_<PoissonClass>;

  std::stringstream poisson_ss;
  poisson_ss << "Poisson";
  if (Dim == TopDim) {
    poisson_ss << TopDim << 'D';
  } else {
    poisson_ss << Dim << TopDim << 'D';
  }

  PyClassPoisson cls_poisson(m, poisson_ss.str().c_str());
  cls_poisson.def(py::init<mesh::System<Dim, TopDim> *>());
  cls_poisson.def("solve", &math::Poisson<Dim, TopDim>::solve);
  cls_poisson.def("grad", &math::Poisson<Dim, TopDim>::grad);
}

template <uint Dim, uint TopDim>
static void declareDifferentialOperators(const py::module &m) {
  using Class = math::DifferentialOperators<Dim, TopDim>;
  using PyClass = py::class_<Class>;

  std::stringstream ss;
  ss << "DifferentialOperators";
  if (Dim == TopDim) {
    ss << TopDim << 'D';
  } else {
    ss << Dim << TopDim << 'D';
  }

  PyClass cls(m, ss.str().c_str());
  cls.def(py::init<mesh::System<Dim, TopDim> *>());
  cls.def("laplace", &Class::laplace);
}

static void declareAttribute(const py::module &m) {
  using AttributeDoubleClass = mesh::Attribute<double>;
  using PyClassAttributeDouble = py::class_<AttributeDoubleClass>;

  PyClassAttributeDouble cls_attribute_double(m, "AttributeDouble");
  cls_attribute_double.def(py::init<mesh::Attribute<double>>())
      .def("__len__", &mesh::Attribute<double>::size)
      .def(
          "__getitem__",
          [](mesh::Attribute<double> *attr, size_t idx) {
            return (*attr)[idx];
          },
          py::is_operator(), rvp::reference_internal)
      .def(
          "__setitem__",
          [](mesh::Attribute<double> *attr, size_t idx, double value) {
            (*attr)[idx] = value;
          },
          py::is_operator(), rvp::reference_internal);
}

template <uint Dim, uint TopDim>
static void declareSystem(const py::module &m) {
  using SystemClass = mesh::System<Dim, TopDim>;
  using PyClassSystem = py::class_<SystemClass>;

  std::stringstream system_ss;
  system_ss << "System";
  if (Dim == TopDim)
    system_ss << TopDim << 'D';
  else
    system_ss << Dim << TopDim << 'D';

  PyClassSystem cls_system(m, system_ss.str().c_str());
  cls_system
      .def("segment",
           py::overload_cast<const std::string &>(&SystemClass::segment),
           rvp::reference_internal)
      .def("segment", py::overload_cast<ID>(&SystemClass::segment),
           rvp::reference_internal)
      .def_property_readonly(
          "mesh", py::overload_cast<>(&SystemClass::mesh, py::const_),
          rvp::reference_internal)
      .def_property_readonly("voronoi", &SystemClass::voronoi,
                             rvp::reference_internal)
      .def("surface", py::overload_cast<ID>(&SystemClass::surface),
           rvp::reference_internal)
      .def("surface",
           py::overload_cast<const std::string &>(&SystemClass::surface),
           rvp::reference_internal)
      .def("interface",
           py::overload_cast<const std::string &, const std::string &>(
               &SystemClass::interface),
           rvp::reference_internal)
      .def("interface", py::overload_cast<ID, ID>(&SystemClass::interface),
           rvp::reference_internal)
      //.def("gradient", &SystemClass::template gradient<double>,
      //     rvp::take_ownership)
      //.def("laplace", &SystemClass::template laplace<double>,
      //     rvp::take_ownership)
      .def("get_raw_address",
           [](SystemClass &foo) { return reinterpret_cast<uint64_t>(&foo); })
      .def("get_or_create_attribute",
           &SystemClass::template getOrCreateAttribute<double>,
           rvp::reference_internal);

  std::stringstream builder_ss;
  builder_ss << "SystemBuilder";
  if (Dim == TopDim) {
    builder_ss << TopDim << 'D';
  } else {
    builder_ss << Dim << TopDim << 'D';
  }
  py::class_<typename mesh::System<Dim, TopDim>::Builder> cls_builder(
      m, builder_ss.str().c_str());
  cls_builder.def(py::init<>())
      .def("add_vertices", &mesh::System<Dim, TopDim>::Builder::addVertices,
           rvp::reference_internal)
      .def_property_readonly("mesh",
                             &mesh::System<Dim, TopDim>::Builder::getMesh,
                             rvp::reference_internal)
      .def("create", &mesh::System<Dim, TopDim>::Builder::create,
           rvp::take_ownership)
      .def(
          "add_segment",
          [](typename mesh::System<Dim, TopDim>::Builder *obj,
             std::string name) { return obj->addSegment(name); },
          rvp::take_ownership);

  if constexpr (TopDim > 0) {
    std::stringstream segment_builder_ss;
    segment_builder_ss << "SegmentBuilder";
    if (Dim == TopDim) {
      segment_builder_ss << TopDim << 'D';
    } else {
      segment_builder_ss << Dim << TopDim << 'D';
    }
    py::class_<
        typename mesh::Segment<Dim, TopDim>::template Builder<Dim, TopDim>>
        cls_segment_builder(m, segment_builder_ss.str().c_str());
    cls_segment_builder
        .def(
            "add_edges",
            [](typename mesh::Segment<Dim, TopDim>::template Builder<
                   Dim, TopDim> *obj,
               const EigenDRef<const MatrixXid> &edges) {
              return obj->addEdges(edges);
            },
            rvp::reference_internal)
        .def("add_segment",
             &mesh::Segment<Dim, TopDim>::template Builder<Dim,
                                                           TopDim>::addSegment,
             rvp::take_ownership)
        .def(
            "add_hole",
            &mesh::Segment<Dim, TopDim>::template Builder<Dim, TopDim>::addHole,
            rvp::take_ownership);
  }
}

PYBIND11_MODULE(pymesh, m) {
  py::class_<mesh::MeshElement>(m, "MeshElement")
      .def_property_readonly("num_vertices", &mesh::MeshElement::getNumVertices)
      .def_property_readonly("id", &mesh::MeshElement::getID)
      .def(
          "__getitem__",
          [](mesh::MeshElement *melm, size_t idx) { return (*melm)[idx]; },
          py::is_operator(), rvp::reference_internal)
      .def("__len__", &mesh::MeshElement::getNumVertices)
      .def("point",
           py::overload_cast<size_t>(&mesh::MeshElement::getPoint, py::const_),
           rvp::reference_internal)
      .def_property_readonly("points", &mesh::MeshElement::getPoints, rvp::move)
      .def_property_readonly("center", &mesh::MeshElement::center, rvp::move);

  py::class_<mesh::MeshElementsProxy>(m, "MeshElementsProxy")
      .def("__len__", &mesh::MeshElementsProxy::size)
      .def(
          "__getitem__",
          [](mesh::MeshElementsProxy *obj, size_t idx) { return (*obj)[idx]; },
          py::is_operator(), rvp::reference_internal)
      .def(
          "__iter__",
          [](mesh::MeshElementsProxy &obj) {
            return py::make_iterator(obj.begin(), obj.end());
          },
          py::keep_alive<
              0, 1>() /* Essential: keep object alive while iterator exists */)
      .def("create", &mesh::MeshElementsProxy::create, rvp::reference_internal)
      .def("add",
           py::overload_cast<const EigenDRef<const MatrixXid> &>(
               &mesh::MeshElementsProxy::add),
           rvp::reference_internal);

  py::class_<mesh::MeshBase>(m, "MeshBase")
      .def_property_readonly("bodies", &mesh::MeshBase::bodies,
                             rvp::reference_internal)
      .def_property_readonly("facets", &mesh::MeshBase::facets,
                             rvp::reference_internal)
      .def_property_readonly("ridges", &mesh::MeshBase::ridges,
                             rvp::reference_internal)
      .def_property_readonly("peaks", &mesh::MeshBase::peaks,
                             rvp::reference_internal)
      .def_property_readonly(
          "pointlist",
          py::overload_cast<>(&mesh::MeshBase::getPointList, py::const_),
          rvp::reference_internal);

  declareSimplex<1, 0>(m);
  declareSimplex<2, 0>(m);
  declareSimplex<3, 0>(m);
  declareSimplex<1, 1>(m);
  declareSimplex<2, 1>(m);
  declareSimplex<3, 1>(m);
  declareSimplex<2, 2>(m);
  declareSimplex<3, 2>(m);
  declareSimplex<3, 3>(m);

  declareMesh0D<1>(m);
  declareMesh0D<2>(m);
  declareMesh0D<3>(m);

  declareMesh1D<1>(m);
  declareMesh1D<2>(m);
  declareMesh1D<3>(m);
  declareMesh2D<2>(m);
  declareMesh2D<3>(m);

  declareSegment<1, 0>(m);
  declareSegment<2, 0>(m);
  declareSegment<3, 0>(m);
  declareSegment<1, 1>(m);
  declareSegment<2, 1>(m);
  declareSegment<3, 1>(m);
  declareSegment<2, 2>(m);
  declareSegment<3, 2>(m);

  declareInterface<1, 1>(m);
  declareInterface<2, 1>(m);
  declareInterface<3, 1>(m);
  declareInterface<2, 2>(m);
  declareInterface<3, 2>(m);

  declareSurface<1, 1>(m);
  declareSurface<2, 1>(m);
  declareSurface<3, 1>(m);
  declareSurface<2, 2>(m);
  declareSurface<3, 2>(m);

  declareSystem<1, 1>(m);
  declareSystem<2, 1>(m);
  declareSystem<2, 2>(m);

  declarePoisson<1, 1>(m);
  declarePoisson<2, 1>(m);
  declarePoisson<2, 2>(m);
  declarePoisson<3, 1>(m);
  declarePoisson<3, 2>(m);
  declarePoisson<3, 3>(m);

  declareDifferentialOperators<1, 1>(m);
  declareDifferentialOperators<2, 1>(m);
  declareDifferentialOperators<2, 2>(m);
  declareDifferentialOperators<3, 1>(m);
  declareDifferentialOperators<3, 2>(m);
  declareDifferentialOperators<3, 3>(m);

  py::enum_<mesh::StorageLocation>(m, "StorageLocation")
      .value("VERTEX", mesh::StorageLocation::VERTEX)
      .value("EDGE", mesh::StorageLocation::EDGE)
      .value("FACE", mesh::StorageLocation::FACE)
      .value("CELL", mesh::StorageLocation::CELL)
      .export_values();

  declareAttribute(m);
}

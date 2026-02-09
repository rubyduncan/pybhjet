#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/buffer_info.h>
#include <pybind11/functional.h>

#include "bhjet_class.hh"
#include "jet_output.hh"

namespace py = pybind11; 

#define SET_ARGS_VEC(classtype, function, type) \
	[](classtype &a, py::array_t<type> array) { \
		py::buffer_info buf = array.request(); \
	    if (buf.ndim != 1) { \
	        throw std::runtime_error("Number of dimensions must be one and/ or array size is wrong"); \
		} \
		a.function(std::vector<type>(array.data(), array.data() + array.size()));}


#define GET_ARGS_VEC(classtype, function, type)        \
	[](classtype &a) {                                 \
		return py::array_t<type>(                      \
			a.function().size(),                       \
			a.function().data());                      \
	},                                                 \
		py::return_value_policy::copy


PYBIND11_MODULE(pybhjet, m){

    py::class_<NumDenPoint>(m, "NumDenPoint")
        .def(py::init<>())
        .def_readonly("momentum", &NumDenPoint::momentum, "Momentum value (p [g cm^-1])")
        .def_readonly("gamma", &NumDenPoint::gamma, "Lorentz factor (g [])")
        .def_readonly("n_p", &NumDenPoint::n_p, "Number density n(p) [# cm^-3 p^-1]")
        .def_readonly("n_g", &NumDenPoint::n_g, "Number density n(g) [# cm^-3 g^-1]");

    py::class_<DataPoint>(m, "DataPoint")
        .def(py::init<>())
        .def_readonly("energy", &DataPoint::energy, "Energy value (nu [Hz])")
        .def_readonly("flux", &DataPoint::flux, "Flux value (mJy)");

    py::class_<JetOutput::JetZoneProperties>(m, "JetZoneProperties")
        .def(py::init<>())
        .def_readonly("jet_bfield", &JetOutput::JetZoneProperties::jet_bfield, "Jet magnetic field")
        .def_readonly("lepton_ndens", &JetOutput::JetZoneProperties::lepton_ndens, "Lepton number density")
        .def_readonly("speed_gamma", &JetOutput::JetZoneProperties::speed_gamma, "Speed gamma")
        .def_readonly("delta", &JetOutput::JetZoneProperties::delta, "Doppler factor delta")
        .def_readonly("tshift", &JetOutput::JetZoneProperties::tshift, "Temperature shift")
        .def_readonly("temp_kev", &JetOutput::JetZoneProperties::temp_kev, "Electron temperature in keV")
        .def_readonly("grid_r", &JetOutput::JetZoneProperties::grid_r, "Grid radius in Rg")
        .def_readonly("delz", &JetOutput::JetZoneProperties::delz, "Zone size in z/Rg")
        .def_readonly("dist_z", &JetOutput::JetZoneProperties::dist_z, "Distance z/Rg")
        .def_readonly("z_delz", &JetOutput::JetZoneProperties::z_delz, "Distance z+delz/Rg")
        .def_readonly("equpar_check", &JetOutput::JetZoneProperties::equpar_check, "Equipartition check")
        .def_readonly("ue_ub", &JetOutput::JetZoneProperties::ue_ub, "Energy density ratio (Ue/Ub)");

    py::class_<JetOutput::JetBaseProperties>(m, "JetBaseProperties")
        .def(py::init<>())
        .def_readonly("pair_content", &JetOutput::JetBaseProperties::pair_content, "Pair content (ne/np)")
        .def_readonly("init_mag", &JetOutput::JetBaseProperties::init_mag, "Initial magnetization")
        .def_readonly("particle_avg_lorentz_factor", &JetOutput::JetBaseProperties::particle_avg_lorentz_factor, "Particle average Lorentz factor")
        .def_readonly("jet_nozzle_end", &JetOutput::JetBaseProperties::jet_nozzle_end, "Jet nozzle ends at")
        .def_readonly("jet_nozzle_optical_depth", &JetOutput::JetBaseProperties::jet_nozzle_optical_depth, "Jet nozzle optical depth");

    py::class_<JetOutput::SpectralProperties>(m, "SpectralProperties")
        .def(py::init<>())
        .def_readonly("disk_lum", &JetOutput::SpectralProperties::disk_lum, "Observed 0.3-5 keV disk luminosity")
        .def_readonly("IC_lum", &JetOutput::SpectralProperties::IC_lum, "Observed 0.3-300 keV Inverse Compton luminosity")
        .def_readonly("xray_lum", &JetOutput::SpectralProperties::xray_lum, "Observed 1-10 keV total luminosity")
        .def_readonly("radio_lum", &JetOutput::SpectralProperties::radio_lum, "Observed 4-6 GHz luminosity")
        .def_readonly("xray_index", &JetOutput::SpectralProperties::xray_index, "X-ray 10-100 keV photon index estimate")
        .def_readonly("radio_index", &JetOutput::SpectralProperties::radio_index, "Radio 10-100 GHz spectral index estimate")
        .def_readonly("jetbase_compactness", &JetOutput::SpectralProperties::jetbase_compactness, "Jet base compactness");

    py::class_<JetOutput::JetProfile>(m, "JetProfile")
        .def(py::init<>())
        .def_readonly("z_rg", &JetOutput::JetProfile::z_rg, "Distance along jet axis (z/Rg)")
        .def_readonly("zone_rg", &JetOutput::JetProfile::zone_rg, "Zone radius (R/Rg)")
        .def_readonly("zone_bfield", &JetOutput::JetProfile::zone_bfield, "Magnetic field in the zone")
        .def_readonly("zone_lepdens", &JetOutput::JetProfile::zone_lepdens, "Lepton number density in the zone")
        .def_readonly("zone_gamma", &JetOutput::JetProfile::zone_gamma, "Lorentz factor in the zone")
        .def_readonly("zone_eltemp", &JetOutput::JetProfile::zone_eltemp, "Electron temperature in the zone");

    // Expose JetOutput class
    py::class_<JetOutput>(m, "JetOutput")
        .def(py::init<>())
        .def_readonly("presyn", &JetOutput::presyn)
        .def_readonly("postsyn", &JetOutput::postsyn)
        .def_readonly("precom", &JetOutput::precom)
        .def_readonly("postcom", &JetOutput::postcom)
        .def_readonly("disk", &JetOutput::disk)
        .def_readonly("bb", &JetOutput::bb)
        .def_readonly("total", &JetOutput::total)
        .def_readonly("cyclosyn_zones", &JetOutput::cyclosyn_zones)
        .def_readonly("compton_zones", &JetOutput::compton_zones)
        .def_readonly("numdens", &JetOutput::numdens)
        .def_readonly("jetprofile", &JetOutput::jetprofile)
        .def_readonly("spectral_properties", &JetOutput::spectral_properties)
        .def_readonly("jet_base_properties", &JetOutput::jet_base_properties)
        .def_readonly("jet_zone_properties", &JetOutput::jet_zone_properties)
        ;

    // Expose BhJetClass - for running 
    py::class_<BhJetClass>(m, "PyBHJet")
        .def(py::init<>(), "Initialize the BHJet model.")
        .def("load_params", &BhJetClass::load_params, "Load parameters from a file.")
        .def("print_parameters", &BhJetClass::print_parameters, "Print all parameters with units.")
        .def("run", &BhJetClass::run, "Run the BHJet model.")
        .def("run_singlezone", &BhJetClass::run_singlezone, "Run the BHJet Single Zone model.")
        .def("get_output", &BhJetClass::get_output, py::return_value_policy::reference, "Retrieve the output from the run.")
        // Expose generic parameter getter and setter
        .def("get_parameter", &BhJetClass::get_parameter, "Get the value of a parameter by name.")
        .def("set_parameter", &BhJetClass::set_parameter, "Set the value of a parameter by name.")
        .def("get_parameter_names", &BhJetClass::get_parameter_names, "Get the names of all parameters.")
        // Implement __getitem__ and __setitem__ for dictionary-like access
        .def("__getitem__", &BhJetClass::get_parameter, "Get the value of a parameter by name.")
        .def("__setitem__", &BhJetClass::set_parameter, "Set the value of a parameter by name.")
        //cutoff function switch
        .def_property("cutoff_type", &BhJetClass::get_cutoff_type, &BhJetClass::set_cutoff_type)
        ;
}
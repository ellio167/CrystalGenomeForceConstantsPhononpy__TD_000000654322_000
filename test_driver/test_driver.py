"""
Force Constants Phononpy OpenKIM Crystal Genome Test Driver
==========================================

This is an example demonstrating usage of the kim-tools package. See https://kim-tools.readthedocs.io for more information.

"""

from kim_tools import CrystalGenomeTestDriver
from kim_tools import get_stoich_reduced_list_from_prototype
from kim_python_utils.ase.core import get_model_energy_cutoff
from numpy.linalg import norm
from math import ceil

class TestDriver(CrystalGenomeTestDriver):
    def _calculate(self, force_influence_distance_fallback: float = 10.0, **kwargs):
        """
        Computes the force constants for a crystal structure

        Args:
            force_influence_distance_fallback
                Fallback value for force influence distance (angstroms)
        """

        cell_replicas = []
        force_constants_matrix = []
        disclaimer = None

        # Find Model's energy influence distance
        energy_influence_distance = 0.0
        for sp1 in self.stoichiometric_species:
            for sp2 in self.stoichiometric_species:
                try:
                    eid = get_model_energy_cutoff(self.kim_model_name, [sp1, sp2])
                except:
                    eid = force_influence_distance_fallback/2.0

                    if eid > energy_influence_distance:
                        energy_influence_distance = eid
        force_influence_distance = 2.0*energy_influence_distance


        # The base class provides self.atoms, an ase.Atoms object representing the initial configuration of the crystal.
        # Use this configuration to evaluate the material property of interest
        original_cell = self.atoms.get_cell()
        num_atoms = len(self.atoms)

        # Compute cell_replicas:
        min_cell_diagonal = min(norm(original_cell[0] - original_cell[1]),
                                norm(original_cell[1] - original_cell[2]),
                                norm(original_cell[2] - original_cell[0]))
        linear_replicas = ceil(force_influence_distance/min_cell_diagonal)
        cell_replicas = [linear_replicas, linear_replicas, linear_replicas]  # Is there a better way to do this?!?

        # Compute force constants
        force_constants_matrix, supercell, supercell_idents = calc_force_constant_matrix(atoms_primitive, cell_replicas, self._calc)

        # Compute band strucuture, DOS, and phonon-stability
        # @@@ TODO @@@ make call to compute post-processed values using force constants
        # @@@ TODO @@@ producing the following output variables"
        # @@@ TODO @@@   wave_vectors
        # @@@ TODO @@@   special_points_indices
        # @@@ TODO @@@   special_points_labels
        # @@@ TODO @@@   response_frequency
        # @@@ TODO @@@   response_eigenvector
        # @@@ TODO @@@   response_frequency_ordered
        # @@@ TODO @@@
        # @@@ TODO @@@   frequency
        # @@@ TODO @@@   density
        # @@@ TODO @@@
        # @@@ TODO @@@   unstable



        # Now it is time to write the output in the format you created in your Property Definition. The base class provides utility methods
        # to facilitate this process.
        # This method initializes the Property Instance and adds the keys common to all Crystal Genome properties.
        # property_name can be the full "property-id" field in your Property Definition, or the "Property Name",
        # which is just the short name after the slash, as used here. You can also specify whether your property
        # includes stress and temperature (no by default), and have the option to specify a disclaimer.

        # force-constants-matrix-crystal property instance
        self._add_property_instance_and_common_crystal_genome_keys("force-constants-matrix-crystal",
                                                                   write_stress=True, write_temp=False, disclaimer=disclaimer)

        # This method adds additional fields to your property instance by specifying the key names you defined
        # in your property definition and providing units if necessary.
        self._add_key_to_current_property_instance("supercell-atom-mapping",supercell_idents,
                                                   units=None)

        # You may also provide a dictionary supplying uncertainty information. It is optional, and normally
        # would not be reported for a deterministic calculation like this, only one involving some statistics,
        # such as molecular dynamics or fitting. There are many possible ways to report uncertainty information,
        # detailed at https://openkim.org/doc/schema/properties-framework/
        uncertainty_info = None
        #
        #uncertainty_info = {
        #    "source-std-uncert-value": ??????)
        #}
        self._add_key_to_current_property_instance("force-constants-matrix",force_constants_matrix,
                                                   units="eV/angstrom^2",uncertainty_info=uncertainty_info)

        # @@@ TODO @@@ write coordiantes file using `original_cell`
        self._add_file_to_current_property_instantce("primitive-coordinates-file", "?????????")

        # @@@ TODO @@@ write supercell file using `supercell`
        self._add_file_to_current_property_instantce("supercell-coordinates-file", "?????????")


        # If your Test Driver reports multiple Property Instances, repeat the process above for each one.

        # phonon-band-structure property instance
        self._add_property_instance_and_common_crystal_genome_keys("phonon-band-structure",
                                                                   write_stress=True, write_temp=False, disclaimer=disclaimer)

        self._add_key_to_current_property_instance("special-points-indices", special_points_indices,
                                                   units=None)
        self._add_key_to_current_property_instance("special-points-labels", special_points_labels,
                                                   units=None)
        self._add_key_to_current_property_instance("response-frequency-ordered", response_frequency_ordered,
                                                   units=None)

        uncertainty_info = None
        self._add_key_to_current_property_instance("wave-vectors", wave_vectors,
                                                   units="?????????",uncertainty_info=uncertainty_info)
        self._add_key_to_current_property_instance("response-frequency", response_frequency,
                                                   units="?????????",uncertainty_info=uncertainty_info)
        self._add_key_to_current_property_instance("response-eigenvector", response_eigenvector,
                                                   units="?????????",uncertainty_info=uncertainty_info)


        # phonon-density-of-states property instance
        self._add_property_instance_and_common_crystal_genome_keys("phonon-density-of-states",
                                                                   write_stress=True, write_temp=False, disclaimer=disclaimer)

        uncertainty_info = None
        self._add_key_to_current_property_instance("frequency", frequency,
                                                   units="?????????",uncertainty_info=uncertainty_info)
        self._add_key_to_current_property_instance("density", density,
                                                   units="?????????",uncertainty_info=uncertainty_info)

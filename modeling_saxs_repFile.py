import IMP
import IMP.saxs
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.saxs
import IMP.pmi.restraints.basic
import argparse
import sys

pdbCA_ = False

saxs_data_ = '../data/glover2_comp2_refined_sx.dat_datgnom.gnom_extrapolated.dat'
outdir_ = 'saxsTests_changeweight/roothierw1000'
topologydir_ = '../data/'

def main(pdbCA, saxs_data, outdir):
    model = IMP.Model()

    #TOPOLOGY FILE
    topologyA = topologydir_ + "topologyA.txt"
    topologyB = topologydir_ +"topologyB.txt"
    topologyAB = topologydir_ +"topologyAB.txt"
    tA = IMP.pmi.topology.TopologyReader(topologyA)
    tB = IMP.pmi.topology.TopologyReader(topologyB)
    tAB = IMP.pmi.topology.TopologyReader(topologyAB)

    bs = IMP.pmi.macros.BuildSystem(model)
    bs.add_state(tA)

    root_hierarchy, dof_obj = bs.execute_macro()

    sele = IMP.atom.Selection(root_hierarchy, resolution=1, molecules=["pnkpA","xrccA","xrccB","lig4"])
    residues = IMP.pmi.tools.OrderedSet(sele.get_selected_particles())

    output_objects = []
    for s in root_hierarchy.get_children():
        for m in s.get_children():
            cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(m,scale=0.1)
            cr.add_to_model()
            output_objects.append(cr)

    ebr = IMP.pmi.restraints.basic.ExternalBarrier(root_hierarchy, radius=220.0, resolution=1)
    ebr.add_to_model()
    output_objects.append(ebr)

    if pdbCA:
        sr = IMP.pmi.restraints.saxs.SAXSRestraint(input_objects=list(residues),
                                                saxs_datafile=saxs_data,
                                                ff_type=IMP.saxs.RESIDUES,
                                                label="Saxs restraint",
                                                weight=10.0,
                                                maxq=0.15)
    
    model.update()
    # Run replica exchange Monte Carlo sampling
    rex=IMP.pmi.macros.ReplicaExchange0(model,
        root_hier=root_hierarchy,                          # pass the root hierarchy

        monte_carlo_sample_objects=dof_obj.get_movers(),   # pass MC movers
        global_output_directory=outdir,
        output_objects=output_objects,
        monte_carlo_temperature=2.5,
        replica_exchange_minimum_temperature=2.5,
        replica_exchange_maximum_temperature=10.0,
        replica_exchange_swap=True,
        num_sample_rounds=8,
        monte_carlo_steps=10,
        number_of_best_scoring_models=10,                    # set >0 to store best PDB files (but this is slow to do online)
        number_of_frames=100,                                # increase number of frames to get better results!
        save_coordinates_mode="25th_score",
        nframes_write_coordinates=1)
    rex.execute_macro()

main(pdbCA_, saxs_data_, outdir_)

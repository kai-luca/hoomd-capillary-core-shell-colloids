// Copyright (c) 2009-2024 The Regents of the University of Michigan.
// Derived from part of HOOMD-blue, released under the BSD 3-Clause License.
//
#ifndef __CAPILLARY_SOFT_SHELL_H__
#define __CAPILLARY_SOFT_SHELL_H__
#include <hoomd/md/PotentialPair.h>
#ifdef ENABLE_HIP
#include <hoomd/md/PotentialPairGPU.h>
#endif //ENABLE_HIP

#ifndef __HIPCC__
#include <string>
#endif

#include "hoomd/HOOMDMath.h"

/*! \file CapillarySoftShell.h
 *     \brief Defines the interaction of two particle attracted through capillary and repulsed via polymers.
 *
 *     We split the repulsive and attractive part in two interactions
 *     so they can be used independently from each other.
 *     */

// need to declare these class methods with __device__ qualifiers when building in nvcc
// // DEVICE is __host__ __device__ when included in nvcc and blank when included into the host
// // compiler
#ifdef __HIPCC__
#define DEVICE __device__
#define HOSTDEVICE __host__ __device__
#else
#define DEVICE
#define HOSTDEVICE
#endif
namespace hoomd
    {
namespace md
    {

class CapillaryInteraction
    {
    public:
        struct param_type
            {
            Scalar cap_charge;

            DEVICE void load_shared(char*& ptr, unsigned int& available_bytes) { }

            HOSTDEVICE void allocate_shared(char*& ptr, unsigned int& available_bytes) const { }
#ifndef __HIPCC__
            param_type() : cap_charge(1) { }
            param_type(pybind11::dict v, bool managed = false)
                {
                    cap_charge = v["cap_charge"].cast<Scalar>();
                }
            pybind11::dict asDict()
                {
                pybind11::dict v;
                v["cap_charge"] = cap_charge;
                return v;
                }
#endif
            }
#if HOOMD_LONGREAL_SIZE == 32
        __attribute__((aligned(8)));
#else
        __attribute__((aligned(16)));
#endif
        //! Construct the CapillaryAtraction evaluator
        /*! \param _rsq Squared particle distance
            \param _rcutsq Squared 
            \param _params
        */
        DEVICE CappilaryInteraction(Scalar _rsq, Scalar _rcutsq, const param_type& _params): rsq(_rsq), rcutsq(_rcutsq), cap_charge(_params.cap_charge)
        {
        }
    // Not based on charges.
    DEVICE static bool needsCharge(){
        return false;
    }
    DEVICE void setCharge(Scalar qi, Scalar qj){}
    DEVICE bool evalForceAndEnergy(Scalar& force_divr, Scalar& pair_eng, bool energy_shift){
        if (rsq<rcutsq){
            Scalar r = fast::sqrt(rsq);
            force_divr = - cap_charge * cap_charge/r;
            pair_eng = - cap_charge*cap_charge*Scalar(log(r));
            if (energy_shift){
                Scalar rcut = fast::sqrt(rcutsq);
                pair_eng -= cap_charge*cap_charge*Scalar(log(r_cut))
            }
            return true;
        }
        else{
            return false;
        }
    }
    DEVICE Scalar evalPressureLRCIntegral(){return 0;}
    DEVICE Scalar evalEnergyLRCIntegral(){return 0;}
        protected:
    Scalar rsq;
    Scalar rcutsq;
    Scalar cap_charge;
    };
class SoftShell
    {
    public:
        struct param_type
            {
            Scalar rcut;
            Scalar strands;

            DEVICE void load_shared(char*& ptr, unsigned int& available_bytes) { }

            HOSTDEVICE void allocate_shared(char*& ptr, unsigned int& available_bytes) const { }
#ifndef __HIPCC__
            param_type() : rcut(1), strands(1){ }
            param_type(pybind11::dict v, bool managed = false)
                {
                    rcut = v["rcut"].cast<Scalar>();
                    strands = v["strands"].cast<Scalar>();
                }
            pybind11::dict asDict()
                {
                pybind11::dict v;
                v["rcut"] = rcut;
                v["strands"] = strands;
                return v;
                }
#endif
            }
#if HOOMD_LONGREAL_SIZE == 32
        __attribute__((aligned(8)));
#else
        __attribute__((aligned(16)));
#endif
    };
}
}

#endif //__CAPILLARY_SOFT_SHELL_H__

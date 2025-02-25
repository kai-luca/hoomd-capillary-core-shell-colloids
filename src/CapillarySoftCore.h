// Copyright (c) 2009-2024 The Regents of the University of Michigan.
// Part of HOOMD-blue, released under the BSD 3-Clause License.

#ifndef __PAIR_EVALUATOR_EXAMPLE_H__
#define __PAIR_EVALUATOR_EXAMPLE_H__

#ifndef __HIPCC__
#include <string>
#endif

#include "hoomd/HOOMDMath.h"

/*! \file CapillarySoftCore.h
    \brief Defines the pair evaluator class for the example potential
*/

// need to declare these class methods with __device__ qualifiers when building in nvcc
// DEVICE is __host__ __device__ when included in nvcc and blank when included into the host
// compiler
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
    //! Define the parameter type used by this pair potential evaluator
    struct param_type
        {
        Scalar q;     //!< Capillary charge

        DEVICE void load_shared(char*& ptr, unsigned int& available_bytes) { }

        HOSTDEVICE void allocate_shared(char*& ptr, unsigned int& available_bytes) const { }

#ifdef ENABLE_HIP
        //! Set CUDA memory hints
        void set_memory_hint() const
            {
            // default implementation does nothing
            }
#endif

#ifndef __HIPCC__
        param_type() : q(1) { }

        param_type(pybind11::dict v, bool managed = false)
            {
            q = v["q"].cast<Scalar>();
            }

        pybind11::dict asDict()
            {
            pybind11::dict v;
            v["q"] = q;
            return v;
            }
#endif
        }
#if HOOMD_LONGREAL_SIZE == 32
        __attribute__((aligned(8)));
#else
        __attribute__((aligned(16)));
#endif

    //! Constructs the pair potential evaluator
    /*! \param _rsq Squared distance between the particles
        \param _rcutsq Squared distance at which the potential goes to 0
        \param _params Per type pair parameters of this potential
    */
    DEVICE CapillaryInteraction(Scalar _rsq, Scalar _rcutsq, const param_type& _params)
        : rsq(_rsq), rcutsq(_rcutsq), q(_params.q)
        {
        }

    //! Example doesn't use charge
    DEVICE static bool needsCharge()
        {
        return false;
        }
    //! Accept the optional charge value
    /*! \param qi Charge of particle i
        \param qj Charge of particle j
    */
    DEVICE void setCharge(Scalar qi, Scalar qj) { }

    //! Evaluate the force and energy
    /*! \param force_divr Output parameter to write the computed force divided by r.
        \param pair_eng Output parameter to write the computed pair energy
        \param energy_shift If true, the potential must be shifted so that
        V(r) is continuous at the cutoff
        \note There is no need to check if rsq < rcutsq in this method.
        Cutoff tests are performed in PotentialPair.

        \return True if they are evaluated or false if they are not because
        we are beyond the cutoff
    */
    DEVICE bool evalForceAndEnergy(Scalar& force_divr, Scalar& pair_eng, bool energy_shift)
        {
        // compute the force divided by r in force_divr
        if (rsq < rcutsq)
            {
            Scalar r = fast::sqrt(rsq);
            Scalar rinv = 1 / r;

            force_divr = -q*q*rinv;

            pair_eng = -q*q*log(r);

            if (energy_shift)
                {
                Scalar rcut = fast::sqrt(rcutsq);
                pair_eng -= q*q*log(rcut);
                }
            return true;
            }
        else
            return false;
        }

    //! Example doesn't eval LRC integrals
    DEVICE Scalar evalPressureLRCIntegral()
        {
        return 0;
        }

    //! Example doesn't eval LRC integrals
    DEVICE Scalar evalEnergyLRCIntegral()
        {
        return 0;
        }

#ifndef __HIPCC__
    //! Get the name of this potential
    /*! \returns The potential name.
     */
    static std::string getName()
        {
        return std::string("example_pair");
        }

    std::string getShapeSpec() const
        {
        throw std::runtime_error("Shape definition not supported for this pair potential.");
        }
#endif

    protected:
    Scalar rsq;    //!< Stored rsq from the constructor
    Scalar rcutsq; //!< Stored rcutsq from the constructor
    Scalar q;      //!< Stored k from the constructor
    };

class SoftShell
    {
    public:
    //! Define the parameter type used by this pair potential evaluator
    struct param_type
        {
        Scalar f;     //!< Capillary charge
        Scalar r_cor;     //!< Capillary charge

        DEVICE void load_shared(char*& ptr, unsigned int& available_bytes) { }

        HOSTDEVICE void allocate_shared(char*& ptr, unsigned int& available_bytes) const { }

#ifdef ENABLE_HIP
        //! Set CUDA memory hints
        void set_memory_hint() const
            {
            // default implementation does nothing
            }
#endif

#ifndef __HIPCC__
        param_type() : f(1),  r_cor(0){ }

        param_type(pybind11::dict v, bool managed = false)
            {
            f = v["f"].cast<Scalar>();
            r_cor = v["r_cor"].cast<Scalar>();
            }

        pybind11::dict asDict()
            {
            pybind11::dict v;
            v["f"] = f;
            v["r_cor"] = r_cor;
            return v;
            }
#endif
        }
#if HOOMD_LONGREAL_SIZE == 32
        __attribute__((aligned(8)));
#else
        __attribute__((aligned(16)));
#endif

    //! Constructs the pair potential evaluator
    /*! \param _rsq Squared distance between the particles
        \param _rcutsq Squared distance at which the potential goes to 0
        \param _params Per type pair parameters of this potential
    */
    DEVICE SoftShell(Scalar _rsq, Scalar _rcutsq, const param_type& _params)
        : rsq(_rsq), rcutsq(_rcutsq), f(_params.f), r_cor(_params.r_cor)
        {
        }

    //! Example doesn't use charge
    DEVICE static bool needsCharge()
        {
        return false;
        }
    //! Accept the optional charge value
    /*! \param qi Charge of particle i
        \param qj Charge of particle j
    */
    DEVICE void setCharge(Scalar qi, Scalar qj) { }

    //! Evaluate the force and energy
    /*! \param force_divr Output parameter to write the computed force divided by r.
        \param pair_eng Output parameter to write the computed pair energy
        \param energy_shift If true, the potential must be shifted so that
        V(r) is continuous at the cutoff
        \note There is no need to check if rsq < rcutsq in this method.
        Cutoff tests are performed in PotentialPair.

        \return True if they are evaluated or false if they are not because
        we are beyond the cutoff
    */
    DEVICE bool evalForceAndEnergy(Scalar& force_divr, Scalar& pair_eng, bool energy_shift)
        {
        // compute the force divided by r in force_divr
        if (rsq < rcutsq)
            {
            Scalar r = fast::sqrt(rsq);
            Scalar rinv = 1 / r;
            Scalar sqrtf = fast::sqrt(f);
            if (r<r_cor){
                pair_eng = 5/18*f*sqrtf*(-log(r/r_cor)-2/(2+sqrtf));
                //force_divr = -q*q*rinv;
            }
            else{
                pair_eng = 5/18*f*sqrtf*(2/(2+sqrtf))*r_cor*rinv*exp(-sqrtf*(r-r_cor)/(2*r_cor));
            }
            if (energy_shift)
                {
                    Scalar rcut=fast::sqrt(rcutsq);
                    if (r<r_cor){
                        pair_eng -= 5/18*f*sqrtf*(-log(rcut/r_cor)-2/(2+sqrtf));
                        //force_divr = -q*q*rinv;
                    }
                    else{
                        pair_eng -= 5/18*f*sqrtf*(2/(2+sqrtf))*r_cor*exp(-sqrtf*(rcut-r_cor)/(2*r_cor)/rcut);
                    }
                }
            return true;
            }
        else
            return false;
        }

    //! Example doesn't eval LRC integrals
    DEVICE Scalar evalPressureLRCIntegral()
        {
        return 0;
        }

    //! Example doesn't eval LRC integrals
    DEVICE Scalar evalEnergyLRCIntegral()
        {
        return 0;
        }

#ifndef __HIPCC__
    //! Get the name of this potential
    /*! \returns The potential name.
     */
    static std::string getName()
        {
        return std::string("example_pair");
        }

    std::string getShapeSpec() const
        {
        throw std::runtime_error("Shape definition not supported for this pair potential.");
        }
#endif

    protected:
    Scalar rsq;    //!< Stored rsq from the constructor
    Scalar rcutsq; //!< Stored rcutsq from the constructor
    Scalar f;      //!< Stored k from the constructor
    Scalar r_cor;      //!< Stored k from the constructor
    };

    } // end namespace md
    } // end namespace hoomd

#endif // __PAIR_EVALUATOR_EXAMPLE_H__

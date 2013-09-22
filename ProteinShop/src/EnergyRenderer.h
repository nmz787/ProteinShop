/***************************************************************************************
*cr								
*cr					Copyright (c) 2004, The Regents of the 
*cr	University of California, through Lawrence Berkeley National Laboratory, 
*cr	Univ. of Calif. at Davis, and Lawrence Livermore National Laboratory 
*cr	(subject to receipt of any required approvals from U.S. Dept. of Energy).  
*cr							All rights reserved.
*cr		
*cr		Please see the accompanying LICENSE file for further information.
*cr
***************************************************************************************/
/***********************************************************************
EnergyRenderer - Protein energy wrapper for class PaletteRenderer.

***********************************************************************/


class EnergyRenderer;


#ifndef ENERGY_RENDERER_INCLUDED
#define ENERGY_RENDERER_INCLUDED


#include <vector>

#include <Geometry/AffineTransformation.h>

#include "AtomClassifier.h"
#include "Buffer.h"
#include "ColorFunction.h"
#include "Globals.h"
#include "VolumeRenderer.h"


/** OpenGL rendering delegate for the volume texturing of molecular force
    fields.  The formulation of the force field is determined by the
    EnergyCalculator plug-in.  One instance of this object is associated with
    each protein for which there is also an EnergyCalculator. */
class EnergyRenderer
{
public:

    /** The coefficient type of the radial basis function that is used to
        generate the volume texture. */
    enum RadiusType
    {
        UNIFORM_RADIUS,         /// Multiplier x 1 angstrom.
        ATOM_RADIUS,            /// Multiplier x atom radius.
        VAN_DER_WAALS_RADIUS    /// Multiplier x Van der Waals charge radius.
    };

    /** The same transformation that is used to position the protein structure
        within the model space must also be used to position the volume texture
        block. */
    typedef Geometry::AffineTransformation<double,3> Transform;

    /** Create an energy rendering delegate in association with a protein.  It
        is an error to break the association after construction.
        @param state
            The ProteinState instance to which this object will assign itself as
            the energyRenderer field.  Both the protein and energyCalculator
            fields of this object must be valid and non-null.  This object
            should be deleted at the same time as all the other fields in
            @a state. */
    EnergyRenderer (ProteinState *state);

    /** Destroy the rendering delegate.  The clearContext() method should be
        called prior to this in order to release OpenGL texture memory. */
    ~EnergyRenderer() {}

    /** Get the name of the atom classifier.
        @return
            The number of the classifier currently assigned. */
    uint classifier() const { return m_classifier; }

    /** Stovepipe method to release OpenGL texture resources prior to object
        destruction.  This method should be called immediately prior to the
        destruction of this object. */
    void clearContext (GLContextData &contextData);

    /** Erase all accumulated energy data.  Accumulation begins anew with the
        next call to updateEnergy(). */
    void clearEnergyAccumulator();

    /** Get the color function associated with a specified classification.
        @param classification
            The classifier range element to get the color function for.
        @return
            The number of the color function, or UINT_MAX is @a classification
            is out of range. */
    uint colorFunction (uint classification) const;

    /** Invoked by the OpenGL widget to render the volume texture.
        @param contextData
            Dictionary storing data for this object that is specific to the
            OpenGL context in which it is being rendered. */
    void glRenderAction (GLContextData &contextData);

    /** Stovepipe method to initialize OpenGL context data structure prior to
        rendering for the first time in a particular context.
        @param contextData
            Dictionary storing data for this object that is specific to the
            OpenGL context in which it is being rendered. */
    void initContext (GLContextData &contextData);

    /** Determine how the normalizing interval is computed by sample().
        @return
            True if sample() automatically computes the normalizing interval
            based on the sample range.  False if it applies whatever normalizing
            interval was previously set. */
    bool isAutoNormalizing() const { return m_autoNormalize; }

    /** Stovepipe method to determine context data initialization.
        @return
            True if initContext() has been called for the specified context,
            false if not. */
    bool isContextInitialized (GLContextData &contextData) const;

    /** Determine the input channel, set on the last sample() call.
        @return
            True if gradient magnitudes are being visualized, false if the
            subset sums of force field terms are being visualized. */
    bool isGradient() const { return m_gradient; }

    /** Get the upper bound of the normalizing interval.
        @return
            The value used in the most recent call to sample(). */
    double maxNormal() const { return m_maxNormal; }

    /** Get the lower bound of the normalizing interval.
        @return
            The value used in the most recent call to sample(). */
    double minNormal() const { return m_minNormal; }

    /** Get the radius multiplier of the basis function.
        @return
            The radius multiplier used in the most recent call to sample(). */
    double radiusMultiplier() const { return m_radiusMultiplier; }

    /** Determine the coefficient type for the radial basis function.
        @return
            The value most recently passed in to sample(). */
    RadiusType radiusType() const { return m_radiusType; }

    /** Sample the atom energies into the texel grid.  No visual results will be
        seen until the texture is rendered.  The sample is based on all energy
        terms that have been accumulated by this object since the last call to
        this function or to clearEnergy().  The original data is then erased,
        and accumulation begins anew.  If no energy updates have been posted
        since the last call, the current EnergyCalculator state is used. */
    void sample();

    /** Set how the normalizing interval is determined by sample().
        @param autoNormalize
            If true, automatically compute the normalizing interval based on the
            sample range.  If false, apply whatever normalizing interval was
            previously set. */
    void setAutoNormalizing (bool autoNormalize);

    /** Assign the atom classifier for this object to use.  Results will not be
        visible until sample() is called again.  If the specified classifier is
        subsequently deleted, a replacement will be selected automatically.  If
        color functions are not specified for part of the new classifier's
        range, they are selected automatically.
        @param classifier
            The number of the classification function to use.  This assignment
            will be accepted if it falls in the range [0,
            AtomClassifier::numClassifiers()) and is not equal to the existing
            value.  If not accepted, this call has no effect. */
    void setClassifier (uint classifier);

    /** Set the color function associated with a specified classification.  If
        the assignment is accepted and sample() has been called at least once,
        the texture is regenerated, but no results will be visible until until
        the texture is rendered again.  If the specified color function is
        subsequently deleted, a replacement will be selected automatically.
        @param classification
            The classifier range element to set the color function for.  If out
            of range, this call is ignored.
        @param colorFunction
            The number of the color function to assign.  This assignment will be
            accepted if it falls in the range [0, ColorFunction::numFunctions())
            and is not equal to the existing value.  If not accepted, this call
            has no effect. */
    void setColorFunction (uint classification, uint colorFunction);

    /** Set the input channel for the next sample() call.  No results will be
        visible until sample() is called again.
        @param gradient
            True to visualize gradient magnitudes, false to visualize the subset
            sums of force field terms. */
    void setGradient (bool gradient);

    /** Set the normalizing range for the sampling procedure.  This is ignored
        if @a autoNormalize == true.  The texture is not regenerated because it
        would make the GUI choppy.
        @param minNormal
            Lower bound of the voxel unit interval.
        @param maxNormal
            Upper bound of the voxel unit interval. */
    void setNormalRange (double minNormal, double maxNormal);

    /** Set the radius multiplier of the basis function.  No changes will be
        visible until sample() is called.
        @param radiusMultiplier
            The radius multiplier to be used by sample().  If less than zero,
            this call has no effect */
    void setRadiusMultiplier (double radiusMultiplier);

    /** Set the coefficient type for the radial basis function.  No changes will
        be visible until sample() is called.
        @param type
            The coefficient type to be used by sample(). */
    void setRadiusType (RadiusType type);

    /** Set the resolution of the texel block.  No results will be visible until
        sample() is called.
        @param resolution
            The resolution to be used by sample().  If less than zero, this call
            has no effect. */
    void setTexelsPerAngstrom (double resolution);

    /** Get the resolution of the texel block.
        @return
            The resolution used in the most recent call to sample(). */
    double texelsPerAngstrom() const { return m_texelsPerAngstrom; }

    /** Catenate a model space transformation to the existing transformation for
        the texture block.  No results will be visible until the texture is
        rendered again.
        @param transform
            The transformation to apply. */
    void transform (const Transform &transform);

    /** Notify this object that the EnergyCalculator has new data. */
    void updateEnergy();

private:

    // VolumeRenderer terminology differs, although it shouldn't
    typedef VolumeRenderer::Voxel Texel;
    typedef ColorFunction::Color Color;

    bool m_autoNormalize;               // store param for resample()
    uint m_classifier;                  // the atom classification function
    std::vector<uint> m_colorNames;     // color functions for atom classes
    std::vector<ColorFunction*> m_colorPtrs;
    Buffer<double> m_energyTerms;       // accumulator for atom energy terms
    bool m_gradient;                    // store param for resample()
    double m_maxNormal, m_minNormal;    // normalizing interval (voxel unit)
    double m_numUpdatesPending;         // number of updateEnergy() calls
    double m_radiusMultiplier;          // store param for resample()
    RadiusType m_radiusType;            // store param for resample()
    VolumeRenderer m_renderer;          // class that does the OpenGL legwork
    Buffer<double,4> m_samples;         // voxel sampling grids
    ProteinState *m_state;              // protein this obj is associated with
    Buffer<Texel,3> m_texels;           // volume texture subimage
    double m_texelsPerAngstrom;         // store param for resample()
    Transform m_transform;              // model space transformation of block

    // configuration file variables
    static bool s_autoSaveGLState;
    static bool s_configLoaded;
    static VolumeRenderer::InterpolationMode s_interpolationMode;
    static uint s_maxBuffer, s_maxTexDim;
    static VolumeRenderer::RenderingMode s_renderingMode;
    static double s_sliceFactor;
    static bool s_textureCaching;
    static VolumeRenderer::TextureFunction s_textureFunction;
    static VolumeRenderer::VoxelAlignment s_voxelAlignment;

    void averageEnergy();
    void buildTexture();
    AtomClassifier *getClassifier();
    void finishUpdate();
    bool isUpdateAvailable() const { return ( m_numUpdatesPending > 0.0 ); }

    static void calcTexParams();
    static void glTexQuery2D (GLint &width, GLint &height);
    static void glTexQuery3D (GLint &width, GLint &height, GLint &depth);
    static void lazyLoadConfig();
    static Texel mapColorToTexel (const Color &color);
};


#endif

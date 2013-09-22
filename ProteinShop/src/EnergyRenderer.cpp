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


#include <cassert>
#include <cmath>
#include <iostream>
using namespace std;

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glext.h>
#else
#include <GL/gl.h>
#include <GL/glext.h>
#endif

#include <Math/Constants.h>

#include "ColorFunction.h"
#include "EnergyRenderer.h"
#include "GLColorOps.h"
#include "GLTransformations.h"
#include "Timer.h"
using namespace MD;


bool EnergyRenderer::s_configLoaded = false;

// allocate storage for the other configuration variables
bool EnergyRenderer::s_autoSaveGLState;
VolumeRenderer::InterpolationMode EnergyRenderer::s_interpolationMode;
uint EnergyRenderer::s_maxBuffer;
uint EnergyRenderer::s_maxTexDim;
VolumeRenderer::RenderingMode EnergyRenderer::s_renderingMode;
double EnergyRenderer::s_sliceFactor;
bool EnergyRenderer::s_textureCaching;
VolumeRenderer::TextureFunction EnergyRenderer::s_textureFunction;
VolumeRenderer::VoxelAlignment EnergyRenderer::s_voxelAlignment;


EnergyRenderer::EnergyRenderer (ProteinState *state) :
    m_autoNormalize (true),
    m_classifier (0),
    m_energyTerms (state->protein->getNumAtoms()),
    m_gradient (true),
    m_maxNormal (0.0),
    m_minNormal (0.0),
    m_numUpdatesPending (0.0),
    m_radiusMultiplier (1.0),
    m_radiusType (VAN_DER_WAALS_RADIUS),
    m_state (state),
    m_texelsPerAngstrom (1.0)
{
    assert ( m_state && m_state->energyCalculator );
    m_state->energyRenderer = this;
    lazyLoadConfig();
    getClassifier();
    m_energyTerms.initElements (0.0);

    // set up the volume renderer
    m_renderer.setVoxelAlignment (s_voxelAlignment);
    m_renderer.setRenderingMode (s_renderingMode);
    m_renderer.setInterpolationMode (s_interpolationMode);
    m_renderer.setTextureFunction (s_textureFunction);
    m_renderer.setSliceFactor (s_sliceFactor);
    m_renderer.setAutosaveGLState (s_autoSaveGLState);
    m_renderer.setTextureCaching (s_textureCaching);
}


void EnergyRenderer::clearContext (GLContextData &contextData)
{
    m_renderer.clearContext (contextData);
}


void EnergyRenderer::clearEnergyAccumulator()
{
    m_numUpdatesPending = 0.0;
    m_energyTerms.initElements (0.0);
}


uint EnergyRenderer::colorFunction (uint classification) const
{
    if ( classification < m_colorNames.size() )
        return m_colorNames[classification];
    else
        return UINT_MAX;
}


void EnergyRenderer::glRenderAction (GLContextData &contextData)
{
    // only render if at least one sample has been taken
    if ( m_state->visualizeEnergy && m_texels.numElements() > 0 )
    {
        glPushMatrix();
        glMultMatrix (m_transform);
        m_renderer.renderBlock (contextData);
        glPopMatrix();
    }
}


void EnergyRenderer::initContext (GLContextData &contextData)
{
    m_renderer.initContext (contextData);
}


bool EnergyRenderer::isContextInitialized (GLContextData &contextData) const
{
    return m_renderer.isContextInitialized (contextData);
}


void EnergyRenderer::sample()
{
    typedef Geometry::Box<double,3> Box;

    static bool firstCall = true;
    if ( firstCall )
    {
        // make these calls when we're sure OpenGL is ready
        firstCall = false;
        calcTexParams();
        
        // rendering mode might change as a result of calcTexParams()
        m_renderer.setRenderingMode (s_renderingMode);
    }
    // don't break the association between renderer and protein
    assert ( m_state->energyRenderer == this );

    // if there are no gradients, ignore the gradient argument
    if ( !m_state->energyCalculator->isGradientSupported() )
        m_gradient = false;

    // try to get some energy terms if none have been accumulated yet
    if ( isUpdateAvailable() )
        averageEnergy();
    else
        updateEnergy();

    // max atom radius = 2.98, max VDW radius = 2.25
    double maxRadius = m_radiusMultiplier;
    if ( m_radiusType == ATOM_RADIUS )
        maxRadius *= 2.98;
    else if ( m_radiusType == VAN_DER_WAALS_RADIUS )
        maxRadius *= 2.25;

    // auto-generate normalizing range if one has not been set
    bool autoNormalize = m_autoNormalize;
    if ( m_minNormal >= m_maxNormal )
        autoNormalize = true;

    // measure the bounding box
    Box box (Box::empty);
    for ( Protein::ConstAtomIterator itr = m_state->protein->atomsBegin();
          itr != m_state->protein->atomsEnd();
          ++itr )
    {
        Position pos = itr->getPosition();
        box.addPoint (pos);
    }
    // pad some extra space around the protein so it doesn't get cut off
    Box::Point corner = box.getVertex (0);
    for ( uint i = 0; i < 3; ++i )
        corner[i] -= maxRadius;
    box.setVertex (0, corner);
    corner = box.getVertex (7);
    for ( uint i = 0; i < 3; ++i )
        corner[i] += maxRadius;
    box.setVertex (7, corner);

    // compute the aspect ratio and resolution of the texel block
    Box::Size aspect = box.getSize();
    int res[3] = {
        int(round(m_texelsPerAngstrom * aspect[0])),
        int(round(m_texelsPerAngstrom * aspect[1])),
        int(round(m_texelsPerAngstrom * aspect[2]))
    };
    // s_maxBuffer is in megabytes, each double is 8 bytes
    AtomClassifier *classifier = getClassifier();
    double testRes = double(classifier->numClasses()) * 8.0 *
                     double(res[0]) * double(res[1]) * double(res[2]);
    double maxBuffer = s_maxBuffer * 1024.0 * 1024.0;
    double correction = 1.0;
    while ( testRes > maxBuffer )
    {
        // keep buffer sizes under the specified limit
        cout << "density " << m_texelsPerAngstrom
             << " exceeds buffer size limit, substituting ";
        m_texelsPerAngstrom *= pow (
            (maxBuffer - correction) / testRes,
            0.33333333333333
        );
        cout << m_texelsPerAngstrom << " texels per angstrom\n";
        res[0] = int (round(m_texelsPerAngstrom * aspect[0]));
        res[1] = int (round(m_texelsPerAngstrom * aspect[1]));
        res[2] = int (round(m_texelsPerAngstrom * aspect[2]));
        testRes = double(classifier->numClasses()) * 8.0 *
                  double(res[0]) * double(res[1]) * double(res[2]);
        correction *= 2.0;
    }
    bool passedTest = false;
    while ( !passedTest )
    {
        // keep munging the dimensions until OpenGL will swallow it
        if ( s_renderingMode == VolumeRenderer::AXIS_ALIGNED )
        {
            // make one query for each axis; all must succeed
            GLint testDim[3][2] = {
                { res[0], res[1] },
                { res[0], res[2] },
                { res[1], res[2] }
            };
            passedTest = true;
            for ( uint i = 0; i < 3 && passedTest; ++i )
            {
                glTexQuery2D (testDim[i][0], testDim[i][1]);
                passedTest = ( testDim[i][0] && testDim[i][1] );
            }
        }
        else
        {
            // only need to make one query
            GLint testDim[3] = { res[0], res[1], res[2] };
            glTexQuery3D (testDim[0], testDim[1], testDim[2]);
            passedTest = ( testDim[0] && testDim[1] && testDim[2] );
        }
        if ( !passedTest )
        {
            cout << "Texture resolution of " << res[0] << " x " << res[1]
                 << " x " << res[2]
                 << " exceeds OpenGL limits,\nreducing density "
                 << m_texelsPerAngstrom << " by 10% to ";
            m_texelsPerAngstrom *= 0.9;
            cout << m_texelsPerAngstrom << " texels per angstrom\n";
            res[0] = int (round(m_texelsPerAngstrom * aspect[0]));
            res[1] = int (round(m_texelsPerAngstrom * aspect[1]));
            res[2] = int (round(m_texelsPerAngstrom * aspect[2]));
        }
        if ( !res[0] || !res[1] || !res[2] )
        {
            cout << "Texture function unsupported, aborting procedure.\n";
            return;
        }
    }
    // get the range to clamp each atom energy into
    double minClamp = m_state->visualizeEnergyMinRange;
    double maxClamp = m_state->visualizeEnergyMaxRange;
    if ( maxClamp <= minClamp )
        maxClamp = minClamp + 100.0;
    uint texelRadius = (uint) round (maxRadius * m_texelsPerAngstrom);

    cout << "---------------------------------------------------------------\n"
         << "resolution    = " << res[0]
         << " x " << res[1]
         << " x " << res[2] << " texels\n"
         << "aspect ratio  = " << aspect[0]
         << " : " << aspect[1]
         << " : " << aspect[2] << " angstroms\n"
         << "density       = " << m_texelsPerAngstrom
         << " texels per angstrom\n"
         << "basal radius  = " << maxRadius << " angstroms (max)\n"
         << "clamp range   = [" << minClamp << ", " << maxClamp << "]\n"
         << "texel radius  = " << texelRadius << " texels (max)\n"
         << "input channel = "
         << (m_gradient ? "gradient norm\n" : "subset sum\n")
         << "radius type   = ";
    switch ( m_radiusType )
    {
        case UNIFORM_RADIUS:
            cout << "uniform\n";
            break;
        case ATOM_RADIUS:
            cout << "proportional to atom radius\n";
            break;
        case VAN_DER_WAALS_RADIUS:
            cout << "proportional to Van der Waals radius\n";
            break;
    }
    // set new buffer dimensions and clear the buffers
    m_samples.setDimSizes (classifier->numClasses(), res[0], res[1], res[2]);
    m_texels.setDimSizes (res[0], res[1], res[2]);
    m_samples.initElements (0.0);
    m_texels.initElements (0);

    // set up all the various loop variables
    corner = box.getVertex (0);
    uint atomIndex = 0;
    double currentRadius = m_radiusMultiplier;
    double invSqr3 = 3.0 / (m_radiusMultiplier * m_radiusMultiplier);
    double invCub2 = 2.0 / (m_radiusMultiplier * m_radiusMultiplier *
                            m_radiusMultiplier);
    double inc[3] = {
        aspect[0] / ((double)res[0]),
        aspect[1] / ((double)res[1]),
        aspect[2] / ((double)res[2])
    };
    if ( autoNormalize )
    {
        // compute the normalizing clamp automatically
        m_minNormal = Math::Constants<double>::max;
        m_maxNormal = Math::Constants<double>::min;
    }
    double minEnergy = Math::Constants<double>::max;
    double maxEnergy = Math::Constants<double>::min;
    double energy, x, y, z, r;
    Position pVox, pAtom;
    uint vox[3], min[3], max[3], idx[3], sIndex, tIndex, atomClass;
    Timer timer;

    // sample the atom energies into the texel grid
    for ( Protein::ConstAtomIterator itr = m_state->protein->atomsBegin();
          itr != m_state->protein->atomsEnd();
          ++itr, ++atomIndex )
    {
        // classify the atom
        atomClass = classifier->classify (*itr, *m_state);
        
        // find the nearest texel to the atom
        pAtom = itr->getPosition();
        for ( uint i = 0; i < 3; ++i )
        {
            // ensure that rounding errors do not put us out of bounds
            vox[i] = (uint) round (
                ((pAtom[i] - corner[i]) / aspect[i]) * res[i]
            );
            if ( vox[i] >= res[i] )
                vox[i] = res[i] - 1;

            // find the bounds of the texel region to sample into
            min[i] = vox[i] - texelRadius;
            max[i] = vox[i] + texelRadius;
            if ( min[i] > vox[i] ) min[i] = 0;
            if ( max[i] > res[i] ) max[i] = res[i];
        }
        // retrieve, clamp, and accumulate the atom energy into the texels
        energy = m_energyTerms[atomIndex];
        if ( energy < minEnergy ) minEnergy = energy;
        else if ( energy > maxEnergy ) maxEnergy = energy;
        if ( energy < minClamp ) energy = minClamp;
        else if ( energy > maxClamp ) energy = maxClamp;

        // figure out the basis function radius to use
        if ( m_radiusType != UNIFORM_RADIUS )
        {
            currentRadius = m_radiusMultiplier;
            if ( m_radiusType == ATOM_RADIUS )
                currentRadius *= itr->getRadius();
            else
                currentRadius*= itr->getVanDerWaalsRadius();
            invSqr3 = 3.0 / (currentRadius * currentRadius);
            invCub2 = 2.0 / (currentRadius * currentRadius * currentRadius);
        }
        // compute the position of the next texel
        pVox[0] = corner[0] + inc[0] * ((double)min[0]);
        for ( idx[0] = min[0]; idx[0] < max[0]; ++idx[0] )
        {
            pVox[1] = corner[1] + inc[1] * ((double)min[1]);
            for ( idx[1] = min[1]; idx[1] < max[1]; ++idx[1] )
            {
                pVox[2] = corner[2] + inc[2] * ((double)min[2]);
                tIndex = m_texels.index (idx[0], idx[1], min[2]);
                sIndex = m_samples.index (atomClass, idx[0], idx[1], min[2]);
                for ( idx[2] = min[2];
                      idx[2] < max[2];
                      ++idx[2], ++tIndex, ++sIndex )
                {
                    // calculate the distance from the texel to the atom
                    x = pVox[0] - pAtom[0];
                    y = pVox[1] - pAtom[1];
                    z = pVox[2] - pAtom[2];
                    r = sqrt (x*x + y*y + z*z);
                    if ( r < currentRadius )
                    {
                        // accumulate the atom energy, flag the texel
                        m_texels[tIndex] = 1;
                        double &sample = m_samples[sIndex];
                        sample += energy * (
                            (invCub2 * r - invSqr3) * r * r + 1.0
                        );
                        if ( autoNormalize )
                        {
                            if ( sample < m_minNormal ) m_minNormal = sample;
                            if ( sample > m_maxNormal ) m_maxNormal = sample;
                        }
                    }
                    pVox[2] += inc[2];
                }
                pVox[1] += inc[1];
            }
            pVox[0] += inc[0];
        }
    }
    timer.elapse();
    cout << "clock time    = " << timer.getTime()
         << " seconds to sample atom energies\n"
         << "energy range  = [" << minEnergy << ", " << maxEnergy << "]\n"
         << "normal range  = [" << m_minNormal << ", " << m_maxNormal << "]\n";

    // quantize and load the texels into the renderer
    buildTexture();
    m_transform = Transform();
    m_renderer.setVoxelBlock (m_texels, res, 0, s_voxelAlignment);
    double extent[3] = { aspect[0], aspect[1], aspect[2] };
    m_renderer.setPosition (corner, extent);
    timer.elapse();
    cout << "clock time    = " << timer.getTime()
         << " seconds to quantize texels\n";

    // clear the buffer for the next set of updates
    clearEnergyAccumulator();
}


void EnergyRenderer::setAutoNormalizing (bool autoNormalize)
{
    m_autoNormalize = autoNormalize;
}


void EnergyRenderer::setClassifier (uint classifier)
{
    if ( classifier != m_classifier )
    {
        // assign classifier and reinitialize
        m_classifier = classifier;
        getClassifier();
    }
}


void EnergyRenderer::setColorFunction (uint classification, uint colorFunction)
{
    AtomClassifier *classifier = getClassifier();
    if ( classification < classifier->numClasses() &&
         m_colorNames[classification] != colorFunction )
    {
        m_colorNames[classification] = colorFunction;
        buildTexture();
    }
}


void EnergyRenderer::setGradient (bool gradient)
{
    m_gradient = gradient;
}


void EnergyRenderer::setNormalRange (double minNormal, double maxNormal)
{
    m_minNormal = minNormal;
    m_maxNormal = maxNormal;
}


void EnergyRenderer::setRadiusMultiplier (double radiusMultiplier)
{
    if ( radiusMultiplier > 0.0 )
        m_radiusMultiplier = radiusMultiplier;
}


void EnergyRenderer::setRadiusType (EnergyRenderer::RadiusType type)
{
    m_radiusType = type;
}


void EnergyRenderer::setTexelsPerAngstrom (double resolution)
{
    if ( resolution > 0.0 )
        m_texelsPerAngstrom = resolution;
}


void EnergyRenderer::transform (const Transform &transform)
{
    // premultiply transformation to match behavior of OpenGL
    m_transform.leftMultiply (transform);
}


void EnergyRenderer::updateEnergy()
{
    if ( m_numUpdatesPending < 1.0 && m_numUpdatesPending > 0.0 )
    {
        // undo the reciprocal and averaging if necessary
        m_numUpdatesPending = 1.0 / m_numUpdatesPending;
        for ( uint i = 0; i < m_energyTerms.numElements(); ++i )
            m_energyTerms[i] *= m_numUpdatesPending;
    }
    // grab the latest data from the calculator
    EnergyCalculator *calculator = m_state->energyCalculator;
    if ( !calculator ) return;
    if ( m_gradient )
    {
        // accumulate the gradient magnitudes
        const double *x = calculator->xGradients();
        const double *y = calculator->yGradients();
        const double *z = calculator->zGradients();
        for ( uint i = 0; i < m_energyTerms.numElements(); ++i )
            m_energyTerms[i] += sqrt (x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
    }
    else
    {
        // accumulate the subset sum of the energy terms
        for ( uint i = 0; i < energyLibrary->getNumEnergyComponents(); ++i )
        {
            if ( calculator->getEnergyComponentState(i) )
            {
                const double *energies = calculator->getAtomEnergies (i);
                if ( !energies ) continue;
                for ( uint j = 0; j < m_energyTerms.numElements(); ++j )
                    m_energyTerms[j] += energies[j];
            }
        }
    }
    // increment the number of updates
    m_numUpdatesPending += 1.0;
}


///////////////////////////////////////////////
/// EnergyRenderer private member functions ///
///////////////////////////////////////////////


/** Average the energy terms that have been accumulated since the last sampling
    pass. */
void EnergyRenderer::averageEnergy()
{
    if ( m_numUpdatesPending > 0.0 )
    {
        if ( m_numUpdatesPending > 1.0 )
            m_numUpdatesPending = 1.0 / m_numUpdatesPending;
        for ( uint i = 0; i < m_energyTerms.numElements(); ++i )
            m_energyTerms[i] *= m_numUpdatesPending;
    }
}


/** Build the texel block from the current samples and flag VolumeRenderer to
    reinitialize texture memory.  Do not initialize any texel that is not
    already nonzero. */
void EnergyRenderer::buildTexture()
{
    // load classifier and ensure that it is safe to execute this procedure
    AtomClassifier *classifier = getClassifier();
    if ( classifier->numClasses() != m_samples.dimSize(0) ||
         m_maxNormal <= m_minNormal )
        return;

    // set up some variables for the quantization loop
    Color vColor, tColor;
    const double intervalCoeff = 1.0 / (m_maxNormal - m_minNormal);
    double voxel, weight, alpha;
    uint j, sampleIndex;
    for ( uint i = 0; i < m_texels.numElements(); ++i )
    {
        // quantize the samples into the texel buffer
        if ( m_texels[i] )
        {
            // only initialize texels that were flagged
            // (avoid mapping empty texels to a nonzero value)
            weight = alpha = 0.0;
            tColor[0] = tColor[1] = tColor[2] = tColor[3] = 0.0;
            for ( j = 0, sampleIndex = i;
                  j < classifier->numClasses();
                  ++j, sampleIndex += m_samples.dimStep(0) )
            {
                // normalize all the sample grids into texel colors
                voxel = (m_samples[sampleIndex] - m_minNormal) * intervalCoeff;
                weight += voxel;
                if ( voxel > alpha ) alpha = voxel;
                m_colorPtrs[j]->map (voxel, vColor);
                vColor *= voxel;
                tColor += vColor;
            }
            // apply the weight and alpha coefficients to get the final texel
            // (weight normalizes components to [0,1], alpha is transparency)
            tColor *= alpha / weight;
            m_texels[i] = mapColorToTexel (tColor);

            // One thing that might speed up the quantization loop is to make
            // m_samples column-major instead of row-major.  That way, memory
            // accesses would be localized.  However, memory accesses within the
            // sampling loop (which takes most of the time) would be delocalized
            // even more than they already are.
        }
    }
    m_renderer.notifyVoxelBlockChanged();
}


/** Get a pointer to the atom classifier and ensure that all is in order with
    the relevant member variables. */
AtomClassifier *EnergyRenderer::getClassifier()
{
    if ( m_classifier >= AtomClassifier::numClassifiers() )
        m_classifier = 0;
    AtomClassifier *classifier = AtomClassifier::get (m_classifier);
    assert ( classifier );
    if ( m_colorNames.size() < classifier->numClasses() )
    {
        // ensure that the color function vectors are consistent
        uint add = classifier->numClasses() - m_colorNames.size();
        uint base = ( m_colorNames.size() ) ? m_colorNames.back() : 5;
        for ( uint i = 0; i < add; ++i )
            m_colorNames.push_back (
                ((i + base) % (ColorFunction::numFunctions() - 1)) + 1
            );
    }
    // truncate the vectors if necessary
    else if ( m_colorNames.size() > classifier->numClasses() )
        m_colorNames.resize (classifier->numClasses());
    m_colorPtrs.resize (classifier->numClasses());
    for ( uint i = 0; i < m_colorNames.size(); ++i )
    {
        // ensure that all the function pointers are available
        m_colorPtrs[i] = ColorFunction::get (m_colorNames[i]);
        assert ( m_colorPtrs[i] );
    }
    return classifier;
}


//////////////////////////////////////////////
/// EnergyRenderer private class functions ///
//////////////////////////////////////////////


/** Compute the maximum OpenGL texture dimensions possible in 2D and 3D, and
    ensure that the selected rendering mode is supported by OpenGL.  If it is
    not, change it. */
void EnergyRenderer::calcTexParams()
{
    uint dim2D = 0, dim3D = 0;

    for ( GLint reqDim = 64; reqDim < 8192; reqDim *= 2 )
    {
        GLint p2D = reqDim, p3D = reqDim;
        glTexQuery2D (p2D, p2D);
        glTexQuery3D (p3D, p3D, p3D);
        if ( p2D > dim2D ) dim2D = p2D;
        if ( p3D > dim3D ) dim3D = p3D;
        if ( !p2D && !p3D )
            break;
    }
    cout << "maximum 2D texture dimension = " << dim2D
         << "\nmaximum 3D texture dimension = " << dim3D << endl;
    if ( s_renderingMode == VolumeRenderer::VIEW_PERPENDICULAR )
    {
        if ( !dim3D )
        {
            cout << "3D texture unsupported, switching to 2D texture\n";
            s_renderingMode = VolumeRenderer::AXIS_ALIGNED;
        }
        else
            s_maxTexDim = dim3D;
    }
    if ( s_renderingMode == VolumeRenderer::AXIS_ALIGNED )
        s_maxTexDim = dim2D;
}


/** Ask OpenGL (make a proxy texture query) if it can handle a 2D texture of the
    specified dimensions, and find out what dimensions will actually be used for
    the texture image. */
void EnergyRenderer::glTexQuery2D (GLint &width, GLint &height)
{
    // make a texture query specific to this class
    glTexImage2D (
        GL_PROXY_TEXTURE_2D,
        0,
        GL_RGBA,
        GLint(ceilingPowerOfTwo(width)),
        GLint(ceilingPowerOfTwo(height)),
        0,
        GL_RGBA,
        GL_UNSIGNED_INT_8_8_8_8,
        0
    );
    glGetTexLevelParameteriv (GL_PROXY_TEXTURE_2D, 0,
                              GL_TEXTURE_WIDTH, &width);
    glGetTexLevelParameteriv (GL_PROXY_TEXTURE_2D, 0,
                              GL_TEXTURE_HEIGHT, &height);
}


/** Ask OpenGL (make a proxy texture query) if it can handle a 3D texture of the
    specified dimensions, and find out what dimensions will actually be used for
    the texture image. */
void EnergyRenderer::glTexQuery3D (GLint &width, GLint &height, GLint &depth)
{
    // make a texture query specific to this class
    glTexImage3D (
        GL_PROXY_TEXTURE_3D,
        0,
        GL_COLOR_INDEX8_EXT,
        GLint(ceilingPowerOfTwo(width)),
        GLint(ceilingPowerOfTwo(height)),
        GLint(ceilingPowerOfTwo(depth)),
        0,
        GL_COLOR_INDEX,
        GL_UNSIGNED_BYTE,
        0
    );
    glGetTexLevelParameteriv (GL_PROXY_TEXTURE_3D, 0,
                              GL_TEXTURE_WIDTH, &width);
    glGetTexLevelParameteriv (GL_PROXY_TEXTURE_3D, 0,
                              GL_TEXTURE_HEIGHT, &height);
    glGetTexLevelParameteriv (GL_PROXY_TEXTURE_3D, 0,
                              GL_TEXTURE_DEPTH, &depth);
}


/** Load configuration variables from the configuration file.  This must be done
    at construction time rather than by a static initializer, because the
    configuration file is not yet loaded when static initializers run. */
void EnergyRenderer::lazyLoadConfig()
{
    if ( !s_configLoaded )
    {
        // can't use static init because config file not yet opened
        s_maxBuffer = int(retrieveValue(
            *configFile, "/EnergyCalculator/maxBuffer", 4096
        ));
        s_autoSaveGLState = retrieveValue (
            *configFile, "/EnergyCalculator/autoSaveGLState", true
        );
        s_sliceFactor = retrieveValue (
            *configFile, "/EnergyCalculator/sliceFactor", 0.707
        );
        s_textureCaching = retrieveValue (
            *configFile, "/EnergyCalculator/textureCaching", true
        );
        s_interpolationMode = VolumeRenderer::InterpolationMode (
            int(retrieveValue(
                *configFile, "/EnergyCalculator/interpolationMode", 1
            ))
        );
        // AXIS_ALIGNED = 2D texture, VIEW_PERPENDICULAR = 3D texture
        s_renderingMode = VolumeRenderer::RenderingMode (
            int(retrieveValue(
                *configFile, "/EnergyCalculator/renderingMode", 1
            ))
        );
        s_textureFunction = VolumeRenderer::TextureFunction (
            int(retrieveValue(
                *configFile, "/EnergyCalculator/textureFunction", 0
            ))
        );
        s_voxelAlignment = VolumeRenderer::VoxelAlignment (
            int(retrieveValue(
                *configFile, "/EnergyCalculator/voxelAlignment", 1
            ))
        );
        s_configLoaded = true;
        static const char *enumNames[4][2] = {
            { "CONSTANT", "LINEAR" },
            { "AXIS_ALIGNED", "VIEW_PERPENDICULAR" },
            { "REPLACE", "MODULATE" },
            { "VERTEX_CENTERED", "CELL_CENTERED" }
        };
        cout << "loaded volume renderer configuration:"
             << "\n\tmax buffer size  = " << s_maxBuffer << " MB"
             << "\n\tautosave GL      = " << (s_autoSaveGLState?"true":"false")
             << "\n\tslice factor     = " << s_sliceFactor
             << "\n\ttexture caching  = " << (s_textureCaching?"true":"false")
             << "\n\tinterpolate mode = " << enumNames[0][s_interpolationMode]
             << "\n\trendering mode   = " << enumNames[1][s_renderingMode]
             << "\n\ttexture function = " << enumNames[2][s_textureFunction]
             << "\n\tvoxel alignment  = " << enumNames[3][s_voxelAlignment]
             << endl;
    }
}


/** Convert a floating-point color vector to a formatted texel of the type used
    by class VolumeRenderer. */
EnergyRenderer::Texel EnergyRenderer::mapColorToTexel (const Color &color)
{
    // assume knowledge of texel format: GL_UNSIGNED_INT_8_8_8_8 RGBA
    uint red   = min (uint(255), uint(color[0] * 256.0));
    uint green = min (uint(255), uint(color[1] * 256.0));
    uint blue  = min (uint(255), uint(color[2] * 256.0));
    uint alpha = min (uint(255), uint(color[3] * 256.0));
    return Texel ((red << 24) | (green << 16) | (blue << 8) | (alpha));
}


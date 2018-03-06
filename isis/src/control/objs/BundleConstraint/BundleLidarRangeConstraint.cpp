#include "BundleLidarRangeConstraint.h"

// Qt Library
#include <QDebug>

// Isis Library
#include "BundleLidarControlPoint.h"
#include "LidarControlPoint.h"
#include "SpicePosition.h"
#include "SpiceRotation.h"

// boost lib
//#include <boost/lexical_cast.hpp>
//#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

using namespace boost::numeric::ublas;

namespace Isis {

  /**
   * Default constructor
   *
   */
  BundleLidarRangeConstraint::BundleLidarRangeConstraint() {
  }


  /**
   * Constructor
   *
   * @param parentObservation parent BundleObservation
   */
  BundleLidarRangeConstraint::
      BundleLidarRangeConstraint(BundleLidarControlPointQsp parentBundleLidarControlPoint) {

    m_parentBundleLidarControlPoint = parentBundleLidarControlPoint;

    // initialize variables
/*
    int nSpkParameters = 0;
    int nCkParameters = 0;

    if (solveSettings->instrumentPositionSolveOption() !=
        BundleObservationSolveSettings::NoPositionFactors) {

      m_numberSpkCoefficients = solveSettings->numberCameraPositionCoefficientsSolved();
      m_numberSpkSegments = solveSettings->numberSpkPolySegments();
      m_numberSpkBoundaries = m_numberSpkSegments - 1;

      nSpkParameters = 3 * m_numberSpkCoefficients;

      // knots contain scaled time
      m_spkKnots = parentObservation->spicePosition()->scaledPolynomialKnots();

      // remove end knots, leaving only knots at segment boundaries
      if (m_spkKnots.size() > 2) {
        m_spkKnots.pop_back();
        m_spkKnots.erase(m_spkKnots.begin());
      }
    }

    if (solveSettings->instrumentPointingSolveOption() !=
        BundleObservationSolveSettings::NoPointingFactors) {

      m_numberCkCoefficients = solveSettings->numberCameraAngleCoefficientsSolved();
      m_numberCkSegments = solveSettings->numberCkPolySegments();
      m_numberCkBoundaries = m_numberCkSegments - 1;

      nCkParameters = 2 * m_numberCkCoefficients;
      if (solveSettings->solveTwist()) {
        nCkParameters += m_numberCkCoefficients;
      }

      // knots contain scaled time
      m_ckKnots = parentObservation->spiceRotation()->scaledPolynomialKnots();

      // remove end knots, leaving only knots at segment boundaries
      if (m_ckKnots.size() > 2) {
        m_ckKnots.pop_back();
        m_ckKnots.erase(m_ckKnots.begin());
      }
    }

    if (nSpkParameters > 0) {
      m_numberConstraintEquations = m_numberSpkBoundaries * (m_numberSpkCoefficients-1) * 3.0;
    }
    if (nCkParameters > 0) {
      if (m_parentObservation->solveSettings()->solveTwist()) {
        m_numberConstraintEquations += m_numberCkBoundaries * (m_numberCkCoefficients-1) * 3.0;
      }
      else {
        m_numberConstraintEquations += m_numberCkBoundaries * (m_numberCkCoefficients-1) * 2.0;
      }
    }

    m_numberParameters = m_numberSpkSegments * nSpkParameters + m_numberCkSegments * nCkParameters;

    m_numberSpkSegmentParameters = nSpkParameters;
    m_numberCkSegmentParameters = nCkParameters;

    constructMatrices();
*/
  }


  /**
   * Destructor
   */
  BundleLidarRangeConstraint::~BundleLidarRangeConstraint() {
  }


  /**
   * Copy constructor
   *
   * Constructs a BundleLidarRangeConstraint from another.
   *
   * @param src Source BundleLidarRangeConstraint to copy.
   */
  BundleLidarRangeConstraint:: BundleLidarRangeConstraint(const BundleLidarRangeConstraint &src) {
    m_parentBundleLidarControlPoint = src.m_parentBundleLidarControlPoint;
/*
    m_numberSegmentParameters = src.m_numberSegmentParameters;
    m_numberParameters = src.m_numberParameters;
    m_numberConstraintEquations = src.m_numberConstraintEquations;

    m_ckKnots = src.m_ckKnots;
    m_spkKnots = src.m_spkKnots;
    m_numberCkCoefficients = src.m_numberCkCoefficients;
    m_numberSpkCoefficients = src.m_numberSpkCoefficients;
    m_numberCkSegments = src.m_numberCkSegments;
    m_numberSpkSegments = src.m_numberSpkSegments;
    m_numberCkBoundaries = src.m_numberCkBoundaries;
    m_numberSpkBoundaries = src.m_numberSpkBoundaries;
    m_numberCkSegmentParameters = src.m_numberCkSegmentParameters;
    m_numberSpkSegmentParameters = src.m_numberSpkSegmentParameters;

    m_designMatrix = src.m_designMatrix;
    m_normalsSpkMatrix = src.m_normalsSpkMatrix;
    m_normalsCkMatrix = src.m_normalsCkMatrix;
    m_rightHandSide = src.m_rightHandSide;
    m_omcVector = src.m_omcVector;
*/
  }


  /**
   * Assignment operator
   *
   * Assigns state of this BundleLidarRangeConstraint from another.
   *
   * @param src Source BundleLidarRangeConstraint to assign state from.
   *
   * @return    BundleLidarRangeConstraint& Returns reference to this
   *            BundleLidarRangeConstraint.
   */
  BundleLidarRangeConstraint &BundleLidarRangeConstraint::
      operator=(const BundleLidarRangeConstraint &src) {

    // Prevent self assignment
    if (this != &src) {
      m_parentBundleLidarControlPoint = src.m_parentBundleLidarControlPoint;
/*
      m_numberSegmentParameters = src.m_numberSegmentParameters;
      m_numberParameters = src.m_numberParameters;
      m_numberConstraintEquations = src.m_numberConstraintEquations;

      m_ckKnots = src.m_ckKnots;
      m_spkKnots = src.m_spkKnots;
      m_numberCkCoefficients = src.m_numberCkCoefficients;
      m_numberSpkCoefficients = src.m_numberSpkCoefficients;
      m_numberCkSegments = src.m_numberCkSegments;
      m_numberSpkSegments = src.m_numberSpkSegments;
      m_numberCkBoundaries = src.m_numberCkBoundaries;
      m_numberSpkBoundaries = src.m_numberSpkBoundaries;
      m_numberCkSegmentParameters = src.m_numberCkSegmentParameters;
      m_numberSpkSegmentParameters = src.m_numberSpkSegmentParameters;

      m_designMatrix = src.m_designMatrix;
      m_normalsSpkMatrix = src.m_normalsSpkMatrix;
      m_normalsCkMatrix = src.m_normalsCkMatrix;
      m_rightHandSide = src.m_rightHandSide;
      m_omcVector = src.m_omcVector;
*/
    }

    return *this;
  }


  /**
   *
   * Compute range constraint between lidar ground point and simultaneous image
   *
   * @author 2018-03-05 Ken Edmundson
   *
   * @param src Source BundleLidarRangeConstraint to assign state from.
   */
  bool BundleLidarRangeConstraint::formRangeConstraint() {

//    bool BundleLidarRangeConstraint::formRangeConstraint(symmetric_matrix<double, upper>&N22,
//                                                         SSparseBlockColumnMatrix& N12,
//                                                         compressed_vector<double>& n1,
//                                                         vector<double>& n2,
//                                                         Camera* camera,
//                                                         const ControlPoint* point,
//                                                         int nImageIndex,
//                                                         double range,
//                                                         double rangeSigma) {
/*
    int i;

    LinearAlgebra::Matrix coeff_range_image(1,m_nNumImagePartials);
    LinearAlgebra::Matrix coeff_range_point3D(1,3);
    LinearAlgebra::Vector coeff_range_RHS(1);

    coeff_range_image.clear();
    coeff_range_point3D.clear();
    coeff_range_RHS.clear();

    std::cout << "image" << std::endl << coeff_range_image << std::endl;
    std::cout << "point" << std::endl << coeff_range_point3D << std::endl;
    std::cout << "rhs" << std::endl << coeff_range_RHS << std::endl;

    // compute partial derivatives for camstation-to-range point condition

    // get ground point in body-fixed coordinates
    double xPoint
        = m_parentBundleLidarControlPoint->rawControlPoint()->GetAdjustedSurfacePoint().GetX().kilometers();
    double yPoint
        = m_parentBundleLidarControlPoint->rawControlPoint()->GetAdjustedSurfacePoint().GetY().kilometers();
    double zPoint
        = m_parentBundleLidarControlPoint->rawControlPoint()->GetAdjustedSurfacePoint().GetZ().kilometers();

    // get spacecraft position in J2000 coordinates
    LinearAlgebra::Vector CameraJ2KXYZ(3);
    CameraJ2KXYZ = camera->instrumentPosition()->Coordinate();
    double xCameraJ2K  = CameraJ2KXYZ[0];
    double yCameraJ2K  = CameraJ2KXYZ[1];
    double zCameraJ2K  = CameraJ2KXYZ[2];

    // get spacecraft position in body-fixed coordinates
    LinearAlgebra::Vector CameraBodyFixedXYZ(3);

    // "InstrumentPosition()->Coordinate()" returns the instrument coordinate in J2000;
    // then the body rotation "ReferenceVector" rotates that into body-fixed coordinates
    // TODO: can we do something so we don't have to do this computation every time? Can the
    //       SpicePosition class store and update (when needed) the body-fixed XYZ coordinates as
    //       well as the J2K coordinates? In other words, just as the camera for the body-fixed or
    //       J2K coordinates directly
    SpiceRotation* bodyRot = camera->bodyRotation();
    CameraBodyFixedXYZ = bodyRot->ReferenceVector(camera->instrumentPosition()->Coordinate());
    double xCamera  = CameraBodyFixedXYZ[0];
    double yCamera  = CameraBodyFixedXYZ[1];
    double zCamera  = CameraBodyFixedXYZ[2];

    // computed distance between spacecraft and point
    double dX = xCamera - xPoint;
    double dY = yCamera - yPoint;
    double dZ = zCamera - zPoint;
    double computed_distance = sqrt(dX*dX+dY*dY+dZ*dZ);

    // observed distance - computed distance
    double observed_computed = range - computed_distance;

    // get matrix that rotates spacecraft from J2000 to body-fixed
    // TODO: can we get this directly as a rotation matrix, not a vector?
    std::vector<double> matrix_Target_to_J2K;
    matrix_Target_to_J2K.resize(9);
    matrix_Target_to_J2K = bodyRot->Matrix();

    double m11 = matrix_Target_to_J2K[0];
    double m12 = matrix_Target_to_J2K[1];
    double m13 = matrix_Target_to_J2K[2];
    double m21 = matrix_Target_to_J2K[3];
    double m22 = matrix_Target_to_J2K[4];
    double m23 = matrix_Target_to_J2K[5];
    double m31 = matrix_Target_to_J2K[6];
    double m32 = matrix_Target_to_J2K[7];
    double m33 = matrix_Target_to_J2K[8];

    // partials w/r to image auxiliaries
    double a1 = m11*xCameraJ2K + m12*yCameraJ2K + m13*zCameraJ2K - xPoint;
    double a2 = m21*xCameraJ2K + m22*yCameraJ2K + m23*zCameraJ2K - yPoint;
    double a3 = m31*xCameraJ2K + m32*yCameraJ2K + m33*zCameraJ2K - zPoint;

    coeff_range_image(0,0) = (m11*a1 + m21*a2 + m31*a3)/computed_distance;
    coeff_range_image(0,1) = (m12*a1 + m22*a2 + m32*a3)/computed_distance;
    coeff_range_image(0,2) = (m13*a1 + m23*a2 + m33*a3)/computed_distance;

    // partials w/r to point
    // TODO: what if we are bundling points in body-fixed XYZ instead of lat, lon, radius?
    //       I think it dramatically simplifies the partial derivatives.
    double lat    = point->GetAdjustedSurfacePoint().GetLatitude().radians();
    double lon    = point->GetAdjustedSurfacePoint().GetLongitude().radians();
    double radius = point->GetAdjustedSurfacePoint().GetLocalRadius().kilometers();

    double sinlat = sin(lat);
    double coslat = cos(lat);
    double sinlon = sin(lon);
    double coslon = cos(lon);

    coeff_range_point3D(0,0) =  radius*(sinlat*coslon*a1 + sinlat*sinlon*a2 - coslat*a3)/computed_distance;
    coeff_range_point3D(0,1) =  radius*(coslat*sinlon*a1 - coslat*coslon*a2)/computed_distance;
    coeff_range_point3D(0,2) =  -(coslat*coslon*a1 + coslat*sinlon*a2 + sinlat*a3)/computed_distance;

    // right hand side
    coeff_range_RHS(0) = observed_computed;

    // multiply coefficients by observation weight
    double dObservationWeight = 1.0/(rangeSigma*0.001); // converting sigma from meters to km
    coeff_range_image   *= dObservationWeight;
    coeff_range_point3D *= dObservationWeight;
    coeff_range_RHS     *= dObservationWeight;

    // form matrices to be added to normal equation auxiliaries
    static vector<double> n1_image(m_nNumImagePartials);
    n1_image.clear();

    // form N11 for the condition partials for image
    static symmetric_matrix<double, upper> N11(m_nNumImagePartials);
    N11.clear();

//  std::cout << "N11" << std::endl << N11 << std::endl;

    std::cout << "image" << std::endl << coeff_range_image << std::endl;
    std::cout << "point" << std::endl << coeff_range_point3D << std::endl;
    std::cout << "rhs" << std::endl << coeff_range_RHS << std::endl;

    N11 = prod(trans(coeff_range_image), coeff_range_image);

    std::cout << "N11" << std::endl << N11 << std::endl;

    int t = m_nNumImagePartials * nImageIndex;

    // insert submatrix at column, row
    m_SparseNormals.InsertMatrixBlock(nImageIndex, nImageIndex, m_nNumImagePartials,
                                      m_nNumImagePartials);

    (*(*m_SparseNormals[nImageIndex])[nImageIndex]) += N11;

    std::cout << (*(*m_SparseNormals[nImageIndex])[nImageIndex]) << std::endl;

    // form N12_Image
    static matrix<double> N12_Image(m_nNumImagePartials, 3);
    N12_Image.clear();

    N12_Image = prod(trans(coeff_range_image), coeff_range_point3D);

    std::cout << "N12_Image" << std::endl << N12_Image << std::endl;

    // insert N12_Image into N12
//  for (i = 0; i < m_nNumImagePartials; i++)
//    for (j = 0; j < 3; j++)
//      N12(i + t, j) += N12_Image(i, j);
    N12.InsertMatrixBlock(nImageIndex, m_nNumImagePartials, 3);
    *N12[nImageIndex] += N12_Image;

//  printf("N12\n");
//  std::cout << N12 << std::endl;

    // form n1
    n1_image = prod(trans(coeff_range_image), coeff_range_RHS);

  std::cout << "n1_image" << std::endl << n1_image << std::endl;

    // insert n1_image into n1
    for (i = 0; i < m_nNumImagePartials; i++)
      n1(i + t) += n1_image(i);

    // form N22
    N22 += prod(trans(coeff_range_point3D), coeff_range_point3D);

  std::cout << "N22" << std::endl << N22 << std::endl;

  std::cout << "n2" << std::endl << n2 << std::endl;

    // form n2
    n2 += prod(trans(coeff_range_point3D), coeff_range_RHS);

  std::cout << "n2" << std::endl << n2 << std::endl;
*/
    return true;
  }


  /**
   * Returns matrix with contribution to position portion of bundle adjustment normal equations from
   * continuity constraints.
   *
   * @return SparseBlockMatrix Matrix with contribution to position portion of bundle adjustment
   *                           normal equations from continuity constraints.
   */
  SparseBlockMatrix &BundleLidarRangeConstraint::normalsSpkMatrix() {
    return m_normalsSpkMatrix;
  }


  /**
   * Returns vector with contribution to bundle adjustment normal equations right hand side from
   * continuity constraints.
   *
   * @return LinearAlgebra::Vector Vector with contribution to bundle adjustment normal equations
   *                               right hand side from continuity constraints.
   */
  LinearAlgebra::Vector
      &BundleLidarRangeConstraint::rightHandSideVector() {
    return m_rightHandSide;
  }


  /**
   * Constructs m_normalsSpkMatrix and m_normalsCkMatrix, m_rightHandSide vector, m_designMatrix,
   * and m_omcVector (or observed - computed vector).
   *
   * @todo need more documentation on technical aspects of approach.
   */
  void BundleLidarRangeConstraint::constructMatrices() {
/*
    if (m_numberConstraintEquations <= 0 )
      return;

    // initialize size of right hand side vector
    // the values in this vector are updated each iteration, but vector size will not change
    m_rightHandSide.resize(m_numberParameters);
    m_rightHandSide.clear();

    // initialize size of design matrix
    // this will not change throughout the bundle adjustment
    m_designMatrix.resize(m_numberConstraintEquations, m_numberParameters, false);
    m_designMatrix.clear();

    // initialize size of "observed - computed" vector
    m_omcVector.resize(m_numberConstraintEquations);

    int designRow=0;

    // spk (position) contribution
    if (m_numberSpkSegments > 1 && m_numberConstraintEquations > 0 && m_numberSpkCoefficients > 1)
      positionContinuity(designRow);

    // ck (pointing) contribution
    if (m_numberCkSegments > 1 && m_numberConstraintEquations > 0 && m_numberCkCoefficients > 1)
      pointingContinuity(designRow);

    int numPositionParameters = m_parentObservation->numberPositionParametersPerSegment();
    int numPointingParameters = m_parentObservation->numberPointingParametersPerSegment();

    // initialize and fill position blocks in m_normalsSpkMatrix
    if (m_numberSpkSegments > 1 && m_numberConstraintEquations > 0 && m_numberSpkCoefficients > 1) {
      m_normalsSpkMatrix.setNumberOfColumns(m_numberSpkSegments);

      for (int i = 0; i < m_numberSpkSegments; i++) {
        m_normalsSpkMatrix.insertMatrixBlock(i, i, numPositionParameters, numPositionParameters);
        LinearAlgebra::Matrix *block = m_normalsSpkMatrix.getBlock(i, i);

        matrix_range<LinearAlgebra::MatrixCompressed>
            mr1 (m_designMatrix, range (0, m_designMatrix.size1()),
                range (i*numPositionParameters, (i+1)*numPositionParameters));

        *block += prod(trans(mr1), mr1);

        if (i > 0) {
          m_normalsSpkMatrix.insertMatrixBlock(i, i-1, numPositionParameters, numPositionParameters);
          block = m_normalsSpkMatrix.getBlock(i, i-1);

          matrix_range<LinearAlgebra::MatrixCompressed>
              mr2 (m_designMatrix, range (0, m_designMatrix.size1()),
                  range ((i-1)*numPositionParameters, i*numPositionParameters));

          *block += prod(trans(mr2), mr1);
        }
      }
    }

    // initialize and fill pointing blocks
    if (m_numberCkSegments > 1 && m_numberConstraintEquations > 0 && m_numberCkCoefficients > 1) {
      m_normalsCkMatrix.setNumberOfColumns(m_numberCkSegments);

      int t = numPositionParameters * m_numberSpkSegments;

      for (int i = 0; i < m_numberCkSegments; i++) {
        m_normalsCkMatrix.insertMatrixBlock(i, i, numPointingParameters, numPointingParameters);
        LinearAlgebra::Matrix *block = m_normalsCkMatrix.getBlock(i, i);

        matrix_range<LinearAlgebra::MatrixCompressed>
            mr1 (m_designMatrix, range (0, m_designMatrix.size1()),
                range (t+i*numPointingParameters, t+(i+1)*+numPointingParameters));

        *block += prod(trans(mr1), mr1);

        if (i > 0) {
          m_normalsCkMatrix.insertMatrixBlock(i, i-1, numPointingParameters, numPointingParameters);
          block = m_normalsCkMatrix.getBlock(i, i-1);

          matrix_range<LinearAlgebra::MatrixCompressed>
              mr2 (m_designMatrix, range (0, m_designMatrix.size1()),
                  range (t+(i-1)*numPointingParameters, t+i*numPointingParameters));

          *block += prod(trans(mr2), mr1);
        }
      }
    }

    // initialize right hand side vector
    updateRightHandSide();
*/
  }


  /**
   * Updates right hand side vector after parameters have been updated at each iteration.
   *
   * @todo need more documentation on technical aspects of approach.
   */
  void BundleLidarRangeConstraint::updateRightHandSide() {
/*
    SpicePosition *position = m_parentObservation->spicePosition();
    SpiceRotation *rotation = m_parentObservation->spiceRotation();

    int designRow = 0;

    bool solveTwist = m_parentObservation->solveSettings()->solveTwist();

    // clear "observed - computed" (omc) and "right-hand side" vectors
//  LinearAlgebra::Vector omcVector(m_numberConstraintEquations);
    m_omcVector.clear();

    m_rightHandSide.clear();

    double t, t2;

    if (m_numberSpkSegments > 1) {
      std::vector<double> seg1CoefX(m_numberSpkCoefficients);
      std::vector<double> seg1CoefY(m_numberSpkCoefficients);
      std::vector<double> seg1CoefZ(m_numberSpkCoefficients);
      std::vector<double> seg2CoefX(m_numberSpkCoefficients);
      std::vector<double> seg2CoefY(m_numberSpkCoefficients);
      std::vector<double> seg2CoefZ(m_numberSpkCoefficients);

      // loop over segment boundaries for zero order in X,Y,Z
      for (int i = 0; i < m_numberSpkBoundaries; i++) {

        position->GetPolynomial(seg1CoefX, seg1CoefY, seg1CoefZ, i);
        position->GetPolynomial(seg2CoefX, seg2CoefY, seg2CoefZ, i+1);

        // TODO: if 1st order, should we be using 3 coefficients here?

        // time (t) and time-squared (t2)
        t = m_spkKnots[i];
        t2 = t*t;

        // 0 order in X
        m_omcVector(designRow) = -seg1CoefX[0] - seg1CoefX[1]*t + seg2CoefX[0] + seg2CoefX[1]*t;
        if (m_numberSpkCoefficients == 3) {
          m_omcVector(designRow) = m_omcVector(designRow) - seg1CoefX[2]*t2 + seg2CoefX[2]*t2;
        }
        designRow++;
        // 0 order in Y
        m_omcVector(designRow) = -seg1CoefY[0] - seg1CoefY[1]*t + seg2CoefY[0] + seg2CoefY[1]*t;
        if (m_numberSpkCoefficients == 3) {
          m_omcVector(designRow) = m_omcVector(designRow) - seg1CoefY[2]*t2 + seg2CoefY[2]*t2;
        }
        designRow++;
        // 0 order in Z
        m_omcVector(designRow) = -seg1CoefZ[0] - seg1CoefZ[1]*t + seg2CoefZ[0] + seg2CoefZ[1]*t;
        if (m_numberSpkCoefficients == 3) {
          m_omcVector(designRow) = m_omcVector(designRow) - seg1CoefZ[2]*t2 + seg2CoefZ[2]*t2;
        }
        designRow++;
      }

      if (m_numberSpkCoefficients == 3) {

        // loop over segment boundaries for 1st order in X,Y,Z
        for (int i = 0; i < m_numberSpkBoundaries; i++) {

          position->GetPolynomial(seg1CoefX, seg1CoefY, seg1CoefZ, i);
          position->GetPolynomial(seg2CoefX, seg2CoefY, seg2CoefZ, i+1);

          // time (t)
          t = m_spkKnots[i];

          // 1st order in X
          m_omcVector(designRow) = -seg1CoefX[1] + seg2CoefX[1];
          if (m_numberSpkCoefficients == 3) {
            m_omcVector(designRow) = m_omcVector(designRow) - 2.0*seg1CoefX[2]*t + 2.0*seg2CoefX[2]*t;
          }
          designRow++;
          // 1st order in Y
          m_omcVector(designRow) = -seg1CoefY[1] + seg2CoefY[1];
          if (m_numberSpkCoefficients == 3) {
            m_omcVector(designRow) = m_omcVector(designRow) - 2.0*seg1CoefY[2]*t + 2.0*seg2CoefY[2]*t;
          }
          designRow++;
          // 1st order in Z
          m_omcVector(designRow) = -seg1CoefZ[1] + seg2CoefZ[1];
          if (m_numberSpkCoefficients == 3) {
            m_omcVector(designRow) = m_omcVector(designRow) - 2.0*seg1CoefZ[2]*t + 2.0*seg2CoefZ[2]*t;
          }
          designRow++;
        }
      }
    }

    if (m_numberCkSegments > 1) {
      std::vector<double> seg1CoefRA(m_numberCkCoefficients);
      std::vector<double> seg1CoefDEC(m_numberCkCoefficients);
      std::vector<double> seg1CoefTWIST(m_numberCkCoefficients);
      std::vector<double> seg2CoefRA(m_numberCkCoefficients);
      std::vector<double> seg2CoefDEC(m_numberCkCoefficients);
      std::vector<double> seg2CoefTWIST(m_numberCkCoefficients);

      // loop over segment boundaries for zero order in RA,DEC,TWIST
      for (int i = 0; i < m_numberCkBoundaries; i++) {

        rotation->GetPolynomial(seg1CoefRA, seg1CoefDEC, seg1CoefTWIST, i);
        rotation->GetPolynomial(seg2CoefRA, seg2CoefDEC, seg2CoefTWIST, i+1);

        // TODO: if 1st order, should we be using 3 coefficients here?

        // time (t) and time-squared (t2)
        t = m_ckKnots[i];
        t2 = t*t;

        // 0 order in RA
        m_omcVector(designRow) = -seg1CoefRA[0] - seg1CoefRA[1]*t + seg2CoefRA[0] + seg2CoefRA[1]*t;
        if (m_numberCkCoefficients == 3) {
          m_omcVector(designRow) = m_omcVector(designRow) - seg1CoefRA[2]*t2 + seg2CoefRA[2]*t2;
        }
        designRow++;
        // 0 order in DEC
        m_omcVector(designRow)
            = -seg1CoefDEC[0] - seg1CoefDEC[1]*t + seg2CoefDEC[0] + seg2CoefDEC[1]*t;
        if (m_numberCkCoefficients == 3) {
          m_omcVector(designRow) = m_omcVector(designRow) - seg1CoefDEC[2]*t2 + seg2CoefDEC[2]*t2;
        }
        designRow++;
        if (solveTwist) {
          // 0 order in TWIST
          m_omcVector(designRow)
              = -seg1CoefTWIST[0] - seg1CoefTWIST[1]*t + seg2CoefTWIST[0] + seg2CoefTWIST[1]*t;
          if (m_numberCkCoefficients == 3) {
            m_omcVector(designRow) = m_omcVector(designRow) - seg1CoefTWIST[2]*t2 + seg2CoefTWIST[2]*t2;
          }
          designRow++;
        }
      }

      if (m_numberCkCoefficients == 3) {

      // loop over segment boundaries for 1st order in RA,DEC,TWIST
        for (int i = 0; i < m_numberCkBoundaries; i++) {

          rotation->GetPolynomial(seg1CoefRA, seg1CoefDEC, seg1CoefTWIST, i);
          rotation->GetPolynomial(seg2CoefRA, seg2CoefDEC, seg2CoefTWIST, i+1);

          // time (t)
          t = m_ckKnots[i];

          // 1st order in RA
          m_omcVector(designRow) = -seg1CoefRA[1] + seg2CoefRA[1];
          if (m_numberCkCoefficients == 3) {
            m_omcVector(designRow) = m_omcVector(designRow) - 2.0*seg1CoefRA[2]*t + 2.0*seg2CoefRA[2]*t;
          }
          designRow++;
          // 1st order in DEC
          m_omcVector(designRow) = -seg1CoefDEC[1] + seg2CoefDEC[1];
          if (m_numberCkCoefficients == 3) {
            m_omcVector(designRow) = m_omcVector(designRow) - 2.0*seg1CoefDEC[2]*t + 2.0*seg2CoefDEC[2]*t;
          }
          designRow++;
          if (solveTwist) {
            // 1st order in TWIST
            m_omcVector(designRow) = -seg1CoefTWIST[1] + seg2CoefTWIST[1];
            if (m_numberCkCoefficients == 3) {
              m_omcVector(designRow)
                  = m_omcVector(designRow) - 2.0*seg1CoefTWIST[2]*t + 2.0*seg2CoefTWIST[2]*t;
            }
            designRow++;
          }
        }
      }
    }

    // NOTE: 1.0e+5 is the square root of the weight applied to the constraint equations
    // the design matrix has been premultipled by the square root of the weight
    // we aren't premultiplying the omcVector so we have access to its raw values to output
    // to the bundleout.txt file as the deltas between 0 and 1st order functions at segment
    // boundaries
    m_rightHandSide = prod(trans(m_designMatrix),m_omcVector*1.0e+5);
*/
  }

  /**
   * Creates and returns formatted QString summarizing continuity constraints for output to
     bundleout.txt file.
   *
   * @return QString Returns formatted QString summarizing continuity constraints for output to
   *                 bundleout.txt file.
   */
  QString BundleLidarRangeConstraint::formatBundleOutputString() {
    QString finalqStr = "";
    QString qStr = "";

/*

    int index = 0;

    if (m_numberSpkBoundaries > 0) {
      qStr = QString("\nContinuity Constraints\n======================\n\n"
                     "Position Segments/Boundaries: %1/%2\n"
                     "         0-order Constraints: %3\n").
                     arg(m_numberSpkSegments).
                     arg(m_numberSpkBoundaries).
                     arg(3*m_numberSpkBoundaries);

      // loop over segment boundaries for 0-order constraints
      for (int i = 0; i < m_numberSpkBoundaries; i++) {
        qStr += QString("            Bndry %1 dX/dY/dZ: %2/%3/%4\n").
                       arg(i+1).
                       arg(m_omcVector(index), 5, 'e', 1).
                       arg(m_omcVector(index+1), 5, 'e', 1).
                       arg(m_omcVector(index+2), 5, 'e', 1);
        index += 3;
      }

      if (m_numberSpkCoefficients > 2) {
        qStr += QString("       1st-order Constraints: %1\n").
                        arg(3*m_numberSpkBoundaries);
        // loop over segment boundaries for 1st-order constraints
        for (int i = 0; i < m_numberSpkBoundaries; i++) {
          qStr += QString("            Bndry %2 dX/dY/dZ: %3/%4/%5\n").
                         arg(i+1).
                         arg(m_omcVector(index), 5, 'e', 1).
                         arg(m_omcVector(index+1), 5, 'e', 1).
                         arg(m_omcVector(index+2), 5, 'e', 1);
          index += 3;
        }
      }
    }

    finalqStr += qStr;

    if (m_numberCkBoundaries > 0) {
      if (!m_parentObservation->solveSettings()->solveTwist()) {
        qStr = QString("\nPointing Segments/Boundaries: %1/%2\n"
                       "         0-order Constraints: %3\n").
                       arg(m_numberCkSegments).
                       arg(m_numberCkBoundaries).
                       arg(2*m_numberCkBoundaries);

        // loop over segment boundaries for 0-order constraints
        for (int i = 0; i < m_numberCkBoundaries; i++) {

          qStr += QString("            Bndry %1 dRa/dDec: %2/%3\n").
                         arg(i+1).
                         arg(m_omcVector(index), 5, 'e', 1).
                         arg(m_omcVector(index+1), 5, 'e', 1);
          index += 2;
        }

        if (m_numberCkCoefficients > 2) {
          qStr += QString("       1st-order Constraints: %1\n").
                          arg(2*m_numberCkBoundaries);
          // loop over segment boundaries for 1st-order constraints
          for (int i = 0; i < m_numberCkBoundaries; i++) {
            qStr += QString("            Bndry %2 dRa/dDec: %3/%4\n").
                           arg(i+1).
                           arg(m_omcVector(index), 5, 'e', 1).
                           arg(m_omcVector(index+1), 5, 'e', 1);
            index += 2;
          }
        }
      }
      else {
        qStr = QString("\nPointing Segments/Boundaries: %1/%2\n"
                       "         0-order Constraints: %3\n").
                       arg(m_numberCkSegments).
                       arg(m_numberCkBoundaries).
                       arg(3*m_numberCkBoundaries);

        // loop over segment boundaries for 0-order constraints
        for (int i = 0; i < m_numberCkBoundaries; i++) {

          qStr += QString("       Bndry %1 dRa/dDec/dTwi: %2/%3/%4\n").
                         arg(i+1).
                         arg(m_omcVector(index), 5, 'e', 1).
                         arg(m_omcVector(index+1), 5, 'e', 1).
                         arg(m_omcVector(index+2), 5, 'e', 1);
          index += 3;
        }

        if (m_numberCkCoefficients > 2) {
          qStr += QString("       1st-order Constraints: %1\n").
                          arg(3*m_numberCkBoundaries);
          // loop over segment boundaries for 1st-order constraints
          for (int i = 0; i < m_numberCkBoundaries; i++) {
            qStr += QString("       Bndry %2 dRa/dDec/dTwi: %3/%4/%5\n").
                           arg(i+1).
                           arg(m_omcVector(index), 5, 'e', 1).
                           arg(m_omcVector(index+1), 5, 'e', 1).
                           arg(m_omcVector(index+2), 5, 'e', 1);
            index += 3;
          }
        }
      }
    }
*/
    finalqStr += qStr;

    return finalqStr;
  }
}







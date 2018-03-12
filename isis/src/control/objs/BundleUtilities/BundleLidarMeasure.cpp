#include "BundleLidarMeasure.h"

// Isis lib
#include "BundleLidarControlPoint.h"
#include "BundleMeasure.h"
#include "BundleObservation.h"
#include "BundleObservationSolveSettings.h"
#include "ControlMeasure.h"
#include "IException.h"

namespace Isis {

  /**
   * Constructor
   *
   * Constructs a BundleLidarMeasure from a ControlMeasure with the passed BundleControlPoint as its
   * parent control point 
   *
   * @param controlMeasure Pointer to the ControlMeasure to store
   * @param bundleControlPoint Pointer to the BundleControlPoint that contains this
   *                           BundleLidarMeasure
   */
  BundleLidarMeasure::BundleLidarMeasure(ControlMeasure *controlMeasure,
                                         BundleLidarControlPoint *point)
    : BundleMeasure(controlMeasure, (BundleControlPoint*) point) {

//    m_controlMeasure = controlMeasure;
//    m_parentControlPoint = bundleControlPoint;
//    m_polyPositionSegmentIndex = 0;
//    m_polyPointingSegmentIndex = 0;
//    m_normalsPositionBlockIndex = -1;
//    m_normalsPointingBlockIndex = -1;
  }


  /**
   * Destructor
   */
  BundleLidarMeasure::~BundleLidarMeasure() {
  }


  /**
   * Copy constructor
   *
   * Constructs a BundleLidarMeasure from another BundleLidarMeasure
   *
   * @param src The source BundleLidarMeasure to copy
   */
//  BundleLidarMeasure::BundleLidarMeasure(const BundleLidarMeasure &src) {
//    m_controlMeasure = src.m_controlMeasure;
//    m_parentControlPoint = src.m_parentControlPoint;
//    m_parentBundleImage = src.m_parentBundleImage;
//    m_parentObservation = src.m_parentObservation;
//    m_polyPositionSegmentIndex = src.m_polyPositionSegmentIndex;
//    m_polyPointingSegmentIndex = src.m_polyPointingSegmentIndex;
//    m_normalsPositionBlockIndex = src.m_normalsPositionBlockIndex;
//    m_normalsPointingBlockIndex = src.m_normalsPointingBlockIndex;
  //}


  /**
   * Assignment operator
   *
   * Assigns the state of this BundleLidarMeasure from another BundleLidarMeasure
   *
   * @param src The source BundleLidarMeasure to assign state from
   *
   * @return BundleLidarMeasure& Returns a reference to this BundleLidarMeasure
   */
  BundleLidarMeasure &BundleLidarMeasure::operator=(const BundleLidarMeasure &src) {
    // Prevent self assignment
    if (this != &src) {
//      m_controlMeasure = src.m_controlMeasure;
//      m_parentControlPoint = src.m_parentControlPoint;
//      m_parentBundleImage = src.m_parentBundleImage;
//      m_parentObservation = src.m_parentObservation;
//      m_polyPositionSegmentIndex = src.m_polyPositionSegmentIndex;
//      m_polyPointingSegmentIndex = src.m_polyPointingSegmentIndex;
//      m_normalsPositionBlockIndex = src.m_normalsPositionBlockIndex;
//      m_normalsPointingBlockIndex = src.m_normalsPointingBlockIndex;
    }

    return *this;
  }

  /**
   * applies lidar range constraint (not sure if this actually makes sense)
   *
   * @param N22 Contribution (3x3) to normal equations matrix for a control point.
   * @param N12 Contribution to normal equations matrix connecting a point and an image.
   * @param n2 Right hand side vector (3x1).
   * @param sparseNormals Reduced normal equations matrix.
   */
  bool BundleLidarMeasure::applyLidarRangeConstraint(LinearAlgebra::MatrixUpperTriangular &N22,
                                                     SparseBlockColumnMatrix &N12,
                                                     LinearAlgebra::Vector &n2,
                                                     LinearAlgebra::VectorCompressed &n1,
                                                     SparseBlockMatrix &sparseNormals) {
    int fred=1;

    // loop over simultaneous measures
//    for (int i = 0; i < m_simultaneousMeasures.size(); i++) {
//      m_lidarRangeConstraints.formConstraint(N22, N12, n2, n1, sparseNormals);
//    }

    return true;
  }
}

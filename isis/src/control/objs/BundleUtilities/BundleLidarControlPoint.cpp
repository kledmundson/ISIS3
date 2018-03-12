#include "BundleLidarControlPoint.h"

// Qt Library
#include <QDebug>

// Isis library
#include "BundleLidarRangeConstraint.h"
#include "ControlMeasure.h"
#include "Latitude.h"
#include "LidarControlPoint.h"
#include "Longitude.h"
#include "SpecialPixel.h"

namespace Isis {


  /**
   * Constructs a BundleLidarControlPoint object from a LidarControlPoint. Only the
   * non-ignored measures are added to the BundleLidarControlPoint.
   *  
   * @param controlPoint Pointer to a ControlPoint that will be used to 
   *                     construct this BundleLidarControlPoint.
   */
  BundleLidarControlPoint::BundleLidarControlPoint(LidarControlPoint *point)
    : BundleControlPoint() {

    // setup vector of BundleLidarMeasures for this control point
    int numMeasures = point->GetNumMeasures();
    for (int i = 0; i < numMeasures; i++) {
      ControlMeasure *controlMeasure = point->GetMeasure(i);
      if (controlMeasure->IsIgnored()) {
        continue;
      }

      addMeasure(controlMeasure);
    }
  }


  /**
   * Creates a BundleMeasure from the given ControlMeasure and appends it to
   * this BundleControlPoint's measure list.
   *
   * @param controlMeasure The ControlMeasure to be converted.
   *
   * @return BundleLidarMeasureQsp A pointer to the new BundleLidarMeasure.
   */
  BundleLidarMeasureQsp BundleLidarControlPoint::addMeasure(ControlMeasure *controlMeasure) {

    BundleLidarMeasureQsp bundleLidarMeasure
        = BundleLidarMeasureQsp( new BundleLidarMeasure(controlMeasure, this) );

    append(bundleLidarMeasure);

    return bundleLidarMeasure;
  }


  /**
   * Copy constructor. Constructs a BundleLidarControlPoint object from an existing
   * BundleLidarControlPoint.
   *
   * TODO: complete
   *  
   * @param src The BundleLidarControlPoint to be copied.
   */
//  BundleLidarControlPoint::BundleLidarControlPoint(const BundleLidarControlPoint &src) {
//    : BundleControlPoint(((ControlPoint) src) {
//    copy(src);
//  }


  /**
   * Destructor for BundleLidarControlPoint.
   */
  BundleLidarControlPoint::~BundleLidarControlPoint() {
  }


  /**
   * Copies given BundleLidarControlPoint to this BundleLidarControlPoint.
   *  
   * @param src The BundleLidarControlPoint to be copied.
   */
  void BundleLidarControlPoint::copy(const BundleLidarControlPoint &src) {

    BundleControlPoint::copy((BundleControlPoint) src);
/*
    // sanity check
    clear();

    m_controlPoint = src.m_controlPoint;

    int numMeasures = src.size();

    for (int i = 0; i < numMeasures; i++)
      append(BundleMeasureQsp( new BundleMeasure(*(src.at(i))) ));

    m_corrections = src.m_corrections;
    m_aprioriSigmas = src.m_aprioriSigmas;
    m_adjustedSigmas = src.m_adjustedSigmas;
    m_weights = src.m_weights;
    m_nicVector = src.m_nicVector;
*/
  }


  /**
   * Adds simultaneous measure to this BundleLidarControlPoint.
   *
   * @param measure QSharedPointer to measure from image acquired simultaneously with lidar point.
   */
  void BundleLidarControlPoint::addSimultaneousMeasure(BundleMeasureQsp measure) {
    m_simultaneousMeasures.append(measure);

    // create new lidar range constraint and add to vector
    BundleLidarRangeConstraintQsp constraintqsp
        = BundleLidarRangeConstraintQsp(new BundleLidarRangeConstraint(measure));
    m_lidarRangeConstraints.append(constraintqsp);
  }


  /**
   * Checks if QVector of any simultaneous measures have the input serial number.
   *
   * @param serialNumberString Serial number to check for.
   */
  bool BundleLidarControlPoint::containsSerialNumber(QString serialNumberString) {
    for (int i = 0; i < m_simultaneousMeasures.size(); i++) {
      if (m_simultaneousMeasures.at(i)->cubeSerialNumber() == serialNumberString) {
        return true;
      }
    }

    return false;
  }


  /**
   * Returns number of image measures simultaneous to lidar observation.
   *
   * @return int Number of image measures simultaneous to lidar observation.
   */
  int BundleLidarControlPoint::numberSimultaneousMeasures() {
    return m_simultaneousMeasures.size();
  }


  /**
   * Returns range of raw lidar control point.
   *
   * @return double Range of raw lidar control point.
   */
  double BundleLidarControlPoint::range() {
    return ((LidarControlPoint*)rawControlPoint())->range();
  }


  /**
   * Returns a priori range sigma of raw lidar control point.
   *
   * @return double a priori range sigma of raw lidar control point.
   */
  double BundleLidarControlPoint::rangeSigma() {
    return ((LidarControlPoint*)rawControlPoint())->sigmaRange();
  }
}
